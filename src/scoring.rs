use crate::config::ZhuntConfig;
use crate::constants::{A, EXP_LIMIT, INT_DBZED, K_RT, SIGMA};
use crate::record::{AntiSynPath, ZScoreOutput, ZScoreRecord};
use crate::sequence::CircularSequence;
use crate::ZhuntError;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use rustc_hash::{FxBuildHasher, FxHashMap};
use std::io;

const PROGRESS_BATCH_SIZE: usize = 1_000;
pub(crate) const SCORE_WORK_BLOCK_POSITIONS: usize = 8_192;

pub fn calculate_zscore<'a>(
    config: &ZhuntConfig,
    min_size: usize,
    max_size: usize,
    input_name: String,
    sequence: &'a CircularSequence,
) -> Result<ZScoreOutput<'a>, ZhuntError> {
    calculate_zscore_with_progress(config, min_size, max_size, input_name, sequence, |_| {})
}

pub fn calculate_zscore_with_progress<'a>(
    config: &ZhuntConfig,
    min_size: usize,
    max_size: usize,
    input_name: String,
    sequence: &'a CircularSequence,
    progress: impl Fn(usize) + Send + Sync,
) -> Result<ZScoreOutput<'a>, ZhuntError> {
    calculate_zscore_with_progress_and_threads(
        config, min_size, max_size, input_name, sequence, None, progress,
    )
}

pub fn calculate_zscore_with_progress_and_threads<'a>(
    config: &ZhuntConfig,
    min_size: usize,
    max_size: usize,
    input_name: String,
    sequence: &'a CircularSequence,
    threads: Option<usize>,
    progress: impl Fn(usize) + Send + Sync,
) -> Result<ZScoreOutput<'a>, ZhuntError> {
    let range = normalized_dinucleotide_range(config, min_size, max_size)?;

    let records = install_scoring_pool(threads, || {
        score_positions_parallel(config, sequence, 0, sequence.len(), range, &progress)
    })?;

    Ok(ZScoreOutput {
        input_name,
        sequence_length: sequence.len(),
        from_dinucleotide: range.from,
        to_dinucleotide: range.to,
        records,
    })
}

pub(crate) fn score_positions_parallel<'a>(
    config: &ZhuntConfig,
    sequence: &'a CircularSequence,
    start: usize,
    end: usize,
    range: DinucleotideRange,
    progress: &(impl Fn(usize) + Send + Sync),
) -> Vec<ZScoreRecord<'a>> {
    debug_assert!(start <= end);
    let positions = end - start;
    if positions == 0 {
        return Vec::new();
    }

    let block_count = positions.div_ceil(SCORE_WORK_BLOCK_POSITIONS);
    let block_records = (0..block_count)
        .into_par_iter()
        .map(|block_index| {
            let block_start = start + block_index * SCORE_WORK_BLOCK_POSITIONS;
            let block_end = (block_start + SCORE_WORK_BLOCK_POSITIONS).min(end);
            score_position_block(config, sequence, block_start, block_end, range, progress)
        })
        .collect::<Vec<_>>();

    let mut records = Vec::with_capacity(positions);
    for block in block_records {
        records.extend(block);
    }
    records
}

pub(crate) fn score_position_block<'a>(
    config: &ZhuntConfig,
    sequence: &'a CircularSequence,
    start: usize,
    end: usize,
    range: DinucleotideRange,
    progress: &(impl Fn(usize) + Send + Sync),
) -> Vec<ZScoreRecord<'a>> {
    debug_assert!(start <= end);

    let mut records = Vec::with_capacity(end - start);
    let mut scratch = ScoringScratch::new(range.to);
    let mut pending_progress = 0;

    for index in start..end {
        records.push(score_position_with_scratch(
            config,
            sequence,
            index,
            range.from,
            range.to,
            &mut scratch,
        ));
        pending_progress += 1;
        if pending_progress == PROGRESS_BATCH_SIZE {
            progress(pending_progress);
            pending_progress = 0;
        }
    }

    if pending_progress > 0 {
        progress(pending_progress);
    }

    records
}

pub(crate) fn install_scoring_pool<R: Send>(
    threads: Option<usize>,
    operation: impl FnOnce() -> R + Send,
) -> Result<R, ZhuntError> {
    let threads = threads.filter(|threads| *threads > 0).unwrap_or_else(|| {
        let cores = std::thread::available_parallelism()
            .map(usize::from)
            .unwrap_or(1);
        if cores >= 8 {
            cores - 1
        } else {
            cores
        }
    });

    let pool = ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|error| ZhuntError::Io(io::Error::other(error)))?;
    Ok(pool.install(operation))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct DinucleotideRange {
    pub(crate) from: usize,
    pub(crate) to: usize,
}

pub(crate) fn normalized_dinucleotide_range(
    config: &ZhuntConfig,
    min_size: usize,
    max_size: usize,
) -> Result<DinucleotideRange, ZhuntError> {
    if min_size == 0 {
        return Err(ZhuntError::InvalidMinimumSize { min_size });
    }

    let to_dinucleotide = max_size.min(config.max_dinucleotides);
    let from_dinucleotide = min_size.min(to_dinucleotide);
    if from_dinucleotide == 0 {
        return Err(ZhuntError::InvalidMinimumSize { min_size });
    }

    Ok(DinucleotideRange {
        from: from_dinucleotide,
        to: to_dinucleotide,
    })
}

#[cfg(test)]
pub(crate) fn score_position<'a>(
    config: &ZhuntConfig,
    sequence: &'a CircularSequence,
    position: usize,
    from_dinucleotide: usize,
    to_dinucleotide: usize,
) -> ZScoreRecord<'a> {
    let mut scratch = ScoringScratch::new(to_dinucleotide);
    score_position_with_scratch(
        config,
        sequence,
        position,
        from_dinucleotide,
        to_dinucleotide,
        &mut scratch,
    )
}

fn score_position_with_scratch<'a>(
    config: &ZhuntConfig,
    sequence: &'a CircularSequence,
    position: usize,
    from_dinucleotide: usize,
    to_dinucleotide: usize,
    scratch: &mut ScoringScratch,
) -> ZScoreRecord<'a> {
    const PI_DEGREES: f64 = 180.0 / std::f64::consts::PI;

    let nucleotides = 2 * to_dinucleotide;
    let window = sequence.window(position, nucleotides);
    assign_bzenergy_indices_into(window.as_ref(), to_dinucleotide, &mut scratch.bzindex);
    best_anti_syn_prefix_paths_into(to_dinucleotide, scratch);

    let mut best: Option<(usize, f64, AntiSynPath)> = None;

    for dinucleotides in from_dinucleotide..=to_dinucleotide {
        let anti_syn_path = scratch.prefixes[dinucleotides - 1];
        materialize_bzenergy(
            anti_syn_path,
            &scratch.bzindex[..dinucleotides],
            &config.exp_dbzed,
            &mut scratch.bzenergy,
        );
        let delta_twist = (A / 2.0) * dinucleotides as f64;
        let delta_linking = {
            let ScoringScratch {
                bzenergy,
                products,
                logcoef,
                exponents,
                ..
            } = scratch;
            find_delta_linking_with_buffers(
                bzenergy,
                &config.bztwist[..dinucleotides],
                delta_twist,
                products,
                logcoef,
                exponents,
            )
        };

        let replace_best = match best.as_ref() {
            Some((_, best_dl, _)) => delta_linking < *best_dl,
            None => true,
        };
        if replace_best {
            best = Some((
                dinucleotides,
                delta_linking,
                AntiSynPath::new(anti_syn_path.mask, anti_syn_path.len),
            ));
        }
    }

    let (best_dinucleotides, best_dl, antisyn) =
        best.expect("from_dinucleotide..=to_dinucleotide is non-empty");
    let best_path = PathState {
        esum: 0,
        mask: antisyn.mask,
        len: antisyn.dinucleotides,
    };
    materialize_bzenergy(
        best_path,
        &scratch.bzindex[..best_dinucleotides],
        &config.exp_dbzed,
        &mut scratch.bzenergy,
    );
    delta_linking_logcoef_into(
        &scratch.bzenergy,
        &mut scratch.products,
        &mut scratch.logcoef,
    );
    let length = 2 * best_dinucleotides;
    let slope = delta_linking_slope_into(
        best_dl,
        &scratch.logcoef,
        &config.bztwist[..best_dinucleotides],
        &mut scratch.exponents,
    )
    .atan()
        * PI_DEGREES;
    let probability = assign_probability(best_dl);
    let sequence_segment = sequence.window(position, length);

    ZScoreRecord {
        start: position + 1,
        length,
        delta_linking: best_dl,
        slope,
        probability,
        sequence: sequence_segment,
        antisyn,
    }
}

struct ScoringScratch {
    bzindex: Vec<usize>,
    prefixes: Vec<PathState>,
    bzenergy: Vec<f64>,
    states: [Vec<PathState>; 2],
    next: [Vec<PathState>; 2],
    seen: [FxHashMap<i64, usize>; 2],
    products: Vec<f64>,
    logcoef: Vec<f64>,
    exponents: Vec<f64>,
}

impl ScoringScratch {
    fn new(max_dinucleotides: usize) -> Self {
        let max_frontier_capacity = anti_syn_frontier_capacity(max_dinucleotides);
        let max_hash_capacity = max_frontier_capacity.saturating_mul(2);
        Self {
            bzindex: Vec::with_capacity(max_dinucleotides),
            prefixes: Vec::with_capacity(max_dinucleotides),
            bzenergy: Vec::with_capacity(max_dinucleotides),
            states: [
                Vec::with_capacity(max_frontier_capacity),
                Vec::with_capacity(max_frontier_capacity),
            ],
            next: [
                Vec::with_capacity(max_frontier_capacity),
                Vec::with_capacity(max_frontier_capacity),
            ],
            seen: [
                FxHashMap::with_capacity_and_hasher(max_hash_capacity, FxBuildHasher),
                FxHashMap::with_capacity_and_hasher(max_hash_capacity, FxBuildHasher),
            ],
            products: Vec::with_capacity(max_dinucleotides),
            logcoef: Vec::with_capacity(max_dinucleotides),
            exponents: Vec::with_capacity(max_dinucleotides),
        }
    }
}

#[cfg(test)]
pub(crate) fn assign_bzenergy_indices(window: &[u8], dinucleotides: usize) -> Vec<usize> {
    let mut bzindex = Vec::with_capacity(dinucleotides);
    assign_bzenergy_indices_into(window, dinucleotides, &mut bzindex);
    bzindex
}

fn assign_bzenergy_indices_into(window: &[u8], dinucleotides: usize, bzindex: &mut Vec<usize>) {
    bzindex.clear();
    bzindex.extend(
        window
            .chunks_exact(2)
            .take(dinucleotides)
            .map(|pair| dinucleotide_index(pair[0], pair[1])),
    );
}

fn dinucleotide_index(first: u8, second: u8) -> usize {
    match (base_number(first), base_number(second)) {
        (Some(first), Some(second)) => first * 4 + second,
        _ => 16,
    }
}

fn base_number(base: u8) -> Option<usize> {
    match base {
        b'a' => Some(0),
        b't' => Some(1),
        b'g' => Some(2),
        b'c' => Some(3),
        b'n' => None,
        _ => None,
    }
}

#[derive(Debug, Clone, PartialEq)]
#[cfg(test)]
pub(crate) struct BestAntiSyn {
    pub(crate) esum: i64,
    pub(crate) bzenergy: Vec<f64>,
    pub(crate) antisyn: AntiSynPath,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct PathState {
    esum: i64,
    mask: u32,
    len: u8,
}

impl PathState {
    fn push_with_esum(&self, state: u8, esum: i64) -> Self {
        Self {
            esum,
            mask: (self.mask << 1) | u32::from(state),
            len: self.len + 1,
        }
    }

    #[cfg(test)]
    fn into_best(self, bzindex: &[usize], exp_dbzed: &[[f64; 17]; 4]) -> BestAntiSyn {
        let mut bzenergy = Vec::with_capacity(usize::from(self.len));
        materialize_bzenergy(self, bzindex, exp_dbzed, &mut bzenergy);

        BestAntiSyn {
            esum: self.esum,
            bzenergy,
            antisyn: AntiSynPath::new(self.mask, self.len),
        }
    }
}

#[cfg(test)]
pub(crate) fn best_anti_syn(bzindex: &[usize], exp_dbzed: &[[f64; 17]; 4]) -> BestAntiSyn {
    best_anti_syn_path(bzindex).into_best(bzindex, exp_dbzed)
}

#[cfg(test)]
fn best_anti_syn_path(bzindex: &[usize]) -> PathState {
    let mut scratch = ScoringScratch::new(bzindex.len());
    scratch.bzindex.extend_from_slice(bzindex);
    best_anti_syn_prefix_paths(bzindex.len(), &mut scratch)
        .pop()
        .expect("at least one anti/syn path is reachable")
}

#[cfg(test)]
fn best_anti_syn_prefix_paths(
    max_dinucleotides: usize,
    scratch: &mut ScoringScratch,
) -> Vec<PathState> {
    best_anti_syn_prefix_paths_into(max_dinucleotides, scratch);
    scratch.prefixes.clone()
}

fn best_anti_syn_prefix_paths_into(max_dinucleotides: usize, scratch: &mut ScoringScratch) {
    let bzindex = &scratch.bzindex;
    assert!(!bzindex.is_empty(), "at least one dinucleotide is required");
    assert!(
        max_dinucleotides <= bzindex.len(),
        "max dinucleotides cannot exceed available energy indexes"
    );

    let first_index = bzindex[0];
    let first_as_esum = INT_DBZED[0][first_index];
    let first_sa_esum = INT_DBZED[1][first_index];
    for state in 0..=1 {
        scratch.states[state].clear();
        scratch.next[state].clear();
        scratch.seen[state].clear();
    }
    scratch.prefixes.clear();

    scratch.states[0].push(PathState {
        esum: first_as_esum,
        mask: 0,
        len: 1,
    });
    scratch.states[1].push(PathState {
        esum: first_sa_esum,
        mask: 1,
        len: 1,
    });
    scratch
        .prefixes
        .push(best_path_from_states(&scratch.states));

    for &dinucleotide_index in &bzindex[1..max_dinucleotides] {
        for state in 0..=1 {
            scratch.next[state].clear();
            scratch.seen[state].clear();
        }

        for (previous_state, paths) in scratch.states.iter().enumerate() {
            for path in paths {
                let as_row = transition_row(previous_state as u8, 0);
                let as_energy = INT_DBZED[as_row][dinucleotide_index];
                let as_esum = path.esum + as_energy;
                insert_distinct_energy_path(
                    &mut scratch.next[0],
                    &mut scratch.seen[0],
                    path.push_with_esum(0, as_esum),
                );

                let sa_row = transition_row(previous_state as u8, 1);
                let sa_esum = path.esum + INT_DBZED[sa_row][dinucleotide_index];
                insert_distinct_energy_path(
                    &mut scratch.next[1],
                    &mut scratch.seen[1],
                    path.push_with_esum(1, sa_esum),
                );
            }
        }

        std::mem::swap(&mut scratch.states, &mut scratch.next);
        scratch
            .prefixes
            .push(best_path_from_states(&scratch.states));
    }
}

fn best_path_from_states(states: &[Vec<PathState>; 2]) -> PathState {
    states
        .iter()
        .flatten()
        .copied()
        .min_by(path_order)
        .expect("at least one anti/syn path is reachable")
}

fn anti_syn_frontier_capacity(max_dinucleotides: usize) -> usize {
    if max_dinucleotides >= usize::BITS as usize - 1 {
        usize::MAX / 4
    } else {
        (1_usize << max_dinucleotides).max(2)
    }
}

fn materialize_bzenergy(
    path: PathState,
    bzindex: &[usize],
    exp_dbzed: &[[f64; 17]; 4],
    bzenergy: &mut Vec<f64>,
) {
    let len = usize::from(path.len);
    debug_assert_eq!(bzindex.len(), len);
    bzenergy.clear();
    if bzenergy.capacity() < len {
        bzenergy.reserve(len - bzenergy.capacity());
    }

    let mut previous_state = None;
    for position in 0..len {
        let shift = len - 1 - position;
        let state = ((path.mask >> shift) & 1) as u8;
        let row = previous_state.map_or(state as usize, |previous| transition_row(previous, state));
        bzenergy.push(exp_dbzed[row][bzindex[position]]);
        previous_state = Some(state);
    }
}

fn insert_distinct_energy_path(
    paths: &mut Vec<PathState>,
    seen: &mut FxHashMap<i64, usize>,
    candidate: PathState,
) {
    if let Some(&index) = seen.get(&candidate.esum) {
        if candidate.mask < paths[index].mask {
            paths[index] = candidate;
        }
    } else {
        seen.insert(candidate.esum, paths.len());
        paths.push(candidate);
    }
}

fn path_order(left: &PathState, right: &PathState) -> std::cmp::Ordering {
    left.esum
        .cmp(&right.esum)
        .then_with(|| left.mask.cmp(&right.mask))
}

pub(crate) fn transition_row(previous_state: u8, current_state: u8) -> usize {
    match (previous_state, current_state) {
        (0, 0) => 0,
        (1, 0) => 3,
        (0, 1) => 2,
        (1, 1) => 1,
        _ => unreachable!("anti/syn states are encoded as 0 or 1"),
    }
}

#[derive(Debug, Clone, PartialEq)]
#[cfg(test)]
pub(crate) struct DeltaResult {
    pub(crate) delta_linking: f64,
    pub(crate) logcoef: Vec<f64>,
}

#[cfg(test)]
pub(crate) fn find_delta_linking(
    best_bzenergy: &[f64],
    bztwist: &[f64],
    delta_twist: f64,
) -> DeltaResult {
    let mut products = Vec::new();
    let mut logcoef = Vec::new();
    let mut exponents = Vec::new();
    let delta_linking = find_delta_linking_with_buffers(
        best_bzenergy,
        bztwist,
        delta_twist,
        &mut products,
        &mut logcoef,
        &mut exponents,
    );

    DeltaResult {
        delta_linking,
        logcoef,
    }
}

fn find_delta_linking_with_buffers(
    best_bzenergy: &[f64],
    bztwist: &[f64],
    delta_twist: f64,
    products: &mut Vec<f64>,
    logcoef: &mut Vec<f64>,
    exponents: &mut Vec<f64>,
) -> f64 {
    delta_linking_logcoef_into(best_bzenergy, products, logcoef);

    linear_search(10.0, 50.0, 0.001, |dl| {
        delta_linking_equation_into(dl, logcoef, bztwist, delta_twist, exponents)
    })
}

fn delta_linking_logcoef_into(
    best_bzenergy: &[f64],
    products: &mut Vec<f64>,
    logcoef: &mut Vec<f64>,
) {
    let dinucleotides = best_bzenergy.len();
    products.clear();
    products.resize(dinucleotides, 1.0);
    logcoef.clear();
    logcoef.resize(dinucleotides, 0.0);

    for i in 0..dinucleotides {
        let mut sum = 0.0;
        for j in 0..(dinucleotides - i) {
            products[j] *= best_bzenergy[i + j];
            sum += products[j];
        }
        logcoef[i] = sum.ln();
    }
}

pub(crate) fn linear_search(
    x1: f64,
    x2: f64,
    tolerance: f64,
    mut func: impl FnMut(f64) -> f64,
) -> f64 {
    let f1 = func(x1);
    let f2 = func(x2);
    if f1 * f2 >= 0.0 {
        return x2;
    }

    let mut dx;
    let mut x;
    if f1 < 0.0 {
        dx = x2 - x1;
        x = x1;
    } else {
        dx = x1 - x2;
        x = x2;
    }

    loop {
        dx *= 0.5;
        let xmid = x + dx;
        let fmid = func(xmid);
        if fmid <= 0.0 {
            x = xmid;
        }
        if dx.abs() <= tolerance {
            return x;
        }
    }
}

fn delta_linking_equation_into(
    dl: f64,
    logcoef: &[f64],
    bztwist: &[f64],
    delta_twist: f64,
    exponents: &mut Vec<f64>,
) -> f64 {
    exponents.clear();
    exponents.extend(
        logcoef
            .iter()
            .zip(bztwist.iter())
            .map(|(logcoef, bztwist)| {
                let z = dl - bztwist;
                logcoef + K_RT * z * z
            }),
    );

    let expmini = exponent_offset(exponents);
    let mut sump = 0.0;
    let mut sumq = 0.0;

    for (exponent, twist) in exponents.iter().zip(bztwist.iter()) {
        let z = (exponent + expmini).exp();
        sumq += z;
        sump += twist * z;
    }
    sumq += (K_RT * dl * dl + SIGMA + expmini).exp();

    delta_twist - sump / sumq
}

fn delta_linking_slope_into(
    dl: f64,
    logcoef: &[f64],
    bztwist: &[f64],
    exponents: &mut Vec<f64>,
) -> f64 {
    exponents.clear();
    exponents.extend(
        logcoef
            .iter()
            .zip(bztwist.iter())
            .map(|(logcoef, bztwist)| {
                let z = dl - bztwist;
                logcoef + K_RT * z * z
            }),
    );

    let expmini = exponent_offset(exponents);
    let x = 2.0 * K_RT;
    let mut sump = 0.0;
    let mut sump1 = 0.0;
    let mut sumq = 0.0;
    let mut sumq1 = 0.0;

    for ((exponent, twist), bztwist_value) in
        exponents.iter().zip(bztwist.iter()).zip(bztwist.iter())
    {
        let z = dl - bztwist_value;
        let y = (exponent + expmini).exp();
        sumq += y;
        sump += twist * y;

        let y1 = y * z * x;
        sumq1 += y1;
        sump1 += twist * y1;
    }

    let y = (K_RT * dl * dl + SIGMA + expmini).exp();
    sumq += y;
    sumq1 += x * dl * y;

    (sump1 - sump * sumq1 / sumq) / sumq
}

fn exponent_offset(exponents: &[f64]) -> f64 {
    let expmini = exponents.iter().copied().fold(0.0_f64, f64::min);
    if expmini < EXP_LIMIT {
        EXP_LIMIT - expmini
    } else {
        0.0
    }
}

pub fn assign_probability(dl: f64) -> f64 {
    const AVERAGE: f64 = 29.653_713_5;
    const STDV: f64 = 2.719_97;
    const SQRT2_INV: f64 = std::f64::consts::FRAC_1_SQRT_2;
    const SQRTPI_INV: f64 = 0.564_189_583_546;

    let mut z = (dl - AVERAGE).abs() / STDV;
    let mut x = z * SQRT2_INV;
    let y = SQRTPI_INV * (-x * x).exp();
    z *= z;
    let mut k = 1.0;
    let mut sum = 0.0;

    loop {
        sum += x;
        k += 2.0;
        x *= z / k;
        if sum + x <= sum {
            break;
        }
    }

    let tail = 0.5 - y * sum;
    if dl > AVERAGE {
        tail
    } else {
        1.0 / tail
    }
}
