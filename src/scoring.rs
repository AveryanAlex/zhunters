use crate::config::ZhuntConfig;
use crate::constants::{A, EXP_LIMIT, INT_DBZED, K_RT, SIGMA};
use crate::record::{AntiSynPath, ZScoreOutput, ZScoreRecord};
use crate::sequence::CircularSequence;
use crate::ZhuntError;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::io;

const PROGRESS_BATCH_SIZE: usize = 1_000;
pub(crate) const SCORE_WORK_BLOCK_POSITIONS: usize = 8_192;

const DELTA_LINKING_MIN: f64 = 10.0;
const DELTA_LINKING_MAX: f64 = 50.0;
const DELTA_LINKING_TOLERANCE: f64 = 0.001;
const DELTA_LINKING_GRID_STEPS: usize = 65_536;
const DELTA_LINKING_GRID_STEP: f64 =
    (DELTA_LINKING_MAX - DELTA_LINKING_MIN) / DELTA_LINKING_GRID_STEPS as f64;
const DELTA_LINKING_NEWTON_STEPS: usize = 6;
const DELTA_LINKING_GRID_ADJUST_LIMIT: usize = 8;

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
        let best_delta_linking = best.as_ref().map(|(_, delta_linking, _)| *delta_linking);
        let Some(delta_linking) = ({
            let ScoringScratch {
                bzenergy,
                products,
                logcoef,
                exponents,
                ..
            } = scratch;
            find_delta_linking_candidate_with_bound(
                bzenergy,
                &config.bztwist[..dinucleotides],
                delta_twist,
                best_delta_linking,
                products,
                logcoef,
                exponents,
            )
        }) else {
            continue;
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
    products: Vec<f64>,
    logcoef: Vec<f64>,
    exponents: Vec<f64>,
}

impl ScoringScratch {
    fn new(max_dinucleotides: usize) -> Self {
        Self {
            bzindex: Vec::with_capacity(max_dinucleotides),
            prefixes: Vec::with_capacity(max_dinucleotides),
            bzenergy: Vec::with_capacity(max_dinucleotides),
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
    scratch.prefixes.clear();

    let mut best_as = PathState {
        esum: INT_DBZED[0][first_index],
        mask: 0,
        len: 1,
    };
    let mut best_sa = PathState {
        esum: INT_DBZED[1][first_index],
        mask: 1,
        len: 1,
    };
    scratch.prefixes.push(best_path(best_as, best_sa));

    // Anti/syn energy is first-order Markov: once the current terminal state is
    // fixed, only that state and the accumulated energy can affect all future
    // transitions. Keep the best prefix ending in AS and the best prefix ending
    // in SA instead of carrying every distinct energy frontier.
    for &dinucleotide_index in &bzindex[1..max_dinucleotides] {
        let next_as = best_path(
            best_as.push_with_esum(0, best_as.esum + INT_DBZED[0][dinucleotide_index]),
            best_sa.push_with_esum(0, best_sa.esum + INT_DBZED[3][dinucleotide_index]),
        );
        let next_sa = best_path(
            best_as.push_with_esum(1, best_as.esum + INT_DBZED[2][dinucleotide_index]),
            best_sa.push_with_esum(1, best_sa.esum + INT_DBZED[1][dinucleotide_index]),
        );

        best_as = next_as;
        best_sa = next_sa;
        scratch.prefixes.push(best_path(best_as, best_sa));
    }
}

fn best_path(left: PathState, right: PathState) -> PathState {
    if path_order(&left, &right).is_le() {
        left
    } else {
        right
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

#[cfg(test)]
fn find_delta_linking_with_buffers(
    best_bzenergy: &[f64],
    bztwist: &[f64],
    delta_twist: f64,
    products: &mut Vec<f64>,
    logcoef: &mut Vec<f64>,
    exponents: &mut Vec<f64>,
) -> f64 {
    find_delta_linking_candidate_with_bound(
        best_bzenergy,
        bztwist,
        delta_twist,
        None,
        products,
        logcoef,
        exponents,
    )
    .expect("unbounded delta-linking search always returns a candidate")
}

fn find_delta_linking_candidate_with_bound(
    best_bzenergy: &[f64],
    bztwist: &[f64],
    delta_twist: f64,
    current_best: Option<f64>,
    products: &mut Vec<f64>,
    logcoef: &mut Vec<f64>,
    exponents: &mut Vec<f64>,
) -> Option<f64> {
    delta_linking_logcoef_into(best_bzenergy, products, logcoef);

    if let Some(best_dl) = current_best {
        // F(dl) = delta_twist - weighted_average_twist(dl) is monotone
        // decreasing. If F(current best) is still positive, this candidate's
        // root must be larger than the current best and therefore cannot win.
        if delta_linking_equation_into(best_dl, logcoef, bztwist, delta_twist, exponents) > 0.0 {
            return None;
        }
    }

    Some(delta_linking_search(
        logcoef,
        bztwist,
        delta_twist,
        exponents,
    ))
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

fn delta_linking_search(
    logcoef: &[f64],
    bztwist: &[f64],
    delta_twist: f64,
    exponents: &mut Vec<f64>,
) -> f64 {
    let f_min =
        delta_linking_equation_into(DELTA_LINKING_MIN, logcoef, bztwist, delta_twist, exponents);
    let f_max =
        delta_linking_equation_into(DELTA_LINKING_MAX, logcoef, bztwist, delta_twist, exponents);
    if f_min * f_max >= 0.0 {
        return DELTA_LINKING_MAX;
    }

    // The Z-HUNT delta-linking equation is monotone decreasing for the normal
    // search interval. Use Newton steps only as a root predictor, then verify
    // the final legacy bisection-grid point with the original equation so the
    // returned delta-linking value stays bit-for-bit compatible with bisection.
    if f_min > 0.0 && f_max < 0.0 {
        let root_estimate =
            approximate_delta_linking_root(f_min, f_max, logcoef, bztwist, delta_twist, exponents);
        if let Some(delta_linking) = snap_delta_linking_to_legacy_grid(
            root_estimate,
            logcoef,
            bztwist,
            delta_twist,
            exponents,
        ) {
            return delta_linking;
        }
    }

    linear_search(
        DELTA_LINKING_MIN,
        DELTA_LINKING_MAX,
        DELTA_LINKING_TOLERANCE,
        |dl| delta_linking_equation_into(dl, logcoef, bztwist, delta_twist, exponents),
    )
}

fn approximate_delta_linking_root(
    f_min: f64,
    f_max: f64,
    logcoef: &[f64],
    bztwist: &[f64],
    delta_twist: f64,
    exponents: &mut Vec<f64>,
) -> f64 {
    let mut lower = DELTA_LINKING_MIN;
    let mut upper = DELTA_LINKING_MAX;
    let mut x = lower + (upper - lower) * f_min / (f_min - f_max);
    if !x.is_finite() || x <= lower || x >= upper {
        x = (lower + upper) * 0.5;
    }

    for _ in 0..DELTA_LINKING_NEWTON_STEPS {
        let (value, derivative) =
            delta_linking_equation_and_derivative_into(x, logcoef, bztwist, delta_twist, exponents);
        if value <= 0.0 {
            upper = x;
        } else {
            lower = x;
        }
        if upper - lower <= DELTA_LINKING_GRID_STEP {
            return (lower + upper) * 0.5;
        }

        let newton = x - value / derivative;
        x = if derivative.is_finite()
            && derivative < 0.0
            && newton.is_finite()
            && newton > lower
            && newton < upper
        {
            newton
        } else {
            (lower + upper) * 0.5
        };
    }

    x
}

fn snap_delta_linking_to_legacy_grid(
    root_estimate: f64,
    logcoef: &[f64],
    bztwist: &[f64],
    delta_twist: f64,
    exponents: &mut Vec<f64>,
) -> Option<f64> {
    if !root_estimate.is_finite() {
        return None;
    }

    let mut index = ((root_estimate - DELTA_LINKING_MIN) / DELTA_LINKING_GRID_STEP).ceil() as isize;
    index = index.clamp(0, DELTA_LINKING_GRID_STEPS as isize);

    for _ in 0..=DELTA_LINKING_GRID_ADJUST_LIMIT {
        let delta_linking = delta_linking_grid_value(index as usize);
        let value =
            delta_linking_equation_into(delta_linking, logcoef, bztwist, delta_twist, exponents);
        if value > 0.0 {
            if index == DELTA_LINKING_GRID_STEPS as isize {
                return None;
            }
            index += 1;
            continue;
        }

        if index == 0 {
            return Some(delta_linking);
        }

        let previous = delta_linking_grid_value(index as usize - 1);
        let previous_value =
            delta_linking_equation_into(previous, logcoef, bztwist, delta_twist, exponents);
        if previous_value <= 0.0 {
            index -= 1;
            continue;
        }

        return Some(delta_linking);
    }

    None
}

fn delta_linking_grid_value(index: usize) -> f64 {
    DELTA_LINKING_MIN + DELTA_LINKING_GRID_STEP * index as f64
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

fn delta_linking_equation_and_derivative_into(
    dl: f64,
    logcoef: &[f64],
    bztwist: &[f64],
    delta_twist: f64,
    exponents: &mut Vec<f64>,
) -> (f64, f64) {
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
    let mut sum_twist_squared = 0.0;

    for (exponent, twist) in exponents.iter().zip(bztwist.iter()) {
        let weight = (exponent + expmini).exp();
        sumq += weight;
        sump += twist * weight;
        sum_twist_squared += twist * twist * weight;
    }
    sumq += (K_RT * dl * dl + SIGMA + expmini).exp();

    let mean = sump / sumq;
    let variance = (sum_twist_squared / sumq - mean * mean).max(0.0);
    (delta_twist - mean, 2.0 * K_RT * variance)
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
