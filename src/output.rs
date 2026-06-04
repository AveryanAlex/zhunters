use crate::config::ZhuntConfig;
use crate::constants::STREAM_CHUNK_POSITIONS;
use crate::record::{ZScoreRecord, ZScoreSummary};
use crate::scoring::{
    install_scoring_pool, normalized_dinucleotide_range, score_position_block, DinucleotideRange,
    SCORE_WORK_BLOCK_POSITIONS,
};
use crate::sequence::CircularSequence;
use crate::ZhuntError;
use std::collections::BTreeMap;
use std::fs;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc;

const WRITER_BUFFERED_CHUNKS: usize = 2;

struct ScoredBlock {
    index: usize,
    bytes: Vec<u8>,
    positions: usize,
}

#[derive(Debug, Clone, Copy)]
struct StreamingPlan<'a> {
    input_name: &'a str,
    sequence: &'a CircularSequence,
    range: DinucleotideRange,
    chunk_size: usize,
    threads: Option<usize>,
}

#[derive(Debug, Clone, Copy)]
pub struct ZScoreRunOptions<'a> {
    pub min_size: usize,
    pub max_size: usize,
    pub input_name: &'a str,
    pub threads: Option<usize>,
}

pub fn write_zscore_file_streaming<P: AsRef<Path>>(
    output_path: P,
    config: &ZhuntConfig,
    min_size: usize,
    max_size: usize,
    input_name: &str,
    sequence: &CircularSequence,
    progress: impl Fn(usize) + Send + Sync,
) -> Result<ZScoreSummary, ZhuntError> {
    write_zscore_file_streaming_with_options(
        output_path,
        config,
        ZScoreRunOptions {
            min_size,
            max_size,
            input_name,
            threads: None,
        },
        sequence,
        progress,
    )
}

pub fn write_zscore_file_streaming_with_options<P: AsRef<Path>>(
    output_path: P,
    config: &ZhuntConfig,
    options: ZScoreRunOptions<'_>,
    sequence: &CircularSequence,
    progress: impl Fn(usize) + Send + Sync,
) -> Result<ZScoreSummary, ZhuntError> {
    let range = normalized_dinucleotide_range(config, options.min_size, options.max_size)?;
    let file = fs::File::create(output_path)?;
    let mut writer = BufWriter::new(file);
    write_zscore_streaming_with_plan(
        &mut writer,
        config,
        StreamingPlan {
            input_name: options.input_name,
            sequence,
            range,
            chunk_size: STREAM_CHUNK_POSITIONS,
            threads: options.threads,
        },
        progress,
    )
}

pub fn write_zscore_streaming<W: Write + Send>(
    writer: &mut W,
    config: &ZhuntConfig,
    min_size: usize,
    max_size: usize,
    input_name: &str,
    sequence: &CircularSequence,
    progress: impl Fn(usize) + Send + Sync,
) -> Result<ZScoreSummary, ZhuntError> {
    write_zscore_streaming_with_options(
        writer,
        config,
        ZScoreRunOptions {
            min_size,
            max_size,
            input_name,
            threads: None,
        },
        sequence,
        progress,
    )
}

pub fn write_zscore_streaming_with_options<W: Write + Send>(
    writer: &mut W,
    config: &ZhuntConfig,
    options: ZScoreRunOptions<'_>,
    sequence: &CircularSequence,
    progress: impl Fn(usize) + Send + Sync,
) -> Result<ZScoreSummary, ZhuntError> {
    let range = normalized_dinucleotide_range(config, options.min_size, options.max_size)?;
    write_zscore_streaming_with_plan(
        writer,
        config,
        StreamingPlan {
            input_name: options.input_name,
            sequence,
            range,
            chunk_size: STREAM_CHUNK_POSITIONS,
            threads: options.threads,
        },
        progress,
    )
}

#[cfg(test)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_zscore_streaming_with_chunk_size<W: Write + Send>(
    writer: &mut W,
    config: &ZhuntConfig,
    range: DinucleotideRange,
    input_name: &str,
    sequence: &CircularSequence,
    chunk_size: usize,
    threads: Option<usize>,
    progress: impl Fn(usize) + Send + Sync,
) -> Result<ZScoreSummary, ZhuntError> {
    write_zscore_streaming_with_plan(
        writer,
        config,
        StreamingPlan {
            input_name,
            sequence,
            range,
            chunk_size,
            threads,
        },
        progress,
    )
}

fn write_zscore_streaming_with_plan<W: Write + Send>(
    writer: &mut W,
    config: &ZhuntConfig,
    plan: StreamingPlan<'_>,
    progress: impl Fn(usize) + Send + Sync,
) -> Result<ZScoreSummary, ZhuntError> {
    install_scoring_pool(plan.threads, || {
        write_zscore_streaming_with_plan_inner(writer, config, plan, progress)
    })?
}

fn write_zscore_streaming_with_plan_inner<W: Write + Send>(
    writer: &mut W,
    config: &ZhuntConfig,
    plan: StreamingPlan<'_>,
    progress: impl Fn(usize) + Send + Sync,
) -> Result<ZScoreSummary, ZhuntError> {
    write_zscore_header(
        writer,
        plan.input_name,
        plan.sequence.len(),
        plan.range.from,
        plan.range.to,
    )?;

    write_score_blocks_streaming(writer, config, plan, &progress)?;

    Ok(ZScoreSummary {
        input_name: plan.input_name.to_owned(),
        sequence_length: plan.sequence.len(),
        from_dinucleotide: plan.range.from,
        to_dinucleotide: plan.range.to,
        records_written: plan.sequence.len(),
    })
}

fn write_score_blocks_streaming<W: Write + Send>(
    writer: &mut W,
    config: &ZhuntConfig,
    plan: StreamingPlan<'_>,
    progress: &(impl Fn(usize) + Send + Sync),
) -> Result<(), ZhuntError> {
    let positions = plan.sequence.len();
    if positions == 0 {
        return Ok(());
    }

    let flush_positions = plan.chunk_size.max(1);
    let block_size = flush_positions.min(SCORE_WORK_BLOCK_POSITIONS);
    let block_count = positions.div_ceil(block_size);
    if rayon::current_num_threads() == 1 {
        return write_score_blocks_sequential(writer, config, plan, progress, block_size);
    }

    let max_in_flight = (WRITER_BUFFERED_CHUNKS * rayon::current_num_threads())
        .max(1)
        .min(block_count);
    let (sender, receiver) = mpsc::sync_channel::<ScoredBlock>(max_in_flight);

    let writer_result = rayon::scope(move |scope| -> io::Result<()> {
        let mut next_to_schedule = 0;
        let mut completed = 0;
        let mut next_to_write = 0;
        let mut buffered = BTreeMap::<usize, ScoredBlock>::new();
        let mut positions_since_flush = 0;
        let mut writer_result = Ok(());

        while next_to_schedule < max_in_flight {
            spawn_score_block(
                scope,
                sender.clone(),
                config,
                plan,
                block_size,
                next_to_schedule,
                progress,
            );
            next_to_schedule += 1;
        }
        while completed < block_count {
            let block = match receiver.recv() {
                Ok(block) => block,
                Err(_) => {
                    writer_result = Err(io::Error::new(
                        io::ErrorKind::BrokenPipe,
                        "scoring worker disconnected before all blocks were written",
                    ));
                    break;
                }
            };
            completed += 1;

            if writer_result.is_ok() {
                buffered.insert(block.index, block);
                debug_assert!(buffered.len() <= max_in_flight);

                while let Some(block) = buffered.remove(&next_to_write) {
                    if let Err(error) = writer.write_all(&block.bytes) {
                        writer_result = Err(error);
                        break;
                    }

                    positions_since_flush += block.positions;
                    next_to_write += 1;
                    if positions_since_flush >= flush_positions || next_to_write == block_count {
                        if let Err(error) = writer.flush() {
                            writer_result = Err(error);
                            break;
                        }
                        positions_since_flush = 0;
                    }
                }
            }

            if next_to_schedule < block_count {
                spawn_score_block(
                    scope,
                    sender.clone(),
                    config,
                    plan,
                    block_size,
                    next_to_schedule,
                    progress,
                );
                next_to_schedule += 1;
            }
        }

        writer_result
    });

    writer_result.map_err(ZhuntError::Io)
}

fn write_score_blocks_sequential<W: Write + Send>(
    writer: &mut W,
    config: &ZhuntConfig,
    plan: StreamingPlan<'_>,
    progress: &(impl Fn(usize) + Send + Sync),
    block_size: usize,
) -> Result<(), ZhuntError> {
    let positions = plan.sequence.len();
    let flush_positions = plan.chunk_size.max(1);
    let block_count = positions.div_ceil(block_size);
    let mut positions_since_flush = 0;

    for block_index in 0..block_count {
        let block = score_and_format_block(config, plan, block_size, block_index, progress);
        writer.write_all(&block.bytes)?;

        positions_since_flush += block.positions;
        if positions_since_flush >= flush_positions || block_index + 1 == block_count {
            writer.flush()?;
            positions_since_flush = 0;
        }
    }

    Ok(())
}

fn spawn_score_block<'scope, 'a: 'scope>(
    scope: &rayon::Scope<'scope>,
    sender: mpsc::SyncSender<ScoredBlock>,
    config: &'scope ZhuntConfig,
    plan: StreamingPlan<'a>,
    block_size: usize,
    block_index: usize,
    progress: &'scope (impl Fn(usize) + Send + Sync),
) {
    scope.spawn(move |_| {
        let block = score_and_format_block(config, plan, block_size, block_index, progress);
        let _ = sender.send(block);
    });
}

fn score_and_format_block(
    config: &ZhuntConfig,
    plan: StreamingPlan<'_>,
    block_size: usize,
    block_index: usize,
    progress: &(impl Fn(usize) + Send + Sync),
) -> ScoredBlock {
    let block_start = block_index * block_size;
    let block_end = (block_start + block_size).min(plan.sequence.len());
    let records = score_position_block(
        config,
        plan.sequence,
        block_start,
        block_end,
        plan.range,
        progress,
    );
    let positions = block_end - block_start;
    let bytes = format_records(records);

    ScoredBlock {
        index: block_index,
        bytes,
        positions,
    }
}

fn format_records(records: Vec<ZScoreRecord<'_>>) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(records.len() * 96);
    for record in records {
        format_record(&record, &mut bytes);
    }
    bytes
}

fn format_record(record: &ZScoreRecord<'_>, bytes: &mut Vec<u8>) {
    push_usize(bytes, record.start);
    bytes.push(b' ');
    push_usize(bytes, record.start + record.length);
    bytes.push(b' ');
    push_usize(bytes, record.length);
    bytes.push(b' ');
    push_fixed_3_width_7(bytes, record.delta_linking);
    bytes.push(b' ');
    push_fixed_3_width_7(bytes, record.slope);
    bytes.push(b' ');
    push_scientific(bytes, record.probability);
    bytes.push(b' ');
    bytes.extend_from_slice(record.sequence.as_ref());
    bytes.extend_from_slice(b"   ");
    push_antisyn(bytes, record.antisyn.mask, record.antisyn.dinucleotides);
    bytes.push(b'\n');
}

fn push_usize(bytes: &mut Vec<u8>, value: usize) {
    push_u64(bytes, value as u64);
}

fn push_u64(bytes: &mut Vec<u8>, mut value: u64) {
    let mut buffer = [0_u8; 20];
    let mut index = buffer.len();
    loop {
        index -= 1;
        buffer[index] = b'0' + (value % 10) as u8;
        value /= 10;
        if value == 0 {
            break;
        }
    }
    bytes.extend_from_slice(&buffer[index..]);
}

fn push_i32_zero_padded_2(bytes: &mut Vec<u8>, value: i32) {
    let magnitude = value.unsigned_abs();
    bytes.push(if value < 0 { b'-' } else { b'+' });
    if magnitude < 10 {
        bytes.push(b'0');
        bytes.push(b'0' + magnitude as u8);
    } else {
        push_u64(bytes, magnitude as u64);
    }
}

fn push_fixed_3_width_7(bytes: &mut Vec<u8>, value: f64) {
    let mut buffer = [0_u8; 32];
    let mut len = 0;

    let negative = value.is_sign_negative();
    let scaled = round_ties_even_to_u64(value.abs() * 1_000.0);
    if negative {
        buffer[len] = b'-';
        len += 1;
    }

    len += write_u64_to_buffer(&mut buffer[len..], scaled / 1_000);
    buffer[len] = b'.';
    len += 1;
    let fraction = scaled % 1_000;
    buffer[len] = b'0' + (fraction / 100) as u8;
    buffer[len + 1] = b'0' + ((fraction / 10) % 10) as u8;
    buffer[len + 2] = b'0' + (fraction % 10) as u8;
    len += 3;

    for _ in len..7 {
        bytes.push(b' ');
    }
    bytes.extend_from_slice(&buffer[..len]);
}

fn push_fixed_6(bytes: &mut Vec<u8>, value: f64) {
    let scaled = round_ties_even_to_u64(value * 1_000_000.0);
    push_u64(bytes, scaled / 1_000_000);
    bytes.push(b'.');
    let fraction = scaled % 1_000_000;
    bytes.push(b'0' + (fraction / 100_000) as u8);
    bytes.push(b'0' + ((fraction / 10_000) % 10) as u8);
    bytes.push(b'0' + ((fraction / 1_000) % 10) as u8);
    bytes.push(b'0' + ((fraction / 100) % 10) as u8);
    bytes.push(b'0' + ((fraction / 10) % 10) as u8);
    bytes.push(b'0' + (fraction % 10) as u8);
}

fn push_scientific(bytes: &mut Vec<u8>, value: f64) {
    if value == 0.0 {
        bytes.extend_from_slice(b"0.000000e+00");
        return;
    }

    if !value.is_finite() {
        write!(bytes, "{value:.6e}").expect("writing formatted records to a Vec cannot fail");
        return;
    }

    if value.is_sign_negative() {
        bytes.push(b'-');
    }
    let value = value.abs();
    let mut exponent = value.log10().floor() as i32;
    let mut mantissa = value / 10_f64.powi(exponent);

    if (mantissa * 1_000_000.0).round() >= 10_000_000.0 {
        mantissa /= 10.0;
        exponent += 1;
    }

    push_fixed_6(bytes, mantissa);
    bytes.push(b'e');
    push_i32_zero_padded_2(bytes, exponent);
}

fn push_antisyn(bytes: &mut Vec<u8>, mask: u32, dinucleotides: u8) {
    let dinucleotides = usize::from(dinucleotides);
    for position in 0..dinucleotides {
        let shift = dinucleotides - 1 - position;
        if (mask >> shift) & 1 == 0 {
            bytes.extend_from_slice(b"AS");
        } else {
            bytes.extend_from_slice(b"SA");
        }
    }
}

fn round_ties_even_to_u64(value: f64) -> u64 {
    let floor = value.floor();
    let fraction = value - floor;
    const TIE_EPSILON: f64 = 1.0e-9;
    let rounded = if fraction > 0.5 + TIE_EPSILON {
        floor + 1.0
    } else if fraction < 0.5 - TIE_EPSILON || (floor as u64).is_multiple_of(2) {
        floor
    } else {
        floor + 1.0
    };
    rounded as u64
}

fn write_u64_to_buffer(buffer: &mut [u8], mut value: u64) -> usize {
    let mut digits = [0_u8; 20];
    let mut index = digits.len();
    loop {
        index -= 1;
        digits[index] = b'0' + (value % 10) as u8;
        value /= 10;
        if value == 0 {
            break;
        }
    }
    let len = digits.len() - index;
    buffer[..len].copy_from_slice(&digits[index..]);
    len
}

pub(crate) fn write_zscore_header<W: Write + ?Sized>(
    writer: &mut W,
    input_name: &str,
    sequence_length: usize,
    from_dinucleotide: usize,
    to_dinucleotide: usize,
) -> io::Result<()> {
    writeln!(
        writer,
        "{} {} {} {}",
        input_name, sequence_length, from_dinucleotide, to_dinucleotide
    )
}

pub fn zscore_output_path(input_name: &str) -> PathBuf {
    PathBuf::from(format!("{input_name}.Z-SCORE"))
}
