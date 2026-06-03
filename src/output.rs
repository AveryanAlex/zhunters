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

struct ScoredBlock<'a> {
    index: usize,
    records: Vec<ZScoreRecord<'a>>,
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
    let (sender, receiver) = mpsc::sync_channel::<ScoredBlock<'_>>(max_in_flight);

    let writer_result = rayon::scope(move |scope| -> io::Result<()> {
        let mut next_to_schedule = 0;
        let mut completed = 0;
        let mut next_to_write = 0;
        let mut buffered = BTreeMap::<usize, Vec<ZScoreRecord>>::new();
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
            let Ok(block) = receiver.recv() else {
                break;
            };
            completed += 1;

            if writer_result.is_ok() {
                buffered.insert(block.index, block.records);
                debug_assert!(buffered.len() <= max_in_flight);

                while let Some(records) = buffered.remove(&next_to_write) {
                    for record in records {
                        if let Err(error) = writeln!(writer, "{record}") {
                            writer_result = Err(error);
                            break;
                        }
                    }
                    if writer_result.is_err() {
                        break;
                    }

                    positions_since_flush += block_len(next_to_write, block_size, positions);
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
        let block_start = block_index * block_size;
        let block_end = (block_start + block_size).min(positions);
        let records = score_position_block(
            config,
            plan.sequence,
            block_start,
            block_end,
            plan.range,
            progress,
        );

        for record in records {
            writeln!(writer, "{record}")?;
        }

        positions_since_flush += block_len(block_index, block_size, positions);
        if positions_since_flush >= flush_positions || block_index + 1 == block_count {
            writer.flush()?;
            positions_since_flush = 0;
        }
    }

    Ok(())
}

fn spawn_score_block<'scope, 'a: 'scope>(
    scope: &rayon::Scope<'scope>,
    sender: mpsc::SyncSender<ScoredBlock<'a>>,
    config: &'scope ZhuntConfig,
    plan: StreamingPlan<'a>,
    block_size: usize,
    block_index: usize,
    progress: &'scope (impl Fn(usize) + Send + Sync),
) {
    scope.spawn(move |_| {
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
        let _ = sender.send(ScoredBlock {
            index: block_index,
            records,
        });
    });
}

fn block_len(block_index: usize, block_size: usize, positions: usize) -> usize {
    let block_start = block_index * block_size;
    (positions - block_start).min(block_size)
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
