use crate::config::ZhuntConfig;
use crate::constants::STREAM_CHUNK_POSITIONS;
use crate::record::{ZScoreRecord, ZScoreSummary};
use crate::scoring::{normalized_dinucleotide_range, score_positions_parallel, DinucleotideRange};
use crate::sequence::CircularSequence;
use crate::ZhuntError;
use std::fs;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc;
use std::thread;

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
    progress: impl Fn(usize) + Sync,
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
    progress: impl Fn(usize) + Sync,
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
    progress: impl Fn(usize) + Sync,
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
    progress: impl Fn(usize) + Sync,
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
    progress: impl Fn(usize) + Sync,
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
    progress: impl Fn(usize) + Sync,
) -> Result<ZScoreSummary, ZhuntError> {
    write_zscore_header(
        writer,
        plan.input_name,
        plan.sequence.len(),
        plan.range.from,
        plan.range.to,
    )?;

    let chunk_size = plan.chunk_size.max(1);
    thread::scope(|scope| -> Result<(), ZhuntError> {
        let (sender, receiver) = mpsc::sync_channel::<Vec<ZScoreRecord>>(0);
        let writer_handle = scope.spawn(move || -> io::Result<()> {
            for records in receiver {
                for record in records {
                    writeln!(writer, "{record}")?;
                }
                writer.flush()?;
            }
            Ok(())
        });

        let mut writer_stopped = false;
        for chunk_start in (0..plan.sequence.len()).step_by(chunk_size) {
            let chunk_end = (chunk_start + chunk_size).min(plan.sequence.len());
            let records = score_positions_parallel(
                config,
                plan.sequence,
                chunk_start,
                chunk_end,
                plan.range,
                plan.threads,
                &progress,
            );

            if sender.send(records).is_err() {
                writer_stopped = true;
                break;
            }
        }
        drop(sender);

        let writer_result = writer_handle
            .join()
            .map_err(|_| ZhuntError::Io(io::Error::other("Z-SCORE writer thread panicked")))?;
        writer_result?;

        if writer_stopped {
            return Err(ZhuntError::Io(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "Z-SCORE writer thread stopped before accepting all score chunks",
            )));
        }

        Ok(())
    })?;

    Ok(ZScoreSummary {
        input_name: plan.input_name.to_owned(),
        sequence_length: plan.sequence.len(),
        from_dinucleotide: plan.range.from,
        to_dinucleotide: plan.range.to,
        records_written: plan.sequence.len(),
    })
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
