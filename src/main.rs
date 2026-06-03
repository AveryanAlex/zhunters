use clap::Parser;
use indicatif::{HumanCount, HumanDuration, ProgressBar, ProgressStyle};
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::process::ExitCode;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;
use zhunters::{
    read_sequence_file, write_zscore_file_streaming_with_options, zscore_output_path,
    ZScoreRunOptions, ZhuntConfig,
};

#[derive(Debug, Parser)]
#[command(
    name = "zhunt",
    about = "Scan DNA sequences for Z-DNA-forming regions",
    override_usage = "zhunt [--threads <threads>] [-o <output>] <windowsize> <minsize> <maxsize> <datafile>"
)]
struct Cli {
    #[arg(long, value_name = "threads")]
    threads: Option<NonZeroUsize>,
    #[arg(short, long, value_name = "output")]
    output: Option<PathBuf>,
    #[arg(value_name = "windowsize")]
    window_size: usize,
    #[arg(value_name = "minsize")]
    min_size: usize,
    #[arg(value_name = "maxsize")]
    max_size: usize,
    #[arg(value_name = "datafile")]
    datafile: String,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("zhunt: {error}");
            ExitCode::FAILURE
        }
    }
}

fn run(cli: Cli) -> Result<(), Box<dyn std::error::Error>> {
    let start = Instant::now();

    let config = ZhuntConfig::new(cli.window_size)?;
    let output_path = cli
        .output
        .clone()
        .unwrap_or_else(|| zscore_output_path(&cli.datafile));

    print_run_summary(&cli, &output_path);

    let sequence = read_sequence_file(&cli.datafile, 2 * config.max_dinucleotides())?;

    println!("✓ Read {} bases", HumanCount(sequence.len() as u64));

    let progress_units = ((sequence.len() as u64) / 1_000).max(1);
    let completed_positions = AtomicU64::new(0);
    let progress = ProgressBar::new(progress_units);
    progress.set_style(progress_style());
    progress.set_message("scoring and streaming positions");

    let summary_result = write_zscore_file_streaming_with_options(
        &output_path,
        &config,
        ZScoreRunOptions {
            min_size: cli.min_size,
            max_size: cli.max_size,
            input_name: &cli.datafile,
            threads: cli.threads.map(NonZeroUsize::get),
        },
        &sequence,
        |positions| {
            let completed = completed_positions.fetch_add(positions as u64, Ordering::Relaxed)
                + positions as u64;
            progress.set_position((completed / 1_000).min(progress_units));
        },
    );
    progress.finish_and_clear();
    let summary = summary_result?;

    println!(
        "✓ Scored {} circular positions ({}..={} dinucleotides)",
        HumanCount(summary.records_written as u64),
        summary.from_dinucleotide,
        summary.to_dinucleotide
    );

    println!(
        "✓ Wrote {} as results became available",
        output_path.display()
    );
    println!("Done in {}", HumanDuration(start.elapsed()));
    Ok(())
}

fn print_run_summary(cli: &Cli, output_path: &std::path::Path) {
    println!("ZHunters scanner");
    println!("────────────────");
    println!("Input      : {}", cli.datafile);
    println!("Output     : {}", output_path.display());
    println!("Window     : {} dinucleotides", cli.window_size);
    if let Some(threads) = cli.threads {
        println!("Threads    : {threads}");
    }
    println!(
        "Size range : {}..={} dinucleotides",
        cli.min_size, cli.max_size
    );
    println!();
}

fn progress_style() -> ProgressStyle {
    ProgressStyle::with_template(
        "[{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}k/{len}k \
         ({percent}%, ETA {eta}) {msg}",
    )
    .expect("progress template is valid")
    .progress_chars("=>-")
}
