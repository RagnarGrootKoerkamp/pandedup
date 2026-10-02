use std::{
    io::{BufWriter, IsTerminal, Write},
    path::PathBuf,
};

use clap::Parser;
use tracing::info;

#[path = "../src/mss.rs"]
mod mss;
#[path = "../src/timing.rs"]
mod timing;

#[derive(Parser)]
struct Args {
    /// Input FASTA/FASTQ file (plain or compressed).
    input: PathBuf,
    /// Output FASTA file. Defaults to input with an .msfa extension.
    #[arg(short, long)]
    output: Option<PathBuf>,
    /// K-mer size.
    #[arg(short, long, default_value_t = 64)]
    k: usize,
}

fn main() {
    tracing_subscriber::fmt()
        .compact()
        .with_target(false)
        .with_writer(std::io::stderr)
        .with_ansi(std::io::stderr().is_terminal())
        .with_timer(tracing_subscriber::fmt::time::ChronoLocal::new(
            "%H:%M:%S".to_string(),
        ))
        .init();
    let args = Args::parse();
    let timing = timing::StageTiming::start();
    info!("Reading input..");
    let (seq, ranges) = packed_seq::PackedSeqVec::from_fastx(&args.input);

    let superstring = if args.k <= 32 {
        mss::masked_superstring::<u64>(args.k, seq, ranges)
    } else {
        mss::masked_superstring::<u128>(args.k, seq, ranges)
    };
    let output = args
        .output
        .unwrap_or_else(|| args.input.with_extension("msfa"));
    let mut writer = BufWriter::new(std::fs::File::create(&output).unwrap());
    writeln!(writer, ">masked-superstring").unwrap();
    writer.write_all(&superstring).unwrap();
    writer.write_all(b"\n").unwrap();
    drop(writer);
    info!(
        "mss: {}, input {} bytes, output {} bases, {} bytes ({})",
        timing.finish(),
        std::fs::metadata(&args.input).unwrap().len(),
        superstring.len(),
        std::fs::metadata(&output).unwrap().len(),
        output.display()
    );
}
