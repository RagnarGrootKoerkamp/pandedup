use std::{
    io::{BufWriter, IsTerminal, Write},
    path::PathBuf,
};

use clap::Parser;

#[path = "../src/mss.rs"]
mod mss;

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
    info!("Reading input..");
    let mut reader = needletail::parse_fastx_file(&args.input).unwrap();
    let mut contigs = Vec::new();
    while let Some(record) = reader.next() {
        contigs.push(record.unwrap().seq().into_owned());
    }

    let superstring = mss::masked_superstring(args.k, &contigs);
    let output = args
        .output
        .unwrap_or_else(|| args.input.with_extension("msfa"));
    let mut writer = BufWriter::new(std::fs::File::create(&output).unwrap());
    writeln!(writer, ">masked-superstring").unwrap();
    writer.write_all(&superstring).unwrap();
    writer.write_all(b"\n").unwrap();
    eprintln!("wrote {} bases to {}", superstring.len(), output.display());
}
