//! Print FASTA contigs shorter than `k` to stdout.
//! Usage: cargo run --example filter_short_contigs -- input.fa 100

use std::{
    io::{self, BufWriter, Write},
    path::PathBuf,
};

use clap::Parser;

#[derive(Parser)]
struct Args {
    /// Input FASTA file.
    input: PathBuf,
    /// Exclude contigs of this length or longer.
    k: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let mut reader = needletail::parse_fastx_file(&args.input)?;
    let mut out = BufWriter::new(io::stdout().lock());

    let mut empty = 0;
    while let Some(record) = reader.next() {
        let record = record?;
        let seq = record.seq();
        if seq.is_empty() {
            empty += 1;
            continue;
        }
        if seq.len() < args.k {
            out.write_all(b">")?;
            out.write_all(record.id())?;
            out.write_all(b"\n")?;
            out.write_all(&seq)?;
            out.write_all(b"\n")?;
        }
    }
    eprintln!("Empty seq: {}", empty);

    Ok(())
}
