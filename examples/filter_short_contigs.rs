//! Print FASTA contigs shorter than `k` to stdout.
//! Usage: cargo run --example filter_short_contigs -- input.fa 100

use std::{
    io::{self, BufWriter, Write},
    path::PathBuf,
};

use clap::Parser;
use helicase::{FastxParser, HelicaseParser, ParserOptions, input::FromFile};

const CONFIG: helicase::Config = ParserOptions::default()
    .ignore_quality()
    .reject_non_actg()
    .config();

#[derive(Parser)]
struct Args {
    /// Input FASTA file.
    input: PathBuf,
    /// Exclude contigs of this length or longer.
    k: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let mut reader = FastxParser::<CONFIG>::from_file(&args.input)?;
    let mut out = BufWriter::new(io::stdout().lock());

    let mut empty = 0;
    while let Some(event) = reader.next() {
        event?;
        let seq = reader.get_dna_string();
        if seq.is_empty() {
            empty += 1;
            continue;
        }
        if seq.len() < args.k {
            out.write_all(b">")?;
            out.write_all(reader.get_header())?;
            out.write_all(b"\n")?;
            out.write_all(seq)?;
            out.write_all(b"\n")?;
        }
    }
    eprintln!("Empty seq: {}", empty);

    Ok(())
}
