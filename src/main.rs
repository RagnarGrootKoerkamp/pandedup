use clap::{Args, Parser, Subcommand, ValueEnum};
use ggcat_api::ExtraElaboration;
use std::{io::IsTerminal, num::NonZeroUsize, path::PathBuf};
use tracing::Level;

#[derive(Parser)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Deduplicate an input archive or FASTA file.
    Dedup(DedupArgs),
    /// Run GGCAT on a FASTA file.
    Ggcat(GgcatArgs),
    /// Build a masked superstring from contigs.
    Mss(MssArgs),
    /// Build a masked superstring by matching unitigs through their graph.
    Matchtigs(MssArgs),
}

#[derive(Args)]
struct DedupArgs {
    /// Input .agc, .tar.gz, or .fa.zst file.
    input: PathBuf,
    /// Output path. Defaults to input.dedup.fa.zst.
    #[arg(short, long)]
    output: Option<PathBuf>,
    /// K-mer size.
    #[arg(short, long, default_value_t = 64)]
    k: usize,
    /// Minimizer phrase window size.
    #[arg(short, long, default_value_t = 100)]
    w: usize,
    /// Number of threads. Defaults to the number of logical cores.
    #[arg(short = 'j', long)]
    threads: Option<usize>,
    /// Skip using the first input as a reference.
    #[arg(long = "no-reference", default_value_t = true, action = clap::ArgAction::SetFalse)]
    reference: bool,
    /// Deduplicate across reverse complements.
    #[arg(long)]
    canonical: bool,
    /// Minimizer length for phrases.
    #[arg(long, default_value_t = 8)]
    mini_k: usize,
}

#[derive(Args)]
struct GgcatArgs {
    /// Input .fa or .fa.zst file.
    input: PathBuf,
    /// Output FASTA file. Defaults to INPUT-kK-MODE.fa.
    #[arg(short, long)]
    output: Option<PathBuf>,
    /// K-mer size.
    #[arg(short, long, default_value_t = 64)]
    k: usize,
    /// Number of threads. Defaults to the number of physical cores.
    #[arg(short = 'j', long)]
    threads: Option<usize>,
    /// GGCAT extra elaboration mode.
    #[arg(long, value_enum, default_value_t = GgcatMode::Simplitigs)]
    mode: GgcatMode,
}

#[derive(Clone, Copy, ValueEnum)]
enum GgcatMode {
    None,
    UnitigLinks,
    GreedyMatchtigs,
    Eulertigs,
    Pathtigs,
    Simplitigs,
    FastEulertigs,
}

impl From<GgcatMode> for ExtraElaboration {
    fn from(mode: GgcatMode) -> Self {
        match mode {
            GgcatMode::None => Self::None,
            GgcatMode::UnitigLinks => Self::UnitigLinks,
            GgcatMode::GreedyMatchtigs => Self::GreedyMatchtigs,
            GgcatMode::Eulertigs => Self::Eulertigs,
            GgcatMode::Pathtigs => Self::Pathtigs,
            GgcatMode::Simplitigs => Self::FastSimplitigs,
            GgcatMode::FastEulertigs => Self::FastEulertigs,
        }
    }
}

#[derive(Args)]
struct MssArgs {
    /// Input FASTA/FASTQ file (plain or compressed).
    input: PathBuf,
    /// Output FASTA file. Defaults to input stem plus -mss.msfa or -greedytigs.msfa.
    #[arg(short, long)]
    output: Option<PathBuf>,
    /// K-mer size.
    #[arg(short, long, default_value_t = 64)]
    k: usize,
    /// Number of threads. Defaults to Rayon's worker count.
    #[arg(short = 'j', long)]
    threads: Option<NonZeroUsize>,
}

fn main() {
    let cli = Cli::parse();
    tracing_subscriber::fmt()
        .compact()
        .with_target(false)
        .with_writer(std::io::stderr)
        .with_ansi(std::io::stderr().is_terminal())
        .with_max_level(Level::DEBUG)
        .with_timer(tracing_subscriber::fmt::time::ChronoLocal::new(
            "%H:%M:%S".to_string(),
        ))
        .init();

    match cli.command {
        Command::Dedup(args) => {
            pandedup::dedup(
                &args.input,
                args.output.as_deref(),
                args.k,
                args.w,
                args.threads,
                args.reference,
                args.canonical,
                args.mini_k,
            );
        }
        Command::Ggcat(args) => {
            pandedup::ggcat::run(&pandedup::ggcat::GgcatConfig {
                input: args.input,
                output: args.output,
                k: args.k,
                threads: args.threads,
                mode: args.mode.into(),
            });
        }
        Command::Mss(args) => {
            pandedup::mss::run(
                &args.input,
                args.output.as_deref(),
                args.k,
                args.threads.map(NonZeroUsize::get),
            );
        }
        Command::Matchtigs(args) => {
            pandedup::matchtigs::run(&args.input, args.output.as_deref(), args.k);
        }
    }
}
