use clap::Parser;
use std::io::IsTerminal;

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

    let args = pandedup::Args::parse();
    pandedup::run(&args);
}
