//! KronoPrime command line entry point.

use std::process::ExitCode;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    ExitCode::from(kronoprime::cli::run(&args))
}
