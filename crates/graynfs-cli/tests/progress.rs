use std::process::Command;

fn cli_binary() -> String {
    env!("CARGO_BIN_EXE_graynfs").to_string()
}

#[test]
fn progress_emits_stderr_lines_during_small_sieve_run() {
    let output = Command::new(cli_binary())
        .args([
            "100160063",
            "--degree",
            "2",
            "--skip-preprocessing",
            "--progress",
        ])
        .output()
        .expect("run graynfs with progress flag");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !stderr.trim().is_empty(),
        "expected progress lines on stderr"
    );
    assert!(
        stderr.contains("progress:") && stderr.contains("stage="),
        "unexpected stderr shape: {stderr}"
    );
}

#[test]
fn json_stdout_stays_pure_with_progress_on() {
    let output = Command::new(cli_binary())
        .args([
            "100160063",
            "--degree",
            "2",
            "--skip-preprocessing",
            "--json",
            "--progress",
        ])
        .output()
        .expect("run graynfs with json and progress flags");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert_eq!(
        stdout.lines().count(),
        1,
        "stdout must be exactly one JSON object line: {stdout}"
    );
    let value: serde_json::Value =
        serde_json::from_str(stdout.trim()).expect("stdout parses as JSON");
    assert!(value.is_object(), "stdout must be a JSON object");
    assert_eq!(value["factors"], serde_json::json!(["10007", "10009"]));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("progress:"),
        "progress must go to stderr only: {stderr}"
    );
}
