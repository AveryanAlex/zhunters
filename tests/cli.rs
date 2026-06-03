use std::fs;
use std::path::PathBuf;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

const REAL_DNA_FIXTURE: &str = "tests/fixtures/nc_043715_128k.fa";
const REAL_DNA_EXPECTED_DIR: &str = "tests/fixtures/expected";
const REAL_DNA_CASES: &[(usize, usize, usize)] = &[
    (1, 1, 1),
    (2, 2, 2),
    (3, 4, 4),
    (4, 6, 8),
    (5, 14, 11),
    (6, 21, 23),
    (7, 38, 38),
    (8, 76, 66),
    (9, 94, 94),
    (10, 133, 133),
    (11, 194, 194),
    (12, 265, 388),
];

#[test]
fn cli_writes_legacy_zscore_file_format() {
    let dir = temp_dir("format");
    let input = dir.join("sample.fa");
    fs::write(&input, b">sample\nACGTNNacgt\n").unwrap();

    let status = Command::new(env!("CARGO_BIN_EXE_zhunt"))
        .arg("4")
        .arg("2")
        .arg("3")
        .arg(input.as_os_str())
        .status()
        .unwrap();
    assert!(status.success());

    let output_path = PathBuf::from(format!("{}.Z-SCORE", input.display()));
    let output = fs::read_to_string(&output_path).unwrap();
    let lines: Vec<&str> = output.lines().collect();

    // Like the C original, input parsing scans the whole file for A/T/G/C/N
    // bytes and does not treat FASTA headers specially; the "a" in "sample"
    // is therefore part of the legacy sequence.
    assert_eq!(lines[0], format!("{} 11 2 3", input.display()));
    assert_eq!(lines.len(), 12);

    for (index, line) in lines.iter().skip(1).enumerate() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        assert_eq!(fields.len(), 8, "line: {line}");
        assert_eq!(fields[0].parse::<usize>().unwrap(), index + 1);
        let start = fields[0].parse::<usize>().unwrap();
        let end = fields[1].parse::<usize>().unwrap();
        let length = fields[2].parse::<usize>().unwrap();
        assert_eq!(end, start + length);
        assert_eq!(fields[6].len(), length);
        assert_eq!(fields[7].len(), length);
        assert!(fields[3].contains('.'));
        assert!(fields[4].contains('.'));
        assert!(fields[5].contains('e'));
    }

    let _ = fs::remove_dir_all(dir);
}

#[test]
fn cli_rejects_files_without_bases() {
    let dir = temp_dir("empty-sequence");
    let input = dir.join("empty.fa");
    fs::write(&input, b">xyz\nXYZ---\n").unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_zhunt"))
        .arg("4")
        .arg("2")
        .arg("3")
        .arg(input.as_os_str())
        .output()
        .unwrap();

    assert!(!output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("input contains no A/T/G/C/N bases"));

    let _ = fs::remove_dir_all(dir);
}

#[test]
fn cli_prints_readable_run_summary() {
    let dir = temp_dir("summary");
    let input = dir.join("sample.txt");
    fs::write(&input, b"ACGTACGT\n").unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_zhunt"))
        .arg("4")
        .arg("1")
        .arg("3")
        .arg(input.as_os_str())
        .output()
        .unwrap();
    assert!(output.status.success());

    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("Z-HUNT 3 scanner"), "{stdout}");
    assert!(
        stdout.contains(&format!("Input      : {}", input.display())),
        "{stdout}"
    );
    assert!(
        stdout.contains("Size range : 1..=3 dinucleotides"),
        "{stdout}"
    );
    assert!(stdout.contains("✓ Read 8 bases"), "{stdout}");
    assert!(stdout.contains("✓ Scored 8 circular positions"), "{stdout}");
    assert!(stdout.contains("✓ Wrote"), "{stdout}");
    assert!(stdout.contains("Done in"), "{stdout}");

    let _ = fs::remove_dir_all(dir);
}

#[test]
fn cli_accepts_optional_threads_argument() {
    let dir = temp_dir("threads");
    let input = dir.join("sample.txt");
    fs::write(&input, b"ACGTACGT\n").unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_zhunt"))
        .arg("--threads")
        .arg("1")
        .arg("4")
        .arg("1")
        .arg("3")
        .arg(input.as_os_str())
        .output()
        .unwrap();
    assert!(output.status.success());

    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("Threads    : 1"), "{stdout}");

    let output_path = PathBuf::from(format!("{}.Z-SCORE", input.display()));
    let output = fs::read_to_string(&output_path).unwrap();
    assert_eq!(output.lines().count(), 9);

    let _ = fs::remove_file(output_path);
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn cli_output_matches_legacy_c_fixture() {
    let dir = temp_dir("legacy-fixture");
    let input = dir.join("sample.txt");
    fs::write(&input, b"ACGTACGTACGT\n").unwrap();

    let status = Command::new(env!("CARGO_BIN_EXE_zhunt"))
        .arg("4")
        .arg("1")
        .arg("3")
        .arg(input.as_os_str())
        .status()
        .unwrap();
    assert!(status.success());

    let output_path = PathBuf::from(format!("{}.Z-SCORE", input.display()));
    let output = fs::read_to_string(&output_path).unwrap();
    let expected = format!(
        "{} 12 1 3\n\
1 7 6  25.133  18.557 2.072370e+01 acgtac   SASASA\n\
2 8 6  24.578  18.943 3.225294e+01 cgtacg   ASASAS\n\
3 9 6  25.133  18.557 2.072370e+01 gtacgt   SASASA\n\
4 10 6  27.130  16.264 5.659386e+00 tacgta   ASASAS\n\
5 11 6  25.133  18.557 2.072370e+01 acgtac   SASASA\n\
6 12 6  24.578  18.943 3.225294e+01 cgtacg   ASASAS\n\
7 13 6  25.133  18.557 2.072370e+01 gtacgt   SASASA\n\
8 14 6  27.130  16.264 5.659386e+00 tacgta   ASASAS\n\
9 15 6  25.133  18.557 2.072370e+01 acgtac   SASASA\n\
10 16 6  24.578  18.943 3.225294e+01 cgtacg   ASASAS\n\
11 17 6  25.133  18.557 2.072370e+01 gtacgt   SASASA\n\
12 18 6  27.130  16.264 5.659386e+00 tacgta   ASASAS\n",
        input.display()
    );
    assert_eq!(output, expected);

    let _ = fs::remove_dir_all(dir);
}

#[test]
fn cli_matches_real_dna_c_ground_truth_cases() {
    let update_tests = std::env::var_os("UPDATE_TESTS").is_some_and(|value| value == "1");

    for &(window_size, min_size, max_size) in REAL_DNA_CASES {
        let input = PathBuf::from(REAL_DNA_FIXTURE);
        let output_path = PathBuf::from(format!("{}.Z-SCORE", input.display()));
        let expected_path = PathBuf::from(REAL_DNA_EXPECTED_DIR)
            .join(format!("zhunt_{window_size}_{min_size}_{max_size}.Z-SCORE"));
        let _ = fs::remove_file(&output_path);

        let output = Command::new(env!("CARGO_BIN_EXE_zhunt"))
            .arg(window_size.to_string())
            .arg(min_size.to_string())
            .arg(max_size.to_string())
            .arg(REAL_DNA_FIXTURE)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "case {window_size} {min_size} {max_size}\nstdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );

        let actual = fs::read(&output_path).unwrap();
        if update_tests {
            fs::write(&expected_path, &actual).unwrap();
        } else {
            let expected = fs::read(&expected_path).unwrap();
            assert_eq!(
                actual,
                expected,
                "case {window_size} {min_size} {max_size} differs from {}; first diff: {}",
                expected_path.display(),
                first_text_diff(&expected, &actual)
            );
        }

        let _ = fs::remove_file(output_path);
    }
}

fn temp_dir(name: &str) -> PathBuf {
    let unique = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("zhunters-{name}-{}-{unique}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn first_text_diff(expected: &[u8], actual: &[u8]) -> String {
    let expected = String::from_utf8_lossy(expected);
    let actual = String::from_utf8_lossy(actual);
    for (line_index, (expected_line, actual_line)) in
        expected.lines().zip(actual.lines()).enumerate()
    {
        if expected_line != actual_line {
            return format!(
                "line {}\nexpected: {expected_line}\nactual:   {actual_line}",
                line_index + 1
            );
        }
    }
    format!(
        "line count/length differs: expected {} bytes, actual {} bytes",
        expected.len(),
        actual.len()
    )
}
