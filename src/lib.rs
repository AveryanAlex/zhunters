mod config;
mod constants;
mod error;
mod output;
mod record;
mod scoring;
mod sequence;

pub use config::ZhuntConfig;
pub use constants::DBZED;
pub use error::ZhuntError;
pub use output::{
    write_zscore_file_streaming, write_zscore_file_streaming_with_options, write_zscore_streaming,
    write_zscore_streaming_with_options, zscore_output_path, ZScoreRunOptions,
};
pub use record::{AntiSynPath, ZScoreOutput, ZScoreRecord, ZScoreSummary};
pub use scoring::{
    assign_probability, calculate_zscore, calculate_zscore_with_progress,
    calculate_zscore_with_progress_and_threads,
};
pub use sequence::{
    parse_sequence_bytes, parse_sequence_reader, read_sequence_file, CircularSequence,
};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{A, INT_DBZED};
    use crate::output::write_zscore_streaming_with_chunk_size;
    use crate::scoring::{
        assign_bzenergy_indices, best_anti_syn, find_delta_linking, linear_search, score_position,
        transition_row, BestAntiSyn, DinucleotideRange,
    };
    use std::io::{self, Write};

    fn bytes_to_ascii(bytes: &[u8]) -> String {
        bytes.iter().map(|byte| *byte as char).collect()
    }

    #[test]
    fn parse_sequence_filters_lowercases_and_wraps() {
        let sequence = parse_sequence_bytes(b">id\nAC GTnxX\n", 5).unwrap();
        assert_eq!(sequence.len(), 5);
        assert_eq!(bytes_to_ascii(&sequence.bases), "acgtn");
        assert_eq!(
            bytes_to_ascii(sequence.window(0, 10).as_ref()),
            "acgtnacgtn"
        );
    }

    #[test]
    fn parse_sequence_reader_streams_small_reads() {
        let reader = OneByteReader::new(b">id\nAC GTnxX\n");
        let sequence = parse_sequence_reader(reader, 5).unwrap();

        assert_eq!(sequence.len(), 5);
        assert_eq!(bytes_to_ascii(&sequence.bases), "acgtn");
        assert_eq!(
            bytes_to_ascii(sequence.window(0, 10).as_ref()),
            "acgtnacgtn"
        );
    }

    #[test]
    fn dinucleotide_indexes_match_legacy_table_order() {
        let sequence = parse_sequence_bytes(b"aattggcccann", 0).unwrap();
        let indexes = assign_bzenergy_indices(&sequence.bases, 6);
        assert_eq!(indexes, vec![0, 5, 10, 15, 12, 16]);
    }

    #[test]
    fn anti_syn_dynamic_programming_matches_exhaustive_search() {
        let config = ZhuntConfig::new(6).unwrap();
        let cases = [
            vec![0],
            vec![14],
            vec![16],
            vec![0, 5],
            vec![10, 11, 14],
            vec![16, 16, 16],
            vec![0, 5, 10, 15],
            vec![1, 6, 11, 12, 16],
            vec![4, 12, 5, 0, 7, 1, 2, 14, 16, 9, 16, 12],
        ];

        for case in cases {
            let dp = best_anti_syn(&case, &config.exp_dbzed);
            let brute = brute_force_best_anti_syn(&case, &config.exp_dbzed);
            assert_eq!(dp.esum, brute.esum, "case {case:?}");
            assert_eq!(dp.antisyn, brute.antisyn, "case {case:?}");
            assert_eq!(dp.bzenergy, brute.bzenergy, "case {case:?}");
        }
    }

    #[test]
    fn linear_search_finds_root() {
        let root = linear_search(0.0, 2.0, 1e-9, |x| x - 1.25);
        assert!((root - 1.25).abs() < 1e-8);
    }

    #[test]
    fn probability_keeps_legacy_tail_ratio_behavior() {
        assert!((assign_probability(29.653_713_5) - 2.0).abs() < 1e-12);
        assert!(assign_probability(20.0) > 1.0);
        assert!(assign_probability(40.0) < 1.0);
    }

    #[test]
    fn calculates_records_in_sequence_order() {
        let config = ZhuntConfig::new(4).unwrap();
        let sequence = parse_sequence_bytes(b"ACGTACGT", 8).unwrap();
        let output = calculate_zscore(&config, 1, 2, "sample".to_owned(), &sequence).unwrap();

        assert_eq!(output.sequence_length, 8);
        assert_eq!(output.from_dinucleotide, 1);
        assert_eq!(output.to_dinucleotide, 2);
        assert_eq!(output.records.len(), 8);
        for (index, record) in output.records.iter().enumerate() {
            assert_eq!(record.start, index + 1);
            assert_eq!(record.length, record.sequence.len());
            assert_eq!(record.length, record.antisyn.len());
        }
    }

    #[test]
    fn progress_callback_reports_completed_positions() {
        use std::sync::atomic::{AtomicUsize, Ordering};

        let config = ZhuntConfig::new(4).unwrap();
        let sequence = parse_sequence_bytes(b"ACGTACGT", 8).unwrap();
        let calls = AtomicUsize::new(0);

        let output = calculate_zscore_with_progress(
            &config,
            1,
            2,
            "sample".to_owned(),
            &sequence,
            |positions| {
                calls.fetch_add(positions, Ordering::Relaxed);
            },
        )
        .unwrap();

        assert_eq!(output.records.len(), sequence.len());
        assert_eq!(calls.load(Ordering::Relaxed), sequence.len());
    }

    #[test]
    fn streaming_writer_matches_collected_output_and_flushes_chunks() {
        use std::sync::atomic::{AtomicUsize, Ordering};

        let config = ZhuntConfig::new(4).unwrap();
        let sequence = parse_sequence_bytes(b"ACGTACGT", 8).unwrap();
        let collected = calculate_zscore(&config, 1, 2, "sample".to_owned(), &sequence).unwrap();
        let mut expected = Vec::new();
        collected.write_legacy(&mut expected).unwrap();

        let calls = AtomicUsize::new(0);
        let test_thread_id = std::thread::current().id();
        let mut writer = RecordingWriter::default();
        let summary = write_zscore_streaming_with_chunk_size(
            &mut writer,
            &config,
            DinucleotideRange { from: 1, to: 2 },
            "sample",
            &sequence,
            2,
            None,
            |positions| {
                calls.fetch_add(positions, Ordering::Relaxed);
            },
        )
        .unwrap();

        assert_eq!(writer.bytes, expected);
        assert_eq!(writer.flushes, 4);
        assert!(
            writer
                .write_thread_ids
                .iter()
                .any(|thread_id| *thread_id != test_thread_id),
            "record chunks should be written by the background writer thread"
        );
        assert_eq!(summary.records_written, sequence.len());
        assert_eq!(summary.from_dinucleotide, 1);
        assert_eq!(summary.to_dinucleotide, 2);
        assert_eq!(calls.load(Ordering::Relaxed), sequence.len());
    }

    #[test]
    fn position_five_case_uses_corrected_exact_scoring() {
        let config = ZhuntConfig::new(12).unwrap();
        let sequence = parse_sequence_bytes(b"tacattaatcatagcgntgnncag", 24).unwrap();
        let window = sequence.window(0, 24);
        let bzindex = assign_bzenergy_indices(window.as_ref(), 12);
        let anti = best_anti_syn(&bzindex, &config.exp_dbzed);
        let delta = find_delta_linking(&anti.bzenergy, &config.bztwist[..12], (A / 2.0) * 12.0);
        assert_eq!(format!("{}", anti.antisyn), "ASASASASASASASASASASSASA");
        assert_eq!(delta.delta_linking, 37.1014404296875);

        let record = score_position(&config, &sequence, 0, 8, 12);
        assert_eq!(
            format!("{record}\n"),
            "1 17 16  34.147   5.873 4.925847e-02 tacattaatcatagcg   ASASASASASASASAS\n"
        );
    }

    fn brute_force_best_anti_syn(bzindex: &[usize], exp_dbzed: &[[f64; 17]; 4]) -> BestAntiSyn {
        #[derive(Default)]
        struct BruteBest {
            esum: i64,
            states: Vec<u8>,
            bzenergy: Vec<f64>,
            initialized: bool,
        }

        fn recurse(
            position: usize,
            bzindex: &[usize],
            exp_dbzed: &[[f64; 17]; 4],
            states: &mut Vec<u8>,
            bzenergy: &mut Vec<f64>,
            esum: i64,
            best: &mut BruteBest,
        ) {
            if position == bzindex.len() {
                if !best.initialized || esum < best.esum {
                    best.esum = esum;
                    best.states = states.clone();
                    best.bzenergy = bzenergy.clone();
                    best.initialized = true;
                }
                return;
            }

            let previous_state = states.last().copied();
            let as_row = previous_state.map_or(0, |state| transition_row(state, 0));
            let as_energy = INT_DBZED[as_row][bzindex[position]];
            let as_esum = esum + as_energy;
            states.push(0);
            bzenergy.push(exp_dbzed[as_row][bzindex[position]]);
            recurse(
                position + 1,
                bzindex,
                exp_dbzed,
                states,
                bzenergy,
                as_esum,
                best,
            );
            states.pop();
            bzenergy.pop();

            let sa_row = previous_state.map_or(1, |state| transition_row(state, 1));
            let sa_esum = esum + INT_DBZED[sa_row][bzindex[position]];
            states.push(1);
            bzenergy.push(exp_dbzed[sa_row][bzindex[position]]);
            recurse(
                position + 1,
                bzindex,
                exp_dbzed,
                states,
                bzenergy,
                sa_esum,
                best,
            );
            states.pop();
            bzenergy.pop();
        }

        let mut best = BruteBest::default();
        recurse(
            0,
            bzindex,
            exp_dbzed,
            &mut Vec::new(),
            &mut Vec::new(),
            0,
            &mut best,
        );

        let mut mask = 0_u32;
        for state in &best.states {
            mask = (mask << 1) | u32::from(*state);
        }

        BestAntiSyn {
            esum: best.esum,
            bzenergy: best.bzenergy,
            antisyn: AntiSynPath::new(mask, best.states.len() as u8),
        }
    }

    #[derive(Default)]
    struct RecordingWriter {
        bytes: Vec<u8>,
        flushes: usize,
        write_thread_ids: Vec<std::thread::ThreadId>,
    }

    impl Write for RecordingWriter {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            self.write_thread_ids.push(std::thread::current().id());
            self.bytes.extend_from_slice(buf);
            Ok(buf.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            self.flushes += 1;
            Ok(())
        }
    }

    struct OneByteReader<'a> {
        bytes: &'a [u8],
        position: usize,
    }

    impl<'a> OneByteReader<'a> {
        fn new(bytes: &'a [u8]) -> Self {
            Self { bytes, position: 0 }
        }
    }

    impl std::io::Read for OneByteReader<'_> {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            if self.position == self.bytes.len() {
                return Ok(0);
            }

            buf[0] = self.bytes[self.position];
            self.position += 1;
            Ok(1)
        }
    }
}
