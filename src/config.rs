use crate::constants::{A, B, DBZED, RT};
use crate::ZhuntError;

/// Precomputed constants for a run of the Z-HUNT scanner.
#[derive(Debug, Clone)]
pub struct ZhuntConfig {
    pub(crate) max_dinucleotides: usize,
    pub(crate) bztwist: Vec<f64>,
    pub(crate) exp_dbzed: [[f64; 17]; 4],
}

impl ZhuntConfig {
    pub fn new(max_dinucleotides: usize) -> Result<Self, ZhuntError> {
        if max_dinucleotides == 0 {
            return Err(ZhuntError::InvalidWindowSize {
                window_size: max_dinucleotides,
            });
        }

        let mut ab = B + B;
        let mut bztwist = Vec::with_capacity(max_dinucleotides);
        for _ in 0..max_dinucleotides {
            ab += A;
            bztwist.push(ab);
        }

        let mut exp_dbzed = [[0.0; 17]; 4];
        for (row, row_values) in DBZED.iter().enumerate() {
            for (column, energy) in row_values.iter().enumerate() {
                exp_dbzed[row][column] = (-energy / RT).exp();
            }
        }

        Ok(Self {
            max_dinucleotides,
            bztwist,
            exp_dbzed,
        })
    }

    pub fn max_dinucleotides(&self) -> usize {
        self.max_dinucleotides
    }
}
