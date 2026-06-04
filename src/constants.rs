/// RT in kcal/mol at 25 °C, using the CODATA 2022 molar gas constant.
pub(crate) const RT: f64 = 0.5924849497028442;
pub(crate) const A: f64 = 2. * (1. / 10.5 + 1. / 12.);
pub(crate) const B: f64 = 0.4;
pub(crate) const K_RT: f64 = -1100. / 4363.;
pub(crate) const SIGMA: f64 = 10. / RT;
pub(crate) const STREAM_CHUNK_POSITIONS: usize = 256 * 1024;

/// Delta BZ Energy of dinucleotide, indexed by anti/syn transition row and
/// dinucleotide index (AA, AT, ..., CC, NN/unknown).
pub const DBZED: [[f64; 17]; 4] = [
    [
        4.40, 6.20, 3.40, 5.20, 2.50, 4.40, 1.40, 3.30, 3.30, 5.20, 2.40, 4.20, 1.40, 3.40, 0.66,
        2.40, 4.26,
    ],
    [
        4.40, 2.50, 3.30, 1.40, 6.20, 4.40, 5.20, 3.40, 3.40, 1.40, 2.40, 0.66, 5.20, 3.30, 4.20,
        2.40, 4.26,
    ],
    [
        6.20, 6.20, 5.20, 5.20, 6.20, 6.20, 5.20, 5.20, 5.20, 5.20, 4.00, 4.00, 5.20, 5.20, 4.00,
        4.00, 4.26,
    ],
    [
        6.20, 6.20, 5.20, 5.20, 6.20, 6.20, 5.20, 5.20, 5.20, 5.20, 4.00, 4.00, 5.20, 5.20, 4.00,
        4.00, 4.26,
    ],
];

/// Integer centi-kcal/mol representation of [`DBZED`]. Using exact integer
/// sums avoids legacy floating-point path-selection artifacts while preserving
/// the same energy table, including the `N`/unknown column.
pub(crate) const INT_DBZED: [[i64; 17]; 4] = [
    [
        440, 620, 340, 520, 250, 440, 140, 330, 330, 520, 240, 420, 140, 340, 66, 240, 426,
    ],
    [
        440, 250, 330, 140, 620, 440, 520, 340, 340, 140, 240, 66, 520, 330, 420, 240, 426,
    ],
    [
        620, 620, 520, 520, 620, 620, 520, 520, 520, 520, 400, 400, 520, 520, 400, 400, 426,
    ],
    [
        620, 620, 520, 520, 620, 620, 520, 520, 520, 520, 400, 400, 520, 520, 400, 400, 426,
    ],
];
