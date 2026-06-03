pub(crate) const RT: f64 = 0.59004; // 0.00198 * 298
pub(crate) const A: f64 = 0.357; // 2 * (1/10.5 + 1/12)
pub(crate) const B: f64 = 0.4;
pub(crate) const K_RT: f64 = -0.252_120_1; // -1100/4363
pub(crate) const SIGMA: f64 = 16.948_003_53; // 10/RT
pub(crate) const EXP_LIMIT: f64 = -600.0;
pub(crate) const PI_DEGREES: f64 = 57.295_779_51; // 180/pi
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
