use std::borrow::Cow;
use std::fmt::{self, Display, Formatter};
use std::io::{self, Write};

#[derive(Debug, Clone, PartialEq)]
pub struct ZScoreRecord<'a> {
    pub start: usize,
    pub length: usize,
    pub delta_linking: f64,
    pub slope: f64,
    pub probability: f64,
    pub sequence: Cow<'a, [u8]>,
    pub antisyn: AntiSynPath,
}

impl Display for ZScoreRecord<'_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} {} {} {:7.3} {:7.3} {} {}   {}",
            self.start,
            self.start + self.length,
            self.length,
            self.delta_linking,
            self.slope,
            Scientific(self.probability),
            AsciiBytes(&self.sequence),
            self.antisyn
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AntiSynPath {
    pub(crate) mask: u32,
    pub(crate) dinucleotides: u8,
}

impl AntiSynPath {
    pub(crate) fn new(mask: u32, dinucleotides: u8) -> Self {
        Self {
            mask,
            dinucleotides,
        }
    }

    pub fn len(self) -> usize {
        2 * usize::from(self.dinucleotides)
    }

    pub fn is_empty(self) -> bool {
        self.dinucleotides == 0
    }
}

impl Display for AntiSynPath {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let dinucleotides = usize::from(self.dinucleotides);
        for position in 0..dinucleotides {
            let shift = dinucleotides - 1 - position;
            if (self.mask >> shift) & 1 == 0 {
                f.write_str("AS")?;
            } else {
                f.write_str("SA")?;
            }
        }
        Ok(())
    }
}

struct AsciiBytes<'a>(&'a [u8]);

impl Display for AsciiBytes<'_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let text = std::str::from_utf8(self.0).map_err(|_| fmt::Error)?;
        f.write_str(text)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ZScoreOutput<'a> {
    pub input_name: String,
    pub sequence_length: usize,
    pub from_dinucleotide: usize,
    pub to_dinucleotide: usize,
    pub records: Vec<ZScoreRecord<'a>>,
}

impl ZScoreOutput<'_> {
    pub fn write_legacy<W: Write>(&self, mut writer: W) -> io::Result<()> {
        crate::output::write_zscore_header(
            &mut writer,
            &self.input_name,
            self.sequence_length,
            self.from_dinucleotide,
            self.to_dinucleotide,
        )?;
        for record in &self.records {
            writeln!(writer, "{record}")?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ZScoreSummary {
    pub input_name: String,
    pub sequence_length: usize,
    pub from_dinucleotide: usize,
    pub to_dinucleotide: usize,
    pub records_written: usize,
}

pub(crate) struct Scientific(pub(crate) f64);

impl Display for Scientific {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write_c_scientific(f, self.0)
    }
}

fn write_c_scientific(f: &mut Formatter<'_>, value: f64) -> fmt::Result {
    if value == 0.0 {
        return write!(f, "0.000000e+00");
    }

    if !value.is_finite() {
        return write!(f, "{value:.6e}");
    }

    let sign = if value.is_sign_negative() { "-" } else { "" };
    let value = value.abs();
    let mut exponent = value.log10().floor() as i32;
    let mut mantissa = value / 10_f64.powi(exponent);

    if (mantissa * 1_000_000.0).round() >= 10_000_000.0 {
        mantissa /= 10.0;
        exponent += 1;
    }

    write!(f, "{sign}{mantissa:.6}e{exponent:+03}")
}

#[cfg(test)]
mod tests {
    use super::Scientific;

    #[test]
    fn c_scientific_uses_two_digit_signed_exponent() {
        assert_eq!(format!("{}", Scientific(1.0)), "1.000000e+00");
        assert_eq!(format!("{}", Scientific(0.0123)), "1.230000e-02");
        assert_eq!(format!("{}", Scientific(123.0)), "1.230000e+02");
    }
}
