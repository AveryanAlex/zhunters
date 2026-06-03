use std::error::Error;
use std::fmt::{self, Display, Formatter};
use std::io;

#[derive(Debug)]
pub enum ZhuntError {
    Io(io::Error),
    InvalidWindowSize { window_size: usize },
    InvalidMinimumSize { min_size: usize },
    EmptySequence,
}

impl Display for ZhuntError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(err) => write!(f, "{err}"),
            Self::InvalidWindowSize { window_size } => {
                write!(
                    f,
                    "window size must be greater than zero, got {window_size}"
                )
            }
            Self::InvalidMinimumSize { min_size } => {
                write!(f, "minimum size must be greater than zero, got {min_size}")
            }
            Self::EmptySequence => write!(f, "input contains no A/T/G/C/N bases"),
        }
    }
}

impl Error for ZhuntError {}

impl From<io::Error> for ZhuntError {
    fn from(value: io::Error) -> Self {
        Self::Io(value)
    }
}
