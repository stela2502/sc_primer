use thiserror::Error;

pub type PrimerResult<T> = Result<T, PrimerError>;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum PrimerError {
    #[error("unknown chemistry: {0}")]
    UnknownChemistry(String),

    #[error("invalid grammar: {0}")]
    InvalidGrammar(String),

    #[error("invalid sequence coordinates: {0}")]
    InvalidCoordinates(String),

    #[error("rhapsody whitelist error: {0}")]
    Rhapsody(String),
}
