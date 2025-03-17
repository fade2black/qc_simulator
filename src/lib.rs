// `rustup override set nightly` to switch to nightly
// and
// `rustup override set stable` to switch back

pub mod integer;
pub mod matrix;
pub mod utils;
pub mod vector;
pub mod vectors;
pub use crate::matrix::Matrix;
pub use crate::vector::Vector;
pub use crate::vectors::*;
