use std::error::Error;

use lattice::integer;
use num_bigint::BigInt;

fn main() -> Result<(), Box<dyn Error>> {
    
    let y = BigInt::from(99544);
    let p = BigInt::from(611837);
    let q = BigInt::from(611953);

    let x = integer::mod_quad_eq(&y, &p, &q);

    if let Some(x) = x {
        println!("{x}");
    } else {
        println!("No solution!");
    }

    Ok(())
}
