use num_bigint::BigInt;
use num_traits::{One, Zero};

#[derive(Debug, PartialEq)]
pub enum LegendreResult {
    One,
    MinusOne,
    Zero,
}
/// Computes `a^{(p-1)/2} mod p`
pub fn legendre(a: &BigInt, p: &BigInt) -> LegendreResult {
    let one = BigInt::from(1u8);
    let pm1 = p - &one;

    let exp = &pm1 >> 1;
    let res = a.modpow(&exp, p);

    if res == one {
        LegendreResult::One
    } else if res == pm1 {
        LegendreResult::MinusOne
    } else {
        LegendreResult::Zero
    }
}

/// Solves the equation x^2 ≡ y (mod pq), for primes p and q.
/// # Safety:
/// The caller must ensure that `p` and `q` are prime.
/// If this condition is not met,
/// the function will result in undefined behavior.
pub fn mod_quad_eq(y: &BigInt, p: &BigInt, q: &BigInt) -> Option<BigInt> {
    let a1 = modsqrt(y, p)?;
    let a2 = modsqrt(y, q)?;
    crt_with2cong(&a1, p, &a2, q)
}

/// Computes x such that x^2 = a (mod p).
/// # Safety:
/// The caller must ensure that modulus is prime.
/// If this condition is not met,
/// the function will result in undefined behavior.
pub fn modsqrt(a: &BigInt, modulus: &BigInt) -> Option<BigInt> {
    tonelli_shanks(a, modulus)
}

fn tonelli_shanks(a: &BigInt, p: &BigInt) -> Option<BigInt> {
    let zero = BigInt::zero();
    let two = BigInt::from(2u8);

    // Step 1: Check if n is a quadratic residue modulo p using the Legendre symbol
    if legendre(a, p) != LegendreResult::One {
        return None; // No solution if n is not a quadratic residue
    }

    // Step 2: Tonelli-Shanks
    // Find q and s such that p - 1 = q * 2^s (with q odd)
    let mut q = p - 1u8;
    let mut s: u64 = 0;
    while (&q % 2u8) == zero {
        q >>= 1;
        s += 1;
    }

    if s == 1 {
        return Some(a.modpow(&((p + 1u8) / 4u8), p));
    }

    // Step 3: Find a non-residue z
    let mut z = BigInt::from(2u8);
    while legendre(&z, p) == LegendreResult::One {
        z += 1u32;
    }

    let mut c = z.modpow(&q, p);
    let mut t = a.modpow(&q, p);
    let mut r = a.modpow(&((q + 1u32) >> 1), p);
    let mut m = s;

    // Step 4: Iteratively apply the Tonelli-Shanks algorithm
    while (&t - 1u8) % p != zero {
        // Find the smallest i such that t^(2^i) ≡ 1 (mod p)
        let mut t2 = t.modpow(&two, p);
        let mut i = 1u64;
        while (&t2 - 1u8) % p != zero {
            t2 = t2.modpow(&two, p);
            i += 1;
        }

        let b = c.modpow(&(BigInt::one() << (m - i - 1)), p);
        r = (r * &b) % p;
        c = b.modpow(&two, p);
        t = (t * &c) % p;
        m = i;
    }

    Some(r)
}

/// Solves the system of two congruences using the Chinese Remainder Theorem
/// x ≡ a1 (mod m1)
/// x ≡ a2 (mod m2)
pub fn crt_with2cong(a1: &BigInt, m1: &BigInt, a2: &BigInt, m2: &BigInt) -> Option<BigInt> {
    // Ensure that m1 and m2 are coprime
    if gcd(m1, m2) != BigInt::one() {
        return None; // No solution
    }

    let m2_inv = m2.modinv(m1)?;
    let m1_inv = m1.modinv(m2)?;

    let part1 = a1 * m2 * m2_inv;
    let part2 = a2 * m1 * m1_inv;

    let result = (part1 + part2) % (m1 * m2);

    Some(result)
}

pub fn gcd(x: &BigInt, y: &BigInt) -> BigInt {
    let mut a = x.clone();
    let mut b = y.clone();

    while !b.is_zero() {
        let rem = &a % &b;
        a = b;
        b = rem;
    }
    a
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gcd() {
        assert_eq!(
            BigInt::from(2u8),
            gcd(&BigInt::from(10u8), &BigInt::from(14u8))
        );
        assert_eq!(BigInt::one(), gcd(&BigInt::from(10u8), &BigInt::from(19u8)));
        assert_eq!(BigInt::one(), gcd(&BigInt::one(), &BigInt::from(8u8)));

        assert_eq!(
            BigInt::from(2u8),
            gcd(
                &BigInt::from(12344389379474u64),
                &BigInt::from(34352353535234u64)
            )
        );
        assert_eq!(
            BigInt::one(),
            gcd(
                &BigInt::from(526734534543534554363244354235u128),
                &BigInt::from(234385345792309030345938593434u128)
            )
        );

        let a = BigInt::parse_bytes("199832844067622864359887407765124902064417551700433100578054818355344417262310911844894521".as_bytes(), 10).expect("Invlaid numner");
        let b = BigInt::parse_bytes(
            "15396585372166771612609511228027291957327038922374351562405570589084162289".as_bytes(),
            10,
        )
        .expect("Invlaid numner");
        assert_eq!(BigInt::from(361371693343444081u64), gcd(&a, &b));
    }

    #[test]
    fn test_legendre() {
        let minus_one = LegendreResult::MinusOne;
        let one = LegendreResult::One;
        let zero = LegendreResult::Zero;

        assert_eq!(one, legendre(&BigInt::from(1u8), &BigInt::from(3u32)));
        assert_eq!(
            minus_one,
            legendre(&BigInt::from(15u32), &BigInt::from(19u32))
        );
        assert_eq!(zero, legendre(&BigInt::from(15u32), &BigInt::from(5u32)));
        assert_eq!(one, legendre(&BigInt::from(30u32), &BigInt::from(127u32)));
        assert_eq!(
            minus_one,
            legendre(&BigInt::from(18u32), &BigInt::from(53u32))
        );
        assert_eq!(zero, legendre(&BigInt::from(331u32), &BigInt::from(331u32)));
        assert_eq!(
            minus_one,
            legendre(&BigInt::from(12345u32), &BigInt::from(331u32))
        );

        assert_eq!(
            minus_one,
            legendre(&BigInt::from(39048230328u64), &BigInt::from(611953u32))
        );
        assert_eq!(
            one,
            legendre(&BigInt::from(390482303325228u64), &BigInt::from(611953u32))
        );
        assert_eq!(
            minus_one,
            legendre(&BigInt::from(390482303325228u64), &BigInt::from(611827u32))
        );

        assert_eq!(
            one,
            legendre(
                &BigInt::from(390482303394388243242342325228u128),
                &BigInt::from(592303u32)
            )
        );
        assert_eq!(
            one,
            legendre(
                &BigInt::from(390482303394388243242342325228u128),
                &BigInt::from(911u32)
            )
        );
        assert_eq!(
            one,
            legendre(
                &BigInt::from(39048230339438824324234232522834532430u128),
                &BigInt::from(1000000005721u64)
            )
        );

        let bignum = BigInt::parse_bytes(
            "1239048230339438824324235937673894389289349482774892734232522834532430345678"
                .as_bytes(),
            10,
        )
        .expect("Invlaid numner");
        assert_eq!(one, legendre(&bignum, &BigInt::from(1000000004239u64)));

        let bignum = BigInt::parse_bytes("39048230339438824324235937673894389289838277378399874923874349482774892734232522834532430".as_bytes(), 10).expect("Invlaid numner");
        assert_eq!(
            minus_one,
            legendre(&bignum, &BigInt::from(1000000000211u64))
        );

        let bignum = BigInt::parse_bytes("543904823033943882432423593767389438928983827999737839987492387434948277489273423252283453243394820".as_bytes(), 10).expect("Invlaid numner");
        assert_eq!(one, legendre(&bignum, &BigInt::from(1000000000163u64)));

        let prime = BigInt::parse_bytes("98696044010893586188344909998761511353136994072407906264133493762200448224192052430017734037185522318240259".as_bytes(), 10).expect("Invlaid numner");
        let bignum = BigInt::parse_bytes(
            "543904823099737839987492387434948277489273423252283453243394820".as_bytes(),
            10,
        )
        .expect("Invlaid numner");
        assert_eq!(minus_one, legendre(&bignum, &prime));

        let bignum = BigInt::parse_bytes(
            "5439048230997378399874923874349482774892734232522834532433948203231".as_bytes(),
            10,
        )
        .expect("Invlaid numner");
        assert_eq!(one, legendre(&bignum, &prime));
    }

    #[test]
    fn test_modsqrt_small_nums() {
        let two = BigInt::from(2u8);
        let zero = BigInt::zero();

        // (a, p) pairs to test solution x^2 = a mod p
        let pairs: [(u128, u128); 18] = [
            (5, 19),
            (17, 19),
            (1241, 7),
            (2141234, 32573),
            (2141237, 32573),
            (13, 103867),
            (1, 99149),
            (99145, 99149),
            (5, 99149),
            (51, 99149),
            (939023, 1000000005721),
            (939023, 1000000004239),
            (111, 1000000004239),
            (1111, 1000000004239),
            (111111, 1000000004239),
            (605655, 1000000002403),
            (390482303394388243242342325228, 592303),
            (39048230339438824324234232522834532430, 1000000005721),
        ];

        for (a, p) in pairs {
            let a = BigInt::from(a);
            let p = BigInt::from(p);
            let res = modsqrt(&a, &p).unwrap();
            assert_eq!((res.modpow(&two, &p) - a) % p, zero);
        }
    }

    #[test]
    fn test_modsqrt_big_nums() {
        let two = BigInt::from(2u8);
        let zero = BigInt::zero();

        // (a, p) pairs to test solution x^2 = a mod p
        let pairs: [(&str, &str); 3] = [
            (
                "1239048230339438824324235937673894389289349482774892734232522834532430345678",
                "1000000004239",
            ),
            (
                "543904823033943882432423593767389438928983827999737839987492387434948277489273423252283453243394820",
                "1000000000163",
            ),
            (
                "5439048230997378399874923874349482774892734232522834532433948203231",
                "98696044010893586188344909998761511353136994072407906264133493762200448224192052430017734037185522318240259",
            ),
        ];

        for (a, p) in pairs {
            let a = BigInt::parse_bytes(a.as_bytes(), 10).expect("Invlaid numner");
            let p = BigInt::parse_bytes(p.as_bytes(), 10).expect("Invlaid numner");
            let res = modsqrt(&a, &p).unwrap();
            assert_eq!((res.modpow(&two, &p) - a) % p, zero);
        }
    }

    #[test]
    fn test_modsqrt_when_no_solutions() {
        let pairs: [(u128, u128); 5] = [(0, 3), (18, 19), (19, 19), (124, 18477), (1846, 1847)];
        for (a, p) in pairs {
            let a = BigInt::from(a);
            let p = BigInt::from(p);
            let res = modsqrt(&a, &p);
            assert_eq!(None, res);
        }

        let pairs: [(&str, &str); 4] = [
            (
                "39048230339438824324235937673894389289838277378399874923874349482774892734232522834532430",
                "1000000000211",
            ),
            (
                "543904823099737839987492387434948277489273423252283453243394820",
                "98696044010893586188344909998761511353136994072407906264133493762200448224192052430017734037185522318240259",
            ),
            (
                "5439048230997378399874923874349482774892734232522834532433948203232",
                "98696044010893586188344909998761511353136994072407906264133493762200448224192052430017734037185522318240259",
            ),
            (
                "98696044010893586188344909998761511353136994072407906264133493762200448224192052430017734037185522318240258",
                "98696044010893586188344909998761511353136994072407906264133493762200448224192052430017734037185522318240259",
            ),
        ];

        for (a, p) in pairs {
            let a = BigInt::parse_bytes(a.as_bytes(), 10).expect("Invlaid numner");
            let p = BigInt::parse_bytes(p.as_bytes(), 10).expect("Invlaid numner");
            let res = modsqrt(&a, &p);
            assert_eq!(None, res);
        }
    }

    #[test]
    fn test_crt() {
        let zero = BigInt::zero();
        let args = [
            ((1, 19),(1, 17)),
            ((0, 23),(1, 3)),
            ((19, 19),(17, 17)),
            ((3, 7), (4, 11)),
            ((123124, 104021), (423412, 104281)),
            ((12, 26249), (31435345, 43)),
            ((1234567, 929), (7654321, 881))
        ];

        for ((a1, m1), (a2, m2)) in args {
            let b = crt_with2cong(
                &BigInt::from(a1),
                &BigInt::from(m1),
                &BigInt::from(a2),
                &BigInt::from(m2),
            )
            .unwrap();

            assert_eq!(zero, (&b - &a1) % (m1));
            assert_eq!(zero, (&b - &a2) % (m2));
        }
    }
}
