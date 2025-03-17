use crate::Vector;

/// Checks whether a set of vectors are orthogonal.
///
/// A set of vectors is called _orthogonal set_
/// if each pair of distinct vectors from the set is orthogonal.
pub fn is_orthogonal(vectors: &Vec<Vector>) -> bool {
    let n = vectors.len();

    for i in 0..(n - 1) {
        for j in (i + 1)..n {
            if !vectors[i].is_orthogonal_to(&vectors[j]) {
                return false;
            }
        }
    }

    true
}

/// Takes a set of linearly independent vectors and transforms them into a set of orthogonal vectors.
/// Notes that the resulting orthogonal vectors from this process span the same subspace
/// as the original set of vectors. 
pub fn gram_schmidt(vectors: &Vec<Vector>) -> Vec<Vector> {
    let mut basis = vec![];
    let zero = Vector::zero(vectors[0].size());

    basis.push(vectors[0].clone());
    for x in vectors[1..].iter() {
        let sum = basis
            .iter()
            .fold(zero.clone(), |acc, b| acc.add(&x.proj_on(&b)));
        basis.push(x.sub(&sum));
    }

    basis
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_orthogonal() {
        let vectors = vec![
            Vector::from_int(&[1, 0, 0]),
            Vector::from_int(&[0, 1, 0]),
            Vector::from_int(&[0, 0, 1]),
        ];

        assert!(is_orthogonal(&vectors));
    }

    #[test]
    fn test_gram_schmidt_trivial() {
        let v1 = Vector::from_int(&[3, 6, 0]);
        let vectors = vec![v1];

        let basis = gram_schmidt(&vectors);
        assert_eq!(basis[0], Vector::from_int(&[3, 6, 0]));
    }

    #[test]
    fn test_gram_schmidt1() {
        let v1 = Vector::from_int(&[3, 6, 0]);
        let v2 = Vector::from_int(&[1, 2, 2]);

        let vectors = vec![v1, v2];

        let basis = gram_schmidt(&vectors);

        assert_eq!(basis[0], Vector::from_int(&[3, 6, 0]));
        assert_eq!(basis[1], Vector::from_int(&[0, 0, 2]));
    }

    #[test]
    fn test_gram_schmidt2() {
        let vectors = vec![
            Vector::from_int(&[1, 0, 0]),
            Vector::from_int(&[0, 1, 0]),
            Vector::from_int(&[0, 0, 1]),
        ];

        let basis = gram_schmidt(&vectors);
        assert_eq!(basis[0], Vector::from_int(&[1, 0, 0]));
        assert_eq!(basis[1], Vector::from_int(&[0, 1, 0]));
        assert_eq!(basis[2], Vector::from_int(&[0, 0, 1]));
    }

    #[test]
    fn test_gram_schmidt3() {
        let vectors = vec![
            Vector::from_int(&[1, 1, 1, 1]),
            Vector::from_int(&[0, 1, 1, 1]),
            Vector::from_int(&[0, 0, 1, 1]),
        ];

        let basis = gram_schmidt(&vectors);
        assert_eq!(basis[0], Vector::from_int(&[1, 1, 1, 1]));
        assert_eq!(basis[1], Vector::new(&[-0.75, 0.25, 0.25, 0.25]));
        assert_eq!(
            basis[2],
            Vector::new(&[0.0, -2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0])
        );
    }
}
