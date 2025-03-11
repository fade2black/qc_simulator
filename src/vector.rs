use crate::utils::*;

#[derive(Debug, Clone)]
pub struct Vector(Vec<f64>);

impl Vector {
    pub fn new(v: &[f64]) -> Self {
        Self(v.to_vec())
    }

    /// Creates a new vector from an array of integer numbers.
    pub fn from_int(arr: &[i64]) -> Self {
        let vec = arr.iter().map(|&x| x as f64).collect();
        Vector(vec)
    }

    /// Returns `i`-th coordinate of the vector.
    /// ```
    /// use lattice::Vector;
    ///
    /// let v = Vector::from_int(&[1, 2]);
    /// assert_eq!(v.get(0), 1.0);
    /// assert_eq!(v.get(1), 2.0);
    /// ```
    pub fn get(&self, i: usize) -> f64 {
        self.0[i]
    }

    /// Creates a zero vector of dim `n`.
    /// ```
    /// use lattice::Vector;
    ///
    /// let v = Vector::zero(2);
    /// assert_eq!(v.get(0), 0.0);
    /// assert_eq!(v.get(1), 0.0);
    /// ```
    pub fn zero(n: usize) -> Self {
        Self(vec![0f64; n])
    }

    /// Returns the number of elements in the vector.
    /// ```
    /// use lattice::Vector;
    ///
    /// let v = Vector::from_int(&[1, 1]);
    /// assert_eq!(v.size(), 2);
    /// ```
    pub fn size(&self) -> usize {
        self.0.len()
    }

    /// Computes dot product of two vectors.
    /// ```
    /// use lattice::Vector;
    ///
    /// let v1 = Vector::from_int(&[2, 1, -2]);
    /// let v2 = Vector::from_int(&[1, 2, 3]);
    /// assert_eq!(v1.dot(&v2), 2.0*1.0 + 1.0*2.0 + (-2.0)*3.0);
    /// ```
    pub fn dot(&self, other: &Self) -> f64 {
        if self.size() != other.size() {
            panic!(
                "vectors of size {} and {} are not compatible for dot product",
                self.size(),
                other.size(),
            );
        }
        self.0.iter().zip(other.0.iter()).map(|(x, y)| x * y).sum()
    }

    pub fn inner_product(&self, other: &Self) -> f64 {
        self.dot(other)
    }

    /// Computes a product of the vector by a scalar `c`.
    /// ```
    /// use lattice::Vector;
    ///
    /// let v1 = Vector::from_int(&[1, -2, 0]);
    /// let v2 = v1.scalar_prod(-2.0);
    /// let v3 = Vector::from_int(&[-2, 4, 0]);
    /// assert_eq!(v2, v3);
    /// ```
    pub fn scalar_prod(&self, c: f64) -> Self {
        Self(self.0.iter().map(|x| c * x).collect())
    }

    /// Adds two vectors.
    /// ```
    /// use lattice::Vector;
    ///
    /// let v1 = Vector::from_int(&[1, 1, -2]);
    /// let v2 = Vector::from_int(&[0, 1,  2]);
    /// let v3 = Vector::from_int(&[1, 2, 0]);
    /// assert_eq!(v1.add(&v2), v3);
    /// ```
    pub fn add(&self, other: &Self) -> Self {
        if self.size() != other.size() {
            panic!(
                "vectors of size {} and {} are not compatible for addition",
                self.size(),
                other.size(),
            );
        }
        Self(
            self.0
                .iter()
                .zip(other.0.iter())
                .map(|(x, y)| x + y)
                .collect(),
        )
    }

    /// Subtracts one vector from another one.
    /// ```
    /// use lattice::Vector;
    ///
    /// let v1 = Vector::from_int(&[1, 1, -2]);
    /// let v2 = Vector::from_int(&[0, 1,  3]);
    /// let v3 = Vector::from_int(&[1, 0, -5]);
    /// assert_eq!(v1.sub(&v2), v3);
    /// ```
    pub fn sub(&self, other: &Self) -> Self {
        if self.size() != other.size() {
            panic!(
                "vectors of size {} and {} are not compatible for subtraction",
                self.size(),
                other.size(),
            );
        }
        Self(
            self.0
                .iter()
                .zip(other.0.iter())
                .map(|(x, y)| x - y)
                .collect(),
        )
    }

    /// Computes the norm (length) of the vector and commonly denoted by `||v||`.
    /// ```
    /// use lattice::Vector;
    ///
    /// let v = Vector::from_int(&[1, 1, -2]);
    /// assert_eq!(v.norm(), 6f64.sqrt());
    /// ```
    pub fn norm(&self) -> f64 {
        let sum: f64 = self.0.iter().map(|&x| x * x).sum();
        (sum as f64).sqrt()
    }

    /// See `norm`.
    pub fn length(&self) -> f64 {
        self.norm()
    }

    /// Checks whether two vectors are orthogonal.
    ///
    /// Recall that two vectors `v` and `u` are said to be orthogonal if `||vu|| = 0`
    /// ```
    /// use lattice::Vector;
    ///
    /// let v1 = Vector::from_int(&[1, -1, 3]);
    /// let v2 = Vector::from_int(&[-1, 2, 1]);
    /// assert!(v1.is_orthogonal_to(&v2));
    /// ```
    pub fn is_orthogonal_to(&self, other: &Self) -> bool {
        if self.size() != other.size() {
            panic!(
                "vectors of size {} and {} are not compatible for orthogonality check",
                self.size(),
                other.size(),
            );
        }

        is_equal(self.dot(&other), 0f64)
    }

    /// Checks whether the vector is a unit vector.
    ///
    /// Recall that a vector is said to be a `unit` vector if its norm (length) is equal to 1.
    /// ```
    /// use lattice::Vector;
    ///
    /// let v = Vector::from_int(&[0, 1, 0]);
    /// assert!(v.is_unit());
    /// ```
    pub fn is_unit(&self) -> bool {
        is_equal(self.norm(), 1.0)
    }

    /// The vector projection on a nonzero vector `u`.
    /// u(v.u/u.u)
    pub fn proj_on(&self, u: &Vector) -> Vector {
        u.scalar_prod(self.dot(&u) / u.dot(u))
    }
}

impl PartialEq for Vector {
    fn eq(&self, other: &Self) -> bool {
        self.0
            .iter()
            .zip(other.0.iter())
            .fold(true, |sum, (x, y)| sum && is_equal(*x, *y))
    }
}

impl Eq for Vector {}

impl From<Vec<i64>> for Vector {
    fn from(arr: Vec<i64>) -> Self {
        let vec = arr.iter().map(|&x| x as f64).collect();
        Vector(vec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_proj1() {
        let v = Vector::new(&[7.0, 6.0]);
        let u = Vector::new(&[4.0, 2.0]);
        let proj = Vector::new(&[8.0, 4.0]);

        assert_eq!(proj, v.proj_on(&u));
    }

    #[test]
    fn test_proj2() {
        let v = Vector::new(&[1.0, 2.0, 3.0]);
        let u = Vector::new(&[3.0, 2.0, 1.0]);
        let proj = Vector::new(&[15.0 / 7.0, 10.0 / 7.0, 5.0 / 7.0]);

        assert_eq!(proj, v.proj_on(&u));
    }

    #[test]
    fn test_proj3() {
        let v = Vector::new(&[-3.0, 2.0, 1.0]);
        let u = Vector::new(&[1.0, 0.0, 0.0]);
        let proj = Vector::new(&[-3.0, 0.0, 0.0]);

        assert_eq!(proj, v.proj_on(&u));
    }
}
