use crate::utils::*;

#[derive(Debug, Clone)]
pub struct Matrix {
    rows: usize,
    cols: usize,
    data: Vec<f64>,
}

impl Matrix {
    pub fn new(rows: usize, cols: usize, data: Vec<f64>) -> Self {
        Self { rows, cols, data }
    }

    pub fn zero(rows: usize, cols: usize) -> Matrix {
        Matrix::new(rows, cols, vec![0.0; rows * cols])
    }

    pub fn identity(rows: usize, cols: usize) -> Matrix {
        let mut mat = Matrix::zero(rows, cols);
        for i in 0..rows {
            mat.set(i, i, 1.0);
        }

        mat
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.cols + j]
    }

    pub fn set(&mut self, i: usize, j: usize, val: f64) {
        self.data[i * self.cols + j] = val;
    }

    pub fn add(&self, rhs: &Self) -> Matrix {
        if self.rows != rhs.rows && self.cols != rhs.cols {
            panic!(
                "inputs {} x {} and {} x {} are not compatible for matrix addition",
                self.rows, self.cols, rhs.rows, rhs.cols,
            );
        }

        let mut data = vec![];
        for i in 0..self.rows {
            for j in 0..self.cols {
                data.push(self.get(i, j) + rhs.get(i, j));
            }
        }

        Self {
            rows: self.rows,
            cols: self.cols,
            data,
        }
    }

    pub fn mul(&self, rhs: &Self) -> Matrix {
        if self.cols != rhs.rows {
            panic!(
                "inputs {} x {} and {} x {} are not compatible for matrix multiplication",
                self.rows, self.cols, rhs.rows, rhs.cols,
            );
        }

        let mut data = vec![];
        for i in 0..self.rows {
            for j in 0..rhs.cols {
                let mut sum = 0f64;
                for k in 0..self.cols {
                    sum += self.get(i, k) * rhs.get(k, j);
                }
                data.push(sum);
            }
        }

        Self {
            rows: self.rows,
            cols: rhs.cols,
            data,
        }
    }

    pub fn transpose(&self) -> Matrix {
        let mut data = vec![];

        for i in 0..self.cols {
            for j in 0..self.rows {
                data.push(self.get(j, i));
            }
        }
        Matrix::new(self.cols, self.rows, data)
    }

    pub fn inverse(&self) -> Matrix {
        if self.rows != self.cols {
            panic!("Matrix must be square.");
        }

        let n = self.rows;
        let mut augmented = Matrix::zero(n, 2 * n);

        for i in 0..n {
            for j in 0..n {
                augmented.set(i, j, self.get(i, j));
                augmented.set(i, j + n, if i == j { 1.0 } else { 0.0 });
            }
        }

        // Gauss elimination
        for pivot_row in 0..n {
            // Find the row with the largest pivot (partial pivoting)
            let mut max_row = pivot_row;
            for row in (pivot_row + 1)..n {
                if augmented.get(row, pivot_row).abs() > augmented.get(max_row, pivot_row).abs() {
                    max_row = row;
                }
            }

            if pivot_row != max_row {
                augmented.swap_rows(pivot_row, max_row);
            }

            let pivot_value = augmented.get(pivot_row, pivot_row);
            if is_equal(pivot_value, 0.0) {
                panic!("Matrix is singular and cannot be inverted.");
            }

            for col in 0..(2 * n) {
                let new_value = augmented.get(pivot_row, col) / pivot_value;
                augmented.set(pivot_row, col, new_value);
            }

            for row in (pivot_row + 1)..n {
                let factor = augmented.get(row, pivot_row);
                for col in 0..(2 * n) {
                    let new_value =
                        augmented.get(row, col) - factor * augmented.get(pivot_row, col);
                    augmented.set(row, col, new_value);
                }
            }
        }

        // Back substitution
        for pivot_row in (0..n).rev() {
            // Eliminate the entries above the pivot
            for row in (0..pivot_row).rev() {
                let factor = augmented.get(row, pivot_row);
                for col in 0..(2 * n) {
                    let new_value =
                        augmented.get(row, col) - factor * augmented.get(pivot_row, col);
                    augmented.set(row, col, new_value);
                }
            }
        }

        // Extract the inverse
        let mut inverse = Matrix::zero(n, n);
        for i in 0..n {
            for j in 0..n {
                inverse.set(i, j, augmented.get(i, j + n));
            }
        }

        inverse
    }

    fn swap_rows(&mut self, row1: usize, row2: usize) {
        for col in 0..self.cols {
            let temp = self.get(row1, col);
            self.set(row1, col, self.get(row2, col));
            self.set(row2, col, temp);
        }
    }

    pub fn determinant(&self) -> f64 {
        if self.rows != self.cols {
            panic!("Matrix must be square.");
        }

        let n = self.rows;
        let mut odd_swaps = false;

        let mut mat = self.clone();

        for pivot_row in 0..n {
            let mut max_row = pivot_row;
            for row in (pivot_row + 1)..n {
                if mat.get(row, pivot_row).abs() > mat.get(max_row, pivot_row).abs() {
                    max_row = row;
                }
            }

            if is_equal(mat.get(max_row, pivot_row), 0.0) {
                return 0.0;
            }

            if pivot_row != max_row {
                mat.swap_rows(pivot_row, max_row);
                odd_swaps = !odd_swaps;
            }

            for row in (pivot_row + 1)..n {
                let factor = mat.get(row, pivot_row) / mat.get(pivot_row, pivot_row);

                for col in pivot_row..n {
                    let new_value = mat.get(row, col) - factor * mat.get(pivot_row, col);
                    mat.set(row, col, new_value);
                }
            }
        }

        let mut det = 1.0;
        for i in 0..n {
            det *= mat.get(i, i);
        }

        if odd_swaps {
            return -det;
        }

        det
    }

    pub fn print(&self) {
        for i in 0..self.rows {
            for j in 0..self.cols {
                print!("{:.4} ", self.get(i, j));
            }
            print!("\n");
        }
    }

    /// Computes three matrices `P`, `L`, `U` for a matrix `A` such that
    /// `PA = LU`. Returns the tripple (`P`, `L`, `U`).
    pub fn lu_decomposition(&self) -> (Matrix, Matrix, Matrix) {
        let mut mat_l = Matrix::identity(self.rows, self.cols);
        let mut mat_p = Matrix::identity(self.rows, self.cols);
        let mut mat_u = self.clone();
        let n = mat_u.rows;

        for k in 0..n {
            println!("ITER #{k}:");
            // Pivoting
            let mut max_row = k;
            for i in k + 1..n {
                if mat_u.get(i, k).abs() > mat_u.get(max_row, k).abs() {
                    max_row = i;
                }
            }
            println!("PIVOT in the column k={k}: {}", mat_u.get(max_row, k));

            // If max_row is not k, swap rows k and max_row in the matrix
            if max_row != k {
                mat_u.swap_rows(k, max_row);
                mat_p.swap_rows(k, max_row);

                for j in 0..k {
                    let temp = mat_l.get(k, j);
                    mat_l.set(k, j, mat_l.get(max_row, j));
                    mat_l.set(max_row, j, temp);
                }
            }

            if is_equal(mat_u.get(k, k), 0.0) {
                panic!("Matrix is singlular an therefore not invertible.");
            }

            for i in k + 1..n {
                let factor = mat_u.get(i, k) / mat_u.get(k, k);
                mat_l.set(i, k, factor);
                for j in k..n {
                    let updated_value = mat_u.get(i, j) - factor * mat_u.get(k, j);
                    mat_u.set(i, j, updated_value);
                }
            }
        }

        (mat_p, mat_l, mat_u)
    }
}

impl PartialEq for Matrix {
    fn eq(&self, rhs: &Self) -> bool {
        if self.rows != rhs.rows || self.cols != rhs.cols {
            return false;
        }

        for i in 0..self.rows {
            for j in 0..self.cols {
                if !is_equal(self.get(i, j), rhs.get(i, j)) {
                    return false;
                }
            }
        }
        true
    }
}

impl Eq for Matrix {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add() {
        let m1 = Matrix::new(2, 3, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let m2 = Matrix::new(2, 3, vec![7.0, 8.0, 9.0, 10.0, 11.0, 12.0]);
        let m3 = Matrix::new(2, 3, vec![8.0, 10.0, 12.0, 14.0, 16.0, 18.0]);

        assert_eq!(m1.add(&m2), m3);
    }

    #[test]
    fn test_mul() {
        let m1 = Matrix::new(2, 3, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let m2 = Matrix::new(3, 1, vec![7.0, 8.0, 9.0]);
        let m3 = Matrix::new(2, 1, vec![50.0, 122.0]);

        assert_eq!(m1.mul(&m2), m3);
    }

    #[test]
    fn test_transpose() {
        let m1 = Matrix::new(2, 3, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let m2 = m1.transpose();

        for i in 0..m1.rows() {
            for j in 0..m1.cols() {
                assert_eq!(m1.get(i, j), m2.get(j, i));
            }
        }
    }

    #[test]
    fn test_identity() {
        let mat = Matrix::identity(3, 3);

        for i in 0..mat.rows() {
            for j in 0..mat.cols() {
                if i == j {
                    assert_eq!(mat.get(i, j), 1.0);
                } else {
                    assert_eq!(mat.get(i, j), 0.0);
                }
            }
        }
    }

    #[test]
    fn test_inverse1() {
        let mat = Matrix::new(1, 1, vec![1.0]);
        let inv = mat.inverse();

        assert_eq!(mat.mul(&inv), Matrix::identity(1, 1));
    }

    #[test]
    fn test_inverse2() {
        let mat = Matrix::new(2, 2, vec![1.0, 2.0, 3.0, 4.0]);
        let inv = mat.inverse();

        assert_eq!(mat.mul(&inv), Matrix::identity(2, 2));
    }

    #[test]
    fn test_inverse3() {
        let mat = Matrix::new(2, 2, vec![1.0, 0.0, 0.0, 1.0]);
        let inv = mat.inverse();

        assert_eq!(mat.mul(&inv), Matrix::identity(2, 2));
    }

    #[test]
    fn test_inverse4() {
        let mat = Matrix::new(3, 3, vec![1.0, -4.0, 2.0, 1.0, 2.0, 3.0, 4.0, 1.0, 5.0]);
        let inv = mat.inverse();

        assert_eq!(mat.mul(&inv), Matrix::identity(3, 3));
    }

    #[test]
    fn test_lu_decomposition1() {
        let mat = Matrix::new(2, 2, vec![1.0, 2.0, 3.0, 4.0]);

        let (p, l, u) = mat.lu_decomposition();
        assert_eq!(p.mul(&mat), l.mul(&u));
    }

    #[test]
    fn test_lu_decomposition2() {
        let mat = Matrix::new(3, 3, vec![4.0, 3.0, 2.0, 3.0, 3.0, 1.0, 2.0, 1.0, 3.0]);

        let l_expected = Matrix::new(
            3,
            3,
            vec![1.0, 0.0, 0.0, 0.75, 1.0, 0.0, 0.5, -2.0 / 3.0, 1.0],
        );

        let u_expected = Matrix::new(
            3,
            3,
            vec![4.0, 3.0, 2.0, 0.0, 0.75, -0.5, 0.0, 0.0, 5.0 / 3.0],
        );

        let (p, l, u) = mat.lu_decomposition();
        assert_eq!(l_expected, l);
        assert_eq!(u_expected, u);
        assert_eq!(p.mul(&mat), l.mul(&u));
    }

    #[test]
    fn test_lu_decomposition3() {
        let data = vec![
            1.2, 8.7, -66.0, 3.0, 333.0, 2.0, 4.0, 53.0, -0.53, -34.0, 22.0, 3.0, -234.0, 3.02,
            1.0, 23.0,
        ];

        let mat = Matrix::new(4, 4, data);
        let (p, l, u) = mat.lu_decomposition();
        assert_eq!(p.mul(&mat), l.mul(&u));
    }

    #[test]
    fn test_determinant() {
        assert!(is_equal(0.0, Matrix::zero(2, 2).determinant()));
        assert!(is_equal(1.0, Matrix::identity(2, 2).determinant()));
        assert!(is_equal(
            -2.0,
            Matrix::new(2, 2, vec![1.0, 2.0, 3.0, 4.0]).determinant()
        ));
        assert!(is_equal(-1.5, Matrix::new(1, 1, vec![-1.5]).determinant()));
        assert!(is_equal(
            0.0,
            Matrix::new(3, 3, vec![1.0, 2.0, 3.0, -1.0, 2.0, -4.0, 3.0, 6.0, 9.0]).determinant()
        ));
        assert!(is_equal(
            0.0,
            Matrix::new(3, 3, vec![1.0, 2.0, 3.0, -2.0, -4.0, -4.0, 3.0, 6.0, 8.0]).determinant()
        ));
        assert!(is_equal(
            0.0,
            Matrix::new(3, 3, vec![0.0, 2.0, 3.0, 0.0, -4.0, -4.0, 0.0, 2.0, 5.0]).determinant()
        ));

        let data = vec![
            1.2, 8.7, -66.0, 3.0, 333.0, 2.0, 4.0, 53.0, -0.53, -34.0, 22.0, 3.0, -234.0, 3.02,
            1.0, 23.0,
        ];
        assert!(is_equal(
            4.1731256061e7,
            Matrix::new(4, 4, data).determinant()
        ));
    }
}
