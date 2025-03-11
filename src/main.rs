use lattice::Matrix;

fn main() {
    let data = vec![
            1.2, 8.7, -66.0, 3.0,
            333.0, 2.0, 4.0, 53.0,
            -0.53, -34.0, 22.0, 3.0,
            -234.0, 3.02, 1.0, 23.0
        ];
       println!("{}", Matrix::new(4, 4, data).determinant());

    // let inv = mat.inverse();
    // inv.print();

    // println!("--------------------");
    // inv.mul(&mat).print();


    // let (l, u) = mat.lu_decomposition();
    // let mat1 = l.mul(&u);
}
