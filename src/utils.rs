pub fn is_equal(x: f64, y: f64) -> bool {
    (x - y).abs() < 1.0E-6
}
