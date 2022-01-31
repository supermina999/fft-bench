use std::f64::consts::PI;
use num::complex::Complex64;

fn rev_bits(a: usize, log: usize) -> usize {
    let mut tmp = a;
    let mut res = 0;
    for _i in 0..log-1 {
        res *= 2;
        res += tmp & 1;
        tmp /= 2;
    }
    res
}

fn log2(a: usize) -> usize {
    let mut tmp = a;
    let mut res = 0;
    while tmp > 0 {
        tmp /= 2;
        res += 1;
    }
    res
}

fn fft(a: &mut Vec<Complex64>, invert: bool) {
    let n = a.len();
    let logn = log2(n);
    for i in 0..n {
        let j = rev_bits(i, logn);
        if i < j {
            let tmp = a[i];
            a[i] = a[j];
            a[j] = tmp;
        }
    }

    let mut ws: Vec<Complex64> = Vec::new();

    for step in 1..logn {
        let len = 1 << step;

        let mut angle = 2.0 * PI / len as f64;
        if !invert {
            angle = -angle;
        }

        ws.clear();
        ws.reserve(len / 2);
        ws.push(Complex64::new(1.0, 0.0));
        let wn = Complex64::new(angle.cos(), angle.sin());
        for i in 0..len/2-1 {
            ws.push(ws[i] * wn);
        }

        for i in (0..n).step_by(len) {
            for j in 0..len/2 {
                let u = a[i + j];
                let v = a[i + j + len / 2] * ws[j];
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
            }
        }
    }
    if invert {
        for i in 0..n {
            a[i] /= n as f64;
        }
    }
}

#[cfg(test)]
mod tests {
    use num::complex::Complex64;
    use crate::fft;

    #[test]
    fn test() {
        let n = 8;
        let mut data: Vec<Complex64> = Vec::with_capacity(n);
        for i in 1..n+1 {
            data.push(Complex64::new(i as f64, 0.0));
        }
        fft(&mut data, false);
        let expected = vec![
            Complex64::new(36.0, 0.0),
            Complex64::new(-4.0, 9.656854),
            Complex64::new(-4.0, 4.0),
            Complex64::new(-4.0, 1.656854),
            Complex64::new(-4.0, 0.0),
            Complex64::new(-4.0, -1.656854),
            Complex64::new(-4.0, -4.0),
            Complex64::new(-4.0, -9.656854),
        ];
        let eps = 0.0001;
        for i in 0..n {
            assert!(data[i].re + eps > expected[i].re && data[i].re - eps < expected[i].re);
            assert!(data[i].im + eps > expected[i].im && data[i].im - eps < expected[i].im);
        }
    }
}

fn main() {
    let n = 1024 * 1024;
    let m = 10;
    let iters = 100;

    let mut data: Vec<Complex64> = Vec::with_capacity(n);
    for i in 0..n {
        data.push(Complex64::new((i % m) as f64, 0.0));
    }

    for i in 0..iters {
        fft(&mut data, i % 2 == 0);
    }
}