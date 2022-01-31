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
            a[j] = a[i];
            a[i] = tmp;
        }
    }

    for step in 1..logn {
        let len = 1 << step;

        let mut angle = 2.0 * PI / len as f64;
        if invert {
            angle = -angle;
        }

        let wn = Complex64::new(angle.cos(), angle.sin());
        for i in (0..n).step_by(len) {
            let mut w = Complex64::new(1.0, 0.0);
            for j in 0..len/2 {
                let u = a[i + j];
                let v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wn;
            }
        }
    }
    if invert {
        for i in 0..n {
            a[i] /= n as f64;
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