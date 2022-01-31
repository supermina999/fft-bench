#include <vector>
#include <algorithm>
#include <cmath>

struct Complex {
    Complex() = default;
    Complex(double re, double im): re(re), im(im) {}

    Complex operator+(Complex other) const {
        return {re + other.re, im + other.im};
    }
    Complex operator-(Complex other) const {
        return {re - other.re, im - other.im};
    }
    Complex operator*(Complex other) const {
        return {re * other.re - im * other.im, re * other.im + im * other.re};
    }
    Complex& operator/=(double val) {
        re /= val;
        im /= val;
        return *this;
    }

    double re = 0.0, im = 0.0;
};

int revBits(int a, int log) {
    int res = 0;
    for(int i = 0; i < log - 1; i++) {
        res *= 2;
        res += a & 1;
        res /= 2;
    }
    return res;
}

int log2(int a) {
    int res = 0;
    while (a > 0) {
        a /= 2;
        res++;
    }
    return res;
}

void fft(std::vector<Complex>& a, bool invert) {
    int n = a.size();
    int logn = log2(n);
    for(int i = 0; i < n; i++) {
        int j = revBits(i, logn);
        if(i < j) {
            std::swap(a[i], a[j]);
        }
    }

    for(int step = 1; step < logn; step++) {
        int len = 1 << step;

        double angle = 2.0 * M_PI / len;
        if(invert) {
            angle = -angle;
        }

        Complex wn(cos(angle), sin(angle));
        for(int i = 0; i < n; i += len) {
            Complex w(1.0, 0.0);
            for(int j = 0; j < len / 2; j++) {
                auto u = a[i + j];
                auto v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w = w * wn;
            }
        }
    }

    if(invert) {
        for(int i = 0; i < n; i++) {
            a[i] /= n;
        }
    }
}

int main() {
    int n = 1024 * 1024;
    int m = 10;
    int iters = 100;

    std::vector<Complex> data;
    data.reserve(n);
    for(int i = 0; i < n; i++) {
        data.emplace_back(i % m, 0.0);
    }

    for(int i = 0; i < iters; i++) {
        fft(data, i % 2 == 0);
    }

    return 0;
}
