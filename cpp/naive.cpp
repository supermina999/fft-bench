#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>

using Complex = std::complex<double>;

int revBits(int a, int log) {
    int res = 0;
    for(int i = 0; i < log - 1; i++) {
        res *= 2;
        res += a & 1;
        a /= 2;
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
        if(!invert) {
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

void test() {
    constexpr int n = 8;
    std::vector<Complex> data;
    data.reserve(n);
    for(int i = 1; i <= n; i++) {
        data.emplace_back(i, 0.0);
    }
    fft(data, false);
    std::vector<Complex> expected = {
        {36, 0},
        {-4,9.656854},
        {-4, 4},
        {-4,1.656854},
        {-4, 0},
        {-4,-1.656854},
        {-4, -4},
        {-4,-9.656854}
    };
    double EPS = 0.0001;
    for(int i = 0; i < n; i++) {
        if(data[i].real() + EPS < expected[i].real() || data[i].real() - EPS > expected[i].real()
            || data[i].imag() + EPS < expected[i].imag() || data[i].imag() - EPS > expected[i].imag()) {
            std::terminate();
        }
    }
}

int main() {
    test();
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
