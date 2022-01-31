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
        a /= 2;
    }
    return res;
}

constexpr int log2(int a) {
    int res = 0;
    while (a > 0) {
        a /= 2;
        res++;
    }
    return res;
}

template<int len>
void fft_impl(std::vector<Complex> &a, bool invert, std::vector<Complex>& ws) {
    double angle = 2.0 * M_PI / len;
    if (!invert) {
        angle = -angle;
    }

    ws.clear();
    ws.resize(len / 2 + 1);
    ws[0] = {1.0, 0.0};
    Complex wn(cos(angle), sin(angle));
    for(int i = 1; i < len / 2; i++) {
        ws[i] = ws[i - 1] * wn;
    }

    for (int i = 0; i < a.size(); i += len) {
        for (int j = 0; j < len / 2; j++) {
            auto u = a[i + j];
            auto v = a[i + j + len / 2] * ws[j];
            a[i + j] = u + v;
            a[i + j + len / 2] = u - v;
        }
    }
}

template<int len, int n>
struct FFT_Impl {
static void run(std::vector<Complex> &a, bool invert, std::vector<Complex>& ws) {
    fft_impl<len>(a, invert, ws);
    FFT_Impl<len * 2, n>::run(a, invert, ws);
}
};

template<int len>
struct FFT_Impl<len, len> {
static void run(std::vector<Complex> &a, bool invert, std::vector<Complex>& ws) {
    fft_impl<len>(a, invert, ws);
}
};

template<int n>
void fft(std::vector<Complex>& a, bool invert) {
    constexpr int logn = log2(n);
    for(int i = 0; i < n; i++) {
        int j = revBits(i, logn);
        if(i < j) {
            std::swap(a[i], a[j]);
        }
    }

    std::vector<Complex> ws;
    FFT_Impl<1, n>::run(a, invert, ws);

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
    fft<n>(data, false);
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
        if(data[i].re + EPS < expected[i].re || data[i].re - EPS > expected[i].re
            || data[i].im + EPS < expected[i].im || data[i].im - EPS > expected[i].im) {
            std::terminate();
        }
    }
}

int main() {
    constexpr int n = 1024 * 1024;
    int m = 10;
    int iters = 100;

    std::vector<Complex> data;
    data.reserve(n);
    for(int i = 0; i < n; i++) {
        data.emplace_back(i % m, 0.0);
    }

    for(int i = 0; i < iters; i++) {
        fft<n>(data, i % 2 == 0);
    }

    return 0;
}
