#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>

using namespace std;

class FlowField {
public:
    FlowField(int n, int m, int t)
        : N(n), M(m), T(t), data(n * m * t, 0.0f) {}

    // Read/write access
    double& operator()(int n, int i, int k) {
        assert(n >= 0 && n < N);
        assert(i >= 0 && i < M);
        assert(k >= 0 && k < T);
        return data[n * (M * T) + i * T + k];
    }

    // Read-only access
    const double& operator()(int n, int i, int k) const {
        assert(n >= 0 && n < N);
        assert(i >= 0 && i < M);
        assert(k >= 0 && k < T);
        return data[n * (M * T) + i * T + k];
    }

private:
    int N, M, T;
    std::vector<double> data;
};
