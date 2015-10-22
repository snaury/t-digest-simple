#include <chrono>
#include <random>
#include <iostream>

#include "merging_tdigest.hpp"

constexpr size_t num_values = 10 * 1000 * 1000;

int main() {
    std::random_device seed;
    std::default_random_engine rng(seed());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    std::vector<double> data;
    std::cout << "Generating random data..." << std::endl;
    for (size_t i = 0; i < num_values; ++i) {
        data.push_back(dist(rng));
    }

    std::cout << "Computing the digest..." << std::endl;
    auto t0 = std::chrono::high_resolution_clock::now();
    MergingTDigest digest(100);
    for (auto value : data) {
        digest.add(value);
    }
    digest.compress();
    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0);
    auto nsperop = dt / (double)num_values;
    std::cout << nsperop.count() << "ns per value" << std::endl;
    std::cout << digest.size() << " centroids" << std::endl;
    std::cout << "p50 = " << digest.quantile(0.5) << std::endl;
    std::cout << "p95 = " << digest.quantile(0.95) << std::endl;
    std::cout << "p99 = " << digest.quantile(0.99) << std::endl;
    return 0;
}
