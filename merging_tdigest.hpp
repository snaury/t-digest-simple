#pragma once

#include <cmath>
#include <cstdint>
#include <vector>
#include <algorithm>

class MergingTDigest {
public:
    struct Centroid {
        double mean;
        int64_t count;

        explicit Centroid(double mean, int64_t count)
          : mean(mean)
          , count(count)
        { }

        Centroid& operator+=(const Centroid& other) {
            count += other.count;
            mean += other.count * (other.mean - mean) / count;
            return *this;
        }

        bool operator<(const Centroid& other) const {
            return mean < other.mean;
        }
    };

private:
    std::vector<Centroid> summary_;
    double epsilon_;
    size_t max_buffer_size_;
    size_t buffer_size_;
    size_t count_;

    static double interpolate(double x, double x1, double y1, double x2, double y2) {
        double d = x2 - x1;
        double a = (x2 - x) / d;
        double b = (x - x1) / d;
        return a * y1 + b * y2;
    }

public:
    MergingTDigest(double compression = 100, size_t max_buffer_size = 2048)
      : epsilon_(1.0 / compression)
      , max_buffer_size_(max_buffer_size)
      , buffer_size_(0)
      , count_(0)
    { }

    void add(double x, int64_t w = 1) {
        add(Centroid(x, w));
    }

    void add(const Centroid& c) {
        if (buffer_size_ >= max_buffer_size_) {
            compress();
        }
        summary_.push_back(c);
        ++buffer_size_;
        ++count_;
    }

    size_t size() {
        compress();
        return summary_.size();
    }

    void compress() {
        if (buffer_size_ > 0) {
            std::stable_sort(summary_.begin(), summary_.end());
            buffer_size_ = 0;
            if (summary_.size() > 3) {
                auto l = summary_.begin();
                auto r = std::next(l);
                int64_t sum = 0;
                while (r != summary_.end()) {
                    // use quantile with the smallest error
                    double ql = (sum + l->count * 0.5) / count_;
                    double err = ql * (1.0 - ql);
                    double qr = (sum + l->count + r->count * 0.5) / count_;
                    double err2 = qr * (1.0 - qr);
                    if (err > err2) {
                        err = err2;
                    }
                    double k = 4 * count_ * err * epsilon_;
                    if (l->count + r->count <= k) {
                        // it is possible to merge left and right
                        *l += *r;
                    } else {
                        // skip to the next pair
                        sum += l->count;
                        if (++l != r) {
                            *l = *r;
                        }
                    }
                    ++r;
                }
                // remove all leftover centroids
                summary_.erase(++l, summary_.end());
            }
        }
    }

    double quantile(double q) {
        compress();
        if (summary_.empty())
            return NAN;
        if (summary_.size() == 1)
            return summary_[0].mean;
        double index = q * (count_ - 1);
        double a = summary_[0].mean;
        double aIndex = (summary_[0].count - 1) * 0.5;
        int64_t sum = summary_[0].count;
        double b = summary_[1].mean;
        double bIndex = sum + (summary_[1].count - 1) * 0.5;
        for (size_t i = 2; i < summary_.size(); ++i) {
            if (index <= bIndex) {
                return interpolate(index, aIndex, a, bIndex, b);
            }
            a = b;
            aIndex = bIndex;
            sum += summary_[i-1].count;
            b = summary_[i].mean;
            bIndex = sum + (summary_[i].count - 1) * 0.5;
        }
        // either there are exactly two centroids or index is past the end
        return interpolate(index, aIndex, a, bIndex, b);
    }
};
