//
// Created by vladim0105 on 03.12.2021.
//

#ifndef T_DIGEST_TDIGESTHISTOGRAM_H
#define T_DIGEST_TDIGESTHISTOGRAM_H

#include <cstdlib>
#include <algorithm>
#include <vector>
#include <cmath>

struct centroid{
    double mean = 0;
    double weight = 0;
};
static inline bool centroid_compare(const centroid a, const centroid b){
    return a.mean < b.mean;
}
class TDigestHistogram {
private:
    int compression = 0;
    // The size of the array holding all the values, faster to use a fixed size array than a dynamic expanding array.
    int capacity = 0;
    int merged_nodes = 0;
    int unmerged_nodes = 0;
    double merged_weight = 0;
    double unmerged_weight = 0;
    double min = __DBL_MAX__;
    double max = __DBL_MIN__;
    std::vector<centroid> nodes;
    static double weighted_average_sorted(double x1, double w1, double x2, double w2) {
        const double x = (x1 * w1 + x2 * w2) / (w1 + w2);
        return std::max(x1, std::min(x, x2));
    }
    static double weighted_average(double x1, double w1, double x2, double w2) {
        if (x1 <= x2) {
            return weighted_average_sorted(x1, w1, x2, w2);
        } else {
            return weighted_average_sorted(x2, w2, x1, w1);
        }
    }
public:
    explicit TDigestHistogram(int compression);
    void add(centroid cent);
    void add(double mean, double weight);
    void add(TDigestHistogram tDigestHistogram);
    int getCount(){
        return merged_nodes+unmerged_nodes;
    }
    double getMin() const{
        return min;
    }
    double getMax() const{
        return max;
    }
    std::vector<centroid> getNodes(){
        return nodes;
    }
    void compress();
    double findQuantile(double q);
};

#endif //T_DIGEST_TDIGESTHISTOGRAM_H
