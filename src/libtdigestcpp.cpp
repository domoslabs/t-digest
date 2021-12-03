//
// Created by vladim0105 on 03.12.2021.
//

#include "../include/libtdigestcpp.h"

TDigestHistogram::TDigestHistogram(int compression) {
    this->compression = compression;
    this->capacity = (6*compression)+10; // No clue to how this relation was found, but some implementations use this.
    this->nodes = std::vector<centroid>(capacity); // Initialize with a size, which it will never have to expand (speed boost).
}

void TDigestHistogram::add(centroid cent) {
    if (getCount() >= capacity){
        compress();
    }
    if(cent.mean < getMin()){
        min = cent.mean;
    }
    if(cent.mean > getMax()){
        max = cent.mean;
    }
    int pos = getCount();
    nodes.at(pos) = cent;
    unmerged_nodes++;
    unmerged_weight += cent.weight;
}

void TDigestHistogram::add(TDigestHistogram tDigestHistogram) {
    compress();
    tDigestHistogram.compress();
    for(int i = 0; i < tDigestHistogram.getNodes().size(); i++){
        add(tDigestHistogram.getNodes().at(i));
    }
}

void TDigestHistogram::add(double mean, double weight) {
    centroid cent;
    cent.mean = mean;
    cent.weight = weight;
    add(cent);
}

void TDigestHistogram::compress() {
    if(unmerged_nodes == 0)
        return;
    std::sort(nodes.begin(), nodes.begin()+getCount(), centroid_compare);
    const double total_weight = unmerged_weight+merged_weight;
    const double denom = 2*M_PI*total_weight*log(total_weight);
    const double normalizer = compression/denom;
    int cur = 0;
    double weight_so_far = 0;
    for(int i = 1; i < getCount(); i++){
        const double proposed_weight = nodes.at(cur).weight+nodes.at(i).weight;
        const double z = proposed_weight * normalizer;
        const double q0 = weight_so_far/total_weight;
        const double q2 = (weight_so_far+proposed_weight)/total_weight;
        const bool should_add = (z <= (q0 * (1 - q0))) && (z <= (q2 * (1 - q2)));
        if (should_add) {
            nodes.at(cur).weight += nodes.at(i).weight;
            const double delta = nodes.at(i).mean - nodes.at(cur).mean;
            const double weighted_delta = (delta * nodes.at(i).weight) / nodes.at(cur).weight;
            nodes.at(cur).mean += weighted_delta;
        } else {
            weight_so_far += nodes.at(cur).weight;
            cur++;
            nodes.at(cur).weight = nodes.at(i).weight;
            nodes.at(cur).mean= nodes.at(i).mean;
        }
        if (cur != i) {
            nodes.at(i).mean = 0.0;
            nodes.at(i).weight = 0.0;
        }
    }
    merged_nodes = cur + 1;
    merged_weight = total_weight;
    unmerged_nodes = 0;
    unmerged_weight = 0;
}

double TDigestHistogram::findQuantile(double q) {
    compress();
    // q should be in [0,1]
    if (q < 0.0 || q > 1.0 || getCount() == 0) {
        return NAN;
    }
    // with one data point, all quantiles lead to Rome
    if (getCount() == 1) {
        return nodes.at(0).mean;
    }

    // if values were stored in a sorted array, index would be the offset we are interested in
    const double index = q * merged_weight;

    // beyond the boundaries, we return min or max
    // usually, the first centroid will have unit weight so this will make it moot
    if (index < 1) {
        return getMin();
    }

    // we know that there are at least two centroids now
    const int n = getCount();

    // if the left centroid has more than one sample, we still know
    // that one sample occurred at min so we can do some interpolation
    const double left_centroid_weight = nodes.at(0).weight;
    const double left_centroid_mean = nodes.at(0).mean;
    if (left_centroid_weight > 1 && index < left_centroid_weight / 2) {
        // there is a single sample at min so we interpolate with less weight
        return getMin() + (index - 1) / (left_centroid_weight / 2 - 1) * (left_centroid_mean - getMin());
    }

    // usually the last centroid will have unit weight so this test will make it moot
    if (index > merged_weight - 1) {
        return getMax();
    }

    // if the right-most centroid has more than one sample, we still know
    // that one sample occurred at max so we can do some interpolation
    const double right_centroid_weight = nodes.at(n-1).weight;
    const double right_centroid_mean = nodes.at(n-1).mean;
    if (right_centroid_weight > 1 && merged_weight - index <= right_centroid_weight / 2) {
        return getMax() - (merged_weight - index - 1) / (right_centroid_weight / 2 - 1) *
                          (getMax() - right_centroid_mean);
    }

    // in between extremes we interpolate between centroids
    double weightSoFar = left_centroid_weight / 2;
    for (int i = 0; i < n - 1; i++) {
        const double node_weight = nodes.at(i).weight;
        const double node_weight_next =  nodes.at(i+1).weight;
        const double node_mean =  nodes.at(i).mean;
        const double node_mean_next = nodes.at(i+1).mean;
        const double dw = (node_weight + node_weight_next) / 2;
        if (weightSoFar + dw > index) {
            // centroids i and i+1 bracket our current point
            // check for unit weight
            double leftUnit = 0;
            if (node_weight == 1) {
                if (index - weightSoFar < 0.5) {
                    // within the singleton's sphere
                    return node_mean;
                } else {
                    leftUnit = 0.5;
                }
            }
            double rightUnit = 0;
            if (node_weight_next == 1) {
                if (weightSoFar + dw - index <= 0.5) {
                    // no interpolation needed near singleton
                    return node_mean_next;
                }
                rightUnit = 0.5;
            }
            const double z1 = index - weightSoFar - leftUnit;
            const double z2 = weightSoFar + dw - index - rightUnit;
            return weighted_average(node_mean, z2, node_mean_next, z1);
        }
        weightSoFar += dw;
    }

    // weightSoFar = totalWeight - weight[n-1]/2 (very nearly)
    // so we interpolate out to max value ever seen
    const double z1 = index - merged_weight - right_centroid_weight / 2.0;
    const double z2 = right_centroid_weight / 2 - z1;
    return weighted_average(right_centroid_mean, z1, getMax(), z2);
}
