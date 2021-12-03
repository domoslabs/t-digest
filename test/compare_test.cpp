#include <iostream>
#include <cassert>
#include "TDigestHistogram.h"
#include "tdigest.h"
#define STREAM_SIZE 1000000

static inline double randMToN(double M, double N)
{
    return M + (rand() / (RAND_MAX / (N - M)));
}


int main()
{

    td_histogram_t *td_ref = td_new(500);
    TDigestHistogram td = TDigestHistogram(500);
    double seeds[STREAM_SIZE];
    for (double & seed : seeds)
    {
        seed = randMToN(0, 100);
    }

    for (double seed : seeds)
    {
        td_add(td_ref, seed, 1);
        td.add(seed, 1);
    }
    assert(td_ref->merged_nodes+td_ref->unmerged_nodes == td.getCount());
    td_compress(td_ref);
    td.compress();
    printf("\n");
    assert(td_ref->merged_nodes+td_ref->unmerged_nodes == td.getCount());
    for (int i = 0; i <= 100; i += 1)
    {
        auto ref_out = td_quantile(td_ref, i / 100.0);
        auto out = td.findQuantile(i / 100.0);
        assert(ref_out == out);
    }
    std::cout << "Test passed!" << std::endl;
}