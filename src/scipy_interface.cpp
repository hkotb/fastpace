#include "scipy_interface.hpp"
#include "scipy/_lib/boost_math/include/boost/math/distributions/binomial.hpp"

// Modified check_dist_and_k function in scipy binomial.hpp to allow k-1 when k=0 (The target is survival function plus pmf, not survival alone)
double calculate_binomial_pmf_plus_survival(int k, int n, double p){
    boost::math::binomial_distribution<> dist(n, p);
    return boost::math::cdf(
        boost::math::complement(dist, k-1)
    );
}
