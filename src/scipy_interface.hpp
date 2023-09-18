#ifndef SCIPY_INTERFACE_HPP
#define SCIPY_INTERFACE_HPP

#ifdef __cplusplus
extern "C" {
#endif

// Function to calculate the binomial PMF plus survival (x >= k) (based on scipy implementation)
double calculate_binomial_pmf_plus_survival(int k, int n, double p);

#ifdef __cplusplus
}
#endif

#endif