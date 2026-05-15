# CCEX_HSP_Breaks
CCEX estimator for heterogeneous spatial panel data with multiple structural breaks and multifactor errors. Implements break date detection via dynamic programming, sequential F-tests for determining break count, and CCEX-IV for regime-specific parameter estimation. Based on Dai, Wang &amp; Zheng (2026).

get_ssr_matrix():
Computes SSR matrix for all possible segments after CCEX filtering.

run_dp():
Dynamic programming to find optimal break partitions for up to 3 breaks.

backtrack():
Traces back break dates from DP output.

postestimate():
Estimates regime-specific parameters using CCEX-IV after breaks are identified.

calc_SupF_Pooled():
Computes SupF test for m breaks vs. 0 breaks.

calc_F_seq_Pooled():
Computes sequential F-test for m vs. m+1 breaks.

CCEX_HSP_Breaks()
Main function: determines number of breaks via sequential tests , then returns break dates and parameter estimates.

Reference:
Dai, Wang &amp; Zheng . Estimation of Heterogeneous Spatial Panel Data Models with Multiple Structural Breaks and a Multifactor Error Structure. Economics Letters, 2016.

