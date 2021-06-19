Changes in sparseFLMM v0.4.1.tar.gz 
* fix bug in computation of FPCs using FAMM approach when N_B is zero

Changes in sparseFLMM v0.4.0.tar.gz (Update date: 2021-01-17):
* add option to allow for nested model structure
* add option to allow skew-symmetric and cylcic smooths


Changes in sparseFLMM v0.3.1.tar.gz (Update date: 2020-10-11):
* depend on higher version of refund


Changes in sparseFLMM v0.2.2.tar.gz (Update date: 2018-03-26):
* replace rBind and cBind by rbind and cbind as they are depricated.

  
  Changes in sparseFLMM v0.2.1.tar.gz (Release date: 2017-09-23):
* fix for triangular covariance estimation; only use products
  \tilde{y}_{ij}\tilde{y}_{i^\prime j^\prime} with 
  t_{ij} \leq t_{i^\prime j^\prime} when i \leq i^\prime and with 
  t_{ij} < t_{i^\prime j^\prime} otherwise. 
