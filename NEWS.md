# sparseFLMM v0.4.2

## changes
- fix Rd cross-references

# sparseFLMM v0.4.1

## changes

- fix bug in computation of FPCs using FAMM approach when N_B is zero

# sparseFLMM v0.4.0

## changes

- add option to allow for nested model structure
- add option to allow skew-symmetric and cylcic smooths


# sparseFLMM v0.3.1

## changes

- depend on higher version of refund


# sparseFLMM v0.2.2

## changes

- replace rBind and cBind by rbind and cbind as they are depricated.

  
# sparseFLMM v0.2.1

## changes

- fix for triangular covariance estimation; only use products
  \tilde{y}_{ij}\tilde{y}_{i^\prime j^\prime} with 
  t_{ij} \leq t_{i^\prime j^\prime} when i \leq i^\prime and with 
  t_{ij} < t_{i^\prime j^\prime} otherwise. 
