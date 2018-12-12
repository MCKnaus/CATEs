# CATEs
Implementation of all estimators that are applied in the Empirical Monte Carlo Study of [Knaus](https://mcknaus.github.io/), [Lechner](https://www.michael-lechner.eu/) and [Strittmatter](http://www.anthonystrittmatter.com/home) (2018). They are based on the packages [grf](https://github.com/grf-labs/grf) and [glmnet](https://github.com/cran/glmnet).

## Example
We have no permission to share the data used in the study. Therefore, the following example uses the observational data generating process of the example of [grf](https://github.com/grf-labs/grf) to illustrate how it works. The function IATEs is a wrapper for the underlying functions and uses all the default settings of the packages. If you want to have more control over the settings, use the respective functions in CATEs_utils.

```R
# Download current version from Github
library(devtools)
install_github(repo="MCKnaus/CATEs")
library(CATEs)

# Generate training sample
n = 4000; p = 20
x_tr = matrix(rnorm(n * p), n, p)
tau_tr = 1 / (1 + exp(-x_tr[, 3]))
d_tr = rbinom(n ,1, 1 / (1 + exp(-x_tr[, 1] - x_tr[, 2])))
y_tr = pmax(x_tr[, 2] + x_tr[, 3], 0) + rowMeans(x_tr[, 4:6]) / 2 + d_tr * tau_tr + rnorm(n)

# Generate validation sample of same size
x_val = matrix(rnorm(n * p), n, p)
tau_val = 1 / (1 + exp(-x_val[, 3]))

# Apply all estimators to the training sample and predict IATEs for validation sample
iates_mat = IATEs(y_tr,d_tr,x_tr,tau_tr,x_val)

# Calculate and print mean MSEs
mMSE = colMeans((iates_mat - tau_val)^2)
names(mMSE) = colnames(iates_mat)
mMSE*1000
```


## References

Knaus, Lechner, Strittmatter (2018). Machine Learning Estimation of Heterogeneous Causal
Effects: Empirical Monte Carlo Evidence, [arXiv](https://arxiv.org/abs/1810.13237)
