# CATEs Estimators

Caller package to implement multiple CATE estimators.

Fork of [CATEs](https://github.com/MCKnaus/CATEs): implementation of
all estimators that are applied in the Empirical Monte Carlo Study of
[Knaus](https://mcknaus.github.io/),
[Lechner](https://www.michael-lechner.eu/) and
[Strittmatter](http://www.anthonystrittmatter.com/home) (2018). They
are based on the packages [grf](https://github.com/grf-labs/grf) and
[glmnet](https://github.com/cran/glmnet).

## Example

```R
library(pacman)
source("R/IATEs.R")
source("R/CATEs_utils.R")

p_load(glmnet, grf, caret)

# %% # Generate training sample
n = 4000; p = 20
x_t = matrix(rnorm(n * p), n, p)
tau_t = 1 / (1 + exp(-x_t[, 3]))
d_t = rbinom(n ,1, 1 / (1 + exp(-x_t[, 1] - x_t[, 2])))
y_t = pmax(x_t[, 2] + x_t[, 3], 0) + rowMeans(x_t[, 4:6]) / 2 + d_t * tau_t + rnorm(n)

# %% # Generate validation sample of same size
x_v     = matrix(rnorm(n * p), n, p)
tau_v   = 1 / (1 + exp(-x_v[, 3]))

# %% without true τ_train (for use in applications)
iates_list = IATEs(y_t, d_t, x_t, x_v)
iates_list  |>  str()
iates_df = as.data.frame(do.call(cbind, iates_list))
iates_df |> str()

# 'data.frame': 4000 obs. of  11 variables:
#  $ RF CMR       : num  0.181 0.389 1.3 1.047 0.419 ...
#  $ CF           : num  0.502 0.381 0.716 0.965 0.532 ...
#  $ RF MOM IPW   : num  0.343 0.244 2.164 0.365 0.534 ...
#  $ RF MOM DR    : num  0.399 0.353 1.066 0.414 0.483 ...
#  $ CF LC        : num  0.415 0.333 0.605 0.475 0.47 ...
#  $ LASSO CMR    : num  0.503 -0.144 1.478 0.624 0.33 ...
#  $ LASSO MOM IPW: num  0.594 0.076 1.156 0.425 0.691 ...
#  $ LASSO MOM DR : num  0.442 0.222 0.844 0.514 0.462 ...
#  $ LASSO MCM    : num  0.6586 0.0472 1.266 0.4265 0.7341 ...
#  $ LASSO MCM EA : num  0.36 0.224 0.686 0.453 0.398 ...
#  $ LASSO RL     : num  0.357 0.301 0.614 0.493 0.38 ...
#

# %% # Calculate and print mean MSEs
mMSE = colMeans((iates_df - tau_v)^2)
names(mMSE) = colnames(iates_df)
round(mMSE*1000, 2)  |> t()

```

## References

Knaus, Lechner, Strittmatter (2020). Machine Learning Estimation of Heterogeneous Causal
Effects: Empirical Monte Carlo Evidence, [arXiv](https://arxiv.org/abs/1810.13237)
