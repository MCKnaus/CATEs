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


