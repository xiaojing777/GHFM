# GHFM Project

This repository contains R scripts for the GHFM project.

## Files
- `HFGM.R`: This R file includes the function of conducting HFGM.
- `HFLM.R`: This R file includes the function of conducting HFLM.
- `Tune_and_pre_clus.R`: This R file includes how to conduct pre-clustering and how to tune parameters to HFGM.

---


## 📘 Example: Simulation and Pre-Clustering

The following code demonstrates how our proposed pre-clustering method generates a finer partition of the true subgroup structure.

### Code
```R

# Load required functions 
# Also, obtain your tuned parameters, e.g., phi
# Set values for:
#   T       : total number of time points
#   n       : sample size
#   sgm     : noise level / difference between true subgroup effects
#   K       : number of pre-clustering groups
#   max.iteration : maximum iterations for the pre-clustering algorithm
#   L_beta  : number of B-spline basis functions for approximating beta_t in the model
#   L_x     : number of B-spline basis functions for representing X
#   grp      : initial value for pre-clustering

beta.basis = create.bspline.basis(rangeval = c(0, T), 
                                  nbasis = L_beta, norder = 4)

# Simulate raw data X as a matrix of dimension (n x T)
x = matrix(rnorm(n * T, mean = 3, sd = 1), n, T)

# Construct a B-spline basis for X and obtain the functional data object.
x.basis = create.bspline.basis(rangeval = c(0, ncol(x)), 
                               nbasis = L_x, norder = 4)
x.fd = Data2fd(y = t(x), basisobj = x.basis)
xt.coef = t(x.fd$coefs)
xt.coef = as.matrix(xt.coef)

# Generate true group-specific beta coefficients with added noise
set.seed(1314)
e1 = rnorm(L, mean = 0, sd = sgm)
set.seed(1315)
e2 = rnorm(L, mean = 0, sd = sgm)

true.b.1 = rep(10, L) + e1
true.b.2 = rep(-10, L) + e2

# Create group-specific beta coefficient matrix for all subjects.
beta.coef = matrix(0, nrow = n, ncol = L_beta)
for (i in 1:(n/2)) {
  beta.coef[i,] = true.b.1
}
for (i in (n/2 + 1):n) {
  beta.coef[i,] = true.b.2
}

# Compute the inner product matrix between X's basis and beta's basis.
Z = inprod(x.fd$basis, beta.basis)

# Generate outcome Y using the functional predictor and the true beta coefficients.
y = matrix(0, nrow = n, ncol = 1)
y.epsilon = rnorm(n)
for (i in 1:n) {
  y[i] = t(xt.coef[i,]) %*% Z %*% beta.coef[i,]
}
y = y + y.epsilon
y = as.matrix(y)

# Apply the pre-clustering method.
pre_grp = beta_clustering(K = K, x, y, beta.basis = beta.basis, 
                          phi, max.iteration, grp)

# Combine subject indices, estimated groups, and true groups for comparison.
group_check = cbind(1:n, pre_grp)
true_grp = c(rep(1, n/2), rep(2, n/2))
group_check = cbind(group_check, true_grp)


---


## 📘 Example: Running HFGM

This code shows how to apply our proposed GHFM method to simulated data

### Code Example
```R

# Load required functions 
# Also, obtain your tuned parameters, e.g., phi, lambda
# Set values for:
#   T       : total number of time points
#   n       : sample size
#   sgm     : noise level / difference between true subgroup effects
#   max.iteration : maximum iterations for the GHFM algorithm
#   L_beta  : number of B-spline basis functions for approximating beta_t in the model
#   L_x     : number of B-spline basis functions for representing X (optional)


library(MASS)
library(fda)


beta.basis = create.bspline.basis(rangeval = c(0,T), nbasis = L_beta, norder = 4)

x = matrix(rnorm(n*(T), mean = 3, sd=1),n,(T))

x.basis = create.bspline.basis(rangeval = c(0, ncol(x)), 
                               nbasis = L_x, norder = 4)
x.fd = Data2fd(y = t(x), basisobj = x.basis)

xt.coef = x.fd$coefs
xt.coef = t(xt.coef)
xt.coef = as.matrix(xt.coef)


beta_means <- c(10, 5, -5, -10)
num_groups <- length(beta_means)
n_per_group <- n / num_groups

beta.coef <- do.call(rbind, lapply(beta_means, function(mu) {
  noise <- rnorm(L_beta, mean = 0, sd = sgm)
  matrix(rep(mu, L_beta) + noise, nrow = n_per_group, ncol = L_beta, byrow = TRUE)
}))

Z = inprod(x.fd$basis,beta.basis)

y = matrix(0, nrow = n, ncol = 1)
y.epsilon = rnorm(n)
for (i in 1:n) {
  y[i] = t(xt.coef[i,])%*%Z%*%beta.coef[i,]
}
y = y+y.epsilon


fuse.tmp = func_fuse_lasso(x.fd, y, beta.basis, b0, lam, phi,
                           max.iteration)

ord.tmp = ord_func_reg(x.fd,y,beta.basis, phi2)$b.hat

sim.b.spl.sq.norm.fuse(beta.coef, beta.basis, fuse.mat)/sim.b.spl.sq.norm.true.beta(beta.coef, beta.basis)
sim.b.spl.sq.norm.ord(beta.coef, beta.basis, ord.tmp)/sim.b.spl.sq.norm.true.beta(beta.coef, beta.basis)





