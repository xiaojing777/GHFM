# GHFM Project

This repository contains R scripts for the GHFM project.

## Files
- `HFGM.R`: This R file includes the function of conducting HFGM.
- `HFLM.R`: This R file includes the function of conducting HFLM.
- `Tune_and_pre_clus.R': This R file includes how to conduct pre-clustering and how to tune parameters to HFGM.


# Example Usage

Here is a simple example to demonstrate how to perform HFGM:

### Code Example
```R
# Load the required functions
# Get the tuned parameters lambda and phi
# Specify any value of L, n, sgm, sgm2


library(MASS)
library(fda)


beta.basis = create.bspline.basis(rangeval = c(0,24), nbasis = L, norder = 4)

D = seq(0, 24, length.out=24)

x = matrix(rnorm(n*(24), mean = 3, sd=1),n,(24))
x = t(x)  
x.fd = Data2fd(argvals=D,y=x)
xt.coef = x.fd$coefs
xt.coef = t(xt.coef)
xt.coef = as.matrix(xt.coef)

x.test = matrix(rnorm((n)*(24), mean = 3, sd = 1), (n), (24))
x.test = t(x.test)
x.test.fd = Data2fd(argvals = D, y=x.test)
x.test.coef = as.matrix(t(x.test.fd$coefs))


e1 = rnorm(L, mean = 0, sd=sgm)
e2 = rnorm(L, mean = 0, sd=sgm)
e3 = rnorm(L, mean = 0, sd=sgm)
e4 = rnorm(L, mean = 0, sd=sgm)

true.b.1 = rep(20,  L)+e1
true.b.2 = rep(6,   L)+e2
true.b.3 = rep(-10, L)+e3
true.b.4 = rep(-40, L)+e4



beta.coef = matrix(0, nrow = n, ncol = L)
for (i in 1:(n/4)) {
  beta.coef[i,] = true.b.1
}

for (i in (n/4+1):(n/2)) {
  beta.coef[i,] = true.b.2
}

for (i in (n/2+1):((3*n)/4)) {
  beta.coef[i,] = true.b.3
}

for (i in ((3*n)/4 + 1):n) {
  beta.coef[i,] = true.b.4
}

Z = inprod(x.fd$basis,beta.basis)

y = matrix(0, nrow = n, ncol = 1)
y.epsilon = rnorm(n)
for (i in 1:n) {
  y[i] = t(xt.coef[i,])%*%Z%*%beta.coef[i,]
}
y = y+y.epsilon

y.test = matrix(0, nrow = (n), ncol = 1)
for (i in 1:(n)) {
  y.test[i] = t(x.test.coef[i,])%*%Z%*%beta.coef[i,]
}
y.test.epsilon = rnorm((n))
y.test = y.test+y.test.epsilon


b0 = matrix(0, nrow = n*L, ncol = 1)
b0[1:((n/4)*L), ] = 20
b0[((n/4)*L+1):((n/2)*L), ] = 6
b0[((n/2)*L+1):(((3*n)/4)*L), ] = -10
b0[(((3*n)/4)*L+1):(n*L), ] = -40
rdm = rnorm(n*L, mean = 0, sd=sgm2)
b0 = b0+rdm


fuse.tmp = func_fuse_lasso(x.fd, y, beta.basis, b0, lam = lambda, phi = phi,
                           max.iteration = 200)
y.fuse.pred = fuse.pred(fuse.tmp)
ord.tmp = ord_func_reg(x.fd,y,beta.basis, phi = 10)$b.hat
#ord.res = rep(ord.tmp,n)
ord.int = ord_func_reg(x.fd,y,beta.basis, phi = 10)$alpha
y.ord.pred = x.test.coef%*%Z%*%ord.tmp + ord.int
dt.ols = cbind(y, t(x))
dt.ols = as.data.frame(dt.ols)
ols = lm(V1~., data = dt.ols, )
dt.ols.test = as.data.frame(cbind(y.test, t(x.test)))
y.ols.fit = predict(ols, dt.ols.test[,-1])
resp.tmp = resp(4,y,x,beta.basis, D)$b.hat
resp.int = resp(4,y,x,beta.basis, D)$alpha.hat
y.resp.pred = matrix(0, nrow = length(y.test), ncol = 1)
for (i in 1:length(y.test)) {
  y.resp.pred[i,] = as.numeric((t(x.test.coef[i,]))%*%Z%*%resp.tmp[i,]) + resp.int[i]
}

sim.b.spl.sq.norm.fuse(beta.coef, beta.basis, fuse.mat)/sim.b.spl.sq.norm.true.beta(beta.coef, beta.basis)
sim.b.spl.sq.norm.ord(beta.coef, beta.basis, ord.tmp)/sim.b.spl.sq.norm.true.beta(beta.coef, beta.basis)
sim.b.spl.sq.norm.fuse(beta.coef, beta.basis, resp.tmp)/sim.b.spl.sq.norm.true.beta(beta.coef, beta.basis)

(pmse.fuse = vecnorm(y.test - y.fuse.pred)/n)
(pmse.ord = vecnorm(y.test - y.ord.pred)/n)
(pmse.ols = vecnorm(y.test - y.ols.fit)/n)
(pmse.resp = vecnorm(y.test - y.resp.pred)/n)
vecnorm(y.test)/n





