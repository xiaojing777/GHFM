library(fda)
library(MASS)

#calculate L2 norm
vecnorm = function(b){
  res = sqrt(sum(b^2))
  return(res)
}

#Homogeneous functional regression
ord_func_reg = function(x.fd, y, beta.basis, phi) {
  # Ensure y is properly formatted as a single-column matrix
  y <- matrix(y, ncol = 1)
  
  # Calculate necessary matrices
  Z0 <- inprod(x.fd$basis, beta.basis)
  Z <- t(x.fd$coefs) %*% Z0
  
  # Regularization matrix with small diagonal adjustment for stability
  R <- eval.penalty(beta.basis, int2Lfd(2))
  diag_adjust <- 0.001 * diag(nrow(R))
  
  # Check matrix dimensions for compatibility
  if(ncol(Z) != nrow(R)) {
    stop("Dimension mismatch between Z and R matrices")
  }
  
  # Calculate coefficients with regularization
  b.hat <- solve(t(Z) %*% Z + phi * R + diag_adjust) %*% t(Z) %*% y
  
  # Calculate intercept
  X.train.26 <- t(x.fd$coefs)
  alpha <- mean(y - X.train.26 %*% Z0 %*% b.hat)
  
  return(list(b.hat = b.hat, alpha = alpha))
}

replicate_rows <- function(mat, k) {
  # Split the matrix into a list of rows
  rows <- split(mat, row(mat))
  
  # Replicate each row k times and recombine them into a new matrix
  new_mat <- do.call(rbind, lapply(rows, function(row) replicate(k, row, simplify = FALSE)))
  
  return(new_mat)
}

#Homogeneous functional regression, return p-value
my_ord_func_reg_p_val = function(x.fd, y, beta.basis, phi){
  y = as.matrix(y)
  Z0 = inprod(x.fd$basis,beta.basis)
  Z = t(x.fd$coefs)%*%Z0
  R = eval.penalty(beta.basis,int2Lfd(2))
  
  b.hat = solve(t(Z)%*%Z+phi*R)%*%t(Z)%*%y
  
  X.train.26 = t(x.fd$coefs)
  alpha = mean(y - X.train.26%*%Z0%*%b.hat)
  
  y.hat = Z %*% b.hat + alpha
  
  y.res = y - y.hat
  SSE1 = sum(y.res^2)
  SSE0 = sum((y-mean(y))^2)
  
  ny = length(y)
  nbasis = beta.basis$nbasis
  (Fratio = ((SSE0-SSE1)/(nbasis+1-1)/(SSE1/(ny-nbasis))))
  
  F.val = qf(0.95, nbasis, ny-nbasis-1)
  
  p.val = 1-pf(Fratio, nbasis-1, ny-nbasis)
  
  res = list(b.hat=b.hat, alpha = alpha,  y.hat = y.hat, p.val = p.val, F.val = F.val)
  return(res)
}

#Clustering wrt beta(t)
beta_clustering = function(K = 100, X, y, beta.basis, phi, max.iteration = 200, 
                           tol = 0.005, grp = NULL, x.basis.num = 15) {
  
  # Input validation
  n <- nrow(X)
  y <- matrix(y, ncol = 1)  # Ensure y is n×1 matrix
  if(nrow(y) != n) stop("X and y have different number of observations")
  
  # Initialize clustering if not provided
  if(is.null(grp)) {
    k2 <- ClusterR::KMeans_rcpp(X, clusters = K, num_init = 5, 
                                max_iters = 100, initializer = 'kmeans++')
    grp <- k2$clusters
  }
  
  # Create basis for functional data
  x.basis <- create.bspline.basis(rangeval = c(0, ncol(X)), 
                                  nbasis = x.basis.num, norder = 4)
  time_points <- seq(0, ncol(X), length.out = ncol(X))
  
  # Convert X to functional data
  X.mat <- as.matrix(X)
  X.fd <- Data2fd(y = t(X.mat), basisobj = x.basis, argvals = time_points)
  X.fd.mat <- t(X.fd$coefs)
  Z0 <- inprod(X.fd$basis, beta.basis)
  
  # Verify matrix dimensions
  stopifnot(ncol(X.fd.mat) == nrow(Z0))  # Should be x.basis.num × beta.basis$nbasis
  
  # Clustering iterations
  it <- 1
  while (it < max.iteration) {
    cat("Pre-clustering Iteration", it, "\n")
    grp.old <- grp
    
    # Calculate cluster-specific coefficients
    b.all <- matrix(NA, nrow = K, ncol = beta.basis$nbasis)
    for (l in 1:K) {
      idx <- which(grp == l)
      if(length(idx) == 0) next  # Skip empty clusters
      
      X.mat.l <- t(X.mat[idx, , drop = FALSE])
      X.fd.l <- Data2fd(y = X.mat.l, basisobj = x.basis, 
                        argvals = time_points)
      y.l <- matrix(y[idx], ncol = 1)
      
      reg_result <- tryCatch(
        ord_func_reg(X.fd.l, y.l, beta.basis, phi),
        error = function(e) NULL
      )
      
      if(!is.null(reg_result)) b.all[l, ] <- reg_result$b.hat
    }
    
    # Remove NA clusters
    valid_clusters <- which(!is.na(b.all[,1]))
    if(length(valid_clusters) == 0) stop("All clusters failed")
    b.all <- b.all[valid_clusters, , drop = FALSE]
    
    # Reassign observations with dimension-safe calculations
    for (i in 1:n) {
      distances <- sapply(1:nrow(b.all), function(j) {
        # Ensure all matrices are properly aligned
        x_vec <- matrix(X.fd.mat[i,], nrow = 1)  # 1 × x.basis.num
        b_vec <- matrix(b.all[j,], ncol = 1)      # beta.basis$nbasis × 1
        
        # Matrix multiplication: (1 × x.basis.num) %*% (x.basis.num × beta.basis$nbasis) %*% (beta.basis$nbasis × 1)
        tmp <- x_vec %*% Z0 %*% b_vec
        (y[i] - tmp[1,1])^2  # Extract scalar value
      })
      grp[i] <- valid_clusters[which.min(distances)]
    }
    
    if(it >= 35 && sum(grp.old != grp) <= tol * n) break
    it <- it + 1
  }
  
  return(grp)
}

#Choose tuning parameters
HFGM_tune = function(x.fd, y, beta.basis, b0, lam, phi, max.iteration, tol=1e-5){
  err.score = matrix(0, nrow=length(lam), ncol=length(phi))
  
  n = length(y)
  L = beta.basis$nbasis
  Z = inprod(x.fd$basis,beta.basis)
  x.coef = x.fd$coefs
  x.coef = t(x.coef)
  
  for (i1 in 1:length(lam)) {
    for (j1 in 1:length(phi)) {
      temp = func_fuse_lasso(x.fd, y, beta.basis, b0, lam=lam[i1], phi=phi[j1], max.iteration=max.iteration, tol=1e-3)
      b = temp[[1]]
      
      y.fit = matrix(0, nrow = n, ncol = 1)
      
      for (k1 in 1:n) {
        x.k = as.matrix(x.coef[k1,])
        x.k = t(x.k)
        b.k = b[(L*k1-L+1):(L*k1) , ]
        y.fit[k1,] = x.k%*%Z%*%b.k
      }
      
      print(y.fit)
      
      resid = y - y.fit
      rss = sum(resid^2)
      #df = n*L
      df = sum(diag(temp$projMat))
      
      err = n*log(rss/n) + log(n)*df
      
      print(err)
      
      err.score[i1,j1] = err
    }
  }
  
  return(err.score)
}

b.func = function(t){
  beta.basis = create.bspline.basis(rangeval = c(0,24), nbasis = L, norder = 4)
  mat <- eval.basis(t,beta.basis)
  res = t(mat)
  return(res)
}

sim.b.spl.sq.norm.ord = function(beta.coef, beta.basis, ord.res){
  Z1 = inprod(beta.basis, beta.basis)
  n = nrow(beta.coef)
  tmp = matrix(0, nrow = n, ncol = 1)
  for (i in 1:n) {
    tmp[i,] = t(beta.coef[i,] - ord.res)%*%Z1%*%(beta.coef[i,] - ord.res)
  }
  res = sqrt(sum(tmp))
  return(res)
}

sim.b.spl.sq.norm.true.beta = function(beta.coef, beta.basis){
  Z1 = inprod(beta.basis, beta.basis)
  n = nrow(beta.coef)
  tmp = matrix(0, nrow = n, ncol = 1)
  for (i in 1:n) {
    tmp[i,] = t(beta.coef[i,])%*%Z1%*%beta.coef[i,]
  }
  res = sqrt(sum(tmp))
  return(res)
}
                              
fuse.pred = function(fuse.tmp){
  fuse.mat = as.matrix(fuse.tmp$b.res)
  fuse.int = fuse.tmp$intercept
  y.fuse.pred = matrix(0, nrow = length(y.test), ncol = 1)
  for (i in 1:length(y.test)) {
    if(length(fuse.int)>1){
      y.fuse.pred[i,] = as.numeric((t(x.test.coef[i,]))%*%Z%*%fuse.mat[i,]) + fuse.int[i]
    }else{
      y.fuse.pred[i,] = as.numeric((t(x.test.coef[i,]))%*%Z%*%fuse.mat[i,]) + fuse.int
      }
  }
  return(y.fuse.pred)
}

sim.b.spl.sq.norm.fuse = function(beta.coef, beta.basis, fuse.res){
  Z1 = inprod(beta.basis, beta.basis)
  n = nrow(beta.coef)
  tmp = matrix(0, nrow = n, ncol = 1)
  for (i in 1:n) {
    tmp[i,] = t(beta.coef[i,]-fuse.res[i,])%*%Z1%*%(beta.coef[i,] - fuse.res[i,])
  }
  res = sqrt(sum(tmp))
  return(res)
}


fuse_group = function(beta.basis, b.mat, tol=0){
  Z = inprod(beta.basis, beta.basis)
  if(tol!=0){
    tau = tol
  }else{
    tmp = b.mat[1,]
    tau = sqrt(t(tmp)%*%Z%*%tmp)*0.9
  }
  
  n = nrow(b.mat)
  grp = 1:n
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      tmp = b.mat[i,] - b.mat[j,]
      difference = sqrt(t(tmp)%*%Z%*%tmp)
      if(difference<=tau){
        grp[j] = grp[i]
        b.mat[j,] = b.mat[i,]
      }
    }
  }
  
  res = list(group = grp, b.hat = b.mat)
  return(res)
}


replicate_rows <- function(mat, k) {
  # Split the matrix into a list of rows
  rows <- split(mat, row(mat))
  
  # Replicate each row k times and recombine them into a new matrix
  new_mat <- do.call(rbind, lapply(rows, function(row) replicate(k, row, simplify = FALSE)))
  
  return(new_mat)
}



# Resp for the case of 3/4 groups
resp = function(clus_num, y, x, beta.basis, D=seq(0, 24, length.out=24)){
  n = length(y)
  L = beta.basis$nbasis
  
  num_clusters <- clus_num
  response_clusters <- kmeans(y, centers = num_clusters)
  cluster_assignments <- response_clusters$cluster
  
  grp1 = which(cluster_assignments==1)
  grp2 = which(cluster_assignments==2)
  grp3 = which(cluster_assignments==3)
  if(clus_num==4){
    grp4 = which(cluster_assignments==4)
  }
  
  # each column of x is a subject
  x1 = x[,grp1]
  x2 = x[,grp2]
  x3 = x[,grp3]
  if(clus_num==4){
    x4 = x[,grp4]
  }
  
  y1 = y[grp1]
  y2 = y[grp2]
  y3 = y[grp3]
  if(clus_num==4){
    y4 = y[grp4]
  }
  
  
  x1.fd = Data2fd(argvals=D,y=x1)
  x2.fd = Data2fd(argvals=D,y=x2)
  x3.fd = Data2fd(argvals=D,y=x3)
  if(clus_num==4){
    x4.fd = Data2fd(argvals=D,y=x4)
  }
  
  
  
  ord.tmp1 = ord_func_reg(x1.fd,y1,beta.basis, phi = 10)$b.hat
  ord.int1 = ord_func_reg(x1.fd,y1,beta.basis, phi = 10)$alpha
  ord.tmp2 = ord_func_reg(x2.fd,y2,beta.basis, phi = 10)$b.hat
  ord.int2 = ord_func_reg(x2.fd,y2,beta.basis, phi = 10)$alpha
  ord.tmp3 = ord_func_reg(x3.fd,y3,beta.basis, phi = 10)$b.hat
  ord.int3 = ord_func_reg(x3.fd,y3,beta.basis, phi = 10)$alpha
  if(clus_num==4){
    ord.tmp4 = ord_func_reg(x4.fd,y4,beta.basis, phi = 10)$b.hat
    ord.int4 = ord_func_reg(x4.fd,y4,beta.basis, phi = 10)$alpha
  }
  
  
  b.hat = matrix(0, nrow = n, ncol = L)
  alpha.hat = matrix(0, nrow = n)
  for (i in grp1) {
    b.hat[i,] = ord.tmp1
    alpha.hat[i,] = ord.int1
  }
  for (i in grp2) {
    b.hat[i,] = ord.tmp2
    alpha.hat[i,] = ord.int2
  }
  for (i in grp3) {
    b.hat[i,] = ord.tmp3
    alpha.hat[i,] = ord.int3
  }
  if(clus_num==4){
    for (i in grp4) {
      b.hat[i,] = ord.tmp4
      alpha.hat[i,] = ord.int4
    }
  }
  
  res = list(b.hat = b.hat, alpha.hat = alpha.hat)
  return(res)
  
}
