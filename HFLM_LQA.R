vecnorm = function(b){
  res = sqrt(sum(b^2))
  return(res)
}

func_fuse_logistic_lasso = function(x.fd, y, beta.basis, b0, lam, phi, eps0, eps, intercept = FALSE, clustering=FALSE, KK=0, label=NULL, max.iteration, tau1, tau2){
  library(Matrix)
  y = as.numeric(y)
  n = length(y)
  if(clustering==FALSE){
    KK = n
    label = 1:n
  }
  L = beta.basis$nbasis
  rng = getbasisrange(beta.basis)
  breaks = c(rng[1],beta.basis$params,rng[2])
  M = length(breaks) - 1
  W = array(0,dim=c(L,L,M))
  for (k in 1:M){
    temp <- inprod(beta.basis,beta.basis,rng=c(breaks[k],breaks[k+1]))
    W[,,k] <- temp
  }
  Z = inprod(x.fd$basis, beta.basis)
  tmp = t(x.fd$coefs)%*%Z
  Big.Gamma = matrix(0, nrow = n, ncol = L*KK)
  for (i in 1:n) {
    Big.Gamma[i , (L*label[i]-L+1):(L*label[i])] = tmp[i,]
  }
  R1 = eval.penalty(beta.basis,int2Lfd(2))
  list2 = NULL
  for (i in 1:KK){
    list2[[i]] = R1
  }
  R = as.matrix(bdiag(list2))
  R0 = diag(KK*L)
  if(intercept==TRUE){
    Big.Gamma = cbind(1, Big.Gamma)
    R = rbind(0, cbind(0, R))
    R0 = rbind(0, cbind(0, R0))
  }
  c0 = sqrt(rng[2]-rng[1])/sqrt(M)/2
  it = 1
  b.l.plus1 = b0
  while (it<=max.iteration) {
    b.l = b.l.plus1
    b.mat = matrix(b.l[(length(b.l)-KK*L+1):length(b.l)], nrow = KK, byrow = TRUE)
    V.l = matrix(0,nrow = KK*L, ncol = KK*L)
    for (i in seq_len(KK-1)) {
      for (j in (i+1):KK) {
        b.ij = b.mat[i,]-b.mat[j,]
        W.sum = matrix(0, nrow = L, ncol = L)
        for (k in 1:M) {
          tmp = as.numeric(t(b.ij)%*%W[,,k]%*%b.ij)
          W.sum = W.sum + W[,,k]/max(sqrt(tmp), eps)
        }
        W.sum = 2*c0*W.sum
        V.l[(i*L-L+1):(i*L) , (i*L-L+1):(i*L)] = V.l[(i*L-L+1):(i*L) , (i*L-L+1):(i*L)] + W.sum
        V.l[(i*L-L+1):(i*L) , (j*L-L+1):(j*L)] = V.l[(i*L-L+1):(i*L) , (j*L-L+1):(j*L)] - W.sum
        V.l[(j*L-L+1):(j*L) , (i*L-L+1):(i*L)] = V.l[(j*L-L+1):(j*L) , (i*L-L+1):(i*L)] - W.sum
        V.l[(j*L-L+1):(j*L) , (j*L-L+1):(j*L)] = V.l[(j*L-L+1):(j*L) , (j*L-L+1):(j*L)] + W.sum
      }
    }
    if(intercept==TRUE){
      V.l = rbind(0, cbind(0, V.l))
    }
    pen.mat = phi*R + eps0*R0 + lam*V.l
    k.grad = 1
    b.grad.k.plus1 = b.l
    while (k.grad<=max.iteration) {
      b.grad.k = b.grad.k.plus1
      p.i = as.vector(1/(1+exp(-Big.Gamma%*%b.grad.k)))
      w.i = pmax(p.i*(1-p.i), 0.001*0.999)
      p.k = -t(Big.Gamma)%*%(y-p.i)/n + 2*pen.mat%*%b.grad.k
      hess.k = t(Big.Gamma)%*%(w.i*Big.Gamma)/n + 2*pen.mat
      b.grad.k.plus1 = b.grad.k - chol2inv(chol(hess.k))%*%p.k
      k.grad = k.grad+1
      if(vecnorm(b.grad.k.plus1-b.grad.k)/max(1, vecnorm(b.grad.k)) < tau1){
        break
      }
    }
    b.l.plus1 = b.grad.k.plus1
    if(vecnorm(b.l.plus1-b.l)/max(1, vecnorm(b.l)) < tau1){
      break
    }
    it = it+1
  }
  b.hat = as.vector(b.l.plus1)
  fuse.intercept = 0
  if(intercept==TRUE){
    fuse.intercept = b.hat[1]
    b.hat = b.hat[-1]
  }
  b.res = matrix(b.hat, nrow = KK, byrow = TRUE)
  grp = cutree(hclust(dist(b.res), method = "single"), h = tau2)
  b.res = as.data.frame(b.res)
  res.final = list(b.hat = b.hat, b.res = b.res, grp = grp, intercept = fuse.intercept)
  return(res.final)
}
