vecnorm = function(b){
  res = sqrt(sum(b^2))
  return(res)
}

func_fuse_logistic_lasso = function(x.fd, y, beta.basis, b0, lam, phi, intercept = FALSE, clustering=FALSE, x.penalty=FALSE, KK=0, label=NULL, max.iteration=100, tol=0.05){
  #Check whether a pre-clustering has been made
  if(clustering==TRUE){
    #sample size: n
    n = ncol(x.fd$coefs)
    
    #number of basis function: L
    L = beta.basis$nbasis
    
    #value range of basis functions: rng
    rng = getbasisrange(beta.basis)
    
    #breaking points in the range of basis functions: breaks
    breaks = c(rng[1],beta.basis$params,rng[2])
    
    #number of sub-intervals: M
    M = length(breaks) - 1
    
    #number of order: norder
    norder = L-M+1
    
    W = array(0,dim=c(L,L,M))
    for (k in 1:M){
      temp <- inprod(beta.basis,beta.basis,rng=c(breaks[k],breaks[k+1]))
      W[,,k] <- temp
    }
    
    
    Z = inprod(x.fd$basis, beta.basis)
    
    
    # Big.Gamma is nL dimension
    tmp = t(x.fd$coefs)%*%Z
    Big.Gamma = matrix(0, nrow = n*L, ncol = 1)
    for (i in 1:n) {
      Big.Gamma[((i-1)*L+1):(i*L), 1] = tmp[i,]
    }
    
    mygrp = label
    #sum.y.E.ga correspond to the first term of gradient
    sum.y.E.ga = matrix(0, nrow = KK*L, ncol = 1)
    for (i in 1:KK) {
      idx.tmp = which(mygrp==i)
      small.sum = matrix(0, nrow = L, ncol = 1)
      for (j in idx.tmp) {
        small.sum = small.sum + as.numeric(y[j])*tmp[j,]
      }
      
      sum.y.E.ga[((i-1)*L+1):(i*L), 1] = small.sum
    }
    
    
    
    R1 = eval.penalty(beta.basis,int2Lfd(2))
    list2 = NULL
    for (i in 1:KK){
      list2[[i]] = R1
    }
    library(Matrix)
    R = as.matrix(bdiag(list2))
    
    
    it = 1
    b.l.plus1 = b0
    
    while (it<max.iteration) {
      
      b.l = b.l.plus1
      
      V.l = matrix(0,nrow = KK*L, ncol = KK*L)
      
      
      ij.start.time = Sys.time()
      
      x.inprod = inprod(x.fd$basis, x.fd$basis)
      x.matt = t(x.fd$coefs)
      
      for (i in 1:KK) {
        for (j in (1:KK)[-i]) {
          
          
          b.i.l = b.l[(i*L-L+1):(i*L)]
          b.j.l = b.l[(j*L-L+1):(j*L)]
          
          tmp1 = as.matrix(x.matt[i,]-x.matt[j,])
          delta_ij = t(tmp1)%*%x.inprod%*% tmp1 + 0.01
          delta_ij = as.numeric(delta_ij)

          
          W.sum = matrix(0, nrow = L, ncol = L)
          
          grp=0
          
          for (k in 1:M) {
            tmp = t(b.i.l-b.j.l)%*%W[,,k]%*%(b.i.l-b.j.l)
            tmp = as.numeric(tmp)
            W.sum = W.sum + W[,,k]/(sqrt(tmp)+0.01)
          }
          
          if(x.penalty==TRUE){
            V.l[(i*L-L+1):(i*L) , (i*L-L+1):(i*L)] = (V.l[(i*L-L+1):(i*L) , (i*L-L+1):(i*L)] + W.sum)/delta_ij
            V.l[(i*L-L+1):(i*L) , (j*L-L+1):(j*L)] = (V.l[(i*L-L+1):(i*L) , (j*L-L+1):(j*L)] + W.sum)/delta_ij
            V.l[(j*L-L+1):(j*L) , (i*L-L+1):(i*L)] = (V.l[(j*L-L+1):(j*L) , (i*L-L+1):(i*L)] + W.sum)/delta_ij
            V.l[(j*L-L+1):(j*L) , (j*L-L+1):(j*L)] = (V.l[(j*L-L+1):(j*L) , (j*L-L+1):(j*L)] + W.sum)/delta_ij
          }else{
            V.l[(i*L-L+1):(i*L) , (i*L-L+1):(i*L)] = V.l[(i*L-L+1):(i*L) , (i*L-L+1):(i*L)] + W.sum
            V.l[(i*L-L+1):(i*L) , (j*L-L+1):(j*L)] = V.l[(i*L-L+1):(i*L) , (j*L-L+1):(j*L)] + W.sum
            V.l[(j*L-L+1):(j*L) , (i*L-L+1):(i*L)] = V.l[(j*L-L+1):(j*L) , (i*L-L+1):(i*L)] + W.sum
            V.l[(j*L-L+1):(j*L) , (j*L-L+1):(j*L)] = V.l[(j*L-L+1):(j*L) , (j*L-L+1):(j*L)] + W.sum
          }
          
          
        }
        
      }
      
      V.l[is.na(V.l)] = 0
      
      ij.end.time = Sys.time()
      ij.time = ij.end.time-ij.start.time
      print(c('ij time:', ij.time))
      print("ij loop finished")
      
      inv.start.time = Sys.time()
      
      k.grad = 0
      
      b.grad.k.plus1 = b.l
      
      while (k.grad<=50) {
        b.grad.k = b.grad.k.plus1
        
        sum.exp.mat = matrix(0, nrow = KK*L, ncol = 1)
        for (i in 1:KK) {
          idx.tmp = which(mygrp==i)
          small.tmp = matrix(0, nrow = L, ncol = 1)
          for (j in idx.tmp) {
            tmp1 = Big.Gamma[((j-1)*L+1):(j*L), 1]
            tmp2 = b.grad.k[((i-1)*L+1):(i*L)]
            tmp3 = Big.Gamma[((j-1)*L+1):(j*L), 1]/as.numeric(1+exp(-t(tmp1)%*%tmp2))
            small.tmp = tmp3 + small.tmp
          }
          
          sum.exp.mat[((i-1)*L+1):(i*L), 1] = small.tmp
        }
        
        p.k = -sum.y.E.ga + sum.exp.mat + 2*(phi*R + lam*V.l)%*%b.grad.k
        
        if(k<10){
          alpha.k = 0.001
        }else{
          alpha.k = 0.00001
        }
        
        b.grad.k.plus1 = b.grad.k - alpha.k*p.k
        k.grad = k.grad+1
        
        if(!is.na(vecnorm(p.k)) && vecnorm(p.k)<=0.01){
          break
        }
      }
      
      b.l.plus1 = b.grad.k.plus1
      
      
      inv.end.time = Sys.time()
      inv.time = inv.end.time - inv.start.time
      print(c("Inverse Time:", inv.time))
      print("Inverse Calculated")
      
      
      it = it+1
      
    }
    
    
    b.hat = b.l.plus1
    
    b.res = matrix(0, nrow = KK, ncol = L)
    
    for (i in 1:KK) {
      b.res[i,] = b.hat[(L*(i-1)+1):(L*(i-1)+L)]
    }
    
    b.res = as.data.frame(b.res)
    
    fuse.intercept = c()
    linshi = as.matrix(b.res)
      X.train.26 = t(x.fd$coefs)
      neiji = inprod(x.fd$basis, beta.basis)
      for (j in 1:KK) {
        tmp1 = which(label==j)
        logit_tmp = 1/(1+exp(-X.train.26[tmp1,]%*%neiji%*%linshi[j,]))
        logit_tmp[which(logit_tmp>0.5)] = 1
        logit_tmp[which(logit_tmp<0.5)] = 0
        tmp2 = as.matrix(y[tmp1] - logit_tmp)
        fuse.intercept[j] = mean(tmp2)
      }
    if(intercept==TRUE){
      fuse.intercept = mean(fuse.intercept)
    }
    
    
    res.final = list(b.hat = b.hat, b.res = b.res, grp = grp, intercept = fuse.intercept)
    
    
    
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #-----------------------------------------------------------------------------------------------------
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  if(clustering==FALSE){
    n = length(y)
    
    L = beta.basis$nbasis
    
    rng = getbasisrange(beta.basis)
    
    breaks = c(rng[1],beta.basis$params,rng[2])
    
    M = length(breaks) - 1
    
    norder = L-M+1
    
    
    W = array(0,dim=c(L,L,M))
    for (k in 1:M){
      temp <- inprod(beta.basis,beta.basis,rng=c(breaks[k],breaks[k+1]))
      W[,,k] <- temp
    }
    
    
    Z = inprod(x.fd$basis,beta.basis)
    
    # Big.Gamma is the summation of all E_i%*%gamma_i
    tmp = t(x.fd$coefs)%*%Z
    Big.Gamma = matrix(0, nrow = n*L, ncol = 1)
    for (i in 1:n) {
      Big.Gamma[((i-1)*L+1):(i*L), 1] = tmp[i,]
    }
    
    
    #sum.y.E.ga correspond to the first term of gradient
    sum.y.E.ga = matrix(0, nrow = n*L, ncol = 1)
    for (i in 1:n) {
      sum.y.E.ga[((i-1)*L+1):(i*L), 1] = y[i]*tmp[i,]
    }
    
    
    
    R1 = eval.penalty(beta.basis,int2Lfd(2))
    list2 = NULL
    for (i in 1:n){
      list2[[i]] = R1
    }
    library(Matrix)
    R = as.matrix(bdiag(list2))
    
    
    it = 1
    b.l.plus1 = b0
    
    grp = 1:n
    
    x.inprod = inprod(x.fd$basis, x.fd$basis)
    x.matt = t(x.fd$coefs)
    
    while (it<max.iteration) {
      
      b.l = b.l.plus1
      
      V.l = matrix(0,nrow = n*L, ncol = n*L)
      
      
      ij.start.time = Sys.time()
      
      x.inprod = inprod(x.fd$basis, x.fd$basis)
      x.matt = t(x.fd$coefs)
      
      
      for (i in 1:n) {
        for (j in (1:n)[-i]) {
          
          
          b.i.l = b.l[(i*L-L+1):(i*L)]
          b.j.l = b.l[(j*L-L+1):(j*L)]
          
          
          
          W.sum = matrix(0, nrow = L, ncol = L)
          
          
          for (k in 1:M) {
            tmp = t(b.i.l-b.j.l)%*%W[,,k]%*%(b.i.l-b.j.l)
            tmp = as.numeric(tmp)
            W.sum = W.sum + W[,,k]/(sqrt(tmp)+0.01)
          }
          
          
          if(x.penalty==TRUE){
            V.l[(i*L-L+1):(i*L) , (i*L-L+1):(i*L)] = (V.l[(i*L-L+1):(i*L) , (i*L-L+1):(i*L)] + W.sum)/delta_ij
            V.l[(i*L-L+1):(i*L) , (j*L-L+1):(j*L)] = (V.l[(i*L-L+1):(i*L) , (j*L-L+1):(j*L)] + W.sum)/delta_ij
            V.l[(j*L-L+1):(j*L) , (i*L-L+1):(i*L)] = (V.l[(j*L-L+1):(j*L) , (i*L-L+1):(i*L)] + W.sum)/delta_ij
            V.l[(j*L-L+1):(j*L) , (j*L-L+1):(j*L)] = (V.l[(j*L-L+1):(j*L) , (j*L-L+1):(j*L)] + W.sum)/delta_ij
          }else{
            V.l[(i*L-L+1):(i*L) , (i*L-L+1):(i*L)] = V.l[(i*L-L+1):(i*L) , (i*L-L+1):(i*L)] + W.sum
            V.l[(i*L-L+1):(i*L) , (j*L-L+1):(j*L)] = V.l[(i*L-L+1):(i*L) , (j*L-L+1):(j*L)] + W.sum
            V.l[(j*L-L+1):(j*L) , (i*L-L+1):(i*L)] = V.l[(j*L-L+1):(j*L) , (i*L-L+1):(i*L)] + W.sum
            V.l[(j*L-L+1):(j*L) , (j*L-L+1):(j*L)] = V.l[(j*L-L+1):(j*L) , (j*L-L+1):(j*L)] + W.sum
          }
          
        }
        
      }
      
      V.l[is.na(V.l)] = 0
      
      ij.end.time = Sys.time()
      ij.time = ij.end.time-ij.start.time
      print(c('ij time:', ij.time))
      print("ij loop finished")
      
      inv.start.time = Sys.time()
      
      k.grad = 0
      
      b.grad.k.plus1 = b.l
      
      while (k.grad<=200) {
        b.grad.k = b.grad.k.plus1
        
        sum.exp.mat = matrix(0, nrow = n*L, ncol = 1)
        for (i in 1:n) {
          tmp1 = Big.Gamma[((i-1)*L+1):(i*L), 1]
          tmp2 = b.grad.k[((i-1)*L+1):(i*L)]
          sum.exp.mat[((i-1)*L+1):(i*L), 1] = Big.Gamma[((i-1)*L+1):(i*L), 1]/as.numeric(1+exp(-t(tmp1)%*%tmp2))
        }
        
        p.k = -sum.y.E.ga + sum.exp.mat + 2*(phi*R + lam*V.l)%*%b.grad.k
        
        if(k<10){
          alpha.k = 0.001
        }else{
          alpha.k = 0.00001
        }
        
        b.grad.k.plus1 = b.grad.k - alpha.k*p.k
        k.grad = k.grad+1
        
        if(!is.na(vecnorm(p.k)) && vecnorm(p.k)<=0.01){
          break
        }
      }
      
      b.l.plus1 = b.grad.k.plus1
      
      
      
      inv.end.time = Sys.time()
      inv.time = inv.end.time - inv.start.time
      print(c("Inverse Time:", inv.time))
      print("Inverse Calculated")
      
      it = it+1
      
    }
    
    
    
    b.hat = b.l.plus1
    
    b.res = matrix(0, nrow = n, ncol = L)
    
    for (i in 1:n) {
      b.res[i,] = b.hat[(L*(i-1)+1):(L*(i-1)+L),]
    }
    
    b.res = as.data.frame(b.res)
    
    fuse.intercept = c()
    if(intercept==TRUE){
      X.train.26 = t(x.fd$coefs)
      for (j in 1:n) {
        linshi = as.matrix(b.res)
        logit_tmp = 1/(1+exp(-X.train.26[j,]%*%Z%*%linshi[j,]))
        tmp = as.numeric(y[j,] - logit_tmp)
        fuse.intercept[j] = mean(tmp)
      }
    }
    
    
    res.final = list(b.hat = b.hat, b.res = b.res, grp = grp, intercept = fuse.intercept)
    
    
  }
  
  
  
  return(res.final)
  
}
