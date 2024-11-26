func_fuse_lasso = function(x.fd, y, beta.basis, b0, lam, phi, intercept = FALSE, clustering=FALSE, x.penalty=FALSE, KK=0, label=NULL, max.iteration=30, tol=0.95){
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
    
    
    S0 = t(x.fd$coefs)%*%Z
    S = matrix(0, nrow = n, ncol = L*KK)
    cnt = 0
    for (k in 1:KK) {
      idx.this.grp = which(label==k)
      for (i in idx.this.grp) {
        tmp = S0[i,]
        cnt = cnt+1
        S[ cnt , (L*k-L+1):(L*k)] = tmp
      }
    }
    
    
    tS.times.S = t(S)%*%S
    tS.times.y = t(S)%*%y
    
    
    
    
    R1 = eval.penalty(beta.basis,int2Lfd(2))
    list2 = NULL
    for (i in 1:KK){
      list2[[i]] = R1
    }
    library(Matrix)
    R = as.matrix(bdiag(list2))
    
    
    it = 1
    b.l.plus1 = b0
    
    grp = 1:KK
    
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
          
          # if(it>15){
          #   if(grp[i] != grp[j]){
          #     if(vecnorm(b.i.l-b.j.l)<=tol){
          #       b.j.l = b.i.l
          #       idx = min(i,j)
          #       grp[i] = idx
          #       grp[j] = idx
          #     }
          #   }
          # }
          
          
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
      
      
      temp.mat = tS.times.S + KK*phi*R + KK*lam*V.l
      #b.l.plus1 = solve(temp.mat)%*%t(S)%*%y
      b.l.plus1 = chol2inv(chol(temp.mat))%*%tS.times.y
      # b.l.plus1 = ginv(temp.mat)%*%t(S)%*%y
      
      
      
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
    if(intercept==TRUE){
      linshi = as.matrix(b.res)
      X.train.26 = t(x.fd$coefs)
      neiji = inprod(x.fd$basis, beta.basis)
      for (j in 1:KK) {
        tmp1 = which(label==j)
        tmp2 = y[tmp1] - X.train.26[tmp1,]%*%neiji%*%linshi[j,]
        tmp2 = as.matrix(tmp2)
        fuse.intercept[j] = mean(tmp2)
      }
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
    
    
    S0 = t(x.fd$coefs)%*%Z
    S = matrix(0, nrow = n, ncol = n*L)
    for (i in 1:n) {
      S[i , (L*i-L+1):(L*i)] = S0[i,]
    }
    
    tS.times.y = t(S)%*%y
    tS.times.S = t(S)%*%S
    
    
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
          
          # if(it>15){
          #   if(grp[i] != grp[j]){
          #     if(vecnorm(b.i.l-b.j.l)<=tol){
          #       b.j.l = b.i.l
          #       idx = min(i,j)
          #       grp[i] = idx
          #       grp[j] = idx
          #     }
          #   }
          # }
          
          
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
      
      temp.mat = tS.times.S + n*phi*R + n*lam*V.l
      # b.l.plus1 = solve(temp.mat)%*%tS.times.y
      b.l.plus1 = chol2inv(chol(temp.mat))%*%t(S)%*%y
      # b.l.plus1 = ginv(temp.mat)%*%t(S)%*%y
      
      
      
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
        tmp = as.numeric(y[j,] - X.train.26[j,]%*%Z%*%linshi[j,])
        fuse.intercept[j] = mean(tmp)
      }
    }
    
    
    res.final = list(b.hat = b.hat, b.res = b.res, grp = grp, intercept = fuse.intercept)
    
    
  }
  
  
  
  return(res.final)
  
}
