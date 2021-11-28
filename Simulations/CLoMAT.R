library(survival)
library(survey)

CLRSKAT.null = function(D,X,S){
  if(is.null(X)){
    g10 = suppressWarnings(clogit(D~1+strata(S)))
  } else{
    g10 = clogit(D~X+strata(S))
  }
  S = as.numeric(levels(S))[S];
  res = g10$res; llk0=tail(g10$loglik, n=1); P0 = matrix(0,max(S),length(D))
  pi0 = D-res; P0[cbind(S,1:length(D))] = pi0; 
  return(list(res=res,pi0=pi0,P0=P0,llk0=llk0,Y=D,X=X,S=S,coef=g10$coef) )
}

CLRBurden = function(obj,G, W.beta) {
  N = dim(G)[2]; maf = colMeans(G)/2
  W = maf^(W.beta[1]-1)*(1-maf)^(W.beta[2]-1);  W = W/sum(W)*N
  X = obj$X; Z = t(t(G)*W); S = obj$S
  ZZ<-apply(Z,1,sum)
  if(is.null(X)){
    Y = obj$Y; llk0 = obj$llk0; rcf = c(0)
    lj = clogit(Y~ZZ+strata(S), init=rcf)
    Qblr = 2*(logLik(lj)-llk0)
    Qbsc = lj$score
  }
  else {
    Y = obj$Y; llk0 = obj$llk0; rcf = c(obj$coef,0)
    lj = clogit(Y~X+ZZ+strata(S), init=rcf)
    Qblr = 2*(logLik(lj)-llk0)
    Qbsc = lj$score
  }
  return(list(plr = pchisq(Qblr,1,lower=F), psc =pchisq(Qbsc,1,lower=F)))
}

CLRMiST = function(D,X,G,S,W.beta) {

  N = dim(G)[2]; maf = colMeans(G)/2
  W = maf^(W.beta[1]-1)*(1-maf)^(W.beta[2]-1);  W = W/sum(W)*N
  Z = t(t(G)*W); ZZ =apply(Z,1,sum)

  #null model
  null.out = CLRSKAT.null(D,X,S);

  #Burden test
  p_burden = CLRBurden(null.out,G, W.beta)$plr;

  #null model with genetic burden
  null.burden = CLRSKAT.null(D,cbind(X,ZZ),S);

  #SKAT-type test
  p_skat = CLRSKAT(null.burden, G, W.beta)$plr;

  #Combine p-values via Fisher's procedure
  return( pchisq(-2*log(p_burden*p_skat), df=4, lower.tail=FALSE)) 
}


CLRSKAT = function(obj,G, W.beta){
  N = dim(G)[2]; maf = colMeans(G)/2
  W = maf^(W.beta[1]-1)*(1-maf)^(W.beta[2]-1);  W = W/sum(W)*N
  pi0 = obj$pi0; P = t(obj$P0); X = obj$X; Z = t(t(G)*W); S = obj$S

  ## Sigma matrix
  if(is.null(X)){
     #ZP = t(Z)%*%P 
      ZP = crossprod(Z,P)  
      VV = crossprod((Z*pi0),Z)-tcrossprod(ZP)
	    #VV = t(Z*pi0)%*%Z-ZP%*%t(ZP)

  }
  else {
     #ZP = t(Z)%*%P 
      ZP = crossprod(Z,P)
	    #XP = t(X)%*%P  
      XP = crossprod(X,P)
      #ZWZ = t(Z*pi0)%*%Z-ZP%*%t(ZP)
      ZWZ = crossprod((Z*pi0),Z)-tcrossprod(ZP)
      #ZWX = t(Z*pi0)%*%X-ZP%*%t(XP)
      ZWX = crossprod((Z*pi0),X)-tcrossprod(ZP,XP)
      #XWX = t(X*pi0)%*%X-XP%*%t(XP)
      XWX = crossprod((X*pi0),X)-tcrossprod(XP)
      VV = ZWZ-ZWX%*%solve(XWX)%*%t(ZWX)
  }

  GL1 = 1/diag(VV)
  Zs = (colSums(obj$res*Z))^2*GL1
  Zs[which(is.na(Zs))] = 0
  Gt = rep(0,N);  idw = 1:N
  ids = which(maf<(10/dim(G)[1]))
  if(length(ids)>0){
    Gt[ids] = Zs[ids]
    idw = (1:N)[-ids]
  }

  print(ids)
  ## LRT
  if(length(idw)>0){
    if(is.null(X)){
        Y = obj$Y; llk0 = obj$llk0; p = 0; rcf = c(0)
        Gt[idw] = suppressWarnings( apply(G[,idw,drop=FALSE], 2, function(Gj){
          lj = clogit(Y~Gj+strata(S), init=rcf)
          #sign(lj$coef[p+1])*sqrt(2*logLik(lj)-2*llk0)
          (2*logLik(lj)-2*llk0)
        }) )
    }
    else {
        Y = obj$Y; llk0 = obj$llk0; p = dim(X)[2]; rcf = c(obj$coef,0)
        Gt[idw] = suppressWarnings( apply(G[,idw,drop=FALSE], 2, function(Gj){
          lj = clogit(Y~X+Gj+strata(S), init=rcf)
          #sign(lj$coef[p+1])*sqrt(2*logLik(lj)-2*llk0)
          (2*logLik(lj)-2*llk0)
        }) )
    }
    if(any(is.na(Gt))){
      ia = which(is.na(Gt))
      Gt[ia] = Zs[ia]
    }
  }
  qgt = Gt/GL1; qzs = Zs/GL1
  qgt[which(is.na(qgt))] = 0; qzs[which(is.na(qzs))] = 0
  Qlr = sum(qgt); Qsc = sum(qzs)
  lam = svd(VV, nu=0,nv=0)$d

  return(list(plr = KAT.pval(Qlr, lam), psc = KAT.pval(Qsc, lam)))
}


### saddlepoint approx: modified from Lumley survey package.
saddle = function(x,lambda){
  d = max(lambda)
  lambda = lambda/d
  x = x/d
  k0 = function(zeta) -sum(log(1-2*zeta*lambda))/2
  kprime0 = function(zeta) sapply(zeta, function(zz) sum(lambda/(1-2*zz*lambda)))
  kpprime0 = function(zeta) 2*sum(lambda^2/(1-2*zeta*lambda)^2)
  n = length(lambda)
  if (any(lambda < 0)) {
    lmin = max(1/(2 * lambda[lambda < 0])) * 0.99999
  } else if (x>sum(lambda)){
    lmin = -0.01
  } else {
    lmin = -length(lambda)/(2*x)
  }
  lmax = min(1/(2*lambda[lambda>0]))*0.99999
  hatzeta = uniroot(function(zeta) kprime0(zeta) - x, lower = lmin, upper = lmax, tol = 1e-08)$root
  w = sign(hatzeta)*sqrt(2*(hatzeta*x-k0(hatzeta)))
  v = hatzeta*sqrt(kpprime0(hatzeta))
  if(abs(hatzeta)<1e-4){
    return(NA)
  } else{
    return( pnorm(w+log(v/w)/w, lower.tail=FALSE) )
  }
}
Sadd.pval = function(Q.all,lambda){
  sad = rep(1,length(Q.all))
  if(sum(Q.all>0)>0){
    sad[Q.all>0] = sapply(Q.all[Q.all>0],saddle,lambda=lambda)
  }
  id = which(is.na(sad))
  if(length(id)>0){
    sad[id] = Liu.pval(Q.all[id], lambda)
  }
  return(sad)
}
### modified Liu method from Lee SKAT-O paper
Liu.pval = function(Q, lambda){
  c1 = rep(0,4); for(i in 1:4){ c1[i] = sum(lambda^i) }
  muQ = c1[1];  sigmaQ = sqrt(2 *c1[2])
  s1 = c1[3]/c1[2]^(3/2);  s2 = c1[4]/c1[2]^2
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2));  d = s1 *a^3 - a^2;  l = a^2 - 2*d
  } else {
    l = 1/s2;  a = sqrt(l);  d = 0
  }
  muX = l+d;  sigmaX = sqrt(2)*a
  
  Q.Norm = (Q - muQ)/sigmaQ
  Q.Norm1 = Q.Norm*sigmaX + muX
  pchisq(Q.Norm1, df = l,ncp=d, lower.tail=FALSE)
}
Liu.qval.mod = function(pval, lambda){
  c1 = rep(0,4)
  c1[1] = sum(lambda); c1[2] = sum(lambda^2)
  c1[3] = sum(lambda^3); c1[4] =sum(lambda^4)
  muQ = c1[1]; sigmaQ = sqrt(2 *c1[2])
  s1 = c1[3]/c1[2]^(3/2); s2 = c1[4]/c1[2]^2
  beta1= sqrt(8)*s1; beta2 = 12*s2; type1 = 0
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2)); d = s1 *a^3 - a^2; l = a^2 - 2*d
  } else {
    type1 = 1; l = 1/s2; a = sqrt(l); d = 0
  }
  muX = l+d; sigmaX = sqrt(2) *a
  df = l
  q.org = qchisq(pval,df=df,lower.tail=FALSE)
  (q.org - df)/sqrt(2*df)*sigmaQ + muQ
}
KAT.pval = function(Q.all, lambda, acc=1e-9,lim=1e6){
  pval = rep(0, length(Q.all))
  i1 = which(is.finite(Q.all))
  for(i in i1){
    tmp = davies(Q.all[i],lambda,acc=acc,lim=lim); pval[i] = tmp$Qq
    if((tmp$ifault>0)|(pval[i]<=0)|(pval[i]>=1)) pval[i] = Sadd.pval(Q.all[i],lambda)
  }
  return(pval)
}