if ( !require(survival,lib.loc=paste(bindir,"/lib/",sep="") ) ) {
  stop("Cannot find survival package");
}

if ( !require(CompQuadForm,lib.loc=paste(bindir,"/lib/",sep="") ) ) {
  stop("Cannot find CompQuadForm package");
}

library(survival)
library(CompQuadForm)

##################################################################
## GENE-LEVEL SKAT TEST UNDER CONDITIONAL LOGISTIC REGRESSION
## INPUT VARIABLES: 
##   n        : total # of individuals
##   genos    : genotype matrix for each gene
##   NS       : number of called samples for each marker
##   AC       : allele count for each marker
##   MAC      : minor allele count for each marker
##   MAF      : minor allele frequency
##   vids     : indices from 1:n after AF/AC threshold
## EXPECTED OUTPUT : list(p, addcols, addnames) for each genos row
##   p        : p-value
##   add      : additional column values
##   cname    : additional column names
##################################################################

## group.clrskat() : CLR-SKAT implementation
## KEY FEATURES : 0/1 collapsing variables ~ rare variants
## TRAITS : BINARY
## RETURNS : PVALUE, LRT_STAT, SCR_PVAL, SCR_STAT
## MISSING VALUE : Imputed as major allele homozygotes

CLRSKAT.null = function(D,X,S){
  if(is.null(X)){
    g10 = suppressWarnings(clogit(D~1+strata(S),method="breslow"))
  } else{
    g10 = clogit(D~X+strata(S),method="breslow") 
  }
  res = g10$res; llk0=tail(g10$loglik, n=1); P0 = matrix(0,max(S),length(D))
  pi0 = D-res; P0[cbind(S,1:length(D))] = pi0; 
  if (!is.null(X)) {
    XP = crossprod(X,t(P0))
    XWX_inv = solve(crossprod((X*pi0),X)-tcrossprod(XP));
  }
  else {XWX_inv=NULL;}
  return(list(res=res,pi0=pi0,P0=P0,llk0=llk0,Y=D,X=X,S=S, coef=g10$coef, XWX_inv=XWX_inv) )
}

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
  f1 =kprime0(lmin) - x;
  f2 =kprime0(lmax) - x;
  if ((sign(f1)+sign(f2))==0) {
    hatzeta = uniroot(function(zeta) kprime0(zeta) - x, lower = lmin, upper = lmax, tol = 1e-08)$root
  } else{
    return(NA)
  }
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

group.clrskat <- function() {
  m <- nrow(genos)
  W.beta = betas;
  W = MAF[vids]^(W.beta[1]-1)*(1-MAF[vids])^(W.beta[2]-1);  W = W/sum(W)*m
  cname <- c("LRT_STAT","SCR_PVAL","SCR_STAT")
  grp = cov[,ncol(cov)]; x = cov[,-c(1,ncol(cov))]; if(dim(as.matrix(x))[2]==0) {x=NULL} 
  print(is.null(x))
  obj = CLRSKAT.null(pheno,x,grp);
  flip <- which(AC>NS)
  genos[flip,] <- 2-genos[flip,]
  genos[is.na(genos)] <- 0  
  if ( qr(genos)$rank < 1 ) {
    m <- 0
  }
  if (m>0) {
    G = t(genos);
    pi0 = obj$pi0; P = t(obj$P0); X = obj$X; Z = t(genos*W); S = obj$S
    if(is.null(X)){
      ZP = crossprod(Z,P) 
      VV = crossprod((Z*pi0),Z)-tcrossprod(ZP)
    }
    else {
      ZP = crossprod(Z,P)
      XP = crossprod(X,P)
      ZWZ = crossprod((Z*pi0),Z)-tcrossprod(ZP)
      ZWX = crossprod((Z*pi0),X)-tcrossprod(ZP,XP)
      XWX_inv=obj$XWX_inv;
      VV = ZWZ-ZWX%*%XWX_inv%*%t(ZWX)
    }
    GL1 = 1/diag(VV)
    Zs = (colSums(obj$res*Z))^2*GL1
    Zs[which(is.na(Zs))] = 0
    Gt = rep(0,m);  idw = 1:m 
    ids = which(MAC[vids]<10)
    if(length(ids)>0){
      Gt[ids] = Zs[ids]
      idw = (1:m)[-ids]
    }
    if(length(idw)>0){
      if(is.null(X)){
          Y = obj$Y; llk0 = obj$llk0; p = 0; rcf = c(0)
          Gt[idw] = suppressWarnings( apply(G[,idw,drop=FALSE], 2, function(Gj){
            lj = tryCatch(logLik(clogit(pheno~Gj+strata(S), init=rcf, method="breslow")), error=function(err) NA);
            (2*lj-2*llk0)
          }) )
      }
      else {
          Y = obj$Y; llk0 = obj$llk0; p = dim(X)[2]; rcf = c(obj$coef,0)
          Gt[idw] = suppressWarnings( apply(G[,idw,drop=FALSE], 2, function(Gj){
            lj = tryCatch(logLik(clogit(pheno~X+Gj+strata(S), init=rcf, method="breslow")), error=function(err) NA);
            (2*lj-2*llk0)
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
    return(list(p=KAT.pval(Qlr, lam),add=c(Qlr,KAT.pval(Qsc, lam),Qsc),cname=cname))
  }
  return(list(p=NA,add=rep(NA,3),cname=cname))  
}
