if ( !require(survival,lib.loc=paste(bindir,"/lib/",sep="") ) ) {
  stop("Cannot find survival package");
}
library(survival)

##################################################################
## GENE-LEVEL BURDEN TEST UNDER CONDITIONAL LOGISTIC REGRESSION
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

## group.clrburden() : CLR-Burden implementation
## KEY FEATURES : 0/1 collapsing variables ~ rare variants
## TRAITS : BINARY
## RETURNS : PVALUE, LRT_STAT, SCR_PVAL, SCR_STAT
## MISSING VALUE : Imputed as major allele homozygotes

CLRSKAT.null = function(D,X,S){
  if(is.null(X)){
    g10 = suppressWarnings(clogit(D~1+strata(S) ,method="breslow"))
  } else{
    g10 = clogit(D~X+strata(S),method="breslow")
  }
  res = g10$res; llk0=tail(g10$loglik, n=1); P0 = matrix(0,max(S),length(D))
  pi0 = D-res; P0[cbind(S,1:length(D))] = pi0; 
  return(list(res=res,pi0=pi0,P0=P0,llk0=llk0,Y=D,X=X,S=S, coef=g10$coef) )
}

group.clrburden <- function() {
  m <- nrow(genos)
  W = MAF[vids]^(betas[1]-1)*(1-MAF[vids])^(betas[2]-1);  W = W/sum(W)*m
  grp = cov[,ncol(cov)]; x = cov[,-c(1,ncol(cov))]; if(dim(as.matrix(x))[2]==0) {x=NULL} 
  obj = CLRSKAT.null(pheno,x,grp);
  flip <- which(AC>NS)
  genos[flip,] <- 2-genos[flip,]
  genos[is.na(genos)] <- 0  
  if ( qr(genos)$rank < 1 ) {
    m <- 0
  }
  ZZ<-rowSums(t(genos*W)); S = obj$S; Y = obj$Y; llk0 = obj$llk0
  cname <- c("LRT_STAT","SCR_PVAL","SCR_STAT")
  if ( var(ZZ) > 0 & m>0) {
    if(is.null(x)){
      rcf = c(0)
      lj = clogit(Y~ZZ+strata(S), init=rcf, method="breslow")
      Qblr = 2*(logLik(lj)-llk0)
      Qbsc = lj$score
    }
    else {
      rcf = c(obj$coef,0)
      lj = clogit(Y~x+ZZ+strata(S), init=rcf, method="breslow")
      Qblr = 2*(logLik(lj)-llk0)
      Qbsc = lj$score
    }
    return(list(p=pchisq(Qblr,1,lower=F),add=c(Qblr,pchisq(Qbsc,1,lower=F),Qbsc),cname=cname))
  }
  return(list(p=NA,add=rep(NA,3),cname=cname))  
}
