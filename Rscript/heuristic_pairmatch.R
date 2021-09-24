effective.sample.size <- function(strata, pheno){
	S <- as.factor(strata);
	ne <- 0;
	for(stra in levels(S)){
		idx <- which(S%in%stra);
		n1 <- length(which(pheno[idx]==1)); # number of cases
		n0 <- length(which(pheno[idx]==0)); # number of controls
		if(n1!=0 & n0!=0) ne <- ne + 2/(1/n1 + 1/n0);
	}
	return(ne);
}

pairmatch.heuristic <- function(D, controls=1, data=NULL, width=Inf){
	D[which(D>=width)] <- NA;
	ids <- order(D, decreasing=FALSE, na.last=NA);  ## most time-consuming step
	cat("Finish ordering distances and start matching ...", "\n");
	m <- controls; 
	n <- dim(D);
	M <- matrix(NA, n[1], m);
	cnt <- rep(0,n[1]);
	rows.rm <- c();
	cols.rm <- c();
	case.id <- rownames(D);
	ctrl.id <- colnames(D);
	rownames(M) <- case.id;
	names(cnt) <- case.id;	
	for(i in 1:length(ids)){
		if(i%%50000 == 0){
			cat(paste("Matching progress: ", round(i/length(ids)*100), "%", sep=""), "\n");
		}
		idx2 <- ceiling(ids[i]/n[1]);
		idx1 <- ids[i]-(idx2-1)*n[1];
		if(idx2 %in% cols.rm || idx1 %in% rows.rm){
			next;
		}else{
			cnt[idx1] <- cnt[idx1]+1;
			M[idx1, cnt[idx1]] <- ctrl.id[idx2];
			cols.rm <- c(cols.rm, idx2); 
			if(cnt[idx1]==m){
				rows.rm <- c(rows.rm, idx1); 
			}
			if(length(rows.rm)==n[1] || length(cols.rm)==n[2]){
				break;
			}
		}
	}
	 partialmatch <- which(cnt<m & cnt>0);
	 if(length(partialmatch)>1){
		 rest.cases <- case.id[partialmatch];
		 rest.ctrls <- M[partialmatch,];
		 rest.ctrls <- rest.ctrls[!is.na(rest.ctrls)];
		 rest.D <- D[rest.cases, rest.ctrls];
		 rest.case.id <- rownames(rest.D);
		 rest.ctrl.id <- colnames(rest.D);		
		 n.potential.ctrls <- apply(!is.na(rest.D),1,sum);
		 rest.case.list <- order(n.potential.ctrls);
		 for(i in rest.case.list){
			 if(n.potential.ctrls[i]>=m){
				 nearest <- order(rest.D[i,])[1:m];
				 if(!is.na(sum(rest.D[i,nearest]))){
					 M[rest.case.id[i],] <- rest.ctrl.id[nearest];
					 rest.D[,nearest] <- NA;
					 cnt[rest.case.id[i]] <- m;
				 }
			 }
		 }
	 }
	cat(paste("Done! ", sum(cnt==m), " 1-to-", m, " matched clusters!", sep=""), "\n");	
	if(!is.null(data)){
		orig.id <- rownames(data);
	}else{
		orig.id <- c(case.id, ctrl.id);	
	}
	strata <- rep(NA,length(orig.id));
	d <- c();
	for(i in which(cnt==m)){
		strata[which(orig.id %in% c(case.id[i],M[i,]))] <- i;
		d <- c(d, D[case.id[i],M[i,]]);
	}
	ne <- sum(cnt==m)*2/(1+1/m);
	cat(paste("Effective Sample Size: ", ne, sep=""), "\n");
	cat(paste("(equivalent number of matched pairs).", sep=""), "\n");
	cat(paste("Summary of within-group distances: ", sep=""), "\n");
	print(summary(d));
	return(list(matched.group=as.factor(strata), effective.sample.size=ne, Distance=unlist(d)));
}


