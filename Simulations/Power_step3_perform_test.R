#----------------------------------------------------------------------------------------------------
# File: Power_step3_perform_test.R
#----------------------------------------------------------------------------------------------------

set.seed(8)

#--Parameters for input------------------------------------------------------------------------------
args = commandArgs(T);

Match = args[1]; #This will be "fullmatch" or integer. Integer means performing matching with 1:n.
casenum = 1000 #sample size of cases (1000)
casetype = args[2] #options: "sharp" or "smooth" 
caliper = "pr95" 
DIR=args[3] #the working directory 
loopindex = as.numeric(args[4]); #This is an index for loop (Set at 1000 in the present simulation design).
#----------------------------------------------------------------------------------------------------


#-----Packages need to be load-----------------------------------------------------------------------
library(CompQuadForm);
library(SKAT);
source(paste(DIR,"/CLoMAT.R",sep="")) 
source(paste(DIR,"/SKATL.R",sep="")) 
#----------------------------------------------------------------------------------------------------


#--Read in data files--------------------------------------------------------------------------------
set.seed(8)
ID = sample(1:10000,casenum,replace=FALSE)

ssd_obj = Open_SSD(paste(DIR,"/data/BiRareMarkerAllIndividualGeneSet1_10000",".SSD",sep=""), paste(DIR,"/data/BiRareMarkerAllIndividualGeneSet1_10000",".Info",sep=""));
#----------------------------------------------------------------------------------------------------


#--Iterations to perform tests-----------------------------------------------------------------------
result = data.frame()
for(loop in 1:loopindex){
	GeneSetID = ID[loop];
	GenoMatrix = Get_Genotypes_SSD(ssd_obj, Set_Index =GeneSetID, is_ID = TRUE);
	wuweights <- function(maf) ifelse(maf>0, dbeta(maf,1,25), 0)

	#--Extract genotypes of cases and controls---------------------------------------------------------
	MatchResultEucli2PC = read.table(paste(DIR,"/power_phenotypes_data/MatchResultEucli2PCForMatch",Match,"cali_",caliper,"_loop_",loop,".txt",sep=""),sep="\t",header=T);

	MatchResultEucli2PC[,"EucliMatch"] = as.factor(MatchResultEucli2PC[,"EucliMatch"]);
	MatchID = match(MatchResultEucli2PC[,"popID"]+1,seq(1,dim(GenoMatrix)[1]));

	Z = GenoMatrix[MatchID,];
	maf = apply(Z,2,function(x){sum(x,na.rm=T)/((length(x)-sum(is.na(x)))*2)} );
	ZeroVarIndex = which(maf==0);
	if(length(ZeroVarIndex)>0){
		Z = Z[,-ZeroVarIndex];
		Z = as.matrix(Z);
	}
	if(dim(Z)[2]==1){
		return(data.frame(CLRSKAT_lrt = -1,
			CLRSKAT_scr = -1,
			CLR_Burden_lrt= -1,
			CLR_Burden_scr= -1,
			SKAT= -1,
			SKAT_lrt= -1,
			Burden = -1,
			CLRMist = -1,
			Mist = -1));	
	}else{
		if(dim(Z)[2]==2){
		Z = Z[,-1];
		maf = sum(Z,na.rm=T)/((length(Z)-sum(is.na(Z)))*2);
		Weight = wuweights(maf)
		Z1 = Z*Weight;
		}else{
		Z = Z[,-1];
		maf = apply(Z,2,function(x){sum(x,na.rm=T)/((length(x)-sum(is.na(x)))*2)} );
		Weight = wuweights(maf)
		Z1 = Z%*%Weight;
		}

	DATA = cbind(MatchResultEucli2PC,Z);
	Z = as.matrix(Z);
	#--------------------------------------------------------------------------------------------------

	#--Fit models under the null-----------------------------------------------------------------------
	Null_skat_2pc_adj = SKAT_Null_Model(Phenotype~PC1+PC2, data=DATA, out_type="D", n.Resampling=0, type.Resampling="bootstrap", Adjustment=FALSE);
	Null_skat_2pc_adj_burden = SKAT_Null_Model(Phenotype~PC1+PC2+Z1, data=DATA, out_type="D", n.Resampling=0, type.Resampling="bootstrap", Adjustment=FALSE);
	SKATL.null <- KAT.null(MatchResultEucli2PC$Phenotype, cbind(MatchResultEucli2PC$PC1, MatchResultEucli2PC$PC2));
	CKAT.null <- CLRSKAT.null(MatchResultEucli2PC$Phenotype,X=NULL,S=MatchResultEucli2PC$EucliMatch);
	#--------------------------------------------------------------------------------------------------

	#----Perform gene-based tests----------------------------------------------------------------------
	Pvalue_clrskat<- unlist(CLRSKAT(obj=CKAT.null,G=Z , W.beta=c(1,25)));
	skat_2pc = SKAT(Z,Null_skat_2pc_adj,weights.beta=c(1,25))$p.value;
	Pvalue_clrburden<- unlist(CLRBurden(obj=CKAT.null,G=Z , W.beta=c(1,25)));
	Burden = SKAT(Z,Null_skat_2pc_adj,r.corr=1,weights.beta=c(1,25))$p.value;
	CLRMist = CLRMiST(D=MatchResultEucli2PC$Phenotype,X=NULL,G=Z,S=MatchResultEucli2PC$EucliMatch,W.beta=c(1,25));
	Mist=pchisq(-2*log(Burden*SKAT(Z,Null_skat_2pc_adj_burden,weights.beta=c(1,25))$p.value), df=4, lower.tail=FALSE);
	SKAT_lrt = SKATL(obj=SKATL.null ,G=Z , W.beta=c(1,25))

	tmp_result = data.frame(CLRSKAT_lrt = Pvalue_clrskat[1],
			CLRSKAT_scr = Pvalue_clrskat[2],
			SKAT=skat_2pc,
			CLR_Burden_lrt=Pvalue_clrburden[1],
			CLR_Burden_scr=Pvalue_clrburden[2],
			Burden = Burden,
			CLRMist = CLRMist,
			Mist = Mist,
			SKAT_lrt=SKAT_lrt);
	}
	result = rbind(result,tmp_result)
}

dir.create(paste(DIR,"/power_output_downsample_",casetype,sep=""))
write.table(result,file=paste(DIR,"/power_output_downsample_",casetype,"/CLR_boot_","caliper_",caliper,"NullModel_match",Match,"c1000_exact_ptem.txt",sep=""),sep="\t",quote=F,col.names=F,row.names=F)
#----------------------------------------------------------------------------------------------------

Close_SSD()
