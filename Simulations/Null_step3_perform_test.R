#----------------------------------------------------------------------------------------------------
# File: Null_step3_perform_test.R
#----------------------------------------------------------------------------------------------------

set.seed(8)

#--Parameters for input------------------------------------------------------------------------------
args = commandArgs(T);

Match = args[1]; #This will be "fullmatch" or integer. Integer means performing matching with 1:n.
loop = as.numeric(args[2]); #This is an index for loop (1:1000).
casetype = args[3] #options: "sharp" or "smooth" 
caliper = "pr95"
DIR=args[4] #the working directory 
#----------------------------------------------------------------------------------------------------

if (loop==1) {dir.create(paste(DIR,"/null_output_downsample_",casetype,sep=""))}

#-----Packages need to be load-------------------------------------------------------------------
library(CompQuadForm);
library(SKAT);
source(paste(DIR,"/CLoMAT.R",sep="")) 
source(paste(DIR,"/SKATL.R",sep="")) 
#----------------------------------------------------------------------------------------------------


#--Extract genotypes of cases and controls-----------------------------------------------------------
MatchResultEucli2PC = read.table(paste(DIR,"/null_phenotypes_data_",casetype,"/MatchResultEucli2PC_",casetype,"_Match",Match,"_cali_",caliper,"_loop_",loop,".txt",sep=""),sep="\t",header=T);

FamInfo = Read_Plink_FAM(paste(DIR,"/data/BiRareMarkerAllIndividual",".fam",sep=""), Is.binary=TRUE, flag1=0);
FamInfo[,c(1,2)] = as.character(as.matrix(FamInfo[,c(1,2)]));

MatchResultEucli2PC[,"popID"] = paste("I",MatchResultEucli2PC[,"popID"],sep="");
MatchResultEucli2PC[,"EucliMatch"] = as.factor(MatchResultEucli2PC[,"EucliMatch"]);
Index = match(FamInfo[,1],MatchResultEucli2PC[,1],nomatch=-1);
IndexSample = which(Index>0);
MatchResultEucli2PC = MatchResultEucli2PC[Index[which(Index>0)],];
#----------------------------------------------------------------------------------------------------


#--Fit models under the null-------------------------------------------------------------------------
Null_skat_2pc_adj = SKAT_Null_Model(Phenotype~PC1+PC2, data=MatchResultEucli2PC, out_type="D", n.Resampling=0, type.Resampling="bootstrap", Adjustment=FALSE);
SKATL.null <- KAT.null(MatchResultEucli2PC$Phenotype, cbind(MatchResultEucli2PC$PC1, MatchResultEucli2PC$PC2));
CKAT.null <- CLRSKAT.null(MatchResultEucli2PC$Phenotype,X=NULL,S=MatchResultEucli2PC$EucliMatch);
#----------------------------------------------------------------------------------------------------


#--Perform gene-based tests--------------------------------------------------------------------------
wuweights <- function(maf) ifelse(maf>0, dbeta(maf,1,25), 0)

ssd_obj = Open_SSD(paste(DIR,"/data/BiRareMarkerAllIndividualGeneSet1_5000",".SSD",sep=""), paste(DIR,"/data/BiRareMarkerAllIndividualGeneSet1_5000",".Info",sep=""));

N.genes <- 5000;
GeneIndex<- seq(1,N.genes,1)

Burden = c();
pskat_2pc_adj=0;
Pvalue_mist= matrix(0,N.genes,1);
Pvalue_clrskat = matrix(0,N.genes,2);
Pvalue_clrburden = matrix(0,N.genes,2);
Pvalue_clrmist= matrix(0,N.genes,1);
pskatl = 0;

counter= 0 
for(i in GeneIndex){
	counter=counter+1
	GenoMatrix = Get_Genotypes_SSD(ssd_obj, Set_Index =i, is_ID = FALSE);
	GenoMatrix = GenoMatrix[IndexSample,];
	maf = apply(GenoMatrix,2,function(x){sum(x,na.rm=T)/((length(x)-sum(is.na(x)))*2)} );
	ZeroVarIndex = which(maf==0);
	if(length(ZeroVarIndex)>0){
		GenoMatrix = GenoMatrix[,-ZeroVarIndex];
	}

	tmp_p = SKAT(Z=GenoMatrix,Null_skat_2pc_adj,weights.beta=c(1,25))$p.value;
	pskat_2pc_adj = c(pskat_2pc_adj,tmp_p);

	tmp_burden_p = SKAT(Z=GenoMatrix,Null_skat_2pc_adj,r.corr=1,weights.beta=c(1,25))$p.value;
	Burden = c(Burden, tmp_burden_p);

	if(dim(GenoMatrix)[2]>1){
	maf = apply(GenoMatrix,2,function(x){sum(x,na.rm=T)/((length(x)-sum(is.na(x)))*2)} );
	if (any(maf>0.5)) {print("ERROR: MAF GREATER THAN HALF!")}
		Weight = wuweights(maf);
		Z = GenoMatrix%*%Weight;
	}else{
		maf = sum(GenoMatrix,na.rm=T)/((length(GenoMatrix)-sum(is.na(GenoMatrix)))*2);
		Weight = wuweights(maf);
		Z = GenoMatrix*Weight;
	}

	DATA = cbind(MatchResultEucli2PC,Z);

	Null_skat_2pc_adj_burden = SKAT_Null_Model(Phenotype~PC1+PC2+Z, data=DATA, out_type="D", n.Resampling=0, type.Resampling="bootstrap", Adjustment=FALSE);
	Pvalue_mist[counter,]  <- pchisq(-2*log(tmp_burden_p*SKAT(GenoMatrix,Null_skat_2pc_adj_burden,weights.beta=c(1,25))$p.value), df=4, lower.tail=FALSE);

	Pvalue_clrskat[counter,]<- unlist(CLRSKAT(obj=CKAT.null,G=GenoMatrix , W.beta=c(1,25)));
	Pvalue_clrburden[counter,] <- unlist(CLRBurden(obj=CKAT.null,G=GenoMatrix , W.beta=c(1,25)));
	Pvalue_clrmist[counter,]   <- CLRMiST(D=MatchResultEucli2PC$Phenotype,X=NULL,G=GenoMatrix,S=MatchResultEucli2PC$EucliMatch,W.beta=c(1,25));

	#SKATL
	pskatl<- c(pskatl, SKATL(obj=SKATL.null ,G=GenoMatrix , W.beta=c(1,25)));
}

pskat_2pc_adj<-pskat_2pc_adj[-1];
pskatl<-pskatl[-1];

write.table(cbind(Pvalue_clrskat,pskat_2pc_adj,Pvalue_clrburden, Burden, Pvalue_clrmist, Pvalue_mist),file=paste(DIR,"/null_output_downsample_",casetype,"/CLR_boot_loop",loop,"caliper_",caliper,"NullModel_match",Match,"c1000_exact_ptem.txt",sep=""),sep="\t",quote=F,col.names=F,row.names=F)

Close_SSD()
#----------------------------------------------------------------------------------------------------










