#-----------------------------------------------------------------------------------------------------
# File: Power_step1_simulate_phenotype.R
#-----------------------------------------------------------------------------------------------------


#--Packages need to be load---------------------------------------------------------------------------
library(SKAT);
#-----------------------------------------------------------------------------------------------------

set.seed(8) 

#--Parameters for input-------------------------------------------------------------------------------
args = commandArgs(T);

CausalPercent = as.numeric(args[1]) #0.2; 0.5
PositivePercent = as.numeric(args[2]) #0.8; 1.0
Constant = as.numeric(args[3]) #0.877; 0.555
casenum = 1000 #sample size of cases (1000)
casetype=args[4] #"smooth" or "sharp"
DIR=args[5] #the working directory 
loopindex = as.numeric(args[6]); #This is an index for loop (Set at 1000 in the present simulation design).
#-----------------------------------------------------------------------------------------------------

dir.create(paste(DIR,"/power_cases_",casetype,sep=""))
dir.create(paste(DIR,"/power_controls_",casetype,sep=""))

#--reference population-------------------------------------------------------------------------------
Pool = read.table(paste(DIR,"/data/M10.20x20.20000.1M.ProPC.coord",sep=""),header=T);
RefPC = read.table(paste(DIR,"/data/M10.20x20.20000.1M.RefPC.coord",sep=""),header=T);
IndexRef = match(RefPC[,1],Pool[,1]);
#-----------------------------------------------------------------------------------------------------

set.seed(8)
ID = sample(1:5000,casenum,replace=FALSE)

maf = read.table(paste(DIR,"/data/AllMarkerAllIndividual.frq",sep=""),header=T);
maf[,"SNP"] = as.character(as.matrix(maf[,"SNP"]));

expit <- function(x) {
  exp(x)/(1+exp(x))
}

logit <- function(x) {
  log(x / (1 - x))
}

Gamma = log((casenum+200)/5000/(1-(casenum+200)/5000));

ssd_obj = Open_SSD(paste(DIR,"/data/BiRareMarkerAllIndividualGeneSet1_5000",".SSD",sep=""), paste(DIR,"/data/BiRareMarkerAllIndividualGeneSet1_5000",".Info",sep=""));

for(loop in 1:loopindex){
GeneSetID = ID[loop];

GenoMatrix = Get_Genotypes_SSD(ssd_obj, Set_Index =GeneSetID, is_ID = TRUE);
IndexRare = match(colnames(GenoMatrix),maf[,"SNP"]);

#--set the beta----------------------------------------------------------------------------------------
	if(length(IndexRare)==0){
	}else{
		if(length(IndexRare)==1){
			RareInfo = maf[IndexRare,];
			RareMatrix = cbind(rep(1,dim(GenoMatrix)[1]),GenoMatrix[,RareInfo[,"SNP"]]);
			Beta = Constant*abs(log10(RareInfo[1,"MAF"]))/2;
			Beta = c(Gamma,Beta);
		}else{
			RareInfo = maf[IndexRare,];
			RareMatrix = cbind(rep(1,dim(GenoMatrix)[1]),GenoMatrix[,RareInfo[,"SNP"]]);
			Numcausal = ceiling(dim(RareInfo)[1]*CausalPercent);
			Beta = rep(0,dim(RareInfo)[1]);
			if(Numcausal==1){
				causalIndex = sample(seq(1,dim(RareInfo)[1]),Numcausal);
				Beta[causalIndex] = Constant*abs(log10(RareInfo[causalIndex,"MAF"]))/2;
				Beta = c(Gamma,Beta);
			}else{
				causalIndex = sample(seq(1,dim(RareInfo)[1]),Numcausal);
				Positive = ceiling(Numcausal*PositivePercent);
				if(Positive==Numcausal){
					Beta[causalIndex] = Constant*abs(log10(RareInfo[causalIndex,"MAF"]))/2;
					Beta = c(Gamma,Beta);
				}else{
					PositiveIndex = sample(causalIndex,Positive);
					Beta[PositiveIndex] = Constant*abs(log10(RareInfo[PositiveIndex,"MAF"]))/2;
					Beta[setdiff(causalIndex,PositiveIndex)] = -Constant*abs(log10(RareInfo[setdiff(causalIndex,PositiveIndex),"MAF"]))/2;
					Beta = c(Gamma,Beta);
				}
			}
		}
  }
	RareMatrix = RareMatrix[,-1]; Beta = Beta[-1]
#-------------------------------------------------------------------------------------------------------

  #--smooth spatial distribution of disease risk--------------------------------------------------------
	if (casetype=="smooth") {
	  probs <- numeric(20000);
	  for(j in 1:400) {
	    tmp = seq((j-1)*50+1,(j-1)*50+50)
	    x = (j-1)%%20; y = (j-1)%/%20
	    probs[tmp] = expit(-0.5-0.03*((y-14)^2+(x-14)^2)+RareMatrix[tmp,]%*%Beta)
	  }
	  Pheno = rbinom(20000, 1, prob=probs)
	}
  #-----------------------------------------------------------------------------------------------------

  #--sharp spatial distribution of disease risk---------------------------------------------------------
  if (casetype=="sharp") {
    probs <- numeric(20000);
    for(j in 1:7){
      tmp1 = (9+j)*20+11;
      tmp2 = tmp1+6;
      for(k in tmp1:tmp2){
        tmp = seq((k-1)*50+1,(k-1)*50+50)
        probs[tmp] = expit(-3+2.6+RareMatrix[tmp,]%*%Beta)
      }
    }
    lp = which(probs==0)
    probs[lp]=expit(-3+RareMatrix[lp,]%*%Beta)
    Pheno = rbinom(20000, 1, prob=probs)
  }
  #-----------------------------------------------------------------------------------------------------
	caseids = setdiff(which(Pheno==1),IndexRef)-1		
	contids = setdiff(which(Pheno==0),IndexRef)-1
	
	write(caseids,file=paste(DIR,"/power_cases_",casetype,"/case",loop,".txt",sep=""))
	write(contids,file=paste(DIR,"/power_controls_",casetype,"/control",loop,".txt",sep=""))
}

Close_SSD()




