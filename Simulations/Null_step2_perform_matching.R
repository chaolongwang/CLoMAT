#-----------------------------------------------------------------------------------------------------
# File: Null_step2_perform_matching.R
#-----------------------------------------------------------------------------------------------------

set.seed(8)

#--Parameters for input-------------------------------------------------------------------------------
args = commandArgs(T);

Match = args[1]; #This will be "fullmatch" or integer. Integer means performing matching with 1:n.
loop = as.numeric(args[2]); #This is an index for loop (1:1000).
casetype = args[3] #options: "sharp" or "smooth" 
casenum = 1000  #options: sample size 
caliper = "pr95"
DIR=args[4] #the working directory 
#-------------------------------------------------------------------------------------------------------


#--Read in files----------------------------------------------------------------------------------------
IndexCase = scan(file= paste(DIR,"/null_cases_",casetype,"/case",loop,".txt",sep=""))
IndexCont = scan(file= paste(DIR,"/null_controls_",casetype,"/control",loop,".txt",sep="")) 

IndexCase = sample(IndexCase, as.numeric(casenum));
Cont_all = IndexCont;

if (loop==1) {dir.create(paste(DIR,"/null_phenotypes_data",sep=""))}
#------------------------------------------------------------------------------------------------------


#--Packages need to be load---------------------------------------------------------------------------
library(optmatch);
source(paste(DIR,"/pairmatch.R",sep="")) 
#---------------------------------------------------------------------------------------------------------


#--run the following code to do control down-sampling--------------------------------------------------
	#----randomly pick 20 controls from each grid (if <20 controls, retain all)--------------------------
	all_controls = IndexCont+1;
	sampled_controls = c();
	control_num=c();
  	for(j in 1:400) {
    tmp = seq((j-1)*50+1,(j-1)*50+50);
   	tmp_controls = tmp[tmp %in% all_controls];
		tmp_controls = sample(tmp_controls, min(length(tmp_controls),20));
    sampled_controls = c(sampled_controls, tmp_controls);
		control_num = c(control_num, length(tmp_controls));
		}
#-------------------------------------------------------------------------------------------------------


#--generate phenotype and combined to top 10 PCs--------------------------------------------------------
Pool = read.table(paste(DIR,"/data/M10.20x20.20000.1M.ProPC.coord",sep=""),header=T);
RefPC = read.table(paste(DIR,"/data/M10.20x20.20000.1M.RefPC.coord",sep=""),header=T);
IndexCont=sampled_controls-1;

#already excluded the references here
IndexRef = match(RefPC[,1],Pool[,1]);
Pool = Pool[-IndexRef,];

Index = match(IndexCase,Pool[,1]);
CasePC = Pool[Index,];
row.names(CasePC) = CasePC[,1];
CasePC = cbind(CasePC,Phenotype = rep(1,dim(CasePC)[1]),grid=IndexCase%/%50);

Index2 = match(IndexCont,Pool[,1]);
PopControl = Pool[Index2,];
row.names(PopControl) = PopControl[,1];
PopControl = cbind(PopControl,Phenotype=rep(0,dim(PopControl)[1]),grid=IndexCont%/%50);

PCA = rbind(CasePC,PopControl);
#------------------------------------------------------------------------------------------------------


#--perform matching------------------------------------------------------------------------------------
if(caliper=="pr95"){
    cali=6.05;
}

if(Match=="fullmatch"){
	Eucli = match_on(Phenotype ~ PC1+PC2, data = PCA,method="euclidean",caliper = cali);
	EucliMatch = fullmatch(Eucli,data=PCA,min.controls=(1/10),max.controls=10,remove.unmatchables=TRUE);
}else{
	Match = as.numeric(Match);
	Eucli = match_on(Phenotype ~ PC1+PC2, data = PCA,method="euclidean");
	matched = pairmatch.heuristic(Eucli, data=PCA, controls=Match, width=cali);
}

if(Match=="fullmatch"){
	collate = data.frame("1:0"=0,"1:1"=0,"1:2"=0,"1:3"=0,"1:4"=0,"1:5+"=0,
			   "2:1"=0,"3:1"=0,"4:1"=0,"5+:1"=0,
			   "0:1"=0,"effective"=0,check.names=F);
	test = summary(EucliMatch)$matched.set.structures;
	for (j in 1: length(test)) {
		collate[1,names(test)[j]] = test[names(test)[j]];	
		}
	collate$effective = summary(EucliMatch)$effective.sample.size;
	collate = cbind(collate,t(quantile(unlist(matched.distances(EucliMatch, Eucli)))));
}else{
	collate = cbind(effective=matched$ne, t(quantile(matched$d)));
	EucliMatch = matched$grp;
}
#------------------------------------------------------------------------------------------------------


#--Combined phenotypes and matched data----------------------------------------------------------------
MatchResultEucli2PC = cbind(PCA,EucliMatch);
IndexNA = which(is.na(EucliMatch)==TRUE);
if(length(IndexNA)>0){
	MatchResultEucli2PC = MatchResultEucli2PC[-IndexNA,];
}
MatchResultEucli2PC = MatchResultEucli2PC[order(MatchResultEucli2PC[,"EucliMatch"],MatchResultEucli2PC[,"Phenotype"]),];
	
Factor = unique(as.character(MatchResultEucli2PC[,"EucliMatch"]));
Factor = cbind(Factor,seq(1,length(Factor)));
Index = match(MatchResultEucli2PC[,"EucliMatch"],Factor[,1]);
MatchResultEucli2PC[,"EucliMatch"] = Factor[Index,2];
write.table(MatchResultEucli2PC,file=paste(DIR,"/null_phenotypes_data/MatchResultEucli2PC_",casetype,"_Match",Match,"_cali_",caliper,"_loop_",loop,".txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F);
#------------------------------------------------------------------------------------------------------


#--randomly sampled 1000 controls without matching (illustration of population stratification)---------
#use 1000 cases and 1000 controls for running unmatched tests
#IndexCont_all = sample(all_controls, 1000);
#PopControl_all = Pool[IndexCont_all,];
#row.names(PopControl_all) = PopControl_all[,1];
#PopControl_all = cbind(PopControl_all,Phenotype=rep(0,dim(PopControl_all)[1]),grid=rep(0,dim(PopControl_all)[1]));
#PCA_unmatched = rbind(CasePC,PopControl_all);
#------------------------------------------------------------------------------------------------------



