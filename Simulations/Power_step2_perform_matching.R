#------------------------------------------------------------------------------------------------------
# File: Power_step2_perform_matching.R
#------------------------------------------------------------------------------------------------------

set.seed(8)

#--Parameters for input--------------------------------------------------------------------------------
args = commandArgs(T);
CausalPercent = as.numeric(args[1]) #0.2; 0.5
PositivePercent = as.numeric(args[2]) #0.8; 1.0
Constant = as.numeric(args[3]) #0.877; 0.555
Match = args[4]; #This will be "fullmatch" or integer. Integer means performing matching with 1:n.
loop = as.numeric(args[5]); #This is an index for loop (1:1000).
casetype = args[6] #options: "sharp" or "smooth" 
casenum = 1000  #sample size of cases (1000)
caliper = "pr95"
DIR=args[7] #the working directory 
#-------------------------------------------------------------------------------------------------------


#--Packages need to be load---------------------------------------------------------------------------
library(optmatch);
source(paste(DIR,"/pairmatch.R",sep="")) 
#---------------------------------------------------------------------------------------------------------

if (loop==1) {dir.create(paste(DIR,"/power_phenotypes_data_",casetype,sep=""))}

#--Read in files----------------------------------------------------------------------------------------
Casepool = scan(file= paste(DIR,"/power_cases_",casetype,"/case_",CausalPercent,"_",PositivePercent,"_",loop,".txt",sep=""))
all_controls = scan(file= paste(DIR,"/power_controls_",casetype,"/control_",CausalPercent,"_",PositivePercent,"_",loop,".txt",sep=""))
IndexCase = sample(Casepool, casenum);

Pool = read.table(paste(DIR,"/data/M10.20x20.20000.1M.ProPC.coord",sep=""),header=T);
RefPC = read.table(paste(DIR,"/data/M10.20x20.20000.1M.RefPC.coord",sep=""),header=T);
IndexRef = match(RefPC[,1],Pool[,1]);
#-------------------------------------------------------------------------------------------------------


#--Run the following code to do control down-sampling---------------------------------------------------
	all_controls = all_controls+1
	#randomly pick 20 controls from each grid (if <20 controls, retain all)
#	if (FALSE) {
#	sampled_controls = c();
#	control_num=c();
#  	for(j in 1:400) {
#      tmp = seq((j-1)*50+1,(j-1)*50+50);
#    	tmp_controls = tmp[tmp %in% all_controls];
#      tmp_controls = sample(tmp_controls, min(length(tmp_controls),20));
#		sampled_controls = c(sampled_controls, tmp_controls);
#		}
#	}
     sampled_controls = sample(all_controls, length(all_controls)*(1/2));
#-------------------------------------------------------------------------------------------------------


#--Generate phenotype and combined to top 10 PCs--------------------------------------------------------
Control = Pool[setdiff(sampled_controls,IndexRef),];
Control = cbind(Control,grid=setdiff(sampled_controls,IndexRef)%/%50)

Index = match(IndexCase,Pool[,1]);
CasePC = Pool[Index,];
row.names(CasePC) = CasePC[,1];
CasePC = cbind(CasePC,Phenotype = rep(1,dim(CasePC)[1]),grid=IndexCase%/%50);
PopControl = cbind(Control,Phenotype=rep(0,dim(Control)[1]));
PCA = rbind(CasePC,PopControl);
#-------------------------------------------------------------------------------------------------------


#--Perform matching-------------------------------------------------------------------------------------
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

N.eff=numeric();
if(Match=="fullmatch"){
	collate = data.frame("1:0"=0,"1:1"=0,"1:2"=0,"1:3"=0,"1:4"=0,"1:5+"=0,
			   "2:1"=0,"3:1"=0,"4:1"=0,"5+:1"=0,
			   "0:1"=0,"effective"=0,check.names=F);
	test = summary(EucliMatch)$matched.set.structures;
	for (j in 1: length(test)) {
		collate[1,names(test)[j]] = test[names(test)[j]];	
		}
	collate$effective = summary(EucliMatch)$effective.sample.size;
	N.eff= summary(EucliMatch)$effective.sample.size;
	collate = cbind(collate,t(quantile(unlist(matched.distances(EucliMatch, Eucli)))));
}else{
	collate = cbind(effective=matched$ne, t(quantile(matched$d)));
	N.eff = matched$ne;
	EucliMatch = matched$grp;
}
#------------------------------------------------------------------------------------------------------


#--Combined phenotypes and matched data----------------------------------------------------------------
MatchResultEucli2PC = cbind(PCA,EucliMatch);
IndexNA = which(is.na(EucliMatch)==TRUE);
if(length(IndexNA)>0){
	MatchResultEucli2PC = MatchResultEucli2PC[-IndexNA,];
}
	
Factor = unique(as.character(MatchResultEucli2PC[,"EucliMatch"]));
Factor = cbind(Factor,seq(1,length(Factor)));
Index = match(MatchResultEucli2PC[,"EucliMatch"],Factor[,1]);
MatchResultEucli2PC[,"EucliMatch"] = Factor[Index,2];
write.table(MatchResultEucli2PC,file=paste(DIR,"/power_phenotypes_data_",casetype,"/MatchResultEucli2PCForMatch",Match,"cali_",caliper,"_",CausalPercent,"_",PositivePercent,"_loop",loop,".txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F);
#------------------------------------------------------------------------------------------------------

