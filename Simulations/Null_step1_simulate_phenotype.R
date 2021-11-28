#-------------------------------------------------------------------------------
# File: Null_step1_simulate_phenotype.R
#-------------------------------------------------------------------------------

set.seed(8)

#--Parameters for input---------------------------------------------------------
args = commandArgs(T);

casetype=args[1] #"smooth" or "sharp"
DIR=args[2] #the working directory 
loopindex = as.numeric(args[3]) #This is an index for loop (Set at 1000 in the present simulation design).
#---------------------------------------------------------------------------------


#--reference population-----------------------------------------------------------
RefPC = read.table(paste(DIR,"/data/M10.20x20.20000.1M.RefPC.coord",sep=""),header=T);
#---------------------------------------------------------------------------------

expit <- function(x) {
  exp(x)/(1+exp(x))
}

logit <- function(x) {
  log(x / (1 - x))
}

probs <- numeric(20000);
grids = numeric(400);

#--smooth spatial distribution of disease risk-----------------------------------
if (casetype=="smooth") {
  for(j in 1:400) {
    tmp = seq((j-1)*50+1,(j-1)*50+50)
    x = (j-1)%%20; y = (j-1)%/%20
    probs[tmp] = expit(-0.5-0.03*((y-14)^2+(x-14)^2))
    grids[j] = expit(-0.5-0.03*((y-14)^2+(x-14)^2))
  }
}
#--------------------------------------------------------------------------------
  
#--sharp spatial distribution of disease risk------------------------------------
if (casetype=="sharp") {
  for(j in 1:7){
    tmp1 = (9+j)*20+11;
    tmp2 = tmp1+6;
    for(i in tmp1:tmp2){
      tmp = seq((i-1)*50+1,(i-1)*50+50)
      grids[i] = expit(-3+2.6)
      probs[tmp] = expit(-3+2.6)
    }
  }
  grids[which(grids==0)]=expit(-3)
  probs[which(probs==0)]=expit(-3)
}
  
dir.create(paste(DIR,"/null_cases_",casetype,sep=""))
dir.create(paste(DIR,"/null_controls_",casetype,sep=""))
for(k in 1:loopindex) {
  y = rbinom(20000, 1, prob=probs)
  #--exclude reference panel ids-------------------------------------------------
  caseids = setdiff(which(y==1)-1,RefPC[,1])		
  contids = setdiff(which(y==0)-1,RefPC[,1])
  write(caseids,file=paste(DIR,"/null_cases_",casetype,"/case",k,".txt",sep="")) #save indivID
  write(contids,file=paste(DIR,"/null_controls_",casetype,"/control",k,".txt",sep=""))
}
#--------------------------------------------------------------------------------
