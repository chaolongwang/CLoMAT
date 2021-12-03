args = commandArgs(T);

caliper = "pr95"; 
CausalPercent = as.numeric(args[1]) #0.2; 0.5
PositivePercent = as.numeric(args[2]) #0.8; 1.0
Match = args[3];  #This will be "fullmatch" or integer. Integer means performing matching with 1:n.
casetype = args[4]; #options: "sharp" or "smooth" 
type = args[5]; #options: "downsample" or "original"
loopindex = as.numeric(args[6]); #This is an index for loop (Set at 1000 in the present simulation design).
DIR = args[7] #the working directory 


SKAT <- c()
CLR_SKAT <- c()
Burden <- c();
CLR_Burden <- c()
Mist <- c() 
CLR_Mist <- c() 

#--read in files
  SKAT_File <- paste(DIR,"/power_output_downsample_",casetype,"/CLR_boot_","caliper_",caliper,"PowerModel_match",Match,"_",CausalPercent,"_",PositivePercent,"c1000_exact_ptem.txt",sep="");
  
  if (file.exists(SKAT_File)) {
    skat_table <- read.table(SKAT_File)
    CLR_Burden = c(CLR_Burden, skat_table$V4);
    CLR_SKAT = c(CLR_SKAT, skat_table$V1);
    SKAT = c(SKAT, skat_table$V3);
    Burden = c(Burden, skat_table$V6);
    CLR_Mist = c(CLR_Mist, skat_table$V7);
    Mist = c(Mist, skat_table$V8);
  }

total1 = cbind(Burden,CLR_Burden,SKAT,CLR_SKAT,Mist,CLR_Mist);

dir.create(paste(DIR,"/plots",sep=""))

#--draw barplot
png(file=paste(DIR,"/plots/Caliper",caliper,"Match",Match,"casetype",casetype,type,"_",CausalPercent,"_",PositivePercent,"barplot.png",sep=""),width=7,height=7,units="in",res=600);

results1<-as.matrix(apply(total1,2,function(x) sum(x<2.5e-6,na.rm = T)/loopindex))
full.res = t(results1)
colours <- c("skyblue","blue","pink","red","yellow","gold")
barplot(full.res,ylim=c(0,1), ylab="Empirical power",beside=T,col=colours,xlab=ifelse(Match %in% c(1,3),paste("1:",Match,sep=""),"Full"),xaxt="n",space=c(0,0))
legend("topleft", c("Burden","CLR-Burden","SKAT","CLR-SKAT","MiST","CLR-MiST"), cex=0.7,bty="n", fill=colours)

dev.off()
