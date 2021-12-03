args = commandArgs(T);

caliper = "pr95"; 
Match = args[1];  #This will be "fullmatch" or integer. Integer means performing matching with 1:n.
casetype = args[2]; #options: "sharp" or "smooth" 
type = args[3]; #options: "downsample" or "original"
loopindex = as.numeric(args[4]); #This is an index for loop (Set at 1000 in the present simulation design).
DIR = args[5] #the working directory 


SKAT <- c()
CLR_SKAT <- c()
Burden <- c();
CLR_Burden <- c()
Mist <- c() 
CLR_Mist <- c() 

#--read in files
for (loop in 1:loopindex) {
  SKAT_File <- paste(DIR,"/null_output_downsample_",casetype,"/CLR_boot_loop",loop,"caliper_",caliper,"NullModel_match",Match,"c1000_exact_ptem.txt",sep="");

  if (file.exists(SKAT_File)) {
    skat_table <- read.table(SKAT_File)
    CLR_Burden = c(CLR_Burden, skat_table$V4);
    CLR_SKAT = c(CLR_SKAT, skat_table$V1);
    SKAT = c(SKAT, skat_table$V3);
    Burden = c(Burden, skat_table$V6);
    CLR_Mist = c(CLR_Mist, skat_table$V7);
    Mist = c(Mist, skat_table$V8);
  }
}

total = cbind(SKAT,CLR_SKAT,Burden,CLR_Burden,Mist,CLR_Mist);

dir.create(paste(DIR,"/plots",sep=""))

#--draw qqplot
png(file=paste(DIR,"/plots/Caliper",caliper,"Match",Match,"casetype",casetype,"type",type,"QQplot.png",sep=""),width=7,height=7,units="in",res=600);
par(mgp=c(3.8,1,0), mai=c(1,1.2,0.1,0.1))
match.loose <- function(a, b) sapply(a, function(x) b[which.min(abs(x-b))])
cols <- c("pink2","red","skyblue","blue","yellow2", "orange");

linpc10p0 = sort(SKAT)
lskat <- qchisq(median(SKAT, na.rm = TRUE),1,lower=F)/qchisq(0.5,1)
n0 = length(linpc10p0)

lobs <- -log10(linpc10p0)
select0 <- 1:n0
funnelx <- seq(log10(1+n0), 0, by=-0.01)
select0 <- unique(match.loose(10^(-funnelx)*(1+n0), select0))
lobs <- lobs[select0]
lexp <- -log10(select0/(1+n0))
funnell<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.975,x,n0+1-x)))
funnelu<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.025,x,n0+1-x)))
plot(lexp, lobs, xlab=expression(Expected ~ ~ -log[10](p)), ylab=expression(Observed ~ ~ -log[10](p)), main="", pch=16, cex=0.5, col=cols[1], xlim=c(0, max(lexp)), ylim=c(0,max(funnelu)), type="n", cex.main=2.5, cex.lab=2, cex.axis=1.8)
polygon(c(funnelx,rev(funnelx)),c(funnell,rev(funnelu)),col="#eeeeee",border="#eeeeee")
lines(c(0,max(lexp)),c(0,max(lexp)))
points(lexp, lobs, pch=16, cex=0.6, col=cols[1])


lobs <- -log10(sort(CLR_SKAT))
lskatl <- qchisq(median(CLR_SKAT, na.rm = TRUE),1,lower=F)/qchisq(0.5,1)
lobs <- lobs[select0]
lexp <- -log10(select0/(1+n0))
points(lexp, lobs, pch=16, cex=0.6, col=cols[2])


lobs <- -log10(sort(Burden))
lclrsc <- qchisq(median(Burden, na.rm = TRUE),1,lower=F)/qchisq(0.5,1)
lobs <- lobs[select0]
lexp <- -log10(select0/(1+n0))
points(lexp, lobs, pch=16, cex=0.6, col=cols[3])

lobs <- -log10(sort(CLR_Burden))
lclrll <- qchisq(median(CLR_Burden, na.rm = TRUE),1,lower=F)/qchisq(0.5,1)
lobs <- lobs[select0]
lexp <- -log10(select0/(1+n0))
points(lexp, lobs, pch=16, cex=0.6, col=cols[4])

lobs <- -log10(sort(Mist))
lclrllm <- qchisq(median(Mist, na.rm = TRUE),1,lower=F)/qchisq(0.5,1)
lobs <- lobs[select0]
lexp <- -log10(select0/(1+n0))
points(lexp, lobs, pch=16, cex=0.6, col=cols[5])

lobs <- -log10(sort(CLR_Mist))
lclrllmist <- qchisq(median(CLR_Mist, na.rm = TRUE),1,lower=F)/qchisq(0.5,1)
lobs <- lobs[select0]
lexp <- -log10(select0/(1+n0))
points(lexp, lobs, pch=16, cex=0.6, col=cols[6])

legend("bottomright", legend=c(as.expression(substitute(paste("SKAT, ", phantom("00000"), lambda[GC]==ll, sep=""), list(ll=format(round(lskat,2), nsmall=2)))),
    as.expression(substitute(paste("CLR-SKAT, ", phantom("0"), lambda[GC]==ll, sep=""), list(ll=format(round(lskatl,2), nsmall=2)))),
    as.expression(substitute(paste("Burden, ",phantom("0000"),lambda[GC]==ll, sep=""), list(ll=format(round(lclrsc,2), nsmall=2)))),
    as.expression(substitute(paste("CLR-Burden, ",lambda[GC]==ll, sep=""), list(ll=format(round(lclrll,2), nsmall=2)))),
    as.expression(substitute(paste("MiST, ", phantom("000000"),lambda[GC]==ll, sep=""), list(ll=format(round(lclrllm,2), nsmall=2)))),
    as.expression(substitute(paste("CLR-MiST, ",phantom("00"),lambda[GC]==ll, sep=""), list(ll=format(round(lclrllmist,2), nsmall=2))))),
    lty=1, lwd=3, col=cols[1:6], bty="n", cex=1.2)

dev.off()
