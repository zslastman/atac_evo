####################################
#       A R G U M E N T S          #
####################################

args          <- commandArgs(trailingOnly=TRUE);
dataFile <- args[1];
plotFile <- args[2];
freqThres<- args[3];
setName  <- args[4];

####################################
#         S E T T I N G S          #
####################################

doSFSplot=FALSE;
# colors
barColorsMain=c("steelblue","gray80");
barColorsHigh=c("steelblue4","gray50");
neutralBarColors=("darkgreen");

genomeWideDivPoly=c(0.00510871850675972090,0.00396358717199217967);
# ylim = c(0,6);

####################################
#     P A R S E    D A T A         #
####################################

table<-read.table(dataFile);
numSamples<-ncol(table)-2;
sfsTable<-matrix(0.0,2,numSamples);
cdfTable<-matrix(0.0,2,numSamples);
# table of neutral and observed SFS
sfsTable[,2:numSamples]<-as.numeric(as.matrix(table[2:3,2:numSamples+2]));
# CDF for neutral and observed SFS
# entry corresponds to fraction of poly with counts <= colIndex
for(sample in 2:numSamples) {
   cdfTable[,sample]<-cdfTable[,sample-1]+as.numeric(as.matrix(table[2:3,sample+2]));
}

#####################################################################################
#     P R O C E S S    A N D    I N T E R P R E T   P O L Y D I V    S T A T S      #
#####################################################################################

# 2X2 poly and div rate matrix (per kbp)
polyDivStats<-matrix(0,2,2);
polyDivStats[]<-as.numeric(as.matrix(table[2:3,3:2])) / rowSums(matrix(as.numeric(as.matrix(table[2:3,1:3])),2,3))
polyDivStats<-polyDivStats*1000;

# high frequency poly rate (neutral and observed)
freqIndex=floor(numSamples*as.numeric(freqThres)*0.01)+1;
highPolyFracs<-polyDivStats[,1]*(1-cdfTable[,freqIndex]);

#1-rho
fracNeutral<-min(1.0,highPolyFracs[1]/highPolyFracs[2]);

# neutral fractions of poly and div
# = fracNeutral * neutral expectation of poly and div counts
neutralPolyDiv<-polyDivStats[2,]*fracNeutral;

Pw=polyDivStats[1,1]-neutralPolyDiv[1];
Dp=polyDivStats[1,2]-neutralPolyDiv[2];

####################################
#             P L O T              #
####################################

pdf(plotFile,height=3,width=5,title="");

# plotting layer 1 - total poly and div
layer1 <- polyDivStats;
# plotting layer 2 - high freq poly
layer2 <- matrix(c(highPolyFracs,0,0),2,2);

ylim=c(0,max(polyDivStats)*1.5);

res=barplot(height=layer1,beside=TRUE,space=c(0.5,1.2),border=NA,col=barColorsMain,ylim=ylim);
# abline(h=1:6,lty=3,col="gray30");
# barplot(height=layer1,beside=TRUE,space=c(0.5,1.5),border=NA,col=barColorsMain,ylim=ylim);
barplot(height=layer2,beside=TRUE,space=c(0.5,1.2),border=NA,col=barColorsHigh,add=TRUE);
segments(x0=res[1,]-0.6,x1=res[2,]+0.6,y0=1000*genomeWideDivPoly,y1=1000*genomeWideDivPoly,col="darkred",lty=3,lwd=2);
numSites=sum(as.numeric(as.matrix(table[2,1:3])))
if(numSites>1000000) {
   numSites=paste(round(numSites/100000)/10,"Mb");
} else {
   numSites=paste(round(numSites/1000),"kb");
}
title(main=paste("Poly & Div Stats for",setName,""),sub=paste("[",numSites,"informative sites ]"),col.main="darkblue",col.sub="gray30",cex.sub=0.8);
# axis labels
mtext(side=1,line=1.5,cex=1.2,at=colMeans(res),text=c("polymorphisms","divergences"));
mtext(side=1,line=-0.7,padj=1,cex=0.7,at=res,col=barColorsHigh,text=c("observed","expected\nby neutral"));
mtext(side=2,line=2.2,cex=1,text=c("rate per kbp"));
# H,L indicators
text(x=res[,1],y=highPolyFracs-ylim[2]*0.05,font=2,cex=0.5,col=c("lightblue",barColorsMain[2]),"H");
if(freqIndex>1) {
   text(x=res[,1],y=highPolyFracs+ylim[2]*0.05,font=2,cex=0.5,col=c("royalblue4",barColorsHigh[2]),"L");
}
text(x=res[2,1]+0.4,pos=4,y=highPolyFracs[2],font=2,cex=0.7,paste(sep='',"f=",freqThres,"%"),col="gray40");
segments(x0=res[2,1]+0.49,x1=res[2,1]+0.56,y0=highPolyFracs[2]-ylim[2]*0.005,y1=highPolyFracs[2]-ylim[2]*0.005,col="gray40");
# horizontal bars for neutral portion of poly and div
segments(x0=res[1,]-0.48, x1=res[1,]+0.48, y0=neutralPolyDiv,y1=neutralPolyDiv, col=neutralBarColors, lwd=2);

# dashed 'shrinkage' indicators for reduction in poly div due to selection
x0s<-res[1,c(1,1,2)]+0.50
x1s<-res[2,c(1,1,2)]-0.50
y0s<-c(highPolyFracs[1],neutralPolyDiv)*1.01;
y1s<-c(highPolyFracs[2],polyDivStats[2,])*0.99;
segments(x0=x0s,x1=x1s,y0=y0s,y1=y1s,lty=2,col="darkgreen",lwd=1.0);

# text indicators for simple estimates
text(x=res[1,1],y=polyDivStats[1,1],pos=3,paste(sep='',"Pw=",round(Pw*10)/10),cex=0.9,col="darkgreen");
text(x=res[1,2],y=polyDivStats[1,2],pos=3,paste(sep='',"Dp=",round(Dp*10)/10),cex=0.9,col="darkgreen");
segments(x0=res[1,],x1=res[1,],y0=polyDivStats[1,]*0.9,y1=polyDivStats[1,]*1.2,col=neutralBarColors,cex=1.5);
print(c(mean(res[,1]),max(polyDivStats[,1])*1000));
text(x=mean(res[,1]),y=max(polyDivStats[,1])*1.15,pos=3,parse(text=paste(sep='',"1-rho~~\"=\"~~",round(fracNeutral*100)/100)),cex=0.9,col="darkgreen");
segments(x0=mean(res[,1]),x1=mean(res[,1]),y0=mean(c(y0s[1],y1s[1])),y1=max(polyDivStats[,1])*1.15,lty=1,col="darkgreen");
dev.off();


####################################
#    SFS PLOTTING  -  SUPPRESSED   #
####################################


if(doSFSplot) {
   pdf(plotFile,height=3,width=6,title="");
   # layout(matrix(c(1,1,1,2),1,4,byrow=TRUE));

   res=barplot(height=sfsTable,beside=TRUE,ylim=c(0,1),space=c(0.0,0.3),border=NA,col=barColors);
   positions=colMeans(res)
   labelPos=c(1,1:(numSamples/10)*10,numSamples-1);
   mtext(side=1,at=positions[labelPos+1],text=labelPos,line=0,cex=0.7,col="gray30");
   mtext(side=1,text="derived allele count",line=1.2);
   mtext(side=2,line=2.2,cex=1,text=c("frequency among poly sites"));
   lines(x=positions[2:ncol(res)],y=cdfTable[1,2:ncol(res)],col=barColors[1]);
   lines(x=positions[2:ncol(res)],y=cdfTable[2,2:ncol(res)],col=barColors[2]);
   mtext(side=3,line=1,at=positions[length(positions)*0.66],col="darkblue",text=paste("Polymorphism and divergence statistics for",setName));
   legend(legend=c("observed","expected by neutral"),fill=barColors,x='right',inset=0.05,border=NA,bty="n",cex=1.2);
   dev.off();
}
