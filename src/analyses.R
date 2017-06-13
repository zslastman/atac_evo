##Does what it says
require(magrittr)
require(ggplot2)
require(dplyr)


fileList = c('4dSites.daf.txt',
             'DNaseSites.daf.txt',
             'neutral_intronicSites.daf.txt',
             'nonsynPolymorphicSites.daf.txt',
             'proxy_neutral_noncodingSites.daf.txt',
             'secondPositionSites.daf.txt',  
              'validSites.daf.txt')

nameList = c('4dSites',
             'DNaseSites',
             'neutral_intronicSites',
             'nonsynPolymorphicSites',
             'proxy_neutral_noncodingSites',
             'secondPositionSites',  
             'validSites')


myData = read.delim(fileList[1])
myData$siteType = nameList[1]
for(i in 2:length(fileList)){
  newData = read.delim(fileList[i])
  newData$siteType = nameList[i]
  myData = rbind(myData,newData)
}

contrast_w_validSites<-function(myData, compName){
myMatrix = matrix(c(myData[myData$siteType==compName,]$rareCount, 
                    myData[myData$siteType==compName,]$commonCount,
                    myData[myData$siteType=="validSites",]$rareCount, 
                    myData[myData$siteType=="validSites",]$commonCount), nrow=2,
                  dimnames = list())
fisherResults<-fisher.test(myMatrix)
return(c(fisherResults$estimate, min= fisherResults$conf.int[1], max = fisherResults$conf.int[2]))
}

for(myName in nameList[1:6]){
  print(myName)
  print(contrast_w_validSites(myData, myName))
}

ggdat = lapply(nameList[1:6],function(myName){
  contrast_w_validSites(myData, myName)
})

pdf('DAF_test.pdf')

ggdat%>%simplify2array%>%t%>%as.data.frame%>%
  rename(odds_ratio=`odds ratio`)%>%
  mutate(set = nameList[1:6])%>%
  {
    qplot(data=.,geom='pointrange',x=set,ymax=max,ylab = 'odds_ratio',ymin=min,y=odds_ratio,ylim=c(0,2.5))+
      coord_flip()+
      geom_hline(yintercept=1)+
      theme_bw()+
      ggtitle('DAF test for alternative sites')
  }

dev.off()