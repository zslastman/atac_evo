library(dplyr)
library(data.table)
library(ggplot2)
library(magrittr)

barwidth=0.8
dodgewidth=0.5

#functiont to stick plots together
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

insightplots <- function(dataname,datafile,outputfile,
  rhomin=0,rhomax=1,ewmin=0,ewmax=NULL,  eamin = 0,eamax = 2
  ){

  #read in data
  ggdata=fread(datafile)
  #filter out organism
  ggdata %<>% filter(!dataID=='organism')
  #arrange by rho
  ggdata  %<>% arrange(rho)

  ggdata$dataID  %<>% factor(.,levels=unique(.))
  # ggdata %<>% group_by(dataID) %>% group_slice(1:3)

  p1 = ggdata %>%
  {
    ggplot(data=.,aes(y=rho,x=dataID)) +
    geom_bar(width=barwidth,stat='identity',aes(fill=I('orange')),position=position_dodge(width=dodgewidth))+
    # scale_y_continuous(limits=c(rhomin,rhomax))+
    geom_errorbar(aes(y=rho,ymin=rho-rho_stderr,ymax=rho+rho_stderr))+
    coord_flip(ylim = c( rhomin, rhomax))+
    theme_bw()
  }


  ggdata %>% {qplot(data=.,geom='bar',stat='indentity',x=dataID,y=rho)}

  names(ggdata)[5] = "E.A."
  names(ggdata)[6] = "E.A._stderr"
  ggdata$E.A. %<>% as.numeric
  ggdata %<>% ungroup
  ggdata  %<>% arrange((E.A.))
  ggdata$dataID  %<>% factor(.,levels=unique(.))
  ggdata$dataID

  p2 = ggdata %>%
  {
    ggplot(data=.,aes(y=E.A.,x=dataID)) +
    geom_bar(width=barwidth,stat='identity',aes(fill=I('purple')),position=position_dodge(width=dodgewidth))+
    # scale_y_continuous(limits=c(-0.1,max(.$E.A.)*1.2))+
    geom_errorbar(aes(y=E.A.,ymin=E.A.-E.A._stderr,ymax=E.A.+E.A._stderr))+
    coord_flip(ylim = c( eamin, eamax))+
    theme_bw()
  }



  names(ggdata)[7] = "E.W."
  names(ggdata)[8] = "E.W._stderr"
  ggdata$E.W. %<>% as.numeric
  ggdata %<>% ungroup
  ggdata  %<>% arrange((E.W.))
  ggdata$dataID  %<>% factor(.,levels=unique(.))
  ggdata$dataID


  p3 = ggdata %>%
  {
    ggplot(data=.,aes(y=E.W.,x=dataID)) +
    geom_bar(width=barwidth,stat='identity',aes(fill=I('darkgreen')),position=position_dodge(width=dodgewidth))+
    # scale_y_continuous(limits=c(0,max(.$E.W.)*1.2))+
    geom_errorbar(aes(y=E.W.,ymin=E.W.-E.W._stderr,ymax=E.W.+E.W._stderr))+
    coord_flip(ylim = c(ewmin,ewmax))+
    theme_bw()
  }


    pdf(outputfile,h=8,w=18)
    # par(mfrow=c(3,1))
    # print(p1)
    # print(p2)
    # print(p3)
    multiplot(p1,p2,p3,cols=3)
    dev.off()
    #file2laptop('tmp.pdf')
}


datafiles  <-  list(
  flanking = "/g/furlong/project/28_B_DNASE/analysis/evolutionary_analyses/INSIGHT/scATAC/combinedResults/10_12h_bedFiles_byTissue_byPeaks_flanking_results_combinedResults.txt",
  fourD = "/g/furlong/project/28_B_DNASE/analysis/evolutionary_analyses/INSIGHT/scATAC/combinedResults/10_12h_bedFiles_byTissue_byPeaks_4d_results_combinedResults.txt"
)

outputfiles  <-  list(
  flanking = "10_12_flanking.tissueResults.insightplot.pdf",
  fourD = "10_12_4d.tissueResults.insightplot.pdf"
)


insightplots(
  datafile = datafiles[[1]],
  outputfile = outputfiles[[1]],
  rhomin = 0.6,rhomax = 0.8,
  eamin = 0,eamax = 2,
  ewmin=10,ewmax=20
)

insightplots(
  datafile = datafiles[[2]],
  outputfile = outputfiles[[2]],
    rhomin = 0.6,rhomax = 0.8,
    eamin = 0,eamax = 2,
  ewmin=10,ewmax=20
)
