#load packages
library(stringr)
library(dplyr)
library(data.table)
library(ggplot2)
library(magrittr)

#get all the data files
datafiles <- '.' %>% list.files(full=TRUE,pattern='bed.*.txt$')
#name them sensibly
tps  <- datafiles %>% basename %>% str_extract('.*?h')
stat  <- datafiles %>% basename %>% str_extract('phastCons|phyloP')
names(datafiles) <- paste0(stat,'_',tps)
#also sort them chronologically
chronorder = tps %>% str_extract('[0-9]+') %>% as.numeric %>% order
datafiles = datafiles[chronorder]

#rad in data and select the source/stat columns
data  <- lapply(datafiles,. %>% fread %>% select(statistic,source))
#make a combined df, using our names as id columns
data %<>% rbindlist(use.names=TRUE,idcol='data_source')
#group by file and source column, get bootstrap confidence limits
conflims <- data %>% group_by(data_source,source) %>% do(mean_cl_boot(.$statistic))

#rename the y column for the plot
colnames(conflims) %<>% str_replace('^y$','statistic')
#and reorder the data srouce column factor
conflims$data_source <- factor(conflims$data_source,levels = names(datafiles))

#generate plot
p1 <- conflims %>% 
  {  
    ggplot(data=.,aes(y=statistic,color=data_source,x=source)) +
    # geom_bar(width=barwidth,stat='identity',aes(fill=I('darkgreen')),position=position_dodge(width=dodgewidth))+
    # scale_y_continuous(limits=c(0,max(.)*1.2))+ 
    geom_pointrange(aes(y=statistic,x=source,ymin=ymin,ymax=ymax))+
    ggtitle('Constraint statistics by Timepoint/algorithm')+
    facet_grid(data_source~.)+
    theme_bw()
  }

#print the plot
pdf('evo_tp_program_barplot.pdf',h=14,w=7)
print(p1)
dev.off()