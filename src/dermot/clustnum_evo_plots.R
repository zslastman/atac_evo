#load packages
library(stringr)
library(dplyr)
library(data.table)
library(ggplot2)
library(magrittr)


#get all the data files
datafiletbl <- data.table(file = config$resfolder %>%
	list.files(full=TRUE,pattern='Tissue')
	)

# datafiles <- config$resfolder %>% list.files(full=TRUE,pattern='phastCons|phyloP')
#name them sensibly
datafiletbl %<>% mutate(
	tp = file %>% basename %>% str_extract('.*?h'),
	#parse the filenames into the algorithm/threshold
	stat = file %>% basename %>% str_extract(regex('(?<=byTissue_).*(?=_combinedResults.txt)')),

	threshold = stat %>% str_extract(regex('([0-9\\.]+|byPeaks)')),
	stat  = str_replace(stat,regex('([0-9\\.]+|byPeaks)_'),''),
	stat  = str_replace(stat,regex('(flanking_results|4d)'),'INSIGHT_\\1'),
	data_source = paste0(stat,'_',tp,'_',threshold),
	data_source = str_replace(data_source,'_NA',''),
	is_insight = stat %>% str_detect('INSIGHT')
)

#also sort them chronologically
datafiletbl %<>% arrange(
	tp %>% str_extract('[0-9]+') %>% as.numeric
)


#read and plot the phyloP and phastcons data
noninsightfiles = datafiletbl %>% filter(!is_insight) %>%{ setNames(.$file,.$data_source)}
#rad in data and select the source/stat columns
phylophast_data  <- lapply(noninsightfiles,. %>% fread %>% select(statistic,source))
#make a combined df, using our data_sources as id columns
phylophast_data %<>% rbindlist(use.names=TRUE,idcol='data_source')
#group by file and source column, get bootstrap confidence limits
conflims <- phylophast_data %>% group_by(data_source,source) %>% do(mean_cl_boot(.$statistic))
conflims %<>%  sleft_join(datafiletbl,by='data_source')

#rename the y column for the plot
colnames(conflims) %<>% str_replace('^y$','statistic')
#and reorder the data srouce column factor
conflims$data_source <- factor(conflims$data_source,levels = names(noninsightfiles))


plotob  <- list()

for(threshold_i in unique(conflims$threshold)){
#generate plot
plotob %<>% append(list( conflims %>% 
	filter(threshold==threshold_i) %>% 
  {  
    ggplot(data=.,aes(y=statistic,color=data_source,x=source)) +
    # geom_bar(width=barwidth,stat='identity',aes(fill=I('darkgreen')),position=position_dodge(width=dodgewidth))+
    # scale_y_continuous(limits=c(0,max(.)*1.2))+ 
    geom_pointrange(aes(y=statistic,x=source,ymin=ymin,ymax=ymax))+
    ggtitle(paste0('Constraint statistics by Timepoint/algorithm, threshold:' ,threshold_i))+
    facet_grid(data_source~.)+
    theme_bw()
  }
 ))
}
#print the plot
pdf(config$plots$phylop_phastcons_plot,h=7,w=7)
plotob
dev.off()

file2laptop(config$plots$phylop_phastcons_plot)


#read and plot the INSIGHT data
insightfiles =  datafiletbl %>% filter(is_insight) %>% {setNames(.$file,.$data_source)}
#rad in data and select the source/stat columns
loadinsight  <- function(x) x %>% fread %>% select(dataID,thres,matches('^rho'),matches('^E\\.A'),matches('^E\\.W'))
insight_data  <- lapply(insightfiles,. %>% loadinsight)
#make a combined df, using our names as id columns
insight_data %<>% rbindlist(use.names=TRUE,idcol='data_source')
#checkw e have unique ids
insight_data %>% group_by(data_source,dataID,thres) %>% tally  %>% .$n %>% eq(1) %>% all
#group by file and source column, get bootstrap confidence limits

#rename the y column for the plot
colnames(conflims) %<>% str_replace('^y$','statistic')
insight_data  %<>% group_by(data_source,dataID,thres)
#reshape the data
insight_data <- rbind(
	insight_data %>% 
	select(E.A.,E.A._stderr) %>% 
	transmute(thres,measure='E.A.',y=E.A.,ymin=E.A.- E.A._stderr,ymax=E.A.+ E.A._stderr),	
	insight_data %>% 
	select(E.W.,E.W._stderr) %>% 
	transmute(thres,measure='E.W.',y=E.W.,ymin=E.W.- E.W._stderr,ymax=E.W.+ E.W._stderr)
)
#checkw e have unique ids
insight_data %>% group_by(data_source,dataID,thres,measure) %>% tally  %>% .$n %>% eq(1) %>% all
#now safely join to get file info - safe join ensures exactly 1 match per line
insight_data %<>%  sleft_join(datafiletbl,by='data_source')

#and reorder the data srouce column factor
insight_data$clustnum <- insight_data$dataID %>% str_extract('\\d+')  %>% as.numeric
insight_data %<>% arrange(clustnum)
insight_data$clustnum %<>% factor(.,levels = unique(.))
#now let's arrange these plots sensibly
insight_data %<>% arrange(measure,tp,threshold)
insight_data$data_source %<>% factor(.,levels=unique(.))
insight_data %>% glimpse

plotob = list()
for(measure_i in c('E.W.','E.A.','rho')){
for(tp_i in unique(insight_data$tp)){
	for(threshold_i in unique(insight_data$threshold)){
		# for(data_source_i in unique(insight_data$data_source)){
			#generate plot
			plotob  %<>% append(list(
				insight_data %>%
				filter(thres %in% unique(thres)[seq(1,length(unique(thres)),by=2)]) %>% 
				mutate(thres = factor(thres)) %>% 
				mutate(statistic=y) %>%
				# filter(data_source == data_source_i) %>%
				filter(tp == tp_i) %>%
				filter(threshold == threshold_i) %>%
				filter(measure==measure_i) %>%
				  {
				    ggplot(data=.,aes(y=statistic,color=thres,x=clustnum)) +
				    # geom_bar(width=barwidth,stat='identity',aes(fill=I('darkgreen')),position=position_dodge(width=dodgewidth))+
				    # scale_y_continuous(limits=c(0,max(.)*1.2))+
				    # geom_pointrange(aes(y=statistic,x=clustnum,ymin=ymin,ymax=ymax))+
				    geom_ribbon(aes(y=statistic,x=clustnum,ymin=ymin,ymax=ymax,fill=thres,group=thres),alpha=.5,color=I('black'))+
				    ggtitle(paste0(measure_i,' for ',tp_i, ' threshold: ',threshold_i))+
				    scale_fill_manual(values = heat.colors(n_distinct(.$thres)))+
					scale_alpha(guide=FALSE)+
				    facet_grid(stat~.,scale='free')+
				    theme_bw()
				  }
			))
		# }
	}
}
}

#
#print the plot
pdf(config$plots$clustcount_insight_plots,h=9,w=7)
plotob
dev.off()
#
file2laptop(config$plots$clustcount_insight_plots)


