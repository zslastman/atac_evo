
#The results in this file apply to two main types of data:

#First, for each time point, we have “bedFiles_per_cluster”. 
#These datasets describe all of the peaks that define (are differentially accessibly within) a given cluster. 
#The identity of a given cluster can be found by looking at “/g/furlong/project/28_B_DNASE/analysis/scATAC/current_tSNE_data/Garfield_Heatmaps/combined_cluster_assignments.tsv”


#load packages
library(stringr)
library(dplyr)
library(data.table)
library(ggplot2)
library(magrittr)

neuralregex = 'Glia|CNS|PNS|Neural|Brain'
#results are split by timepoint, bytissue|by cluster, 4d|flank|flank+dnase
cluster_info_tbl <- config$cluster_assignment_file %>% 
	fread %>% 
	mutate(time=time %>% str_replace('\\-','_')) %>% 
	mutate(time=time %>% str_replace('hr','h'))

config$resfolder %>% list.files
#get all the data files
clustsfiletbl <- data.table(file = config$resfolder %>%
	list.files(full=TRUE,pattern='per_cluster')
	)
# datafiles <- config$resfolder %>% list.files(full=TRUE,pattern='phastCons|phyloP')
#name them sensibly
clustsfiletbl %<>% mutate(tp = file %>% basename %>% str_extract('.*?h'))
	#parse the filenames into the algorithm/threshold
clustsfiletbl %<>% mutate(
	stat = file %>% basename %>% 
	str_extract(regex('(?<=per_cluster_).*(?=_combinedResults.txt)'))
)
clustsfiletbl %<>% mutate(
	data_source = paste0(stat,'_',tp),
	is_insight = stat %>% str_detect('phastCons|phyloP') %>% not
)
#also sort them chronologically
clustsfiletbl %<>% arrange(
	tp %>% str_extract('[0-9]+') %>% as.numeric
)




#load up daf results for the clusters
dafresults <- config$resfolder %>%
	list.files(full=TRUE,patt='bedFiles_per_cluster.*sample_')  %>%
	lapply( fread)  %>%
	rbindlist(use.names=TRUE)
#unique ID col
stopifnot(dafresults$dataID %>% duplicated %>% not)

# #load up data on peak-> cluster assignments - not necessary with new files
# peakclustassignments <- config$peak_clust_bed_glob %>%
# 	Sys.glob %>%
# 	setNames(.,.) %>%
# 	lapply(fread) %>%
# 	rbindlist(use.names=TRUE,idcol='file') %>%
# 	mutate(tp = file %>% basename %>% str_extract('.*?h')) %>%
# 	mutate(clust = file %>% basename %>% str_extract(regex('(?<=clust_)\\d+'))) %>% 
# 	set_colnames(w('file chromosome bed_start bed_end feature_name score strand tp cluster'))


# -------------------------------------------------------------------------------
# --------Plotting hte PHyloPPHastcons results
# -------------------------------------------------------------------------------
#read and plot the phyloP and phastcons data
noninsightfiles = clustsfiletbl %>% filter(!is_insight) %>%{ setNames(.$file,.$data_source)}
#rad in data and select the source/stat columns
phylophast_data  <- lapply(noninsightfiles,. %>% fread)
#make a combined df, using our data_sources as id columns
phylophast_data %<>% rbindlist(use.names=TRUE,idcol='data_source')
phylophast_data %<>% select(-score)
# #get info on file
# phylophast_data %<>%  sleft_join(clustsfiletbl,by='data_source')
#now get cluster assignments
# stopifnot(phylophast_data$tp%in%peakclustassignments$tp)
# phylophast_data %<>% sleft_join(allow.dup=TRUE,by=c('tp','chromosome','bed_start','bed_end','strand'),
# 	peakclustassignments %>%
# 		select(chromosome,bed_start,bed_end,strand,tp,cluster)
# 	)


#group by file and source column, get bootstrap confidence limits
conflims <- phylophast_data %>% group_by(data_source,source) %>% do(mean_cl_boot(.$statistic))
#rename the y column for the plot
colnames(conflims) %<>% str_replace('^y$','statistic')
#now join in the info on the clusters
conflims %<>% mutate(time = str_extract(data_source,'\\d+_\\d+h'))

conflims %<>% tbl_df %>%
	rename(Cluster=source) %>% 
	sleft_join(
		allow.duplic=TRUE,
		cluster_info_tbl %>% 
			select(time,assignment,CatName,Cluster),
		by=c('Cluster','time')
	)
#remove the collision clusters
conflims %<>% filter(!str_detect(assignment,'ollision'))	
#add is_neural column
conflims  %<>% mutate(is_neural = str_detect(assignment,neuralregex))
conflims %<>% arrange(desc(is_neural))
#and reorder the data srouce column factor
conflims$data_source <- factor(conflims$data_source,levels = names(noninsightfiles))
# conflims %<>% arrange(conflims$Cluster %>% str_replace('clust_','') %>% as.numeric)
conflims$Cluster %<>% factor(.,levels=unique(.))


plotob  <- list()
#generate plot
plotob %<>% append(list( conflims %>% 
  {  
    ggplot(data=.,aes(y=statistic,color=is_neural,x=Cluster)) +
    # geom_bar(width=barwidth,stat='identity',aes(fill=I('darkgreen')),position=position_dodge(width=dodgewidth))+
    # scale_y_continuous(limits=c(0,max(.)*1.2))+ 
    geom_pointrange(aes(y=statistic,x=Cluster,ymin=ymin,ymax=ymax))+
    ggtitle(paste0('Constraint statistics for clusters' ))+
    facet_grid(data_source~.,scale='free')+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45,hjust=1))
  }
 ))

#
#print the plot
pdf(config$plots$clusts_phylop_phastcons_plot,h=12,w=7)
print(plotob)
dev.off()
#
file2laptop(config$plots$clusts_phylop_phastcons_plot)




# -------------------------------------------------------------------------------
# --------Same but using only neural/non neural as the grouping, instead of clusters
# -------------------------------------------------------------------------------
phylophast_data
#now join in the info on the clusters
phylophast_data %<>% mutate(time = str_extract(data_source,'\\d+_\\d+h'))

phylophast_data %<>% tbl_df %>%
	rename(Cluster=source) %>% 
	sleft_join(
		allow.duplic=TRUE,
		cluster_info_tbl %>% 
			select(time,assignment,CatName,Cluster),
		by=c('Cluster','time')
	)
#remove the collision clusters
phylophast_data %<>% filter(!str_detect(assignment,'ollision'))	
#now add in the is_nueral column
phylophast_data  %<>% mutate(is_neural = str_detect(assignment,neuralregex))
phylophast_data %<>% arrange(desc(is_neural))
#group by file and source column, get bootstrap confidence limits
conflims <- phylophast_data %>% distinct(chromosome,bed_start,bed_end,is_neural,.keep_all=TRUE) %>%  group_by(data_source,is_neural) %>% do(mean_cl_boot(.$statistic))
#rename the y column for the plot
colnames(conflims) %<>% str_replace('^y$','statistic')
#and reorder the data srouce column factor
conflims$data_source <- factor(conflims$data_source,levels = names(noninsightfiles))
# conflims %<>% arrange(conflims$Cluster %>% str_replace('clust_','') %>% as.numeric)
conflims$is_neural <- conflims$is_neural %>% ifelse('Neural','Non_Neural')


plotob  <- list()
#generate plot
plotob %<>% append(list( conflims %>% 
  {  
    ggplot(data=.,aes(y=statistic,color=is_neural,x=is_neural)) +
    # geom_bar(width=barwidth,stat='identity',aes(fill=I('darkgreen')),position=position_dodge(width=dodgewidth))+
    # scale_y_continuous(limits=c(0,max(.)*1.2))+ 
    geom_pointrange(aes(y=statistic,x=is_neural,ymin=ymin,ymax=ymax))+
    ggtitle(paste0('Constraint statistics for clusters' ))+
    facet_grid(data_source~.,scale='free')+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45,hjust=1))
  }
 ))
#
#print the plot
pdf(config$plots$isneural_phylop_phastcons_plot,h=12,w=7)
print(plotob)
dev.off()
#
file2laptop(config$plots$isneural_phylop_phastcons_plot)



# -------------------------------------------------------------------------------
# --------Plot insight data
# -------------------------------------------------------------------------------
#read and plot the INSIGHT data
insightfiles =  clustsfiletbl %>% filter(is_insight) %>% filter(!file %>% str_detect('sample')) %>% {setNames(.$file,.$data_source)}
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
insight_data %<>% tbl_df
insight_data  %<>% group_by(data_source,dataID,thres)
#reshape the data
insight_data <- rbind(
	insight_data %>% 
	select(rho,rho_stderr) %>% 
	transmute(measure='rho',y=rho,ymin=rho- rho_stderr,ymax=rho+ rho_stderr),		
	insight_data %>% 
	select(E.A.,E.A._stderr) %>% 
	transmute(measure='E.A.',y=E.A.,ymin=E.A.- E.A._stderr,ymax=E.A.+ E.A._stderr),	
	insight_data %>% 
	select(E.W.,E.W._stderr) %>% 
	transmute(measure='E.W.',y=E.W.,ymin=E.W.- E.W._stderr,ymax=E.W.+ E.W._stderr)
)
#checkw e have unique ids
insight_data %>% group_by(data_source,dataID,thres,measure) %>% tally  %>% .$n %>% eq(1) %>% all
#now safely join to get file info - safe join ensures exactly 1 match per line
insight_data %<>%  sleft_join(clustsfiletbl,by='data_source')

insight_data %<>% ungroup %>% mutate(dataID = str_replace(dataID,'\\d+_\\d+h_',''))

stopifnot(insight_data$dataID %in% cluster_info_tbl$Cluster)
stopifnot(insight_data$tp %in% cluster_info_tbl$time)


#now join in the info on the clusters
insight_data  %<>%  tbl_df %>%	rename(Cluster=dataID,time=tp) 


insight_data  %<>%
	sleft_join(
		allow.duplic=TRUE,
		cluster_info_tbl %>% 
			select(time,assignment,CatName,Cluster),
		by=c('Cluster','time')
	)

#remove the collision clusters
phylophast_data %<>% filter(!str_detect(assignment,'ollision'))	
# insight_data %<>% mutate(is_neural = str_detect(assignment,'Glia|CNS|PNS|Neural|Brain'))
insight_data %<>% mutate(is_neural = str_detect(assignment,neuralregex))
insight_data %<>% arrange(is_neural)
#checkw e have unique ids
stopifnot(insight_data %>% group_by(data_source,Cluster,thres,measure) %>% tally  %>% .$n %>% eq(1) %>% all)

#and reorder the data srouce column factor
insight_data$clustnum <- insight_data$Cluster %>% str_extract('\\d+$')  %>% as.numeric
insight_data %<>% arrange(clustnum)
insight_data$clustnum %<>% factor(.,levels = unique(.))
#now let's arrange these plots sensibly
insight_data %<>% arrange(measure,time)
insight_data$data_source %<>% factor(.,levels=unique(.))

insight_data %>% filter(!is_neural)  %>% distinct(time,Cluster)
insight_data %>% filter(!is_neural)  %>% filter(time==time[1]) %>% distinct(Cluster)
insight_data %>% filter(!is_neural)  %>% filter(time==time[1]) %>% distinct(assignment)
insight_data %>% filter(is_neural)  %>% filter(time==time[1]) %>% distinct(assignment)

sortedclusters = unique(insight_data$Cluster)
sortedclusters %<>% .[order(as.numeric(str_extract(.,'\\d+$')))]

thres_i = 0.25
plotob = list()
measures = c('rho','E.A.')
for(measure_i in measures){
for(time_i in unique(insight_data$time)){
		# for(data_source_i in unique(insight_data$data_source)){
			#generate plot
			plotob  %<>% append(list(
				insight_data %>%
				ungroup() %>% 
				# filter(thres %in% unique(thres)[seq(1,length(unique(thres)),by=2)]) %>% 
				filter(thres == thres_i) %>% 
				# mutate(thres = factor(thres)) %>% 
				mutate(statistic=y) %>%
				filter(time == time_i) %>%
				# filter(data_source == data_source[1]) %>%
				filter(measure==measure_i) %>%
				  {
				    ggplot(data=.,aes(y=statistic,color=is_neural,x=Cluster)) +
				    # geom_bar(width=barwidth,stat='identity',aes(fill=I('darkgreen')),position=position_dodge(width=dodgewidth))+
				    # scale_y_continuous(limits=c(0,max(.)*1.2))+
				    # geom_pointrange(aes(y=statistic,x=clustnum,ymin=ymin,ymax=ymax))+
				    geom_pointrange(aes(y=statistic,x=Cluster,ymin=ymin,ymax=ymax,fill=is_neural))+
				    ggtitle(paste0(measure_i,' for ',time_i,' MAF threshold: ', thres_i))+
					scale_alpha(guide=FALSE)+
					scale_x_discrete(limits = sortedclusters,breaks = sortedclusters)+
					scale_y_continuous(name=measure_i)+
				    facet_grid(data_source~.,scale='free')+
				    theme_bw()+
				    theme(axis.text.x = element_text(angle=45,hjust=1))
				  }
			))
		# }
}
}
#
#print the plot
pdf(config$plots$clustcount_insight_plots,h=6,w=7)
plotob
dev.off()
#
file2laptop(config$plots$clustcount_insight_plots)







# list.files(atac_res_folder,pattern= 'bedFiles_per_cluster',full=TRUE) %>%
# 	lapply(fread) %>%
# 	stackdfs(file) %>%
# 	mutate(
# 		time = str_extract(file,'+h'),
# 		neutral = str_extract(file,'+h'),

# 		)

# #now we have the insight results per cluster, per timepoint, and per flanking site type





# #object that define the actual clusters used.
# clustergrs  <-  list.files(atac_res_folder,pattern= 'bedFiles_per_cluster' ) %>%
# 	lapply(rtracklayer::import) 



# #We want to:
# #Address the hypothesis that one particular cluster or tissue is driving the
# #clusternumber -> constraint results.
# #We'll do this simply by plotting the insight,phyloP, and phastCons results per cluster first.


# #get all the data files
# phastcons_files <- atac_res_folder %>% list.files(full=TRUE,pattern='(phastCons)')
# phylop_files <- atac_res_folder %>% list.files(full=TRUE,pattern='(phyloP)')

# #name them sensibly
# tps  <- datafiles %>% basename %>% str_extract('.*?h')










