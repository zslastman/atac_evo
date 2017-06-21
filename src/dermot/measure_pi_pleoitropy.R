# -------------------------------------------------------------------------------
# --------This script is going to measure Pi as a function of 
# --------enhancer pleoitropy, i.e., see if overall levels of polymorphism change
# --------as a funciton of the number of clusters an ehancer is in
# -------------------------------------------------------------------------------
library(Hmisc)


stopifnot(is(clusters,"GRanges"))
stopifnot(clustnum %in% mcols(clusters))


sitefile <- "/g/furlong/garfield/projects/INSIGHT_analyses/valid_sites_file/validSites.txt"



sites <- config$data$validSitesgff3


#calculate pi for each cluster
calculate_pi  <- function(gr,sites){
	#get overlapping variants
	ovsites  <- subsetByOverlaps(sites,gr)
	#pi is simply the fraction of total basepairs with variants
	binconf(sum(width(ovsites) , sum(width(gr)))
}

clusts %>% 
	split(.,.$clustnum) %>% 
	map(  calculate_pi, sites = sites)

#Produce plots of pi as a function of cluster number
#just do it for the peaks for now.



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