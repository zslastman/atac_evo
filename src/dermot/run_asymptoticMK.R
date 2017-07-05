# -------------------------------------------------------------------------------
# --------Run the asymptotik MK test on a bunhc of sites
# -------------------------------------------------------------------------------
require(data.table)
require(stringr)
require(dplyr)
require(devtools)
require(GenomicRanges)
require (magrittr)
require (tidyverse)
require (yaml)
install_or_load('nls2')
ANCESTRAL_P_THRESH <- 0.8
config = yaml.load_file(config_file)
source('./src/asymptoticMK/asymptoticMK_local.R')

insight_polysights  <-  rtracklayer::import(config$data$validSitesgff)
divergences_gr <-  import(config$data$dviergencesbed)

neutral_divergences  <- import(format='gff3',config$data$neutral_divergences)
neutral_polymorphisms  <- read_tsv(config$data$neutral_polymorphisms)

#this function runs the asymptotic
#this program will run asymptotic

frame_data(
~x,~y,
1,
1,
2,
0,
1,
0,
10,
3
) %>%
	mutate(y = 1:nrow(.)) %>%
	mutate(.cumzero = cumsum( x != 0 ) ) %>%
	group_by(.cumzero) %>%
	summarize(funs(sum)) %>%
	select(- .cumzero)


asymptotic_mk_test <- function(
	testregion,
	divergences_gr,
	sites_gr,
	neutral_divergences,
	neutral_polymorphisms,
	xlow = 0.1,
	xhigh = 0.9,
	output = "table"
	){
	#first test the inputs
	stopifnot(is(testregion,'GRanges'))
	stopifnot(is(sites_gr,'GRanges'))
	stopifnot('daf' %>% is_in(mcols(sites_gr) %>% colnames ))
	stopifnot('daf' %>% is_in(colnames(neutral_polymorphisms)))
	#
	#get divergence in test regions
	d  <- length(unique(subsetByOverlaps(divergences_gr,testregion)))
	#
	#pull out the sites for our test region
	testsites <- subsetByOverlaps(sites_gr,testregion) %>% unique
	#
	#Now figure out which neutral blocks we care about for our test region
	neut_blocks <- subsetByOverlaps(neutral_divergences,testregion)
	d0 <- sum(as.numeric(neut_blocks$L))
	ctlsites <- filter(neutral_polymorphisms,block %in% neut_blocks$blockName)
	#now figure out how many bins to segregate the sites into
	#
	bins <- c(testsites$daf,ctlsites$daf)  %>% table
	#
	p = table(testsites$daf)[names(bins)] %>% replace(.,is.na(.),0) %>% as.vector
	p0 = table(ctlsites$daf)[names(bins)] %>% replace(.,is.na(.),0) %>% as.vector
	#
	#
	#get the actual daf frequencies in each bin (rightmost)
	f  = names(bins) %>% as.numeric
	#
	#define df for the function, make sure p0 is never 0
	df <- data.frame(f,p,p0,row.names=NULL)
	df <-
		df %>%
		#by merging with the row above where it is
		group_by(cgroup = cumsum(p0!=0)) %>%
		summarise(f = min(f),p=sum(p),p0 = sum(p0)) %>%
		select(-cgroup)
	#run asymptotic mk
	out <-
		asymptoticMK(d0=d0, d=d, df=df,
			xlow=xlow, xhigh=xhigh, output=output
		) %>%
	#also just add polymorphism stats
	mutate(
		K = length(testsites$daf)/sum(width(testregion)) ,
		K0 = sum(as.numeric(neut_blocks$K))/sum(width(neut_blocks)),
		D=d, D0=d0,
		lenTest = sum(width(testregion)),
		lenBlocks = sum(width(neut_blocks))
	)
	out <- out %>% select()
	map_df(out,~ ifelse(is.factor(.),as.character(.),.) )
}

stop('stopping here')
# sitefile = config$data$validSites
message("Loading and processing polymorphism data")






# -------------------------------------------------------------------------------
# --------Now run the asymptotic MK test
# -------------------------------------------------------------------------------

#read in all our clusters
clustergrs <-
	c(
		Sys.glob("data/scATAC/*_bedFiles/*/byPeaks/*.bed"),
		Sys.glob("data/scATAC/*/per_cluster/*.bed")
	)	 %>%
	setNames(map(.,import),.)



blockgrs <-
	"/g/furlong/garfield/projects/INSIGHT_analyses/precomputed_blockData/5kb_withOverlaps_updated/"  %>%
	"/g/furlong/garfield/projects/INSIGHT_analyses/precomputed_blockData/5kb_withOverlaps_updated/"  %>%
	list.files(full=TRUE) %>%
	str_match('(chr.*?):(\\d+)\\-(\\d+)') %>%
	as_tibble %>%
	.[-1] %>%
	{GRanges(.[[1]],IRanges(start=as.numeric(.[[2]]),as.numeric(.[[3]])))}

#usage
asym_mk_output <-
	clustergrs %>%
	map(safely(asymptotic_mk_test),
		sites_gr=insight_polysights,
		divergences_gr = divergences_gr,
		neutral_divergences,
		neutral_polymorphisms
	)
#check for erros
stopifnot(asym_mk_output %>% map(2) %>% map_lgl(is,'error') %>% sum %>% `==`(0))
#take only results
asym_mk_output  %<>% map(1)

#now parse relevant info out of our files
metadata <- 
	names(asym_mk_output) %>% 
	str_match(regex('(?<=\\/)(\\d+_\\d+h).*?(\\d+)\\.bed')) %>% 
	as.data.frame %>% 
	set_colnames(c('V1','tp','clustnum')) %>% 
	mutate(clustnumtype = ifelse(str_detect(V1,'count'),'pleoitropy','Cluster_ID')) %>% 
	select(clustnumtype,tp,clustnum)
#and bind them all together with the metadata
asym_mk_output <- 
	asym_mk_output %>% 
	bind_rows %>% 
	assert_that(~nrow(.x) == nrow(metadata)) %>% 
	bind_cols(metadata)
asym_mk_output$clustnum %<>% as.character %>% as.numeric %>% factor
#order the tps
orderedtps <- unique(asym_mk_output$tp)[asym_mk_output$tp %>% unique %>%  str_extract('[0-9]+') %>% as.numeric %>% order]
asym_mk_output$tp %<>% factor(.,levels=orderedtps)


#now, get the asym alpha data with the confidence intervals
mkdata <- 
	asym_mk_output %>% 
	select(alpha_asymptotic,CI_low,CI_high,tp,clustnumtype,clustnum)
#and the other stats
	
mkpolystats <- asym_mk_output %>% 
	select(alpha_original,K,K0,D,D0,lenTest,lenBlocks,tp,clustnumtype,clustnum) %>% 
	mutate(p_div = D/lenTest , p0_div = D0/lenBlocks) %>% 
	gather(measure,measureval,-tp,-clustnumtype,-clustnum)

#list of other stats to plot
otherstatlist <- unique(mkpolystats$measure)

#then we'll run this on each of our grange blocks

#first just run it on all the clusters

#then run it on the extragenic ones


#for each time point (panel)
#for each clustnum (x axis)
#we'll have different plots showing different stats
#
# -------------------------------------------------------------------------------
# --------
# -------------------------------------------------------------------------------

mkpolystats
mkpolystats %<>% select(measure,measureval,clustnumtype,clustnum,tp)
#first do the other stats
plotob = list()
for(measure_i in otherstatlist){
for(clustnumtype_i in unique(mkpolystats$clustnumtype)){
	#generate plot
	plotob  %<>% append(list(
		mkpolystats %>%
		filter(measure==measure_i) %>%
		  {
		    ggplot(data=.,aes(y=measureval,x=clustnum)) +
		    geom_point(aes(y=measureval,x=clustnum),size = I(3),alpha=.5,color=I('black'))+
		    ggtitle(paste0(measure_i))+
		    scale_y_continuous(name = measure_i)+
		    scale_x_discrete(paste0(clustnumtype_i))+
			scale_alpha(guide=FALSE)+
		    facet_grid(tp~.,scale='free')+
		    theme_bw()
		  }
	))
}
}
stopifnot(plotob %>% has_length)
stopifnot(plotob %>% map_lgl(~ "ggplot" %in% class(.)) )
#print the plot
pdf(config$plots$mk_stat_plots,h=9,w=7)
plotob
dev.off()
config$plots$mk_stat_plots %>% file.info %$% size
#
file2laptop(config$plots$mk_stat_plots)
#


#And then the asym_mkresults
mkdata %<>% select(alpha_asymptotic,clustnumtype,clustnum,CI_low,CI_high,tp)
plotob = list()
for(clustnumtype_i in unique(mkpolystats$clustnumtype)){
	#generate plot
	plotob  %<>% append(list(
		mkdata %>%
		  {
		    ggplot(data=.,aes(y=alpha_asymptotic,x=clustnum)) +
		    geom_ribbon(aes(ymax = CI_high, ymin = CI_low,x=clustnum),alpha=.5,color=I('black'))+
		    ggtitle(paste0("alpha_asymptotic"))+
		    scale_x_discrete(name = paste0(clustnumtype_i))+
			scale_y_continuous(name = "alpha_asymptotic",limits = c(0,1))+
			scale_alpha(guide=FALSE)+
		    facet_grid(tp~.,scale='free')+
		    theme_bw()
		  }
	))
}
#print the plot
pdf(config$plots$asym_mk_plots,h=9,w=7)
plotob
dev.off()
#
file2laptop(config$plots$asym_mk_plots)
config$plots$asym_mk_plots %>% file.info %$% size




#the issue here is that deciding on the appropriate number of bins is tough, cos you
#kind of need to just look at the sites in test and block before you can do it...
#actually not sure if the number of bins matters...

# source('/g/furlong/project/28_B_DNASE/analysis/evolutionary_analyses/INSIGHT/src/dermot/run_asymptoticMK.R')