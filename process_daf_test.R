# -------------------------------------------------------------------------------
# --------This script will carry out some basic tests on our peaks as a function
# --------of pleoitropy, to see what
# -------------------------------------------------------------------------------
RARETHRESH <- 0.15
datafiletbl   <- read_tsv(config$results$input_metadata)
datafiletbl %>% glimpse
datafiletbl$threshold %>% unique

statfiles <-
	datafiletbl %>%
	filter(threshold=='byPeaks') %>%
	filter(stat %>% str_detect('sample')) %>%
	assert_that(has_rows) %>%
	.$file

polystats <-
	statfiles %>%
	setNames(map(.,read_tsv),.) %>%
	bind_rows(.id = "file") %>%
	mutate(file = basename(file)) %>% 
	mutate(tp = str_extract(file,regex('\\d+_\\d+h'))) %>%
	select(-file)

allsitestats <- 
	"/g/furlong/project/28_B_DNASE/analysis/evolutionary_analyses/INSIGHT/scATAC/daf_results" %>% 
	list.files(full=TRUE) %>% 
	str_select('validSites') %>% 
	read_tsv


polystats %<>% mutate(fractionDivergent = fixedDivergence/totalSites)
polystats %<>% mutate(ratio_rareToCommon = rareCount/commonCount)
polystats %<>% mutate(fractionNonAncestral = nonAncestralPolymorphisms/totalPolymorphisms)
polystats %<>% mutate(fracPoly = totalPolymorphisms/totalSites)

daftest <- function(rareCount,commonCount,...){
	# browser()
	#run a fisher test
	 test <- c(
	 		rareCount,
	 		commonCount+rareCount,
	 		allsitestats$rareCount,
	 		allsitestats$commonCount+allsitestats$rareCount
	 	)  %>%
	 	matrix(ncol=2) %>%
	 	fisher.test

	 #extract stuff from test and return
	 test %>%
	 	.[c('p.value','conf.int','estimate')]%>%
	 	flatten  %>%
	 	setNames(c('p.value','Lower','Upper','PointEst')) %>% 
	 	as.data.frame
}

#nwo runa  fisher test on our rare/common derived allele counts
polystats %<>% 
	mutate(fish=map2(rareCount,commonCount,daftest)) %>% 
	unnest

#numeric clustnum column
polystats %<>% mutate(
	clustnum = dataID %>% str_replace('count\\.','') %>% as.numeric
)

polystats$clustnum <- polystats$clustnum  %>% as.numeric
polystats %<>% arrange(clustnum)
polystats$clustnum %<>% factor(.,levels = unique(.))


# -------------------------------------------------------------------------------
# --------Repeat DAF test myself from site data
# -------------------------------------------------------------------------------
#read in all our clusters
clustergrs <-
	Sys.glob("data/scATAC/*_bedFiles/byTissue/byPeaks/*.bed") %>%
	setNames(map(.,import),.)

# #file info now in the names of our ranges, process tiepoint and
# #clustnum info out of it
# mcols(clustergrs)   <-
# 	clustergrs@ranges@NAMES %>%
# 	str_match(regex('(?<=\\/)(\\d+_\\d+h).*?\\.(\\d+)\\.bed')) %>%
# 	.[,-1] %>%
# 	set_colnames(c('tp','clustnum'))
# #and erase those names
# clustergrs@ranges@NAMES  <- NULL
# mcols(clustergrs)   %<>%  as.data.frame %>% map(as.character)

#this function will test whether derived alleles in a region are more likely to
#be rare than amongst sites in general.
daftest <- function(region,polysites,threshold = RARETHRESH){
	#test our inputs
	polysites %>%
		assert_that(
			~ is(.x,"GRanges"),
			~ length(.x) > 500,
			~ c("daf") %in% colnames(mcols(.x))
		)
	region %>%
		assert_that(
			~ is(.x,"GRanges")
			# ~ sum(width((.x)) > 1e3
		)
	#extract derived allele frequencies from the sites object
	regiondafs <- polysites %>%
		subsetByOverlaps(region) %>% .$daf %>% as.numeric
	nonregiondafs <- polysites %>%
		subsetByNegOverlaps(region) %>% .$daf  %>% as.numeric
	 #logical vectors saying which dafs are below threshold
	 regiondafs <- regiondafs < threshold
	 nonregiondafs <- nonregiondafs < threshold
	 #run a fisher test
	 test <- c(table(nonregiondafs),table(regiondafs))  %>%
	 	matrix(ncol=2) %>%
	 	fisher.test
	 #extract stuff from test and return
	 test %>%
	 	.[c('p.value','conf.int','estimate')]%>%
	 	flatten_dbl  %>%
	 	setNames(c('p.value','Lower','Upper','PointEst'))
}

# regions <- clustergrs %>% split(list(.$tp)) %>% as.list %>% map(~split(.x,.x$clustnum))
# flatnames <- regions %>% at_depth(1,names)  %>% map2(names(.),.,paste,sep='__') %>% flatten_chr

#peform daf tests for all files
dafresults <- clustergrs %>% map(safely(daftest),polysites)
#make sure it all worked
daferrors  <- dafresults  %>% map('error')
stopifnot(daferrors%>%map_lgl(is.null)%>%all)

#also calculate pi


#
dafresulttbl  <-
	dafresults  %>%
	map('result') %>%
	map_df(bind_rows)

dafresulttbl %<>% cbind(
	names(clustergrs) %>%
	str_match(regex('(?<=\\/)(\\d+_\\d+h).*?\\.(\\d+)\\.bed')) %>%
	.[,-1] %>%
	set_colnames(c('tp','clustnum'))
)

dafresulttbl$clustnum <- dafresulttbl$clustnum  %>% as.numeric
dafresulttbl %<>% arrange(clustnum)
dafresulttbl$clustnum %<>% factor(.,levels = unique(.))

# -------------------------------------------------------------------------------
# --------now plot DAF results
# -------------------------------------------------------------------------------


plotob = list()
#generate plot
plotob  %<>% append(list(
	polystats %>%
	# dafresulttbl %>%
	mutate(DAF = PointEst) %>%
	  {
	    ggplot(data=.,aes(y=PointEst,x=clustnum)) +
	    # geom_bar(width=barwidth,stat='identity',aes(fill=I('darkgreen')),position=position_dodge(width=dodgewidth))+
	    scale_y_continuous(name = 'DAF test 95% CI')+
	    scale_x_discrete(name = 'Number of Tissue clusters with Peak')+
	    # geom_pointrange(aes(y=statistic,x=clustnum,ymin=ymin,ymax=ymax))+
	    geom_ribbon(aes(y=PointEst,x=clustnum,ymin=Lower,ymax=Upper),color=I('green'))+
	    ggtitle(paste0("DAF test Odds Ratio, For scATAC peaks, as a function\n# of tissues with expression (clustnum)"))+
	    facet_grid(tp~.,scale='free')+
	    theme_bw()
		}
	))

#
#print the plot
pdf(config$plots$clustcount_daf_plots,h=9,w=7)
plotob
dev.off()
#
file2laptop(config$plots$clustcount_daf_plots)

