# -------------------------------------------------------------------------------
# --------This script will carry out some basic tests on our peaks as a function
# --------of pleoitropy, to see what
# -------------------------------------------------------------------------------
#threshold at which variants are counted as rare
RARETHRESH <- 0.15

#metadata on the various input files
datafiletbl   <- read_tsv(config$results$input_metadata)
datafiletbl %>% glimpse
datafiletbl$threshold %>% unique

#from correct files
statfiles <-
	datafiletbl %>%
	filter(threshold=='byPeaks') %>%
	filter(stat %>% str_detect('sample')) %>%
	assert_that(has_rows) %>%
	.$file
#read in polymorphism data
polystats <-
	statfiles %>%
	setNames(map(.,read_tsv),.) %>%
	bind_rows(.id = "file") %>%
	mutate(file = basename(file)) %>%
	mutate(tp = str_extract(file,regex('\\d+_\\d+h'))) %>%
	select(-file)

#also read in data on the background - all valid or 4d sites
allsitefiles <-
	"/g/furlong/project/28_B_DNASE/analysis/evolutionary_analyses/INSIGHT/scATAC/daf_results" %>%
	list.files(full=TRUE) %>%
	str_select('validSites|4dSites')
#read in
allsitestats <-
	allsitefiles %>%
	setNames(.,.) %>%
	map(read_tsv)  %>%
	bind_rows(.id='file')

#calculate various other stats
polystats %<>% mutate(fractionDivergent = fixedDivergence/totalSites)
polystats %<>% mutate(ratio_rareToCommon = rareCount/commonCount)
polystats %<>% mutate(fractionNonAncestral = nonAncestralPolymorphisms/totalPolymorphisms)
polystats %<>% mutate(fracPoly = totalPolymorphisms/totalSites)

#define a function that carries out the DAF test with a fisher test
daftest <- function(rareCount,commonCount,allrare,allcommon){
	#run a fisher test
	 test <- c(
	 		rareCount,
	 		commonCount+rareCount,
	 		allrare,
	 		allcommon+allrare
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

#split up the background data
allsitestats4d  <- allsitestats %>% filter(file %>% str_detect('4d'))
allsitestats  <- allsitestats %>% filter(! file %>% str_detect('4d'))
#use this split data to partially apply our test functon
allvaliddaftest <- partial(daftest,
	allrare=allsitestats$rareCount,
	allcommon=allsitestats$commonCount
)
fourDvaliddaftest <- partial(daftest,
	allrare=allsitestats4d$rareCount,
	allcommon=allsitestats4d$commonCount
)

#now apply the daftest rowwise to our data
polystats4d <- polystats %>%
	mutate(fish=map2(rareCount,commonCount,fourDvaliddaftest)) %>%
	unnest
#nwo runa  fisher test on our rare/common derived allele counts
polystats %<>%
	mutate(fish=map2(rareCount,commonCount,allvaliddaftest)) %>%
	unnest
#combine the results from both backgrounds in a single frame
polystats <- rbind(
	polystats %>% mutate(background = 'validsites'),
	polystats4d %>% mutate(background = 'fourDsites')
)
#numeric clustnum column, correctly ordered
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
	    ggplot(data=.,aes(y=PointEst,x=clustnum,color = background)) +
	    # geom_bar(width=barwidth,stat='identity',aes(fill=I('darkgreen')),position=position_dodge(width=dodgewidth))+
	    scale_y_continuous(name = 'DAF test 95% CI')+
	    scale_x_discrete(name = 'Number of Tissue clusters with Peak')+
	    # geom_pointrange(aes(y=statistic,x=clustnum,ymin=ymin,ymax=ymax))+
	    geom_ribbon(aes(y=PointEst,x=clustnum,ymin=Lower,ymax=Upper,color=background,fill=background))+
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

