################################################################################
##########Working with compressed genome files
################################################################################
make_file_tabix<-function(file){
 require(tools)
  filenoex = tools::file_path_sans_ext(file)
  file_ext= file_ext(file)
  filegz = paste0(filenoex,'.sort.',file_ext,'.gz')
  fileindex = paste0(filenoex,'.sort.',file_ext,'.gz.tbi')
  system(sprintf('less %s | sort -k 1,1 -k 2,2n | bgzip > %s',file,filegz))
  system(sprintf('tabix -p bed %s',filegz))
  stopifnot(file.exists(fileindex))
  return(filegz)
}


access_tabix <- function(tabixedfile,chrom,start,end) {
  require(tools)
  stopifnot(file.exists(paste0(tabixedfile,'.tbi')))
  stopifnot(file.exists(paste0(tabixedfile)))
  indexstring = paste0(chrom,':',start,'-',end)
  tryCatch({
      fread(sprintf('tabix %s %s ',tabixedfile,indexstring))
  },error=function(e){
      NULL
  })
}

file.copy(sitefile,'data/validsites_insight.txt')
make_file_tabix('data/validsites_insight.txt')



# -------------------------------------------------------------------------------
# --------Run the asymptotik MK test on a bunhc of sites
# -------------------------------------------------------------------------------
require(data.table)
require(stringr)
require(dplyr)
require(devtools)
require(GenomicRanges)
require (magrittr)
ANCESTRAL_P_THRESH <- 0.8
sitefile <- "/g/furlong/garfield/projects/INSIGHT_analyses/valid_sites_file/validSites.txt"


#this function runs the asymptotic

#this program will run asymptotic


asymptotic_mk_test <- function(
	testregion,
	ctlregion,
	sites_gr,
	binsize = 25,
	xlow = 0.1,
	xhigh = 0.9,
	output = "table"
	){
	#first test the inputs
	stopifnot(is(testregion,'GRanges'))
	stopifnot(is(ctlregion,'GRanges'))
	stopifnot(is(sites_gr,'GRanges'))
	stopifnot('daf' %>% is_in(mcols(sites_gr) %>% colnames ))

	#pull out the sites for our test region
	testsites <- subsetByOverlaps(sites_gr,testregion)
	stopifnot(testsites %>% length %>% `>`(binsize) )

	#and for our neutral region
	ctlsites <-  subsetByOverlaps(sites_gr,ctlregion)
	stopifnot(ctlsites %>% length %>% `>`(binsize) )

	#get total polymorphism
	d0  <- length(ctlsites)
	d  <- length(testsites)
	#now figure out how many bins to segregate the sites into
	#to do this we can
	#our target bin size is 50, so we want to pick the minimum
	#then divide our test and ctl into l/(b*2) pieces so it's
	#always b or above elements per bin
	#this code will ensure that there are binsize sites in each bin
	#though it won't care about the proportions in them.
	#	
	stopifnot(min(d0,d) > binsize )
	cutnum  <- floor(min(d0,d)/(binsize*2))
	#
	bins <- c(testsites$daf,ctlsites$daf) %>% cut_number(cutnum)
	#get the 
	p = table(bins[1:length(testsites)]) %>% as.vector
	p0 = table(bins[length(testsites)+1:length(bins)])%>% as.vector
	#get the actual daf frequencies in each bin (rightmost)
	f  =
		bins %>% levels %>%
		(regex('(?<=,)\\-?[0-9\\.]+')) %>%
		as.numeric
	#run asymptotic mk
	asymptoticMK(d0=d0, d=d, 
		df=data.frame(f,p,p0,row.names=NULL),
		xlow=xlow, xhigh=xhigh, output=output, ...
	)

}

stop('stopping here')
# sitefile = config$data$validSites
message("Loading and processing polymorphism data")


insight_polysights <-
	sitefile %>%
	sprintf(
		fmt = "cat %s  | grep %s | cut -f 1,3,7,8,9,10",
		"\\\tP\\\t"
	)%>%
	fread %>%
	set_colnames(c('seqnames','start','maj_ancest_prob',
		'min_ancest_prob','n_maj','n_min')
	)
#get rid of sites where both may be derived
insight_polysights %<>% filter(
		pmax(maj_ancest_prob,min_ancest_prob) > ANCESTRAL_P_THRESH
	)

#decide which one is derived, get it's frequency
insight_polysights  %<>%
	tbl_df %>% 
	mutate(maf = n_min / (n_maj + n_min)) %>%
	mutate(
		daf = ifelse(
			maj_ancest_prob > ANCESTRAL_P_THRESH,
			maf,
			1 - maf
		)
	)

#now make a granges object
insight_polysights %<>%
	with (  {
		GRanges(
			seqnames = seqnames,
			IRanges(start = start,w=1),
			daf = daf,
			maf = maf
		)
	} )



# -------------------------------------------------------------------------------
# --------Now run the asymptotic MK test
# -------------------------------------------------------------------------------

"scATAC/2_4h_bedFiles/per_cluster/"

testregion = GRanges('chr2L',IRanges(1,w=10e6))
ctlregion = GRanges('chr2R',IRanges(1,w=10e6))

#usage
asym_mk_output <- asymptotic_mk_test(
	testregion,
	ctlregion,
	sites_gr=insight_polysights
)


#then we'll run this on each of our grange blocks

#first just run it on all the clusters

#then run it on the extragenic ones




#the issue here is that deciding on the appropriate number of bins is tough, cos you
#kind of need to just look at the sites in test and block before you can do it...
#actually not sure if the number of bins matters...

# source('/g/furlong/project/28_B_DNASE/analysis/evolutionary_analyses/INSIGHT/src/dermot/run_asymptoticMK.R')