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

# insight_polysights  <-  rtracklayer::import(config$data$validSitesgff)

#this function runs the asymptotic
#this program will run asymptotic


asymptotic_mk_test <- function(
	testregion,
	ctlregion,
	sites_gr,
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
	#	
	bins <- c(testsites$daf,ctlsites$daf)  %>% table

	p = table(testsites$daf)[names(bins)] %>% replace(.,is.na(.),0)
	p0 = table(ctlsites$daf)[names(bins)] %>% replace(.,is.na(.),0)

	#get the actual daf frequencies in each bin (rightmost)
	f  = names(bins) %>% as.numeric
	#run asymptotic mk
	asymptoticMK(d0=d0, d=d, 
		df=data.frame(f,p,p0,row.names=NULL),
		xlow=xlow, xhigh=xhigh, output=output, ...
	)

}

stop('stopping here')
# sitefile = config$data$validSites
message("Loading and processing polymorphism data")






# -------------------------------------------------------------------------------
# --------Now run the asymptotic MK test
# -------------------------------------------------------------------------------
import('')
Sys.glob("./scATAC/*bedFiles/*byTissue*/*byPeaks*")
Sys.glob("./scATAC/*bedFiles/byTissue/byPeaks/*.bed")

#usage
asym_mk_output <- asymptotic_mk_test(
	testregion,
	ctlregion,
	sites_gr=insight_polysights
	) %>% 
	select(model,a,b,c,CI_low,CI_high,alpha_original)


#then we'll run this on each of our grange blocks

#first just run it on all the clusters

#then run it on the extragenic ones




#the issue here is that deciding on the appropriate number of bins is tough, cos you
#kind of need to just look at the sites in test and block before you can do it...
#actually not sure if the number of bins matters...

# source('/g/furlong/project/28_B_DNASE/analysis/evolutionary_analyses/INSIGHT/src/dermot/run_asymptoticMK.R')