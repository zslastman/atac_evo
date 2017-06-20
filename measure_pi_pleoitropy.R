# -------------------------------------------------------------------------------
# --------This script is going to measure Pi as a function of 
# --------enhancer pleoitropy, i.e., see if overall levels of polymorphism change
# --------as a funciton of the number of clusters an ehancer is in
# -------------------------------------------------------------------------------



stopifnot(is(clusters,"GRanges"))
stopifnot(clustnum %in% mcols(clusters))


sitefile <- "/g/furlong/garfield/projects/INSIGHT_analyses/valid_sites_file/validSites.txt"



sites <- config$data$validSites %>% fread


#calculate pi for each cluster
calculate_pi  <- function(gr,sites){
	#get overlapping variants
	ovsites  <- subsetByOverlaps(sites,gr)
	#pi is simply the fraction of total basepairs with variants
	sum(width(ovsites) / sum(width(gr))
}

clusts$pi  <-  calculate_pi(clusts,sites)

#Produce plots of pi as a function of cluster number
#just do it for the peaks for now.