ANCESTRAL_P_THRESH <- 0.8
LOW_ANCESTRAL_P_THRESH <- 0.2
# -------------------------------------------------------------------------------
# --------Read in valid sites used by insight, getting daf where possible
# -------------------------------------------------------------------------------
sitefile <-  config$data$validSites
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

#mark sites where both may be derived
insight_polysights %<>% mutate(
		both_derived = pmax(maj_ancest_prob,min_ancest_prob) < ANCESTRAL_P_THRESH
	)


#decide which one is derived, get it's frequency
insight_polysights  %<>%
	tbl_df %>% 
	#get minor allele frequency
	mutate(maf = n_min / (n_maj + n_min)) %>%
	mutate(
		#get daf
		daf = ifelse(
			maj_ancest_prob > ANCESTRAL_P_THRESH,
			maf,
			1 - maf
		),
		#NA when both are derived
		daf = ifelse( 
			pmax(maj_ancest_prob,min_ancest_prob) < ANCESTRAL_P_THRESH,NA,daf)
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

insight_polysights %>% export(config$data$validSitesgff)

stopifnot(config$data$validSitesgff %>% import %>% is("GRanges"))

#also git divsites
divergences <- 
	sitefile %>%
	sprintf(
		fmt = "cat %s | grep %s | awk '($7 < %s)' |	cut -f 1,3",
		"\\\tM\\\t",
		# fmt = "cat %s  | head| grep %s | cut -f 1,3,7,8,9,10",
		LOW_ANCESTRAL_P_THRESH
	)%>%
	fread

divergences %<>% {GRanges(.$V1,IRanges(as.numeric(.$V2),w=1))}

divergences %>% export(config$data$dviergencesbed)