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
	with ({
		GRanges(
			seqnames = seqnames,
			IRanges(start = start,w=1),
			daf = daf,
			maf = maf
		)
	})

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
divergences %>% export(config$data$divergencesbed)


config$data$insight_ctl_flank_folder <- "/g/furlong/garfield/projects/INSIGHT_analyses/precomputed_blockData/50kb_withOverlaps_final"

#read in table of stats
neutral_divergences <-
	config$data$insight_ctl_flank_folder %>%
	list.files(full=TRUE,pattern="blockData_as_table") %>%
	fread
#as a Granges object
neutral_divergences <- 
	neutral_divergences[[1]] %>%
	str_match(regex('(.*):(\\d+)\\-(\\d+)')) %>%
	.[,-1] %>%
	{GRanges(.[,1],IRanges(as.numeric(.[,2]),as.numeric(.[,3])))} %>% 
	{mcols(.) <- neutral_divergences;.}

#we also need to get the list of polymorphisms for each block.
#since this might be redundant we'll need to just get all the relevant 
#polymorphisms, and list their blocks as well.
neutral_polymorphisms <- 	
	config$data$insight_ctl_flank_folder %>%
	paste0("/neutralSites/") %>% 
	list.files(full=TRUE,pattern="neutral.txt")  %>% 
	setNames(.,basename(.)) %>% 
	map(read_tsv,skip=1,col_names=FALSE,col_types = cols()) %>% 
	bind_rows(.id='block')

neutral_polymorphisms  %<>% 
	set_colnames(
		c('block','type','site','poly','maj_ancest_prob',
		'min_ancest_prob','n_maj','n_min')
	) %>% 
	select(site,maj_ancest_prob,n_min,n_maj,block) %>% 
	mutate(maf = n_min / (n_maj + n_min)) %>%
	mutate(daf = ifelse(
			maj_ancest_prob > ANCESTRAL_P_THRESH,
			maf,
			1 - maf
			)) %>% 
	select(block,site,daf)  %>% 
	mutate(block = str_replace(block,'.neutral.*',''))

neutral_divergences %>% write_tsv(config$data$neutral_divergences)
neutral_polymorphisms %>% write_tsv(config$data$neutral_polymorphisms)

