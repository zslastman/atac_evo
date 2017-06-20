# -------------------------------------------------------------------------------
# --------This wil read in the locations of the ATAC seq peaks and
# --------output gff 3 files, with timeoint and info on which clusters
# --------each is in 
# -------------------------------------------------------------------------------


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


#read read in the location o four clusters
noninsightfiles = clustsfiletbl %>% filter(data_source %>% str_detect('phylo')) %>%{ setNames(.$file,.$data_source)}
#rad in data and select the source/stat columns
phylophast_data  <- lapply(noninsightfiles,fread)
#make a combined df, using our data_sources as id columns
phylophast_data %<>% rbindlist(use.names=TRUE,idcol='data_source')

phylophast_data   %<>% 
	mutate(time = str_extract(data_source,'\\d+_\\d+h')) %>% 
	select(seqnames=chromosome,start=bed_start,end=bed_end,strand=strand,source,time)
phylophast_data %<>% distinct

phylophast_data %<>% DT2GR
