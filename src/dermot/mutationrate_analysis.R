
sharpdata <- fread(config$data$sharp_mutdata)
keightdata <- fread(config$data$keightly_mutdata)
sharpdata$Location %<>% str_replace_all(', \\d+','')

#load up our two sets of mutations as GRanges.
keightdata <- tidyGR(GRanges(keightdata$Chrom,IRanges(keightdata$Location,width=1)))
sharpdata <- tidyGR(GRanges(sharpdata$Chrom,IRanges(as.numeric(sharpdata$Location),width=1)))


# -------------------------------------------------------------------------------
# --------Testing the relationshiop between promoter shape and mutation rate
# -------------------------------------------------------------------------------
#Now our question is first of all.... 
#do we have different mutation rates in broad vs narrow regions
cagewinds <- import("data/Supple_Tables_for_paper/Schor_TableS1.gff3")
#filter out the ultrabroad ones (usually intragenic)
cagewinds$isnarrow  <- as.numeric(cagewinds$shape.index) > -1


cagewinds$hasmut_keight <- overlapsAny(cagewinds,keightdata)
cagewinds$hasmut_sharp <- overlapsAny(cagewinds,sharpdata)


with(cagewinds,table(hasmut_keight,isnarrow))
with(cagewinds,fisher.test(hasmut_keight,isnarrow))

cagewindschr2 <- cagewinds %>% subset(seqnames%in%c('chr2R','chr2L'))



with(
	cagewindschr2),
	table(hasmut_sharp,isnarrow)
)

with(
	cagewindschr2),
	fisher.test(hasmut_sharp,isnarrow)
)

with(
	cagewindschr2),
	glm(hasmut_sharp~as.numeric(shape.index))
) %>% summary

#also let's check total muts, though possibly not needed
bmutnum = length(sharpdata %>% subsetByOverlaps(cagewindschr2 %>% subset(!isnarrow))
nmutnum = length(sharpdata %>% subsetByOverlaps(cagewindschr2 %>% subset(isnarrow)))

c(bmutnum,
	cagewindschr2 %>% subset(!isnarrow) %>% width %>% sum,
	nmutnum,
	cagewindschr2 %>% subset(isnarrow) %>% width %>% sum
	) %>% matrix(ncol=2) %>%  fisher.test




# -------------------------------------------------------------------------------
# --------Testing to see if there's a relationship between cluster number and mut rate
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


#read and plot the phyloP and phastcons data
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

#now split the regions by time/source
#and get confidence intervals on the odds of them having mutations
keightovconflims <- 
	phylophast_data %>%
	split(.,paste0(.$source,'__',.$time)) %>%
	lapply(. %>% 
	{
		ov = overlapsAny(.,keightdata)
		data.frame(binconf(sum(ov),length(ov)),time=.$time[1],Cluster=.$source[1])
	}
	)

keightovconflims %<>% rbindlist

plotkeight  <- 
	keightovconflims %>% 
	{
		qplot(main='keightly mutation rates',data=.,y=PointEst,ymax=Upper,ymin=Lower,geom='linerange',x=Cluster)+
		facet_grid(time~.)+
		theme(axis.text.x = element_text(angle=45,hjust=1))

	}

sharpovconflims <- 
	phylophast_data %>%
	subset(seqnames%in%c('chr2R','chr2L')) %>% 
	split(.,paste0(.$source,'__',.$time)) %>%
	lapply(. %>% 
	{
		ov = overlapsAny(.,sharpdata)
		data.frame(binconf(sum(ov),length(ov)),time=.$time[1],Cluster=.$source[1])
	}
	)

sharpovconflims %<>% rbindlist


plotsharp  <- 
	sharpovconflims %>% 
	{
		qplot(main='sharp mutation rates (chr2)',data=.,y=PointEst,ymax=Upper,ymin=Lower,geom='linerange',x=Cluster)+
		facet_grid(time~.)+
		theme(axis.text.x = element_text(angle=45,hjust=1))
	}

pdf(config$plots$mut_cluster_plot)
print(plotkeight)
print(plotsharp)
dev.off()
