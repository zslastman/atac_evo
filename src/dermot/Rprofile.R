# -------------------------------------------------------------------------------
# --------This is the boilerplate R code for the eRNA project that gets run for
# --------all of my (dermot's) code. Hopefully we'll have a single one for the project
# --------so far.
# -------------------------------------------------------------------------------

config_file <- "atac_evo_config.yaml"
#set the location of our libraries
stopifnot(file.exists("./env"))
stopifnot(file.exists("./data"))
stopifnot(file.exists("./src"))

#load the conda installed biocinstaller
library(BiocInstaller)

#this function loads packages and installs them if they don't exist
install_or_load <- function(libname, cranurl = "http://cran.rstudio.com/"){
  package_missing <- ! require(
        warn.conflicts = FALSE,
        character.only = TRUE,
        libname
  )
  if (package_missing){
      #the biocLite function (from the bioconductor website)
      #will install from CRAN if a package isn't on biocnductor
      biocLite(
        suppressUpdates = FALSE,
        suppressAutoUpdate = TRUE,
        libname
      )
      require(libname, warn.conflicts = FALSE, character.only = TRUE)
  }
  else{
    require(libname, warn.conflicts = FALSE, character.only = TRUE)
  }
}




install_or_load('data.table')
install_or_load('GenomicRanges')
install_or_load('rtracklayer')
install_or_load('dplyr')
install_or_load('yaml')
install_or_load('XML')
install_or_load('tidyverse')
install_or_load('BSgenome.Dmelanogaster.UCSC.dm3')
install_or_load('stringr')
install_or_load('checkmate')
install_or_load('tidyverse')



# #clear objects
# rm(list=ls())

#load the config file as a list
config = yaml.load_file(config_file)

desc=dplyr::desc
select=dplyr::select
summarize=dplyr::summarize
summarise=dplyr::summarise
equals=magrittr::equals
slice=dplyr::slice
mutate=dplyr::mutate
import=rtracklayer::import
matches=dplyr::matches
#define some variables for GRanges objects
chrs.keep<-seqnames(Dmelanogaster)[!grepl('U|M', (seqnames(Dmelanogaster)))]
si<-seqinfo(Dmelanogaster)[chrs.keep]


################################################################################
##########Functions
################################################################################
#shorthand functions
not=  magrittr::not
eq =  magrittr::equals
gt = is_greater_than = `>`
lt = is_less_than = `<`
rename = dplyr::rename 
revalue= plyr::revalue



shape.index<-function(v){
  v<-v[v!=0]
  s<-sum(v)
  2+sum((v/s)*log2(v/s))
}

# Carries out Views on a GRange object ------------------------------------
GRViews<-function(rle,gr){
  stopifnot(is(rle,'RleList'))
  stopifnot(is(gr,'GenomicRanges'))

  gr = sort(gr) #sort
  
  v=Views(rle[seqlevels(gr)],as( gr[,NULL],'RangesList') )[seqlevels(gr)] #get views
  
  v=v[sapply(v,function(x){length(x)!=0})] #get a vector, leaving out empty chromosomes
  
  return(v)
}

mergeGR2DT = function(mergegr){
  grcols = mergegr%>%vapply(is,TRUE,'GRanges')%>%which
  mergedt=cbind(
    mergegr[[grcols[[1]]]]%>%GR2DT,
    mergegr[[grcols[[2]]]]%>%GR2DT
  )
  cnames = colnames(mergedt)
  cnames[duplicated(cnames)]%<>%paste0('.1')
  setnames(mergedt,cnames)
  mergedt
}

subsetByNegOverlaps<-function(gr,gr2,...){
  isov = countOverlaps(gr,gr2,...)>0
  gr[!isov]
}

GR2DT = function(gr){
  if(is.data.table(gr)) {
    warning('already a data table')
    return(gr)
  }
  #all columns must be vectors
  for(i in seq_len(ncol(mcols(gr)))){
    mcols(gr)[[i]] = as.vector(mcols(gr)[[i]])
  }

  dt = as.data.frame(gr,row.names=NULL,optional=FALSE)%>%as.data.table
  dt$strand= as.character(dt$strand)
  # setkey(dt,seqnames,strand,start,end)

  stopifnot( all(colnames( mcols(gr)) %in% colnames(dt) )  )
  dt
}

getStrandedviewlists <- function(gr,SRLs){
  require(quietly=TRUE,Matrix)
  #grange object and list of simple RLElists goes in
  #array comes out with dimensions [n-granges,GrangeWidth,tracknumber]
  #original order of granges  
  origorder = match(gr,sort(gr))
  #sort as this is how views will output the data
  gr = sort(gr)
  #needs to be a list so I can iterate over it (otherwise would iterate over chrs)
  if(is(SRLs,'SimpleRleList')) SRLs = list(track=SRLs)
  #iterate over the tracks
    SRL=SRLs[[1]]
  viewtracklist=lapply(SRLs,function(SRL){#use mclapply to iterate over the SRLSs
    #if the track has two subtracks for each strand
    if( identical(names(SRL)[1:3],c('pos','neg','both')) ) {#if we have pos and neg run seperately on them
      cat('.')
      ispos  = as.logical(strand(gr)=='+')#Splitting our granges by strand
      isneg  = as.logical(strand(gr)=='-')
      isstar = as.logical(strand(gr)=='*')
      #strands are okay
      stopifnot(all( ispos | isneg | isstar))
      #now fill in appropriate rows of the matrix, with the views from the appropriate pos or neg track
      viewlist=list()
      if( any(ispos) ){  viewlist[ ispos ] =  as.list( unlist (  IRanges::viewApply(GRViews(gr=gr[ ispos ],SRL$pos),as.numeric,simplify=FALSE) ))   }
      if( any(isstar) ){  viewlist[ isstar ] =  as.list( unlist (GRViews(gr=gr[ isstar ],SRL$both) %>% viewApply(as.numeric,simplify=FALSE) ))  }
      if( any(isneg) ){  viewlist[ isneg ] =   as.list( unlist (GRViews(gr=gr[ isneg ],SRL$neg) %>% viewApply(as.numeric,simplify=FALSE) ))  }
      #we need to flip the columns of the  matrix for negatives
      viewlist[isneg] = lapply(viewlist[isneg],rev)
      #return in the original order
      return(viewlist[origorder])
    }else{#Or if not stranded....
      cat('.')
      #get the views as a matrix
     viewlist =  as.list( unlist (GRViews(gr=gr,SRL) %>% viewApply(as.numeric,simplify=FALSE) )) 
     isneg  = as.logical(strand(gr)=='-')
      #now reverse the columns that are negative stranded
     viewlist[isneg] = lapply(viewlist[isneg],rev)
      #return in the original order
     return(viewlist[origorder]) 
    }
  })
  stopifnot(not(vapply(viewtracklist,is.error,TRUE)))
   #
  lapply(seq_along(gr),function(i){
    mat = sapply(viewtracklist,'[[',i)
    Matrix(sparse=TRUE,mat)
  }) 
}

add_ov_width <- function(dt,colname = 'ov_width'){ 
  stopifnot(c('end','end.1','start','start.1')%in% colnames(dt))
  stopifnot( ! 'ov_width' %in% colnames(dt) )
  dt[[colname]] <- pmin(dt$end,dt$end.1) - pmax(dt$start,dt$start.1) 
  dt
}

add_rangeid <- function(gr,colname){
  mcols(gr)[[colname]]=paste0(seqnames(gr),'_',start(gr),'_',end(gr))
  gr
}
DT2GR = function(dt,seqinf=si,checksi=TRUE){

  if(is(dt,'GenomicRanges')) {
    warning('already a GRanges Object')
    return(dt)
  }


  stopifnot(c('seqnames','start')%in%colnames(dt))
  stopifnot(c('width')%in%colnames(dt)|c('end')%in%colnames(dt))
  if(checksi){stopifnot(dt[['seqnames']] %in% seqlevels(seqinf))
  }else{seqinf=NULL}
  
  hasend=FALSE
  haswidth=FALSE

  if('end' %in% colnames(dt) ){
    stopifnot (dt[['end']] %>% `>`(0) %>%all)
    hasend=TRUE
  }
  if('width' %in% colnames(dt) ){
    stopifnot (dt[['width']] %>% `>`(0) %>%all)
    haswidth=TRUE
  }
  
  stopifnot(dt[['start']] %>% is.numeric)
  stopifnot(hasend|haswidth )
  
  if(haswidth & ! hasend ){
    dt[['end']]  = dt[['start']]+dt[['width']]-1 
  } 
  if(hasend ){

  } 

  #strand
  if(! 'strand' %in% colnames(dt)){
    dt[['strand']] ='*'
  }

  stopifnot(dt[['strand']] %in% c('+','-','*'))
  



  mdatcols = colnames(dt) %>% setdiff(c('seqnames','start','width','strand','end')) 
  #create basic granges object
  if(checksi){
    gr=GRanges(dt[['seqnames']],IRanges(dt[['start']],dt[['end']]),strand=dt[['strand']],seqinfo=seqinf)
  }else{    gr=GRanges(dt[['seqnames']],IRanges(dt[['start']],dt[['end']]),strand=dt[['strand']])}

  #add in the metadata if we need to
  if(length(mdatcols)){
    if(is.data.table(dt)){ mcols(gr) = dt[,mdatcols,with=FALSE]%>%as("DataFrame")
    }else{ mcols(gr) = dt[,mdatcols]%>%as("DataFrame")}
  }

    stopifnot(all(colnames(dt) %in% c(colnames(mcols(gr)),'seqnames','start','end','width' ,'strand')))

  gr
}


nameFilt<-function(x,pat,neg=F,...){
  nms = names(x)
  inm = grep(pat,nms,...)
  if (neg) inm = inm * -1
  x[inm]
}

gr2srl <- function(gr,weight='score',...){
  list(
    pos =  gr%>%subset(strand=='+')%>%coverage(weight=weight,...),
    neg =  gr%>%subset(strand=='-')%>%coverage(weight=weight,...)
  )
}

trimRanges<-function(gr,filt=FALSE){
  #expand ranges within limits of the chrs
  start(gr)%<>%pmax(1)
  chrends = seqlengths(gr)[seqnames(gr)%>%as.character]
  if(filt==TRUE) {
    gr=gr[end(gr)<chrends]
    gr=gr[start(gr)>1]
  }
  if(!any(is.na(seqlengths(gr)))) end(gr)%<>%pmin(chrends)
  gr
}

tidyGR<-function(gr,sinfo=si){
   #add chr if needed
   if(!grepl('chr',seqlevels(gr))[[1]]){
     seqlevels(gr)=paste0('chr',seqlevels(gr))
   }
   #get the chromosomes in the object we want to keep
   chrs2keep = intersect(seqlevels(gr),si@seqnames)
   #now drop the chromosomes we don't want (M or U)
   gr=keepSeqlevels(gr,chrs2keep)
   #modify the seqlevels slot
   seqlevels(gr)=si@seqnames
   #modify the seqinfo slot
   seqinfo(gr)=si
   #and finally trim the ranges so they're within the limits now defined for the chromosomes
   gr = gr [ start(gr) <= seqlengths(gr) [as.character(seqnames(gr)) ] ]
   gr=trimRanges(gr)
   return(gr)
 }
nomcols <- function(x){ x %>% .[,NULL] }



####This kludgy function let's me source things with sublime text 
#When the source code is on my laptop and I"m working on the server
#it takes in an anddress that looks like i'ts on the mac, rsyncs
#the file from the mac to the server, then sources that.
source <- function(file,...){
  require(stringr)
  if(file %>% str_detect('(/Volumes/)|(Users/harnett/)')){
    message('syncing local files to server')
    #and change file to a path on the server
    file=str_replace(file,'Volumes','g')
    file=str_replace(file,'^.*/src','./src')
    #
  }
  if(any(grepl('#\\s*SPINFILE',readLines(file,5)))){
    require(quietly=TRUE,knitr)
    message('spinning file')
    # look for a 'reports' folder and put the report in there
    oldwd = getwd()
    newwd = file.path(oldwd,'analysis/Reports')
    stopifnot(file.exists(newwd))
    on.exit(setwd(oldwd))
    setwd(newwd)
    # spin it!
    outfile = knitr::spin(file,FALSE)
    knit2html(outfile)
    # rmarkdown::render(outfile)
    # system("/usr/lib/rstudio-server/bin/pandoc/pandoc +RTS -K512m -RTS ggvis_knit_issue.md --to html --from markdown+autolink_bare_uris+ascii_identifiers+tex_math_single_backslash-implicit_figures --output ggvis_knit_issue.html --smart --email-obfuscation none --self-contained --standalone --section-divs --template /g/software/linux/pack/r-3.2.2/centos-6/lib64/R/library/rmarkdown/rmd/h/default.html --variable 'theme:bootstrap' --mathjax --variable 'mathjax-url:https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --no-highlight --variable highlightjs=/g/software/linux/pack/r-3.2.2/centos-6/lib64/R/library/rmarkdown/rmd/h/highlight ")
  }
  else{
    message(paste0('sourcing file: ',file))
    base::source(file,...)
  }
}


is.error <- function(x) inherits(x, "try-error")


updateDisplay<-function(n=NA){
  try(rep(100,{dev.off()}),silent=TRUE)

  if(is.na(n)) localhost=readLines('~/display.txt')[[1]]
  if(!is.na(n)) localhost = paste0('localhost:',n,'.0')

  X11.options(display=localhost)

  system(paste0("export DISPLAY=",localhost))

  message(paste0('display set to ',localhost))
}
setRecover <- function(){ options(error=recover)}
recoverOnce <- function(){ options(error=function(x){setNoRecover();recover()})}
setNoRecover <- function(){ options(error=NULL)}



stackdfs = function(dflist,varname='Lvar'){
  for(i in seq_along(dflist)){ 
    stopifnot(! varname %in% colnames(dflist[[i]]))
  }
  
  stopifnot(!is.null(names(dflist)))
  
  out=
    dflist%>%
    lapply(as.data.frame,stringsAsFactors=FALSE)%>%
    mapply(.,XXXvarname=names(.),FUN=cbind,SIMPLIFY=F)%>%

  rbindlist(fill=TRUE,use.names=TRUE)

  setnames(out,'XXXvarname',varname)

  out
}

#safe left join, can specifi
sleft_join<-function(dt1,dt2,by,allow.missing=FALSE,allow.duplicates=FALSE,...){
  stopifnot(!'tmpid' %in% colnames(dt1))
  # stopifnot(!'tmpid' %in% colnames(dt2))
  dt1[['tmpid']] =  seq_len(nrow(dt1))
  # dt2[['tmpid']] =  seq_len(nrow(dt2))
  newcols = setminus(colnames(dt2),colnames(dt1))
  #there must be new cols
  stopifnot(length(newcols)!=0)
  #
  out=left_join(dt1,dt2,by=by,...)
  #
  if(!allow.missing) {
    if(nrow(anti_join(dt1%>%select(one_of(by)),dt2%>%select(one_of(by)),by=by))!=0){
     stop('Some elements missing from dt2')
   }
  }
  if(!allow.duplicates){
    if( any(duplicated(out$tmpid))){
      stop('Some elements from dt1 duplicated in dt2')
    }
  }
  #  
  out =out %>% select(-starts_with('tmpid'))
  dt1 =dt1 %>% select(-starts_with('tmpid'))
  #
  out
}

#
setminus<-function(x,y){
  setdiff(x,intersect(x,y))
}

#' Select Strings, like vgrep
#' 
#' Intro Sentence.
#' 
#' Paragraph explaining more
#' 
#' @param str string vector to filter
#' @param pattern the pattern (can be a regex(pat) pattern)
#' @param inv invert selection?
#' @examples str_select(letters[1:10],'a|d')
#' @return subset of str
#' @export 

str_select <- function(str,pattern,invert=FALSE){
  if(invert == TRUE){
    lfunc = not
  }else if (invert == FALSE){
   lfunc = identity
  }else{
    stop('invert should be TRUE or false')
  }
  str[ lfunc( str_detect(str, pattern) ) ]  
}


w  <- function(str,wstring='[^\n ]+'){
  stopifnot(is.character(str))
  out=str_extract_all(str,wstring)
  if(length(out)==1 ) out = out[[1]]
  out
}

#' reduce
#' 
#' reduce a gr, duplicate the rows, add in the old metadata.
#' 
#' @param gr a granges object
#' @examples reduce_with_mcols(gr)
#' 
#' @return a granges object 
#' @export 

reduce_with_mcols=function(gr,...){
  #create reduced gr with reverse mapping
  gr_red =  reduce(gr,...,with.revmap=TRUE)
  #expand it out
  expanded_gr2 <- rep(gr_red, elementLengths(gr_red$revmap))
  #add in the mcols of the appopriate rows from the original entry
  mcols(expanded_gr2) %<>% cbind( mcols(gr)[unlist(gr_red$revmap),])
  #we now have the reduced gr, with each reduced entry dduplicated
  #once for every range in it, and the mcols attached.
  expanded_gr2
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


list2tab = function(l,namecol='val'){
  stopifnot(is.list(l))
  stopifnot(!is.null(names(l)))
  lens = vapply(l,length,666)
  val=unlist(l,rec=FALSE,use.names=FALSE)
  if(length(val)==0) return(NULL)
  totlen = length(val)
  val = unlist(val)
  stopifnot(length(val)==totlen)
  stopifnot(is.vector(val))
  out = data.table(name=rep(names(l),times=lens),val = val)
  
  #set the colnames with input colname
  newcolnames = colnames(out)
  newcolnames[[2]] = namecol
  setnames(out,newcolnames)

  out
}



#this function puts alphabetic labels at the corner of multipanel figures.
labelplots<-function(rownum,colnum,labsize=30,xoffset=0.02,yoffset=0.02,labels=toupper(letters)){
  stopifnot( ! (rownum)%%1 , length(rownum)==1)
  stopifnot( ! (colnum)%%1 , length(colnum)==1)
  #positions
  xpos = seq(from=0,by = 1/colnum,length.out = colnum) %>% add(xoffset)
  ypos = seq(from=1,by = - 1/rownum,length.out = rownum) %>% subtract(yoffset)
  #labels
  #label counter
  k=0
  for(i in 1:rownum){
    for(j in 1:colnum){
      k = k+1
      lab = labels[k]
      grid.text(lab,gp=gpar(fontsize=labsize),y =  unit(ypos[i],"npc") ,x = unit(xpos[j]  ,"npc"))
    }
  }
}
# labelplots(2,3)

lognz<-function(x){x[x!=0]<-log10(x[x!=0]);x }



#' countOverlapsBam
#' 
#' count a Grange objects' overlaps with a bam file.
#' 
#' 
#' @param gr - the GRanges object
#' @param bam - the bam file
#' @examples counts <- countOverlapsBam(genes_gr,cagebam,chromsizes)
#' @return a vector of overlap counts
#' @export 
countOverlapsBam <- function(gr,bamfile){
  stopifnot(gr %>% is("GenomicRanges"))
  stopifnot(gr %>% length %>% gt(0))
  stopifnot(file.exists(bamfile))
  stopifnot(file.exists(bamfile))
  stopifnot(str_detect(bamfile,'\\.bam$'))
  stopifnot(system("which bedtools",ignore.stderr=TRUE)==0)  #temp dir to be deleted on exit
  #create a temporary directory for the operation
  #name bed and write it to temp folder
  gr$name <- seq_along(gr)
  gr <- sort(gr)
  gr$score <- 666
  #write the ranges as a bed file
  bedfilename <- basename(tempfile(fileext='.bed'))
  on.exit(file.remove(bedfilename))
  export(gr[,c('name','score')],bedfilename)
  #now run bedtools count overlaps
  cmd  <- sprintf("intersectBed -a %s -b %s -c | cut -f 7",
    bedfilename,bamfile)
  ovnums <- system(cmd,intern=TRUE)
  ovnums %<>% as.numeric
  #return int the correct order
  ovnums[order(gr$name)]
}

getStrandedExprMat <- function(gr,SRLs,viewfunc=viewSums,parfunc=function(...){mclapply(mc.cores=10,...)}){
  #function that takes in some SimpleRleLists (e.g. RNAseq data)
  #which have a 'positive' and 'negative' sublist structure
  #and then gets uses the GRViews function to get the strand specific
  #score for each of the input regions in each of the SRLs 
  stopifnot( is (gr,'GRanges'))
  if(is(SRLs,'SimpleRleList')) SRLs = list(SRLs)
  stopifnot( is (SRLs,'list'))
    origorder = match(gr,sort(gr))
    gr = sort(gr)
    gstrand = strand(gr)%>%as.character
    SRLs[[1]]->SRL
    out=parfunc(SRLs,function(SRL){#use mclapply to iterate over the SRLSs
      if( identical(names(SRL)[1:2],c('pos','neg')) ) {#if we have pos and neg run seperately on them
        cat('.')
        # stopifnot( seqlevels(gr) %in% SRL[[1]]%>%names)
        v=rep(0,length(gr))
        ispos  =     gstrand=='+'#Splitting our granges by strand
        isneg  =     gstrand=='-'
        isstar =     gstrand=='*'
        expect_true(all( ispos | isneg | isstar))
        if(any(isneg)){  v[ isneg ] = unlist(viewfunc(GRViews(SRL$neg,gr[ isneg ]))) }
        if(any(ispos)){  v[ ispos ] = unlist(viewfunc(GRViews(SRL$pos,gr[ ispos ] ))) }
        if(any(isstar)){ v[ isstar ] = unlist(viewfunc(GRViews(SRL$pos,gr[ isstar ] ))) +
        unlist(viewfunc(GRViews(SRL$neg,gr[ isstar ] ))) }
        expect_is(v,c('numeric','vector'))
        v[origorder]
      }else{#else just run GRViews on it
        # stopifnot( seqlevels(gr) %in% SRL%>%names)
        unlist(viewfunc(GRViews(SRL,gr%>%tidyGR)))[origorder]
      }
  })
  out%<>%simplify2array
  out
}





# bam2coverage ------------------------------------------------------------
#Read in a bam file and return the coverage as and RleList object - no normalization
bam2coverage<-function(bam,doshift=F,doresize=F,stranded=F,fragment.length=NA,use.spp=F,
  is.FirstMateRead=NA,is.ProperPair=NA,use.size=F,maxsize=NA,maxfraglength=300,resizeSide='start'){
  library(Rsamtools)
  cat(paste0('attempting to read bam file  ',bam,'  \n'))
  #scan the bam file
  reads <- scanBam(
    file = bam,
    param = ScanBamParam(
      what = c('pos', 'rname', 'strand', 'qwidth','isize'),
      flag = scanBamFlag(isUnmappedQuery = FALSE,isFirstMateRead=is.FirstMateRead,isProperPair=is.ProperPair)
    )
  )[[1]]

  message('..bam read..')
  #fix the names - necessary for some bams,by adding 'chr'
  if(!length(reads[['rname']][grepl(reads[['rname']],pattern='^2R$')])==0 ){
    reads[['rname']] <- Rle(paste0('chr',reads[['rname']]))
  }


  #and fixing the name of the mitochondrial genome
   if(!length(reads[['rname']][grepl(reads[['rname']],pattern='itochond')])==0 ){
     reads[['rname']][grepl(reads[['rname']],pattern='itochond')]<-  'chrM'
   }

  #make sure all the seqnames are okay
  if(!all(unique(reads[['rname']]) %in% names(seqinfo(Dmelanogaster)))){
    message(paste0('these chrs arent right',unique(reads[['rname']])))
    stop()
  }
  # stopifnot(length(unique(reads[['qwidth']]))==1)
  #don't keep chrM or chrU
  reads.keep<-as.logical(reads[['rname']] %in% chrs.keep)
  #reads.keep<-as.logical(reads[['rname']] %in% 'chr2R')

  # keep isize forlater
  isize=reads$isize[reads.keep]

  #now as a GRange object
  reads<-GRanges(reads[['rname']][reads.keep],
             IRanges(reads[['pos']][reads.keep],
                     width=reads[['qwidth']][reads.keep]),
             strand=Rle(as.character(reads[['strand']][reads.keep])),
             seqinfo=si
  )


  if(!is.na(maxsize)){
    reads=reads[isize<maxsize]
    isize=isize[isize<maxsize]
  }
  #coverageplot(peaks.pos$chr10[wpeaks[1]], peaks.neg$chr10[wpeaks[1]])
  
  
  #this calculates the cross correlation betwen strands to estimate the fragment length
  #fragment.length<-estimate.mean.fraglen(g[seqnames(g)=='chr2R'],method='correlation')
  
  #we'll just us ccf, and only on chr2L to save time (won't differ for other chrs)
  if((doshift || doresize)&(is.na(fragment.length))){
    message(paste0('...calculating fragment length....'))
    fragment.length<-which.max(
      ccf(plot=F,
                                 as.vector(coverage(reads[strand(reads)=='+' & seqnames(reads) == 'chr2R'])[['chr2R']])[1000000:5000000],
                                 as.vector(coverage(reads[strand(reads)=='-' & seqnames(reads) == 'chr2R'])[['chr2R']])[1000000:5000000],
                                 lag.max=maxfraglength)$acf
      )
    fragment.length<-abs(fragment.length)
    message(paste0('.fragment length ',fragment.length))
    
    stopifnot(fragment.length>50)
    
  }
  
  ispos=strand(reads)=='+'

  #shift reads
  if(doshift==T){
    ifpos=strand(reads)=='+'
    #move the negative ones towards the start, stopping at one
    reads[!ispos]<-shift(reads[!ifpos],-pmin(round(fragment.length*0.5),start(reads[!ifpos])-1))
    reads[ispos]<-shift(reads[ifpos],pmin(round(fragment.length*0.5),seqlengths(reads[ifpos])[as.character(seqnames(reads[ifpos]))]-end(reads[ifpos])))
  }

  if(use.size){
    doresize=T
    fragment.length=abs(isize)
  }

  #do resize
  if(doshift==F & doresize==T){
    message(paste0('..reads length resized... '))
    reads<-resize(reads,fragment.length,resizeSide)
  } 

  #and return the coverage
  if(stranded==F){return(coverage(reads))
  }else{
    return(#or the stranded coverage
      list(pos=coverage(subset(reads,strand=='+')),
         neg=coverage(subset(reads,strand=='-'))
      )
    )
  }  
}

#' Load a single saved object
#' 
#' Import a single object saved with 'save' and return its value.
#' 
#' @param f a file
#' @examples e.g. - an example
#' this 
#' @return returns 
#' @export 

load_obj <- function(f)
{
    env <- new.env()
    nm <- load(f, env)[1]
    env[[nm]]
}



srl2gr = function(srl){
  if(is(srl,'SimpleRleList')) {
      srl%>%as("GRanges")%>%{seqinfo(.)=si;.}%>%subset(score!=0)%>%setstrand('*')
  }else if (is(srl$pos,'SimpleRleList') && is(srl$neg,'SimpleRleList') ){
    c(
      srl$pos%>%as("GRanges")%>%{seqinfo(.)=si;.}%>%subset(score!=0)%>%setstrand('+'),
      srl$neg%>%as("GRanges")%>%{seqinfo(.)=si;.}%>%subset(score!=0)%>%setstrand('-')
    )
  }else {
    stop('bad input, need a simple rle list!')
  }
}

# cov=  G(10,score=10,seqinfo=si)%>%coverage(weight='score')
# GRViews(cov,G(1))

gr2srl <- function(gr,weight='score',...){
  list(
    pos =  gr%>%subset(strand=='+')%>%coverage(weight=weight,...),
    neg =  gr%>%subset(strand=='-')%>%coverage(weight=weight,...)
  )
}

setstrand <- function(gr,chr){strand(gr)=chr;gr}

flipstrand<-function(gr){
  s=strand(gr)%>%as.character
  strand(gr) = ifelse(s=='-','+','-')
  gr
}


#' get_system_name
#' 
#' gets the name of the machine your on.
#' 
#' @return returns name of machine
#' @export 

get_system_name <- function(){
  system(intern=TRUE,'uname -a')
}
# get_system_name()


#' gets the name of the sshng machine from a server
#' 
#' This function is just a convience function when working ont he file server.
#' 
#' @param username - 
#' @examples e.g. - an example
#' plots2Laptop() 
#' @return Nothing - this acts on files 
#' @export 

get_laptop_name <- function(){
  username=system('echo $USER',int=TRUE)
  #get the name of the laptop
    # laptopname = system('who | grep "mac-"',int=TRUE)
    laptopname = system('who --lookup',int=TRUE)
    if(! length(laptopname)>0) stop('no machines found')
    laptopname =base::grep(laptopname,value=TRUE,pattern=username)
    if(! length(laptopname)>0) stop('no machines matching username')
    laptopname =  str_extract(laptopname,perl('(?<=\\d\\d:\\d\\d\\s\\().*(?=.)'))
    laptopname =base::grep(laptopname,value=TRUE,pattern='pts/',invert=TRUE)

    if(length(laptopname)>1){
      for(i in seq_along(laptopname)){
        isvalidIP = !system(paste0('ssh -q ' ,username,'@',laptopname[i], ' exit' ))
        if(isvalidIP) break
      }
      stopifnot(isvalidIP)
      laptopname=laptopname[i]
    }
    laptopname#shoudl select the valid one if multiple exist
#   
}

#' Copies files to my laptop
#' 
#' This function is just a convience function when working ont he file server.
#' 
#' @param tfile - The plot to copy 
#' @param dest - the location the files will go to 
#' @examples e.g. - an example
#' plots2Laptop() 
#' @return Nothing - this acts on files 
#' @export 



#transfer plots to laptop
file2laptop = function(  
  tfile ,
  folder = '~/Desktop/'
){

  stopifnot(file.exists(tfile))
  localmachine=get_laptop_name()
  myMachineFolder = paste0(localmachine,':',folder)
    
  ifNotTransferred= system(paste0('scp -r ' , tfile , ' ' , myMachineFolder))
  if(ifNotTransferred) warning(paste0('Transfer to ',localmachine, ' failed'))
  invisible(NULL)

}



#' group_slice
#' 
#' get N groups.
#' #' 
#' @param dt - a grouped tbl
#' @param v - a vector for slicing
#' @examples e.g. - an example
#' this 
#' @return returns 
#' @export 

group_slice<-function(dt,v){
  stopifnot('data.frame'%in%class(dt))
  stopifnot(v>0,(v%%1)==0)
  if(n_groups(dt)==1) warning('Only 1 group to slice')
  #get the grouping variables
  groups=
    unique(dt%>%select())%>%
    group_by()%>%
    slice(v)
  out=inner_join(dt,groups)
  out
}

# source('/g/furlong/project/28_B_DNASE/analysis/evolutionary_analyses/INSIGHT/src/dermot/Rprofile.R')