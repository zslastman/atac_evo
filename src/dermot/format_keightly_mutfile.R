inputfilename <- "data/keightly_et_al_2009_mutdata.txt"
formatfilename <- "data/keightly_et_al_2009_mutdata_formatted.txt"
mutdata <- fread(inputfilename)



mutdata = mutdata[1,,]

ischrom  <- . %>% is_in(c('2L','2R','Y','4','3L','3R','X'))

ischrominds <- ischrom(mutdata)

linelengths  <- ischrominds %>% cumsum %>% Rle %>% .@lengths


#function which, given a vector of integers and a single large vector,
#splits the large vector into chunks of the sizes specified. Can also do rles
chunkvect <- function(x, chunks, as.rle=F) {
  stopifnot(sum(chunks)==length(x))
  chunks=cumsum(chunks)
  chunkstarts=c(1,(chunks[-length(chunks)]+1))
  if(as.rle==T){func=Rle
  }else{func=identity}
  out=lapply(seq_along(chunks),function(n){
    func(x[chunkstarts[n]:chunks[n]])
  })
  names(out)=names(chunks)
  return(out)
}
#split the data by it's line number
lines <- chunkvect(mutdata %>% as.character,linelengths)
#figure out which lines are 8 entries and which are 9
isshortinds <- lines %>% vapply(length,666) %>% eq(8)
#no add in th emissing col so all lines are the same
lines[isshortinds] %<>% lapply(. %>% append(.,'n'))
#now conver the thing into a matrix and make it a data frame
linedf <- lines %>% simplify2array %>% t %>% as.data.frame(log=FALSE,header=TRUE)
#fix the automatic assignment of T to TRUE by as.data.frame
for(ind in 3:6){levels(linedf[[ind]]) %<>% str_replace('TRUE','T') }
#and now remove the first line and adda  character
colnames(linedf) = linedf[1,,drop=TRUE] %>% unlist %>% as.character
linedf %<>% .[-1,]
#finally write it
linedf %>% write.table(quote=FALSE,row.names=FALSE,col.names=TRUE,formatfilename)
