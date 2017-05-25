ab1Files2SequenceDF <- function(file.directory, basecall.ratio = 0.33, file.type = ".ab1"){
  
  # check for and install packages if missing----
  if( !require(plyr) ) {install.packages("plyr")}
  if( !require(sangerseqR) ) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("sangerseqR")
  }
  library(plyr)
  library(sangerseqR)
  
  # make directory and file variables
  directory <- file.directory
  
  # clean file info from sequence files names----
  generateSeqFilesInfo <- function(directory, file.type = ".ab1"){
    file.path <- list.files(directory, recursive = T, full.names = T)
    file.path <- file.path[grep(pattern = file.type,x = file.path)]
    
    file.info <- t(as.data.frame(
      strsplit(x = sub(".ab1","",x = file.path),"_"),
      stringsAsFactors = F
    ))
    file.info <- cbind.data.frame(file.path,file.info[,c(3:4)],stringsAsFactors = F)
    
    rownames(file.info) <- NULL
    colnames(file.info) <- c("file.path","sample.name.full","primer.name.full")
    file.info
  }
  
  # extract sequences given .ab1 file ----
  extractSequences <- function(sanger.file, basecall.ratio = 0.33){
    # make seqs objects
    sanger.seq <- readsangerseq(sanger.file)
    basecall.seq <- makeBaseCalls(sanger.seq,ratio = basecall.ratio)
    # put seqs into data frame row
    row.out <- cbind.data.frame(
      as.character(sanger.file),
      as.character(sanger.seq@primarySeq),
      as.character(basecall.seq@primarySeq),
      as.character(basecall.seq@secondarySeq),
      stringsAsFactors = F
    )
    names(row.out) <- c(
      "file.path",
      "seq.consensus",
      "seq.basecall.primary",
      "seq.basecall.secondary"
    )
    print(row.out)
  }
  
  # make dfs and put together ----
  file.info <- generateSeqFilesInfo(directory)  # remember file extention
  seqs <- ldply(
    file.info$file.path, 
    extractSequences, 
    basecall.ratio = 0.33
  )
  df <- join(file.info,seqs,by = "file.path")
  
  
  # output data frame ----
  df
}