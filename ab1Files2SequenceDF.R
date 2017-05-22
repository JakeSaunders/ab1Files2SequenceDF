
# load packages----
source("https://bioconductor.org/biocLite.R")
biocLite("sangerseqR")
install.packages("tidyr")
install.packages("reshape2")
biocLite("DECIPHER")

# check libraries
library(sangerseqR)
library(tidyr)
library(reshape2)
library(DECIPHER)
library(plyr)

# set wd
setwd("C:/Users/saundecj/Google Drive/@R/Sanger")

# clean file info from sequence files names ----

library(plyr)

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

library(sangerseqR)

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

library(plyr)

file.info <- generateSeqFilesInfo("Sanger/seqfiles/")  # remember file extention

seqs <- ldply(
  file.info$file.path, 
  extractSequences, 
  basecall.ratio = 0.33
)

test <- join(file.info,seqs,by = "file.path")

save(test,file = "sequencedataframe.Rdata")
write.csv(test,file = "sangerseqorder2017-05.csv")

#### TODO: Just need to put everything in a wrapper and run on files

###ideas for later
# lenght in basepairs in data frame

#################################################################
#### save this for later, use when aligning like samples ########
# files.info <- cbind.data.frame(
#   files.info,
#   colsplit(files.info$V3,pattern = "-",names = c("sample.name","sample.number")),
#   colsplit(files.info$V4,pattern = "-",names = c("primer","direction")),
#   file.name
# )
#
#
# files.info$direction <- sub("for|1-for|3-for",replacement = "F",
#                             sub("rev|r-rev|4-rev",replacement = "R",files.info$direction)
# )
#################################################################

##### note #####
# AlignSeqs
# AlignTranslation
# DECIPHER package
