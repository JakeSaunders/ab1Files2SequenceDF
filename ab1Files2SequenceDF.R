
# load packages
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

# clean data from sequence files names
file.name <- list.files("seqfiles/",recursive = T)
file.name <- file.name[grep(pattern = ".ab1",x = file.name)]

file.info <- sub(".ab1","",x = file.name)
files.info <- as.data.frame(t(as.data.frame(strsplit(x = file.info, "_"))))

files.info <- cbind.data.frame(
  files.info,
  colsplit(files.info$V3,pattern = "-",names = c("sample.name","sample.number")),
  colsplit(files.info$V4,pattern = "-",names = c("primer","direction")),
  file.name
)
rownames(files.info) <- NULL
files.info$direction <- sub("for|1-for|3-for",replacement = "F",
                        sub("rev|r-rev|4-rev",replacement = "R",files.info$direction)
  )
rm(file.info,file.name)
write.csv(files.info, file="SangerSeqRunInfo.csv")

# select degen worm pcr sanger seq files
files.info <- files.info[43:82,6:10]

##### maybe start here
sanger.file1 <- paste0("seqfiles/",files.info$file.name)
sanger.file2 <- paste0("seqfiles/",files.info$file.name[2])

extractSequences <- function(sanger.file){
  # make seqs objects
  sanger.seq <- readsangerseq(sanger.file)
  basecall.seq <- makeBaseCalls(sanger.seq)
  # put seqs into data frame row
  row.out <- cbind.data.frame(
    as.character(sanger.file),
    as.character(sanger.seq@primarySeq),
    as.character(basecall.seq@primarySeq),
    as.character(basecall.seq@secondarySeq),
    stringsAsFactors = F
  )
  names(row.out) <- c(
    "source.file",
    "seq.consensus",
    "seq.basecall.primary",
    "seq.basecall.secondary"
  )
  print(row.out)
}

sanger.files <- paste0("seqfiles/",files.info$file.name)
files.seq <- ldply(sanger.files,extractSequences)

##### note #####
# AlignSeqs
# AlignTranslation
# DECIPHER package
