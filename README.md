### Returns dataframe containing information from .ab1 files from Sanger sequencing.

#### Description

ab1Files2SequenceDF is a function which reads .ab1 files from Sanger Sequencing into data frame containing columns for file path, sample info, primer info, consensus sequence and a primary and secondary basecall sequence produced from the Bioconductor package sanderseqR's makeBaseCalls function.

#### Usage
ab1Files2SequenceDF(file.directory= "...", basecall.ratio = 0.33, file.type = ".ab1", recursive.dir = TRUE)

#### Arguments:

file.directory 
The directory containing the sequencing files

basecall.ratio = 0.33 
A variable passed to the makeBaseCalls function. see ?makeBaseCalls for details.

file.type = ".ab1" 
File extention of sequencing files.

recursive.dir = TRUE 
Are sequencing files stored in subdirectory.

### Examples:
ab1Files2SequenceDF(file.directory="sangerSeq/")

# if files are in subdirectories 
ab1Files2SequenceDF(file.directory="sangerSeq/", recursive.dir = TRUE)

#### Author(s)

CJ "Jake" Saunders
