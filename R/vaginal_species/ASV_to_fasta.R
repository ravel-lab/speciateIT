##
## Extracting the header (ASV sequences) from a ASV file and crreating the correponding FASTA file
##


asvs.dir <- "~/projects/ASV_files/data/"

##
## CONTRA
##
file <- paste0(asvs.dir,"CONTRA_all_runs_dada2_abundance_table.rda")
##save(pt, file=file)
load(file)
dim(pt) # 4607 20053

ASVs <- colnames(pt)

file <- paste0(asvs.dir,"CONTRA_ASVs.fa")
cmd <- paste0("rm -f ",file)
system(cmd)

asv.i <- 1
for ( asv in ASVs )
{
    s <- paste0(">ASV",asv.i,"\n")
    asv.i <- asv.i + 1
    cat(s, file=file, append=TRUE)
    cat(asv, file=file, append=TRUE)
    cat("\n", file=file, append=TRUE)
}
file
## "~/projects/ASV_files/data/CONTRA_ASVs.fa"

## https://stackoverflow.com/questions/61374573/how-can-i-eliminate-duplicated-sequences-in-fasta-file

## https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#CDHITEST

## pushd ~/projects/ASV_files/data/
## cd-hit-est -i CONTRA_ASVs.fa -o CONTRA_ASVs_nr.fa -c 1.00 -n 10 -d 0 -M 16000 -T 8

## % grep '>' CONTRA_ASVs.fa | wc -l
##    20053

## % grep '>' CONTRA_ASVs_nr.fa | wc -l
##    18546



## Rebuilding V3V4 models
