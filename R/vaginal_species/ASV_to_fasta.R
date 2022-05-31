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


## clstr_to_clstr2.pl -i CONTRA_ASVs_nr.fasta.clstr -o CONTRA_ASVs_nr.clstr2


##
## DOUCHING
##
file <- paste0(asvs.dir,"Douching_all_runs_dada2_abundance_table.rda")
##save(pt, file=file)
load(file)
dim(pt) #  973 21695

ASVs <- colnames(pt)

file <- paste0(asvs.dir,"DOUCHING_ASVs.fa")
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
## "~/projects/ASV_files/data/DOUCHING_ASVs.fa"

## pushd ~/projects/ASV_files/data/
## cd-hit-est -i DOUCHING_ASVs.fa -o DOUCHING_ASVs_nr.fa -c 1.00 -n 10 -d 0 -M 16000 -T 8

## % grep '>' DOUCHING_ASVs.fa | wc -l
## 21695

## % grep '>' DOUCHING_ASVs_nr.fa | wc -l
## 16071



##
## GALE
##
file <- paste0(asvs.dir,"GALE_all_runs_dada2_abundance_table.rda")
##save(pt, file=file)
load(file)
dim(pt) #  1275 4617

ASVs <- colnames(pt)

file <- paste0(asvs.dir,"GALE_ASVs.fa")
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
## "~/projects/ASV_files/data/GALE_ASVs.fa"

## pushd ~/projects/ASV_files/data/
## cd-hit-est -i GALE_ASVs.fa -o GALE_ASVs_nr.fa -c 1.00 -n 10 -d 0 -M 16000 -T 8
## 4616  finished       4201  clusters


##
## HMP
##
file <- paste0(asvs.dir,"HMP_all_runs_dada2_abundance_table.rda")
##save(pt, file=file)
load(file)
dim(pt) #   3715 15753

ASVs <- colnames(pt)

file <- paste0(asvs.dir,"HMP_ASVs.fa")
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
## "~/projects/ASV_files/data/HMP_ASVs.fa"

## pushd ~/projects/ASV_files/data/
## cd-hit-est -i HMP_ASVs.fa -o HMP_ASVs_nr.fa -c 1.00 -n 10 -d 0 -M 16000 -T 8
## 15752  finished      14228  clusters



##
## LSVF
##
file <- paste0(asvs.dir,"LSVF_all_runs_dada2_abundance_table.rda")
##save(pt, file=file)
load(file)
dim(pt) # 1819 7701

ASVs <- colnames(pt)

file <- paste0(asvs.dir,"LSVF_ASVs.fa")
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
## "~/projects/ASV_files/data/LSVF_ASVs.fa"

## pushd ~/projects/ASV_files/data/
## cd-hit-est -i LSVF_ASVs.fa -o LSVF_ASVs_nr.fa -c 1.00 -n 10 -d 0 -M 16000 -T 8
## 7700  finished       7141  clusters



##
## V400
##
file <- paste0(asvs.dir,"V400_all_runs_dada2_abundance_table.rda")
##save(pt, file=file)
load(file)
dim(pt) # 337 2831

ASVs <- colnames(pt)

file <- paste0(asvs.dir,"V400_ASVs.fa")
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
## "~/projects/ASV_files/data/V400_ASVs.fa"

## pushd ~/projects/ASV_files/data/
## cd-hit-est -i V400_ASVs.fa -o V400_ASVs_nr.fa -c 1.00 -n 10 -d 0 -M 16000 -T 8
## 2830  finished       2708  clusters
