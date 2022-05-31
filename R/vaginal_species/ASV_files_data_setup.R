##
## Processing ASV files
##

##in ~/.Rprofile
##library(DT)
##dt.show <- function(x){ datatable(x, filter = 'top', options = list(pageLength=100, autoWidth = TRUE))}

library(data.table)
## data.table 1.14.3 IN DEVELOPMENT built 2022-05-29 12:17:01 UTC using 8 threads (see ?getDTthreads).  Latest news: r-datatable.com

get.pt <- function(file, has.tx.row=TRUE)
{
    ct <- fread(file, fill=TRUE)
    print(dim(ct))

    ct <- as.data.frame(ct)

    if ( has.tx.row )
    {
        ct <- ct[-2,] ## removing taxon row
    }

    rownames(ct) <- ct[,1]
    ct <- ct[,-1]

    colnames(ct) <- ct[1,]
    ct <- ct[-1,]

    for ( i in seq(ncol(ct)) )
    {
        ct[,i] <- as.numeric(ct[,i])
        ct[is.na(ct[,i]),i] <- 0
    }

    ct <- as.matrix(ct)
    rs <- rowSums(ct, na.rm=TRUE)
    ##pt <- ct/rs

    list(rs=rs,
         ct=ct)
}


asvs.dir <- "~/projects/ASV_files/data/"


##
## CONTRA
##
file <- paste0(asvs.dir,"CONTRA_all_runs_dada2_abundance_table_w_taxa.csv")
##x <- read.csv(file, row.names=1)
r <- get.pt(file)

myHist(log10(r$rs))

rs <- log10(r$rs)

dim(ct)
ct <- r$ct[rs>=4,]
dim(ct)
ct[1:2,1:2]

pt <- ct / rowSums(ct)

pt[1:2,1:2]

sum(pt[1,])

file <- paste0(asvs.dir,"CONTRA_all_runs_dada2_abundance_table.rda")
save(pt, file=file)


##
## DARE
##
## file <- paste0(asvs.dir,"DARE_all_runs_dada2_abundance_table.csv")
## r <- get.pt(file, has.tx.row=FALSE)

## file <- paste0(asvs.dir,"DARE_all_runs_dada2_abundance_table.csv")
## save(r$pt, file=file)


##
## Douching
##
file <- paste0(asvs.dir,"Douching_all_runs_dada2_abundance_table_w_taxa.csv")
r <- get.pt(file)

rs <- log10(r$rs)
myHist(rs)

thld <- 4
sum(rs>=thld)

dim(r$ct)
ct <- r$ct[rs>=thld,]
dim(ct)

ct[1:2,1:2]

pt <- ct / rowSums(ct)

pt[1:2,1:2]

sum(pt[1,])


file <- paste0(asvs.dir,"Douching_all_runs_dada2_abundance_table.rda")
save(pt, file=file)


##
## GALE
##
file <- paste0(asvs.dir,"GALE_all_runs_dada2_abundance_table_w_taxa.csv")
r <- get.pt(file)

rs <- log10(r$rs)
myHist(rs)

thld <- 4
sum(rs>=thld)

dim(r$ct)
ct <- r$ct[rs>=thld,]
dim(ct)

ct[1:2,1:2]

pt <- ct / rowSums(ct)

pt[1:2,1:2]

sum(pt[1,])


file <- paste0(asvs.dir,"GALE_all_runs_dada2_abundance_table.rda")
save(pt, file=file)


##
## HMP
##
file <- paste0(asvs.dir,"HMP_all_runs_dada2_abundance_table_w_taxa.csv")
r <- get.pt(file)
## dim: 4888 15754

rs <- log10(r$rs)
myHist(rs)

thld <- 4
sum(rs>=thld)

dim(r$ct)
ct <- r$ct[rs>=thld,]
dim(ct)
## 3715 15753


ct[1:2,1:2]

pt <- ct / rowSums(ct)

pt[1:2,1:2]

sum(pt[1,])


file <- paste0(asvs.dir,"HMP_all_runs_dada2_abundance_table.rda")
save(pt, file=file)



##
## HMP ???
##
file <- paste0(asvs.dir,"hmp_asv_merge_s2.csv")
## r <- get.pt(file)
## dim: 18568  6122 ???



##
## LSVF
##
file <- paste0(asvs.dir,"LSVF_all_runs_dada2_abundance_table_w_taxa.csv")
r <- get.pt(file)
## dim: 2285 7702

rs <- log10(r$rs)
myHist(rs)

thld <- 4
sum(rs>=thld)

dim(r$ct)
ct <- r$ct[rs>=thld,]
dim(ct)

ct[1:2,1:2]

pt <- ct / rowSums(ct)

pt[1:2,1:2]

sum(pt[1,])


file <- paste0(asvs.dir,"LSVF_all_runs_dada2_abundance_table.rda")
save(pt, file=file)


##
## CHARM
##
file <- paste0(asvs.dir,"MCHRM_all_runs_dada2_abundance_table_w_taxa.csv")
r <- get.pt(file)
## dim: 1098 2596
## Error in `.rowNamesDF<-`(x, value = value) :
##   duplicate 'row.names' are not allowed
## In addition: Warning message:
## non-unique value when setting 'row.names': ‘CHARM004.visit1.V2.215737’

##
## V400
##
file <- paste0(asvs.dir,"V400_all_runs_dada2_abundance_table_w_taxa.csv")
r <- get.pt(file)
## dim: 433 2832

rs <- log10(r$rs)
myHist(rs)

thld <- 4
sum(rs>=thld)

dim(r$ct)
ct <- r$ct[rs>=thld,]
dim(ct)

ct[1:2,1:2]

pt <- ct / rowSums(ct)

pt[1:2,1:2]

sum(pt[1,])


file <- paste0(asvs.dir,"V400_all_runs_dada2_abundance_table.rda")
save(pt, file=file)
