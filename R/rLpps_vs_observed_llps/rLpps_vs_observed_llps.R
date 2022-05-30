##
## Comparing the distribution of random lpp's with observed lpp's
##

## src/est_error_thlds.cc was compiled with OUTPUT_LPPS_TO_FILE flag set to 1
## which resulted in log10(pp)'s of the given node ref seq's and its siblings
## being written to a file for all species

## in ~/devel/speciateIT/data/V3V4_vaginal_Mar2021

## % est_error_thlds -v -d V3V4_vag_mcDir -c 0.9
## lpp's produced by est_error_thlds were saved in V3V4_vag_mcDir/Lpps


## % sp_model_seq_lpps -v -d V3V4_vag_mcDir -s 1000
## lpp's of random model sequences produced by sp_model_seq_lpps were saved in V3V4_vag_mcDir/rLpps


file1 <- "~/devel/speciateIT/data/V3V4_vaginal_Mar2021/V3V4_vag_mcDir/Lpps/BVAB1.csv"
lpps <- read.csv(file1, header=FALSE)
colnames(lpps) <- c("seqID","tx","lpp")
head(lpps)

(tx.freq <- table(lpps$tx))
## BVAB1    g_Blautia      g_Dorea g_Howardella   g_Moryella
##     7         1820          355           12           22

(taxa <- names(tx.freq))
## "BVAB1"        "g_Blautia"    "g_Dorea"      "g_Howardella" "g_Moryella"



file2 <- "~/devel/speciateIT/data/V3V4_vaginal_Mar2021/V3V4_vag_mcDir/rLpps/BVAB1.csv"
rlpps <- read.csv(file2, header=FALSE, row.names=1)

setequal(rownames(rlpps),taxa)

(ref.tx <- rownames(rlpps)[1])

r <- plot.dlpp(lpps, ref.tx)

rr <- plot.drlpp(rlpps, ref.tx, lpps)

plot(rr$dlpp[[ref.tx]], las=1)
idx <- lpps$tx==ref.tx
x <- lpps$lpp[idx]
points(x, rep(0, length(x)), pch="|", cex=0.75, col='red')


pics.dir <- "~/devel/speciateIT/R/rLpps_vs_observed_llps/pics/"

file <- paste0(pics.dir, ref.tx,".pdf")
pdf(file, width=12, height=6)
op <- par(mfrow=c(1,2), mar=c(4, 3.75, 1.75, 0.75), mgp=c(2.25,0.5,0),tcl = -0.3)
r <- plot.dlpp(lpps, ref.tx)
rr <- plot.drlpp(rlpps, ref.tx, lpps)
par(op)
dev.off()
system(paste0("open ",file))
file


##
## Generating such plots for all species
##

dir <- "~/devel/speciateIT/data/V3V4_vaginal_Mar2021/V3V4_vag_mcDir/Lpps"
files <- dir(dir, pattern =".csv")

spp <- gsub(".csv","", files)

pdf.files <- c()
for ( sp in spp )
{
    cat("\r",sp)
    file1 <- paste0("~/devel/speciateIT/data/V3V4_vaginal_Mar2021/V3V4_vag_mcDir/Lpps/",sp,".csv")
    lpps <- read.csv(file1, header=FALSE)
    colnames(lpps) <- c("seqID","tx","lpp")
    ##
    file2 <- paste0("~/devel/speciateIT/data/V3V4_vaginal_Mar2021/V3V4_vag_mcDir/rLpps/",sp,".csv")
    if ( !file.exists(file2) )
    {
        next
    }
    rlpps <- read.csv(file2, header=FALSE, row.names=1)
    ##ref.tx <- rownames(rlpps)[1]
    ref.tx <- sp
    file <- paste0(pics.dir, ref.tx,".pdf")
    pdf(file, width=12, height=6)
    op <- par(mfrow=c(1,2), mar=c(4, 3.75, 1.75, 0.75), mgp=c(2.25,0.5,0),tcl = -0.3)
    r <- plot.dlpp(lpps, ref.tx)
    rr <- plot.drlpp(rlpps, ref.tx, lpps)
    par(op)
    dev.off()
    ##system(paste0("open ",file))
    pdf.files <- c(pdf.files, file)
}
