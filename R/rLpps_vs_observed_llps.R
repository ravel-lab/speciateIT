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

plot.dlpp(lpps, ref.tx)


plot.dlpp <- function(lpps, ref.tx)
{
    tx.freq <- table(lpps$tx)
    taxa <- names(tx.freq)
    dlpp <- list()
    x <- c()
    y <- c()
    for ( tx in taxa )
    {
        idx <- lpps$tx==tx
        lpp <- lpps$lpp[idx]
        dlpp[[tx]] <- density(lpp)
        x <- c(x, dlpp[[tx]]$x)
        y <- c(y, dlpp[[tx]]$y)
    }

    x.lim <- range(x)
    y.lim <- range(y)

    plot(dlpp[[tx]]$x, dlpp[[tx]]$y, type='n', las=1, ylab="Density", xlab="log10(pp)",
         xlim=x.lim, ylim=y.lim, main=ref.tx)
    i <- 1
    for ( tx in taxa )
    {
        lines(dlpp[[tx]]$x, dlpp[[tx]]$y, col=i)
        i <- i + 1
    }
    legend("topleft", legend=taxa, lty=1, col=seq(taxa), inset=0.05)

    invisible(list(dlpp=dlpp, tx.freq=tx.freq))
}


plot.drlpp(rlpps, ref.tx)

plot.drlpp <- function(rlpps, ref.tx)
{
    taxa <- rownames(rlpps)
    dlpp <- list()
    x <- c()
    y <- c()
    for ( tx in taxa )
    {
        lpp <- as.numeric(rlpps[tx,])
        dlpp[[tx]] <- density(lpp)
        x <- c(x, dlpp[[tx]]$x)
        y <- c(y, dlpp[[tx]]$y)
    }

    x.lim <- range(x)
    y.lim <- range(y)

    plot(dlpp[[tx]]$x, dlpp[[tx]]$y, type='n', las=1, ylab="Density", xlab="log10(pp)",
         xlim=x.lim, ylim=y.lim, main=ref.tx)
    i <- 1
    for ( tx in taxa )
    {
        lines(dlpp[[tx]]$x, dlpp[[tx]]$y, col=i)
        i <- i + 1
    }
    legend("topleft", legend=taxa, lty=1, col=seq(taxa), inset=0.05)

    invisible(list(dlpp=dlpp, tx.freq=tx.freq))
}
