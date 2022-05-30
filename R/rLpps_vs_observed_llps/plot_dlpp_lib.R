plot.dlpp <- function(lpps, ref.tx)
{
    tx.freq <- table(lpps$tx)
    taxa <- names(tx.freq)
    ext.taxa <- c()
    for ( tx in taxa )
    {
        ext.taxa <- c(ext.taxa, paste0(tx," (",tx.freq[tx][[1]],")"))
    }

    dlpp <- list()
    x <- c()
    y <- c()
    for ( tx in taxa )
    {
        idx <- lpps$tx==tx
        if ( sum(idx) > 2 )
        {
            lpp <- lpps$lpp[idx]
            dlpp[[tx]] <- density(lpp)
            x <- c(x, dlpp[[tx]]$x)
            y <- c(y, dlpp[[tx]]$y)
        }
    }

    x.lim <- range(x)
    y.lim <- range(y)
    tx <- names(dlpp)[1]

    plot(dlpp[[tx]]$x, dlpp[[tx]]$y, type='n', las=1, ylab="Density", xlab="log10(pp)",
         xlim=x.lim, ylim=y.lim, main=ref.tx)
    i <- 1
    for ( tx in names(dlpp) )
    {
        lines(dlpp[[tx]]$x, dlpp[[tx]]$y, col=i)
        i <- i + 1
    }
    legend("topleft", legend=ext.taxa, lty=1, col=seq(taxa), inset=0.05)

    invisible(list(dlpp=dlpp, tx.freq=tx.freq))
}

plot.drlpp <- function(rlpps, ref.tx, lpps=NULL, pt.cex=0.75, pt.pch="|")
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

    if ( !is.null(lpps) )
    {
        tx.freq <- table(lpps$tx)
        ext.taxa <- c()
        for ( tx in taxa )
        {
            ext.taxa <- c(ext.taxa, paste0(tx," (",tx.freq[tx][[1]],")"))
        }
    }

    plot(dlpp[[tx]]$x, dlpp[[tx]]$y, type='n', las=1, ylab="Density", xlab="log10(pp)",
         xlim=x.lim, ylim=y.lim, main=ref.tx)
    i <- 1
    for ( tx in taxa )
    {
        lines(dlpp[[tx]]$x, dlpp[[tx]]$y, col=i)
        if ( !is.null(lpps) )
        {
            idx <- lpps$tx==tx
            x <- lpps$lpp[idx]
            points(x, rep(0, length(x)), pch=pt.pch, cex=pt.cex, col=i)
        }
        i <- i + 1
    }

    if ( !is.null(lpps) )
    {
        legend("topleft", legend=ext.taxa, lty=1, col=seq(taxa), inset=0.05)
    } else {
        legend("topleft", legend=taxa, lty=1, col=seq(taxa), inset=0.05)
    }

    invisible(list(dlpp=dlpp, tx.freq=tx.freq))
}
