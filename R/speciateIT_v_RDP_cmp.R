##
## Comparing CV results
##

setwd("/Users/pgajer/devel/speciateIT/data")

file <- "speciateIT_v_RDP_misassigned_taxa.csv"
x <- read.csv(file)

head(x)

idx <- x$tool=="RDP" & x$region=="V4"
rdp.v4 <- x[idx,]
rdp.v4 <- rdp.v4[,-c(3,4)]

idx <- x$tool=="speciateIT" & x$region=="V4"
sIT.v4 <- x[idx,]
sIT.v4 <- sIT.v4[,-c(3,4)]

rdp.sIT.v4 <- c()
for ( id in unique(x$correct_lineage_sp) )
{
    idx <- rdp.v4$correct_lineage_sp==id
    r <- sum(rdp.v4[idx,"nSeqs_misassigned"])
    idx <- sIT.v4$correct_lineage_sp==id
    s <- sum(sIT.v4[idx,"nSeqs_misassigned"])
    rdp.sIT.v4 <- rbind(rdp.sIT.v4, c(id, r, s))
}
rdp.sIT.v4 <- as.data.frame(rdp.sIT.v4)
rownames(rdp.sIT.v4) <- rdp.sIT.v4[,1]
rdp.sIT.v4 <- rdp.sIT.v4[,-1]
colnames(rdp.sIT.v4) <- c("RDP","sIT")
rdp.sIT.v4[,1] <- as.numeric(rdp.sIT.v4[,1])
rdp.sIT.v4[,2] <- as.numeric(rdp.sIT.v4[,2])


plot(log(rdp.sIT.v4[,1]), log(rdp.sIT.v4[,2]))


##
##
##
file <- "speciateIT_v_RDP_10foldCV_scores.csv"
y <- read.csv(file)

idx <- y$region=="V4" & y$tool=="RDP"
rdp.v4 <- y[idx,]

idx <- y$region=="V3V4" & y$tool=="RDP"
rdp.v3v4 <- y[idx,]

idx <- y$region=="V1V3" & y$tool=="RDP"
rdp.v1v3 <- y[idx,]


idx <- y$region=="V4" & y$tool=="speciateIT"
sIT.v4 <- y[idx,]

idx <- y$region=="V3V4" & y$tool=="speciateIT"
sIT.v3v4 <- y[idx,]

idx <- y$region=="V1V3" & y$tool=="speciateIT"
sIT.v1v3 <- y[idx,]



spp <- unique(y$correct_lineage_sp)
cmp.v4 <- matrix(nrow=length(spp), ncol=5)
colnames(cmp.v4) <- c("n", "ncRDP", "pcRDP", "ncsIT", "pcsIT")
rownames(cmp.v4) <- spp
cmp.v4 <- as.data.frame(cmp.v4)
for ( sp in spp )
{
    idx <- rdp.v4$correct_lineage_sp==sp
    nr <- sum(idx)
    nr.cor <- sum(rdp.v4$score_cor[idx])
    pr.cor <- 100*nr.cor / nr
    idx <- sIT.v4$correct_lineage_sp==sp
    ns <- sum(idx)
    ns.cor <- sum(sIT.v4$score_cor[idx])
    ps.cor <- 100*ns.cor / ns
    if ( nr != ns ) {
        stop("nr != ns")
    }
    cmp.v4[sp,] <- c(nr, nr.cor, pr.cor, ns.cor, ps.cor)
}



plot(cmp.v4$pcRDP, cmp.v4$pcsIT)

plot(cmp.v4$pcRDP-cmp.v4$pcsIT)

myHist(cmp.v4$pcRDP)

myHist(cmp.v4$pcsIT)

myHist(log2(cmp.v4$pcsIT/cmp.v4$pcRDP))
mode.1D(log2(cmp.v4$pcsIT/cmp.v4$pcRDP)) #  -0.01106072
median(log2(cmp.v4$pcsIT/cmp.v4$pcRDP), na.rm=TRUE) #  -0.05658353

100*sum(cmp.v4$ncRDP)/sum(cmp.v4$n) # 75.50902
100*sum(cmp.v4$ncsIT)/sum(cmp.v4$n) # 73.83403


##
## V1V3
##

spp <- unique(y$correct_lineage_sp)
cmp.v1v3 <- matrix(nrow=length(spp), ncol=5)
colnames(cmp.v1v3) <- c("n", "ncRDP", "pcRDP", "ncsIT", "pcsIT")
rownames(cmp.v1v3) <- spp
cmp.v1v3 <- as.data.frame(cmp.v1v3)
for ( sp in spp )
{
    idx <- rdp.v1v3$correct_lineage_sp==sp
    nr <- sum(idx)
    nr.cor <- sum(rdp.v1v3$score_cor[idx])
    pr.cor <- 100*nr.cor / nr
    idx <- sIT.v1v3$correct_lineage_sp==sp
    ns <- sum(idx)
    ns.cor <- sum(sIT.v1v3$score_cor[idx])
    ps.cor <- 100*ns.cor / ns
    if ( nr != ns ) {
        stop("nr != ns")
    }
    cmp.v1v3[sp,] <- c(nr, nr.cor, pr.cor, ns.cor, ps.cor)
}

plot(cmp.v1v3$pcRDP, cmp.v1v3$pcsIT)
plot(cmp.v1v3$pcRDP-cmp.v1v3$pcsIT)

myHist(cmp.v1v3$pcRDP)

myHist(cmp.v1v3$pcsIT)

myHist(log2(cmp.v1v3$pcsIT/cmp.v1v3$pcRDP))
mode.1D(log2(cmp.v1v3$pcsIT/cmp.v1v3$pcRDP)) # 0.004612865
median(log2(cmp.v1v3$pcsIT/cmp.v1v3$pcRDP), na.rm=TRUE) # 0

100*sum(cmp.v1v3$ncRDP)/sum(cmp.v1v3$n) # 90.2355
100*sum(cmp.v1v3$ncsIT)/sum(cmp.v1v3$n) # 87.82411


##
## V3V4
##

spp <- unique(y$correct_lineage_sp)
cmp.v3v4 <- matrix(nrow=length(spp), ncol=5)
colnames(cmp.v3v4) <- c("n", "ncRDP", "pcRDP", "ncsIT", "pcsIT")
rownames(cmp.v3v4) <- spp
cmp.v3v4 <- as.data.frame(cmp.v3v4)
for ( sp in spp )
{
    idx <- rdp.v3v4$correct_lineage_sp==sp
    nr <- sum(idx)
    nr.cor <- sum(rdp.v3v4$score_cor[idx])
    pr.cor <- 100*nr.cor / nr
    idx <- sIT.v3v4$correct_lineage_sp==sp
    ns <- sum(idx)
    ns.cor <- sum(sIT.v3v4$score_cor[idx])
    ps.cor <- 100*ns.cor / ns
    if ( nr != ns ) {
        stop("nr != ns")
    }
    cmp.v3v4[sp,] <- c(nr, nr.cor, pr.cor, ns.cor, ps.cor)
}

plot(cmp.v3v4$pcRDP, cmp.v3v4$pcsIT)
plot(cmp.v3v4$pcRDP-cmp.v3v4$pcsIT)

myHist(cmp.v3v4$pcRDP)
myHist(cmp.v3v4$pcsIT)

myHist(log2(cmp.v3v4$pcsIT/cmp.v3v4$pcRDP))
mode.1D(log2(cmp.v3v4$pcsIT/cmp.v3v4$pcRDP)) #  -0.002846234
median(log2(cmp.v3v4$pcsIT/cmp.v3v4$pcRDP), na.rm=TRUE) # -0.06453129

100*sum(cmp.v3v4$ncRDP)/sum(cmp.v3v4$n) # 84.93189
100*sum(cmp.v3v4$ncsIT)/sum(cmp.v3v4$n) # 79.09445
