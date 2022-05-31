##
## Generating 3d embeddings of ASVs to identify ASV abundance profile clusters
##

library(uwot)     # umap()
library(rgl)

asvs.dir <- "~/projects/ASV_files/data/"

##
## CONTRA
##
file <- paste0(asvs.dir,"CONTRA_all_runs_dada2_abundance_table.rda")
##save(pt, file=file)
load(file)
dim(pt) # 4607 20053

pt.3d <- umap(t(pt),
             n_components = 3,
             metric="cos",
             min_dist = 0,
             n_neighbors = 15,
             verbose = T)

rownames(pt.3d) <- colnames(pt)

file <- paste0(asvs.dir,"CONTRA_all_runs_dada2_pt_3d.rda")
save(pt.3d, file=file)

plot3d(pt.3d)

##
## Douching
##
file <- paste0(asvs.dir,"Douching_all_runs_dada2_abundance_table.rda")
## save(pt, file=file)
load(file)
dim(pt) #  973 21695

pt.3d <- umap(t(pt),
             n_components = 3,
             metric="cos",
             min_dist = 0,
             n_neighbors = 15,
             verbose = T)

rownames(pt.3d) <- colnames(pt)

file <- paste0(asvs.dir,"Douching_all_runs_dada2_pt_3d.rda")
save(pt.3d, file=file)

plot3d(pt.3d)

myHist(pt.3d[pt.3d[,2]< 50,2])

idx <- pt.3d[,3]> -50 & pt.3d[,2]<50
pt <- pt[,idx]

idx <- pt.3d[,3]<50 & pt.3d[,1]<50 & pt.3d[,2]> -50
plot3d(pt.3d[idx,])
pt <- pt[,idx]


idx <- pt.3d[,3]> -50 & pt.3d[,1]<50
plot3d(pt.3d[idx,])

pt <- pt[,idx]
