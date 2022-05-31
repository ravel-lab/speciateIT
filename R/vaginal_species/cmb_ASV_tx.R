##
## Taxonomy of ASV's from different 16S rRNA projects
##

## This is an attempt to generate a list of species detected in vaginal environment.
## The issue is that many of these tx classifications are misclassifications

##
## Combining *_ASVs_sIT/taxonRank.txt files
##

## % classify -v -i ~/projects/ASV_files/data/CONTRA_ASVs_nr.fa -d sIT_models_V3V4 -o CONTRA_ASVs_sIT
## % classify -v -i ~/projects/ASV_files/data/DOUCHING_ASVs_nr.fa -d sIT_models_V3V4 -o DOUCHING_ASVs_sIT
## % classify -v -i ~/projects/ASV_files/data/GALE_ASVs_nr.fa -d sIT_models_V3V4 -o GALE_ASVs_sIT
## % classify -v -i ~/projects/ASV_files/data/HMP_ASVs_nr.fa -d sIT_models_V3V4 -o HMP_ASVs_sIT
## % classify -v -i ~/projects/ASV_files/data/LSVF_ASVs_nr.fa -d sIT_models_V3V4 -o LSVF_ASVs_sIT
## % classify -v -i ~/projects/ASV_files/data/V400_ASVs_nr.fa -d sIT_models_V3V4 -o V400_ASVs_sIT

## % cat CONTRA_ASVs_sIT/taxonRank.txt DOUCHING_ASVs_sIT/taxonRank.txt GALE_ASVs_sIT/taxonRank.txt HMP_ASVs_sIT/taxonRank.txt LSVF_ASVs_sIT/taxonRank.txt V400_ASVs_sIT/taxonRank.txt > ASVs_taxonRank.txt



##
## Extracting unique taxons from ASVs_taxonRank.txt
##

# % head ASVs_taxonRank.txt
##
## Desulfurispira_natronophila	1109	5.97972608648765
## Chrysiogenes_arsenatis	789	4.25428663862828
## Ruminococcus_gnavus	576	3.10579100614688
## Dictyoglomus_thermophilum	568	3.06265501995039

file <- "~/devel/speciateIT/data/ASVs_taxonRank.txt"
txs <- read.table(file, header=FALSE)[,1]

txs <- unique(txs)
length(txs) # 2686

cbind(sort(txs))
##
##    [1,] "Abiotrophia_sp_Abiotrophia_defectiva"
##    [2,] "Abyssivirga_alkaniphila"
##    [3,] "Acetanaerobacterium_sp"
##    [4,] "Acetitomaculum_ruminis"
##    [5,] "Acetivibrio_ethanolgignens"
##    [6,] "Acetivibrio_sp"
##    [7,] "Acetobacter_indonesiensis_Acetobacter_cibinongensis_Acetobacter_thailandicus"
##    [8,] "Acetobacter_orientalis"
##    [9,] "Acetomicrobium_faecale"


sp <- "Acetobacter_indonesiensis_Acetobacter_cibinongensis_Acetobacter_thailandicus"
s <- strsplit(sp, "_")[[1]]
length(s)

uq.spp <- c()
for ( sp in txs )
{
    s <- strsplit(sp, "_")[[1]]
    if ( length(s) <= 2 ){
        uq.spp <- c(uq.spp, sp)
    } else {
        ii <- seq(1, length(s)-1, by=2)
        for ( i in ii )
        {
            sp <- paste0(s[i],"_",s[i+1])
            uq.spp <- c(uq.spp, sp)
        }
    }
}

uq.spp <- unique(uq.spp)
length(uq.spp) # 3518


file <- "/Users/pgajer/devel/speciateIT/post_merge_reference_files/V3V4_spp_new.lineage"
reference_sequences/V3V4_trimmed_noEuks_nr_Complete.fa
post_merge_reference_files/V3V4_spp_new.tx -o sIT_models_V3V4
