

file <- "/Users/pgajer/devel/speciateIT/data/V3V4_vaginal_Mar2021/V3V4_vag_mcDir/error_thlds.txt"
x <- read.table(file)
nrow(x) # 1450
head(x)
##                  V1        V2
## 1             Taxon Threshold
## 2  c_Actinobacteria  0.989005
## 3   o_Bacteroidales  0.989005
## 4       g_Chlamydia  0.989005
## 5      p_Firmicutes  0.989005
## 6 o_Fusobacteriales  0.989005
