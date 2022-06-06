##
## Generating an example of a CSV file holding a numerical data for readTable test
##

nr <- 3
nc <- 4
tbl <- matrix(runif(nr*nc), nrow=nr, ncol=nc)
rownames(tbl) <- paste0("r",seq(nr))

file <- "~/devel/speciateIT/data/test_data/sample_3x4_num_tbl_with_rownames_and_no_header.csv"
write.table(tbl, file, quote=FALSE, col.names=FALSE, sep=",")
