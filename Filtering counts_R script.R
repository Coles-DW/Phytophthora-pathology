setwd("")
install_course_packages.R
install.packages("usethis")
expdata <- read.csv("Unique counts_all samples_05102020.csv")
str(expdata)

mat <- as.matrix(expdata[, -1])
rownames(mat) <- expdata[, 1]
str(mat)

# zeros any cells (gene in a sample) where read count is less than a threshold, here 10
# then convert cells above that threshold to '1', resulting in a binary matrix
mat[mat < 10] <- 0
mat[mat > 0] <- 1

# the number of genes that have expression above the threshold in at least X replicates
sum(rowSums(mat) >= 1)
sum(rowSums(mat) >= 3)

# remove any genes that do not meet the threshold in at least X replicates
expdata.fil <- expdata[rowSums(mat) >= 1, ]

write.csv(expdata.fil, 'filtered_unique_counts.csv', row.names=F)
write.table(expdata.fil, 'filteredcounts1.csv', row.names=F, sep=';')