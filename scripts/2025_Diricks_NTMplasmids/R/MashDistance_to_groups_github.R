library(polysat)
#Input (distance matrix, e.g. output from mash) 
#Note, make sure there is no "-" in your file names, because it will convert it in the row header to . but not in column header; "_" is OK
#NOTE! Mash distance table can contain e_05 or e_06,... scientific notations, these are not treated correctly in this script! Change them beforehand to 0 or change all cells to "number"!
#Note, it turns down 0.054 to 0.05 and then they are included in group! Sometimes does not fit with phylogeny, always check!
#Note: If mash distance is 0, mashtree sometimes outputs an empty cell in distance matrix!! Check this!!
File <- "Mash_distance"
matrix <- read.table (File, sep="\t", header=TRUE)
View(matrix)
rn <- matrix[,1]
matrix_g <- matrix [,-1]
rownames(matrix_g) <- rn
#matrix_g <- as.numeric(matrix_g)
matrix_numeric <- as.matrix(matrix_g)
View(matrix_numeric)

groupings_0.05 <- assignClones(matrix_numeric, threshold=0.05)
#Output (csv files with names of samples included in matrix and their grouping)
write.table(groupings_0.05,"Grouping_0.05.csv",sep=",")
