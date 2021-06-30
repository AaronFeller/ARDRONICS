## ARDRA-Analyzer ##############################################################
##
## Analyzes ARDRA base-pair lengths from samples and compares the results
## to a library of known samples.
## Input: CSV files of sample and library
## Output: CSV files of comparisons and graphical representation of UPGMA
## Current Version: 1.6 automation for sample file only

## Old Versions:
## version 1.5 new binning with less weight error
## version 1.4 first successful iteration
## Version 1.3 jaccard_UPGMA_formatting
## Version 1.2 rework_for_distance_matrix
## Version 1.1 looping_work
## Version 1.0 beginning_work

## Required libraries:
#install.packages("ade4")
#install.packages("ape")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DECIPHER")
#install.packages("phylogram")


Ardronics <- function () {

## Libraries:
library(ade4)
library(ape)
library(DECIPHER)
library(phylogram)

## SETTINGS: ###################################################################
{
work_dir <- readline(prompt = "Enter the location of data file (home/R): ")
samp_file <- readline(prompt = "Enter the name of the data file: ")
## Plus/Minus Variable
plus_minus <- readline(prompt = "Enter the variance (standard is 20): ")
## Minimum Length
min_bp <- readline(prompt = "Enter the minimum band size (in bp): ")
## Number of Enzyme Cuts
enz_count <- readline(prompt = "Enter the number of enzymes: ")

## Set working Directory #######################################################
setwd(work_dir)

## Uploading Library Data ######################################################
# NOT USING FOR RUNNING CODE ONLY ON SAMPLES
# lib_data = read.csv(lib_file, header=TRUE)

## Uploading Sample Data #######################################################
sample_data = read.csv(samp_file, header=TRUE)
# FOR TESTING sample_data = "NULL"

#For only sample data below, combining samples and library in next line
merge_data <- sample_data

## Combine two files
# merge_data <- merge.data.frame(lib_data, 
#                                sample_data, 
#                                all.x = T, 
#                                all.y = T, 
#                                no.dups = T)

## Change data point under minimum length to NA ################################
merge_data[merge_data <= min_bp] <- NA
}


###FULL START###
{

## Initialize Library ##########################################################
{
jaccard_data <- matrix(data = NA, 
                       nrow = nrow(merge_data)+2, 
                       ncol = 1500*enz_count +1)
"Bands" -> jaccard_data[1,1]
"Enzyme" -> jaccard_data[2,1]
}
## Bring in organism name from merged data #####################################
for (i in 1:nrow(merge_data)) {
  toString(merge_data[i, 1]) -> jaccard_data[i+2,1]
}

## Bring in band lengths/enzymes for Jaccard index #############################

## TEST THIS ## merge_data[2] <- tolower(merge_data[2])0
ux <- unique(merge_data[2])
toString(ux[2,1])


#MAKE THIS GET THE ENZYME ^^^

#first 1-1500
for (i in 1:nrow(ux)){
print(toString(ux[i,1]))}
  
for (k in 1:1500) {
  jaccard_data[1, k+1] <- k
  jaccard_data[2, k+1] <- toString(ux[1,1])
}

for (l in 1:1500) {
  jaccard_data[1, l+1501] <- l
  jaccard_data[2, l+1501] <- toString(ux[2,1])
}

# THIS IS WHERE THE EDITS NEED TO BE
# for (k in 1:nrow(merge_data)) {
#   for (l in 3:ncol(merge_data))  {
#     if (!is.na(merge_data[k,l]))    { 
#             jaccard_data[1, l-1+(k-1)*(ncol(merge_data)-2)] <- 
#         toString(merge_data[k, l])
#       jaccard_data[2, l-1+(k-1)*(ncol(merge_data)-2)] <- 
#         toString(merge_data[k, 2])
#     }
#   }
# }

## Clean-Up ####################################################################
{
## Delete NA Columns
jaccard_data <- jaccard_data[, ! apply(jaccard_data, 2, 
                                       function(x)all(is.na(x)))]
## Delete Duplicate Rows
jaccard_data <- unique(jaccard_data[])
## Delete Duplicate Data Points
jaccard_data <- unique(jaccard_data, MARGIN = 2)
## Fill with 0s
jaccard_data[is.na(jaccard_data)] <- 0
}
## Compare band lengths and create binary ######################################
for (i in 1:nrow(merge_data)) {
  for (j in 3:ncol(merge_data))  { 
    if (is.na(merge_data[i, j])){next}
    for (k in 2:ncol(jaccard_data))    {
      
    ### Extra If function for band size inserted here! ###

      if (tolower(toString(merge_data[i , 2])) == tolower(toString(jaccard_data[2 , k]))
        &&
        as.numeric(merge_data[i, j]) <= (as.numeric(jaccard_data[1 , k])
                                         +plus_minus)
        &&
        as.numeric(merge_data[i, j]) >= (as.numeric(jaccard_data[1 , k])
                                         -plus_minus)) { 
        jaccard_data[trunc(3+((i-1)/enz_count)), k] <- 1}
    }
  }
}

## Create File for Jaccard Analysis ############################################
{
jaccard_matrix <- jaccard_data
jaccard_matrix <- jaccard_matrix[-1,]                         #removes a row
jaccard_matrix <- jaccard_matrix[-1,]                         #removes a row
rownames(jaccard_matrix) <- jaccard_matrix[,1]                #names the rows
jaccard_matrix <- jaccard_matrix[,-1]                         #removes a column
jaccard_matrix <- as.data.frame.matrix(jaccard_matrix)        #change to df
jaccard_matrix[] <- lapply(jaccard_matrix, function(x){if(is.factor(x)) 
  as.numeric(as.character(x)) else x})                        #converts to numbs                        
}
## Jaccard Analysis from R package 'ade4' ######################################
{
# Method 1 is Jaccard index (1901) S3 coefficient of Gower & Legendre
m <- jaccard_matrix
m1 <- jaccard_matrix
# For Methods 1-10 see readme for package ade4  
d <- dist.binary(m1, 
                 method = 1, 
                 diag = FALSE, 
                 upper = FALSE)  
}

## Apply hierarchical clustering ###############################################
{
hc <- hclust(dist(d))

## Plot the Dendrogram #########################################################
jaccard_result <- plot(hc, 
### ADD TITLE HERE                       main = 
                       labels=m$ID, 
                       hang = 0, 
                       cex = .55,)
#jaccard_hcd <- plot(hcd, labels=m$ID, cex = .55, type = "triangle")
}

###FULL STOP###
}
  
#end program
  
}

## Export CSV of the Distance Matrix ###########################################
hcd = as.dendrogram(hc)
dist_matrix_as_df <- as.data.frame(as.matrix(d))
write.csv(dist_matrix_as_df,"/home/aaron/R/ARDRA-Analyzer/dissm1.csv", 
          row.names = TRUE)

dissm1 <- dist_matrix_as_df

# ## THIS IS IT! ###############################################################
hcd = as.dendrogram(hc)
WriteDendrogram(hcd, file = "testingDendroFileSave")
tets1 <- ReadDendrogram("testingDendroFileSave")
y <- as.phylo(tets1)


hins_dend <- as.dendrogram(hins_UPGMA)
WriteDendrogram(hins_dend, file = "testing2_20210602")
tets2 <- ReadDendrogram("testing2_20210602")
z <- as.phylo(tets2)
 
 
all.equal.phylo(y, z)
comparePhylo(y, z, plot = T)

yz <- rbind(y, z)
CADM.global(yz, 2, 115)

Y <- cophenetic(y)
Z <- cophenetic(z)

rownames(Y) <- sub(".", "", rownames(Y))
colnames(Y) <- sub(".", "", colnames(Y))


Zfix <- rownames(Y)
Z <- Z[Zfix, Zfix]
YZ <- rbind(Y, Z)
CADM.global(YZ, 2, 115, nperm=99, make.sym=TRUE, weights=NULL)
# IT WORKS


# 
# ## CADM.global() ## This might be the comparison! ############################
# 
# 
# a <- rtree(5)
# b <- rtree(5)
# A <- cophenetic(a)
# B <- cophenetic(b)
# x <- rownames(A)
# B <- B[x, x]
# M <- rbind(A, B)
# CADM.global(M, 2, 5)
# 
