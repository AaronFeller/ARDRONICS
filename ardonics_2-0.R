################################################################################
## ARDRA-Analyzer ##############################################################
## SETTINGS: ONLY change this part! ############################################
## When you are ready to run the analysis, press ctrl + shift + enter and wait.
################################################################################

# Set your working directory here:
work_dir = "/home/aaron/R/ARDRA-Analyzer"

# Enter your sample file here:
samp_file = "20211109_Digest_Plate1_All.csv"

# Enter your library file here:
lib_file = "20210505_lib_clip_HaeIII.csv"

## Enter the +/- bin size (10 will be for bands on either side making it a size 20 bin)
plus_minus <- 20

## Enter the minimum cutoff length for bands in base pairs
min_bp <- 100

## Enter the number of restriction enzymes used in the analysis (This setup is for single digests only)
enz_count <- 1


################################################################################
## ARDRA-Analyzer ##############################################################
## Analyzes ARDRA base-pair lengths from samples and compares the results
## to a library of known samples.
## Input: CSV files of sample and library
## Output: CSV files of comparisons and graphical representation of UPGMA
## Current Version: 1.5 Updated binning code for reduced weight error

## Required libraries:
# install.packages("ade4")
# install.packages("ape")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DECIPHER")
# install.packages("phylogram")
# install.packages("tibble")
# install.packages("dplyr")



## WORKING CODE BELOW ## DO NOT CHANGE #########################################

## Load Libraries:
library(ade4)
library(ape)
library(DECIPHER)
library(phylogram)
library(tibble)
library(dplyr)

################################################################################
## Set working Directory #######################################################
{
  setwd(work_dir)
  ## Uploading Library Data ######################################################
  lib_data = read.csv(lib_file, header=TRUE)
  ## Uploading Sample Data #######################################################
  sample_data = read.csv(samp_file, header=TRUE)
  # FOR TESTING sample_data = "NULL"
  ## Combine two files
  
  for (i in 1:(ncol(lib_data)-ncol(sample_data))) {
    sample_data <- sample_data %>%
      add_column(new_col = NA)
  }
  
  colnames(sample_data) <- colnames(lib_data)
  merge_data <- rbind(lib_data, sample_data) 
  
  
  # merge_data <- rbind(lib_data,
  #   sample_data,
  #   all.x = T,
  #   all.y = T,
  #   no.dups = T)
  
  ## Change data point under minimum length to NA ################################
  merge_data[merge_data <= min_bp] <- NA
}
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
  ux <- unique(merge_data[2])
  #toString(ux[2,1])
  #first 1-1500
  #for (i in 1:nrow(ux)){
  #  print(toString(ux[i,1]))}

for (enznum in 1:length(ux)){
  for (k in 1:1500) {
    jaccard_data[1, (k+1)+(1500*(enznum-1))] <- k
    jaccard_data[2, (k+1)+(1500*(enznum-1))] <- toString(ux[enznum,1])
  }
}
# #first loop to setup enzyme names for jaccard matrix
#   for (k in 1:1500) {
#     jaccard_data[1, k+1] <- k
#     jaccard_data[2, k+1] <- toString(ux[1,1])
#   }
# 
# #Second loop for second enzyme name for jaccard matrix
# if (length(ux) == 2) {  
#     for (l in 1:1500) {
#       jaccard_data[1, l+1501] <- l
#       jaccard_data[2, l+1501] <- toString(ux[2,1])
#     }
#   } else {break}
#   
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
for (l in 3:nrow(jaccard_data)){ # 3:28 (all samples in jaccard_data)
  for (i in 1:nrow(merge_data)) # 1:52 (all samples in merge_data)
    #If the Sample name doesn't match the row, SKIP
  { if (jaccard_data[l,1] != merge_data[i,1]){next}
    for (j in 3:ncol(merge_data))  # 3:12 (all band lengths in merge_data)
      
      #If the datapoint is NA, SKIP
    { if (is.na(merge_data[i, j])){next}
      for (k in 2:ncol(jaccard_data)){ # 2:3001 (all blanks in jaccard_data)
        
        #If enzyme in merge_data DOES NOT match enzyme in jaccard_data: SKIP
        if (tolower(toString(merge_data[i,2])) != tolower(toString(jaccard_data[2,k]))){next}
        
        #If the band size plus the variance is less than or equal to the target...
        if (merge_data[i,j] <= (as.numeric(jaccard_data[1,k])+plus_minus)
            
            #And if the band size minus variance is greater than or equal to the target...
            && as.numeric(merge_data[i,j]) >= (as.numeric(jaccard_data[1,k])-plus_minus))
        {
          #set the position in the jaccard matrix equal to 1
          jaccard_data[l,k] <- 1
        }
      }
    }
  }
}

#as.numeric(merge_data[1,3])
#as.numeric(jaccard_data[1,3001])
# make sure jaccard data goes from 3 to # of samples
#nrow(merge_data) -- 52


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

hc <- hclust(dist(d))

## Plot the Dendrogram #########################################################
jaccard_result <- plot(hc, 
                       labels=m$ID, 
                       hang = 0, 
                       cex = .55,)

