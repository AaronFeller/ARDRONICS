library(seqRFLP)
library(tidyverse)
#library(reshape2)

## Set your primers
forward_primer <- "AGAGTTTGATCMTGGCTCAG"
reverse_primer <- "GGTTACCTTGTTACGACTT"

## Reverse compliment of both forward and reverse primers.
forward_primer_rc <- revComp(forward_primer)
reverse_primer_rc <- revComp(reverse_primer)

## Import FASTA file
fasta_file <- read.fasta("/home/aaron/R/ARDRA-Analyzer/20210505_lib.fa")



## Clip FASTA file for primers
clipped_fasta_file <- clipprobe(fasta_file,
                               forward_primer,
                               reverse_primer_rc,
                               tol = 1,
                               clipped.only = TRUE)

 ## Data for clipped library
 clipped_length <- length(gnames.fas(clipped_fasta_file))
 clipped_length
 length(gnames.fas(fasta_file))
 setdiff(gnames.fas(fasta_file), gnames.fas(clipped_fasta_file))

write(clipped_fasta_file,
          "/home/aaron/R/ARDRA-Analyzer/representative_library_clipped.fa")



################################################################################

## This imports the enzyme data
data(enzdata)
## Naming for single cut enzymes you'd like to use
single_enznames <- c("HaeIII", "HhaI")

## Multi-Cut for entire library ################################################
full_data <- matrix(data = NA, nrow = 0, ncol = 27)
colnames(full_data) <- c("Sample", "enznames", "1", "2", "3", "4", "5", "6", 
                         "7", "8", "9", "10", "11", "12", "13", "14", "15", 
                         "16", "17", "18", "19", "20", "21", "22", "23", "24", 
                         "25")

for (i in 1:length(single_enznames)) {
  frag_data <- frag.dat(fasta_file, 
                        single_enznames[i], 
                        enzdata = enzdata)
  frag_data <- frag_data %>% select(1,4)
  frag_data <- separate(frag_data, 
                        col = "fragment_Length", 
                        into = c("1", "2", "3", "4", "5", "6", "7", "8", "9", 
                        "10", "11", "12", "13", "14", "15", "16", "17", "18", 
                        "19", "20", "21", "22", "23", "24", "25"), 
                        sep = ",")
  
  frag_data <- tibble::rownames_to_column(frag_data, "Sample")
  full_data <- rbind(full_data, frag_data)
  }

## Double Cutting! #############################################################
## Double Cut Enzyme Names
double_enznames <- c("HaeIII", "HhaI")

for (i in 1:(length(clipped_fasta_file)/2))   {
  double_cutter <- enzCut(clipped_fasta_file[i*2], 
                          double_enznames, 
                          enzdata = enzdata)
  double_cutter_result <- data.frame(lapply(double_cutter[3], 
                                            function(x) t(data.frame(x))))
  double_cutter_result <- add_column(double_cutter_result, 
                                     paste(double_enznames[1],"&",double_enznames[2]), 
                                     .before = 1)
  double_cutter_result <- add_column(double_cutter_result, 
                                     clipped_fasta_file[(i*2)-1], 
                                     .before = 1)
  double_cutter_result
  add_cols <- matrix(data = NA, 
                     nrow = 1, 
                     ncol = (27-ncol(double_cutter_result)))
  double_cutter_result <- cbind(double_cutter_result, 
                                add_cols)
  colnames(double_cutter_result) <- c("Sample", "enznames", "1", "2", "3", "4", "5", "6", 
                           "7", "8", "9", "10", "11", "12", "13", "14", "15", 
                           "16", "17", "18", "19", "20", "21", "22", "23", "24", 
                           "25")
  full_data <- rbind(full_data, 
                     double_cutter_result)
  }
  




write.csv(full_data, 
          "/home/aaron/R/ARDRA-Analyzer/20210505_lib_clip.csv", 
          row.names = FALSE)

  

#plotenz(clipped_fasta_file, 
#        single_enznames, 
#        enzdata, 
#        side = TRUE, 
#        type = c("RFLP", "TRFLP"), 
#        Terminal = c("T5", "T3")
#        )