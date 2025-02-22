################################################################################
#                                                                              #
# MothurTools Functions Source Code                                            #
#                                                                              #
################################################################################
#                                                                              #
# Written by: Mario Muscarella                                                 #
#                                                                              #
# Last update: 2019/01/04 by N. Wisnoski to drop reshape and work with stringr
#  
#                                                                              #
################################################################################
#                                                                              #
# Notes: This code provides functions to be used in the analysis of   #
#        16S rRNA sequence data post mothur anlaysis                           #
#                                                                              #
# Issues: Slow performance reading in OTU tables (common R issue)              #
#                                                                              #
# Recent Changes:                                                              #
#                                                                              #
# Future Changes (To-Do List):                                                 #
#         1. Design functions to work with shared files in memory              #
#         2. Add warnings                                                      #
#                                                                              #
################################################################################
library("stringr")
#library("data.table")

# Import OTU Site-by-Species Matrix
read.otu <- function(shared = " ", cutoff = "0.03"){
  matrix <- data.table::fread(shared, header=T, fill=TRUE, sep="\t")
  matrix.cutoff <- subset(matrix, matrix$label == cutoff)
  matrix.out    <- as.matrix(matrix.cutoff[1:dim(matrix.cutoff)[1],
                                           4:(3+mean(matrix.cutoff$numOtus))])
  row.names(matrix.out) <- matrix.cutoff$Group
  return(matrix.out)
}

# Import Taxonomy Information
read.tax <- function(taxonomy = " ", format = "rdp"){
  tax_raw <- read.delim(taxonomy)                 # load genus-level data
  if (format == "rdp"){
    tax <- str_split_fixed(tax_raw[,3], pattern = "\\;", 6)
    colnames(tax) <- c("Domain","Phylum","Class","Order","Family","Genus")
    tax <- cbind.data.frame(OTU = tax_raw[,1], tax)
    for (i in 2:7){
      tax[,i] <- gsub("\\(.*$", "", tax[,i])
    }
  } else {
    stop("This funciton currently only works for RDP taxonomy")
  }
  return(tax)
}