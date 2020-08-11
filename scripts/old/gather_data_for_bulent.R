####################################
# Gather Tufts TB data for Bulent.
# DMG
# Mar. 20 , 2019
####################################
library(rhdf5)
library(tidyverse)

setwd("/media/dmgatti/hdb/projects/TB/")

# Read in the QTL mapping input data.
load("data/phenotypes/GB_Tufts_mapping_input.Rdata")
markers = snps[rownames(snps) %in% dimnames(probs)[[3]],]
stopifnot(rownames(markers) == dimnames(probs)[[3]])

####################################
# Write out genoprobs in HDF5 format.
h5filename = "results/bulent/beamer_do_haplotype.h5"
h5createFile(h5filename)

chrs = distinct(markers, chr) %>% pull(chr)

for(chrom in chrs) {
  
  print(str_c("Writing chr ", chrom))
  chr_grp = str_c("chr", chrom)
  h5createGroup(h5filename, chr_grp)
  mkrs = filter(markers, chr == chrom)
  chr_probs = probs[,,mkrs$marker]
  h5write(chr_probs,           h5filename, str_c(chr_grp, "/probs"))
  h5write(rownames(chr_probs), h5filename, str_c(chr_grp, "/samples"))
  h5write(colnames(chr_probs), h5filename, str_c(chr_grp, "/founders"))
  h5write(mkrs,                h5filename, str_c(chr_grp, "/markers"))
  
} # for(chr)

# Check the work.
h5ls(h5filename)

# Close all HDF5 files.
h5closeAll()

####################################
# Write out QTL mapping results.





