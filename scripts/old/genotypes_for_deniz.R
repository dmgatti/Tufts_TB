# Gather genotype data to send to Deniz @ RPI.
# Daniel Gatti
# July 21, 2019
##############################################
library(argyle)

base_dir = '/media/dmgatti/hdb/projects/TB/data/genotypes'
snp_file = "http://csbio.unc.edu/MUGA/snps.gigamuga.Rdata"

# Loads in 'snps'.
load(url(snp_file))

geno_dirs = list.dirs(base_dir, full.names = TRUE)[-1]

result = NULL

for(d in geno_dirs) {

  print(d)  
  geno = read.beadstudio(prefix = "", snps = snps, in.path = d)

  gt = as.data.frame(geno)
  if(is.null(result)) {
    result = gt
  } else {
    result = cbind(result, gt[,-(1:6), drop = FALSE])
  } # else
    
} # for(d)

write.csv(result, file = paste0(base_dir, "/tufts_genotypes.csv"), row.names = FALSE)
