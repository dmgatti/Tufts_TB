################################################################################
# Association mapping permutations.
# Using simulated phenotype data with MAF = 2, sample size = 600 and 
# effect size = 0.5. We run 1,000 simulated phenotypes.
################################################################################
options(stringsAsFactors = F)
library(DOQTL)

setwd("/hpcdata/dgatti/PowerSim/")

# Load in the DO founder probs.
#load("DO.gen8.1129.samples.founder.probs.Rdata")

# Load in the MUGA SNPs.
#load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))
#muga_snps = muga_snps[!is.na(muga_snps[,3]),]
#snps = muga_snps
#snps = snps[snps[,1] %in% dimnames(probs)[[3]],]

snp.file = "/hpcdata/dgatti/SNP/cc.snps.NCBI38.txt.gz"

maf = 2
s = 600
e = 0.5
print(paste("maf", maf, "sample", s, "effect", e))

# This loads pheno and qtl.
#load(file = paste("assoc.pheno.qtl.maf", maf, "sample.size", s ,"effect.size", 
#     e, "Rdata", sep = "."))

#set.seed(234238907)

#pheno.perms = matrix(0, length(pheno[[1]]), length(pheno), dimnames = 
#              list(names(pheno[[1]]), 1:1000))
#for(i in 1:length(pheno)) {
#  pheno.perms[,i] = pheno[[i]]
#} # for(i)

#probs = probs[dimnames(probs)[[1]] %in% rownames(pheno.perms),,]

setwd("assoc")

#save(pheno.perms, probs, snps, snp.file, file = "perm.data.Rdata")

for(i in 1:50) {
  print(date())
  load(file = "perm.data.Rdata")
  perms = assoc.map.perms(pheno = pheno.perms, pheno.col = i, probs = probs, 
        snps = snps, model = "additive",
        scan = "one", output = "p-value", snp.file = snp.file,
        nperm = 1000)
  save(perms, file = paste("perms", i, "Rdata", sep = "."))
  rm(list = ls()[ls() != "i"])
  gc()
  print(date())
} # for(i)






