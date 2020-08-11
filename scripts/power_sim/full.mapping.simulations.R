################################################################################
# Full model mapping simulations.
# Using simulated phenotype data with MAF = 2, sample size = 600 and 
# effect size = 0.5. We run 1,000 simulated phenotypes.
################################################################################
options(stringsAsFactors = F)
library(DOQTL)

setwd("/hpcdata/dgatti/PowerSim/")

# Load in the DO founder probs.
load("DO.gen8.probs.Rdata")

# Load in the kinship matrix.
load("DO.gen8.1129.samples.kinship.Rdata")
rownames(K) = make.names(rownames(K))
colnames(K) = rownames(K)
K = K[rownames(K) %in% dimnames(probs)[[1]], colnames(K) %in% dimnames(probs)[[1]]]

stopifnot(all(rownames(K) == dimnames(probs)[[1]]))
stopifnot(all(colnames(K) == dimnames(probs)[[1]]))

# Load in the MUGA SNPs.
load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))
snps = muga_snps
snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
snps = snps[snps[,2] %in% 1:19,]
probs = probs[,,dimnames(probs)[[3]] %in% snps[,1]]

stopifnot(all(snps[,1] == dimnames(probs)[[3]]))

maf = 2
s = 100
e = 0.5
print(paste("maf", maf, "sample", s, "effect", e))

# This loads pheno and qtl.
load(file = paste("pheno.qtl.maf", maf, "sample.size", s ,"effect.size", 
     e, "Rdata", sep = "."))
pheno = lapply(pheno, function(z) { names(z) = make.names(names(z)); z })

# Because the LLP samples hae been dropped, we have to sample from the
# rest of the samples in the 100 sample size group.

setwd("full")
for(i in 1:length(pheno)) {
  print(paste("pheno", i))
  t1 = Sys.time()
  pheno[[i]] = pheno[[i]][names(pheno[[1]]) %in% dimnames(probs)[[1]]]
  pheno[[i]] = pheno[[i]][1:600]
  keep = which(dimnames(probs)[[1]] %in% names(pheno[[i]]))
  ph = matrix(pheno[[i]])
  rownames(ph) = names(pheno[[i]])
  res = scanone(pheno = ph, pheno.col = 1, probs = probs[keep,,],
        K = K[keep, keep], snps = snps)
  save(res, file = paste("pheno", i, ".QTL.Rdata", sep = ""))
  print(Sys.time() - t1)
} # for(i)

