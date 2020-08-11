################################################################################
# Association mapping simulations.
# Using simulated phenotype data with MAF = 2, sample size = 600 and 
# effect size = 0.5. We run 1,000 simulated phenotypes.
################################################################################
options(stringsAsFactors = F)
library(DOQTL)

# Load in the DO founder probs.
load("DO.gen8.1129.samples.founder.probs.Rdata")

# Load in the kinship matrix.
load("DO.gen8.1129.samples.kinship.Rdata")
colnames(K) = rownames(K)
all(rownames(K) == dimnames(probs)[[1]])
all(colnames(K) == dimnames(probs)[[1]])

# Load in the MUGA SNPs.
load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))
snps = muga_snps
snps = snps[snps[,1] %in% dimnames(probs)[[3]],]

snp.file = "/hpcdata/dgatti/SNP/cc.snps.NCBI38.txt.gz"

maf = 2
s = 600
e = 0.5
print(paste("maf", maf, "sample", s, "effect", e))

# This loads pheno and qtl.
load(file = paste("assoc.pheno.qtl.maf", maf, "sample.size", s ,"effect.size", 
     e, "Rdata", sep = "."))

setwd("assoc")
for(i in 1:length(pheno)) {
  print(paste("pheno", i))
  t1 = Sys.time()
  keep = which(dimnames(probs)[[1]] %in% names(pheno[[i]]))
  ph = matrix(pheno[[i]])
  rownames(ph) = names(pheno[[i]])
#  for(c in 1:19) {
    c = which(qtl[[i]]$sim)
    res = assoc.map(pheno = ph, pheno.col = 1, probs = probs[keep,,],
          K = K[keep, keep], chr = c, start = 0, end = 200e6, snps = snps, 
          snp.file = snp.file)
    save(res, file = paste("pheno", i, ".chr", c, ".QTL.Rdata", sep = ""))
#  } # for(c)
  print(Sys.time() - t1)
} # for(i)

