# Read in the simulated phenotypes and measure Type I error using 
# association mapping.
setwd("/hpcdata/dgatti/PowerSim/")
library(DOQTL)

# Load in the phenotypes.
pheno = read.csv("null.phenotypes.for.typeI.error.csv")
rownames(pheno) = pheno[,1]
pheno = as.matrix(pheno[,-1])

# Load in the founder probs and kinship.
load("DO.gen8.1129.samples.founder.probs.Rdata")
load("DO.gen8.1129.samples.kinship.Rdata")

# Load in MUGA SNPs and keep only Chr 1:19.
load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))
muga_snps = muga_snps[muga_snps[,2] %in% 1:19,]
probs = probs[,,dimnames(probs)[[3]] %in% muga_snps[,1]]
muga_snps = muga_snps[muga_snps[,1] %in% dimnames(probs)[[3]],]
probs = probs[,,match(muga_snps[,1], dimnames(probs)[[3]])]

# Load in the permutations for significance thresholds.
load("~/PowerSim/assoc/assoc.perms.Rdata")
thr = quantile(perms, probs = 0.05)

# Set up the SNP file to use.
snp.file = "/hpcdata/dgatti/SNP/cc.snps.NCBI38.txt.gz"

# Pick 600 samples and map them.
setwd("assoc_null")
TypeI.error = rep(0, 1000)
for(i in 1:1000) {
  print(paste("perm", i))
  keep = sample(x = 1:nrow(pheno), size = 600)
  t1 = proc.time()[3]

  for(c in 1:19) {
    a = assoc.map(pheno = pheno[keep,], pheno.col = i, probs = probs[keep,,],
        K = K[keep,keep], snps = muga_snps, chr = c, start = 0, end = 200,
        output = "p-value", snp.file = snp.file)
    TypeI.error[i] = as.numeric(TypeI.error[i] | any(a[,12] < thr))
  } # for(c)

  print(proc.time()[3] - t1)
} # for(i)
write(TypeI.error, file = "assoc.null.typeI.error.txt", sep = "\n", 
      ncolumns = 1)
