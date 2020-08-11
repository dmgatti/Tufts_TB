# Run power simulations with null phenotypes.
library(DOQTL)
library(regress)
setwd("/hpcdata/dgatti/PowerSim/")
set.seed(1193857299)
args = commandArgs(trailingOnly = TRUE)
n = as.numeric(args[1]) # Number of samples.

# Load in the genotype and kinship data.
load("Petr/DO.gen8.1129.samples.founder.probs.Rdata")
load("Petr/DO.gen8.1129.samples.kinship.LOCO.Rdata")

# Load in SNPs.
load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))
snps = muga_snps[muga_snps[,1] %in% dimnames(probs)[[3]],]
snps[,2] = as.numeric(as.character(snps[,2]))

# Generate the phenotypes.
pheno = matrix(rnorm(1000 * n), nrow = n, ncol = 1000,
        dimnames = list(sample(rownames(probs), n), 1:1000))

# Perform QTL mapping using the LOCO method.
# Save only the LOD scores.
permlod = matrix(0, dim(probs)[3], 1000, dimnames =
          list(dimnames(probs)[[3]], 1:1000))
chrrng = table(snps[,2])
chrrng = c(0, cumsum(chrrng))
samples = match(rownames(pheno), rownames(probs))

for(i in 1:ncol(pheno)) {

  for(c in 1:19)  {

    snprng = (chrrng[c]+1):chrrng[c+1]

    mod = regress(pheno[,i] ~ 1, ~K[[c]][samples, samples], pos = TRUE)
    cormat = mod$sigma[[1]] * K[[c]][samples, samples] + 
             mod$sigma[[1]] * diag(1, nrow(pheno))

    tmp = scanone.eqtl(expr = pheno[,i,drop=F],
          probs = probs[samples,,snprng],
          K = cormat, snps = snps[snprng,])

    permlod[snprng,i] = tmp[,1]

  } # for(c)

} # for(i)

save(permlod, file = paste("loco/null.lod.ss", n, ".Rdata", sep = ""))

# Harvest the QTL and get the Type I error.
#setwd("/hpcdata/dgatti/PowerSim/")

# Load in permutations to get threstholds.
#load("qtlrel/perms.ss1000.Rdata")
#thr = quantile(apply(perms, 2, max), 0.95)

#for(i in n) {
#  load(file = paste("loco/null.lod.ss", i, ".Rdata", sep = ""))
#  max.lod = apply(permlod, 2, max)
#  print(mean(max.lod >= thr))
#} # for(i)

