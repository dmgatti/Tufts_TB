# Run power simulations with simualted phenotypes.
library(DOQTL)
library(regress)
setwd("/hpcdata/dgatti/PowerSim/")
args = commandArgs(trailingOnly = TRUE)

# Load in the genotype and kinship data.
load("Petr/DO.gen8.1129.samples.founder.probs.Rdata")
load("Petr/DO.gen8.1129.samples.kinship.single.Rdata")

# Load in the phenotype and simulated QTL data.
load(paste("Petr/",args[1], sep = ""))

# Load in SNPs.
load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))
snps = muga_snps[muga_snps[,1] %in% dimnames(probs)[[3]],]
snps[,2] = as.numeric(as.character(snps[,2]))

# Perform QTL mapping using the LOCO method.
# Save only the LOD scores.
permlod = matrix(0, dim(probs)[3], length(pheno), dimnames =
          list(dimnames(probs)[[3]], 1:length(pheno)))

for(i in 1:length(pheno)) {

  print(i)

  ph = as.matrix(pheno[[i]])
  samples = match(rownames(ph), rownames(probs))

  mod = regress(ph ~ 1, ~K[samples, samples], pos = TRUE)
  cormat = mod$sigma[[1]] * K[samples, samples] + 
           mod$sigma[[2]] * diag(1, length(ph))

  tmp = scanone.eqtl(expr = ph, probs = probs[samples,,],
          K = cormat, snps = snps)

  permlod[,i] = tmp[,1]

} # for(i)

save(permlod, file = sub("pheno\\.qtl", "qtlrel/lod", args[1]))



