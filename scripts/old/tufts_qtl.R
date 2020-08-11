################################################################################
# Mapping to TB traits in Gillian Beamer's DO data.
# Daniel Gatti
# dan.gatti@jax.org
# Nov. 21, 2016
################################################################################
library(DOQTL)
setwd("/hpcdata/dgatti/Tufts/")

# Load the file generated in tufts_gather_data.R.
load(file = "data/GB_Tufts_mapping_input.Rdata")

setwd("QTL")

for(i in 1:ncol(pheno.rz)) {

  t1 = proc.time()[3]
  print(paste(i, "of", ncol(pheno.rz)))
  pheno.name = colnames(pheno.rz)[i]

  qtl = scanone(pheno = pheno.rz, pheno.col = i, probs = probs, K = K,
        addcovar = covar, snps = snps)
  saveRDS(qtl, file = paste0(pheno.name, "_QTL.rds"))

  png(paste0(pheno.name, "_QTL.png"), width = 1000, height = 800, res = 128)
  plot(qtl, sig.thr = 7.2, main = pheno.name)
  dev.off()

  for(j in 1:19) {

    png(paste0(pheno.name, "_coef_chr", j,".png"), width = 1000, height = 800, 
        res = 128)
    coefplot(qtl, chr = j, main = pheno.name)
    dev.off()

  } # for(j)

  print(proc.time()[3] - t1)

} # for(i)


pheno.name = "Survival"
qtl = readRDS(file = paste0(pheno.name, "_QTL.rds"))
max.snp = qtl$lod$A[which.max(qtl$lod$A[,7]),]
pxg.plot(pheno = pheno.rz, pheno.col = pheno.name, probs = probs, 
        snp.id = max.snp[1,1], snps = snps)
