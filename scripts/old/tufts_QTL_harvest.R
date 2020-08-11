################################################################################
# Harvest the TB traits in Gillian Beamer's DO data.
# Daniel Gatti
# dan.gatti@jax.org
# Nov. 21, 2016
################################################################################
library(DOQTL)
setwd("/hpcdata/dgatti/Tufts/")

snp.file = "/hpcdata/shared/sanger/current_snps/REL-1505/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"

# Load the file generated in tufts_gather_data.R.
load(file = "data/GB_Tufts_mapping_input.Rdata")

perms = readRDS("/hpcdata/gac/projects/Svenson_DO850/QTL/perms/all_perms.rds")
thr = quantile(perms[,1], probs = c(0.95, 0.8))

setwd("QTL")

result = data.frame(phenotype = "", marker = "", chr = "",
         pos = 0, lod = 0, perc.var = 0, proximal = 0, distal = 0,
         stringsAsFactors = FALSE)

for(i in 1:ncol(pheno.rz)) {

  t1 = proc.time()[3]
  print(paste(i, "of", ncol(pheno.rz)))
  pheno.name = colnames(pheno.rz)[i]

  qtl = readRDS(file = paste0(pheno.name, "_QTL.rds"))

  png(paste0(pheno.name, "_QTL.png"), width = 1000, height = 800, res = 128)
  plot(qtl, sig.thr = c(7.2, 6), sig.col = c("red", "goldenrod"), 
       main = pheno.name)
  dev.off()

  for(j in 1:19) {

    png(paste0(pheno.name, "_coef_chr", j,".png"), width = 1000, height = 800, 
        res = 128)
    coefplot(qtl, chr = j, main = pheno.name)
    dev.off()

  } # for(j)

  png(paste0(pheno.name, "_coef_chrX.png"), width = 1000, height = 800, 
      res = 128)
  coefplot(qtl, chr = "X", sex = "F", main = pheno.name)
  dev.off()

  # Use a threshold of 6 to make association plots.
  loci = qtl$lod$A[qtl$lod$A$lod > 6,]
  chr = unique(loci[,2])

  for(j in chr) {

    interval = bayesint(qtl, chr = j)
    result = rbind(result, c(pheno.name, unlist(interval[2,c(1:3,7,5)]),
             interval[1,3], interval[3,3]))

    assoc = assoc.map(pheno = pheno.rz, pheno.col = i, probs = probs, K = K[[j]],
            addcovar = covar, snps = snps, chr = j, start = interval[2,3] - 2,
            end = interval[2,3] + 2, snp.file = snp.file)
    saveRDS(assoc, file = paste0(pheno.name, "_assoc_chr", j,".rds"))

    png(paste0(pheno.name, "_assoc_chr", j,".png"), width = 2400, height = 1800, 
        res = 200)
    tmp = assoc.plot(assoc, thr = max(assoc[,12]) - 0.5, show.sdps = TRUE)
    dev.off()

  } # for(j)

  print(proc.time()[3] - t1)

} # for(i)

write.csv(result, file = "GLB_Tufts_QTL_summary.csv", row.names = FALSE,
          quote = FALSE)

