################################################################################
# Association mapping harvest.
# Using simulated phenotype data with MAF = 2, sample size = 600 and 
# effect size = 0.5. We run 1,000 simulated phenotypes.
################################################################################
options(stringsAsFactors = F)
library(DOQTL)

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

setwd("/hpcdata/dgatti/PowerSim/assoc")
qtl.files = dir(pattern = "QTL.Rdata$")
perm.files = dir(pattern = "^perms")

qn = as.numeric(gsub("^pheno|\\.chr[0-9]+\\.QTL.Rdata$", "", qtl.files))
pn = as.numeric(gsub("^perms\\.|\\.Rdata$", "", perm.files))

qtl.files = qtl.files[order(qn)]
perm.files = perm.files[order(pn)]
qn = sort(qn)
pn = sort(pn)
qtl.files = qtl.files[qn %in% pn]
perm.files = perm.files[pn %in% qn]
stopifnot(all(qn == pn))

found = rep(0, length(qtl.files))

for(i in 1:length(pheno)) {
  print(paste("pheno", i))
  t1 = proc.time()[3]
  load(qtl.files[i])
  load(perm.files[i])
  thr = quantile(perms, 0.95)

  res = res[,c(1,2,12)]
  res = res[res[,3] >= thr,]
  if(nrow(res) > 0) {
    found[i] = any(abs(res[,2] * 1e-6 - qtl[[i]]$loc[qtl[[i]]$sim]) < 5)
  } # if(nrow(res) > 0)
  
  print(proc.time()[3] - t1)
} # for(i)

save(found, file = "../assoc.power.results.txt")

