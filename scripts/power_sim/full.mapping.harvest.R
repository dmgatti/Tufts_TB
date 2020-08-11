setwd("/hpcdata/dgatti/PowerSim/full_sims")
load("../pheno.qtl.maf.2.sample.size.600.effect.size.0.5.Rdata")

ph.f = dir(pattern = "QTL.Rdata")
prm.f = dir(pattern = "perms.Rdata")

ph.num = gsub("^pheno|\\.QTL.Rdata$", "", ph.f)
prm.num = gsub("^pheno|\\.perms.Rdata$", "", prm.f)
stopifnot(all(ph.num == prm.num))

detected = rep(F, length(ph.num))

for(i in 1:length(ph.f)) {
  print(i)
  load(ph.f[i])
  load(prm.f[i])
  q = qtl[[i]][qtl[[i]]$sim,]
  res = res$lod$A[res$lod$A[,2] == q[1,1],]
  thr = quantile(perms, 0.95)
  res = res[res[,7] >= thr,]
  detected[i] = any(abs(res[,3] - q[1,2]) < 5)
} # for(i)

write(detected, file = "full.sims.results.txt", ncolumns = 1, sep = "\n")

