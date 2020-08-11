setwd("/hpcdata/dgatti/PowerSim/full_null")

ph.f = dir(pattern = "QTL.Rdata")
prm.f = dir(pattern = "perms.Rdata")

ph.num = gsub("^null|\\.QTL.Rdata$", "", ph.f)
prm.num = gsub("^null|\\.perms.Rdata$", "", prm.f)
stopifnot(all(ph.num == prm.num))

typeI.error = matrix(0, 101, length(ph.num))

for(i in 1:length(ph.f)) {
  load(ph.f[i])
  load(prm.f[i])
  thr = quantile(perms, 0:100/100)
  gt = outer(res$lod$A[,7], thr, ">")
  typeI.error[,i] = as.numeric(colSums(gt) > 0)
} # for(i)

write.table(typeI.error, file = "full.null.sims.results.txt", sep = "\t")


