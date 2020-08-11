################################################################################
# Mapping to TB traits in Gillian Beamer's DO data.
# Daniel Gatti
# dan.gatti@jax.org
# Nov. 21, 2016
################################################################################
library(DOQTL)
library(survival)
setwd("/hpcdata/dgatti/Tufts/")

# Load the file generated in tufts_gather_data.R.
load(file = "data/GB_Tufts_mapping_input.Rdata")

setwd("QTL")

  # Keep only the mice from the first set.
  pheno = pheno[pheno$Survival <= 35,]
  samples = intersect(rownames(pheno), rownames(probs))
  pheno = pheno[samples,]
  probs = probs[samples,,]
  covar = covar[samples,]
  for(i in 1:length(K)) {
    K[[i]] = K[[i]][samples, samples]
  } # for(i)

  # Make a survival object.
  surv.obj = Surv(pheno$Survival, pheno$Survival < 35)

  i = 19

  t1 = proc.time()[3]
  print(paste(i, "of", ncol(pheno.rz)))
  pheno.name = colnames(pheno.rz)[i]

  # Calculate the null logLik.
  dose = covar[,2]
  null.mod = coxph(surv.obj ~ dose)
  null.ll = logLik(null.mod)

  ll = rep(0, dim(probs)[3])
  coef = matrix(0, dim(probs)[3], 8, dimnames = list(dimnames(probs)[[3]],
         LETTERS[1:8]))
  for(s in 1:dim(probs)[3]) {

    print(s)
    x[,-1] = probs[,-1,s]
    mod = coxph(surv.obj ~ dose + probs[,-1,s])
    ll[s] = logLik(mod) 
    coef[s,] = coef(mod)  

  } # for(s)

  lrs = ll - null.ll
  lod = lrs / (2 * log(10))
  coef[,2:8] = coef[,2:8] + coef[,1]
  coef = coef - rowMeans(coef)
  qtl = cbind(snps[rownames(coef),1:3], lod, coef)

#  qtl = scanone(pheno = pheno.rz, pheno.col = i, probs = probs, K = K,
#        addcovar = covar, snps = snps)

  saveRDS(qtl, file = paste0(pheno.name, "_QTL_survival_model.rds"))

  pdf("Survival_coef.pdf", width = 8, height = 7)
  for(j in 1:19) {

    ss = which(qtl[,2] == j)
    layout(matrix(1:2, 2, 1))
    par(plt = c(c(0.08, 0.99, 0.05, 0.9)))
    plot(qtl[ss,3], qtl[ss,4], type = "l", lwd = 2, las = 1, xlab = "",
         ylab = "LOD", xaxt = "n")
    mtext(side = 3, line = 0.5, text = paste("Survival:", "Chr", j))
    ylim = range(qtl[,5:ncol(qtl)], na.rm = T)
    par(plt = c(c(0.08, 0.99, 0.1, 0.99)))
    plot(qtl[ss,3], qtl[ss,5], type = "l", lwd = 2, col = do.colors[1,3],
         las = 1, xlab = "", ylab = "LOD", ylim = c(-2,2))
    for(k in 2:8) {
      points(qtl[ss,3], qtl[ss,4+k], type = "l", lwd = 2, col = do.colors[k,3])
    } # for(k)

  } # for(j)
  dev.off()

  print(proc.time()[3] - t1)




