setwd("/hpcdata/dgatti/PowerSim/")

###
# Harvest QTLRel results.
files = dir(path = "qtlrel", pattern = "^lod.*Rdata$", full.names = T)
load(files[1])

# Get MUGA SNP locations.
load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))
snps = muga_snps[muga_snps[,1] %in% rownames(permlod),]
snps[,2] = as.numeric(as.character(snps[,2]))

# Load in the permutations and get the 0.05 threshold.
load("qtlrel/perms.ss1000.Rdata")
thr = quantile(apply(perms, 2, max), 0.95)

# Create the results data.frame.
results = data.frame(maf    = rep(0, length(files)),
                     effect = rep(0, length(files)),
                     sample = rep(0, length(files)),
                     pvar   = rep(0, length(files)),
                     power  = rep(0, length(files)))
# Go through each file and determine if we detected the QTL in the
# correct location and record the percent variance explained based 
# on the LOD score.
for(i in 1:length(files)) {

  print(i)

  # Load in the simulated QTL locations for this maf, smaple and effect
  # size.
  load(sub("^qtlrel/lod", "Petr/pheno.qtl", files[i]))

  # Load in the simulation results.
  load(files[i])
  spl = strsplit(files[i], split = "\\.")[[1]]
  results$maf[i] = spl[3]
  results$effect[i] = spl[10]
  results$sample[i] = spl[6]

  pvar  = rep(NA, ncol(permlod))
  detected = rep(F, ncol(permlod))
  n = as.numeric(results$sample[i])

  for(j in 1:ncol(permlod)) {
    # Get the simulated QTL location.
    sim.qtl = qtl[[j]][qtl[[j]]$sim,]
    # Some QTL were simulated at NA locations.
    if(!is.na(sim.qtl[1,2])) {
      # Keep LOD scores above the threshold on the simulated chromosome.
      rng = which(snps[,2] == sim.qtl[1,1])
      lod = permlod[rng,j]
      # Record the % variance explained at the simulated QTL.
      lod.at.qtl = lod[which.min(abs(sim.qtl[1,2] - snps[rng,3]))]
      pvar[j] = 1.0 - exp(lod.at.qtl * (2 * log(10) / -n))
      # Keep LODs above the threshold.
      lod = lod[lod > thr]
      # Record whether the QTL was detected.
      if(length(lod) > 0) {
        locs = snps[match(names(lod), snps[,1]),3]
        detected[j] = any(abs(locs - sim.qtl[1,2]) <= 5.0, na.rm = T)
      } # if(length(lod) > 0)
    } # if(!is.na(sim.qtl[1,2])
  } # for(j)

  results$pvar[i]  = mean(pvar, na.rm = T)
  results$power[i] = mean(detected, na.rm = T)

} # for(i)

save(results, file = "power.sim.qtlrel.results.Rdata")


###
# Harvest LOCO results.
files = dir(path = "loco", pattern = "^lod.*Rdata$", full.names = T)
load(files[1])

# Get MUGA SNP locations.
load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))
snps = muga_snps[muga_snps[,1] %in% rownames(permlod),]
snps[,2] = as.numeric(as.character(snps[,2]))

# Load in the permutations and get the 0.05 threshold.
load("qtlrel/perms.ss1000.Rdata")
thr = quantile(apply(perms, 2, max), 0.95)

# Create the results data.frame.
results2 = data.frame(maf   = rep(0, length(files)),
                     effect = rep(0, length(files)),
                     sample = rep(0, length(files)),
                     pvar   = rep(0, length(files)),
                     power  = rep(0, length(files)))

# Go through each file and determine if we detected the QTL in the
# correct location and record the percent variance explained based
# on the LOD score.
for(i in 1:length(files)) {

  print(i)

  # Load in the simulated QTL locations for this maf, smaple and effect
  # size.
  load(sub("^loco/lod", "Petr/pheno.qtl", files[i]))

  # Load in the simulation results.
  load(files[i])
  spl = strsplit(files[i], split = "\\.")[[1]]
  results2$maf[i] = spl[3]
  results2$effect[i] = spl[10]
  results2$sample[i] = spl[6]

  pvar  = rep(NA, ncol(permlod))
  detected = rep(F, ncol(permlod))
  n = as.numeric(results2$sample[i])
  
  for(j in 1:ncol(permlod)) {
    # Get the simulated QTL location.
    sim.qtl = qtl[[j]][qtl[[j]]$sim,]
    # Some QTL were simulated at NA locations.
    if(!is.na(sim.qtl[1,2])) {
      # Keep LOD scores above the threshold on the simulated chromosome.
      rng = which(snps[,2] == sim.qtl[1,1])
      lod = permlod[rng,j]
      # Record the % variance explained at the simulated QTL.
      lod.at.qtl = lod[which.min(abs(sim.qtl[1,2] - snps[rng,3]))]
      pvar[j] = 1.0 - exp(lod.at.qtl * (2 * log(10) / -n))
      # Keep LODs above the threshold.
      lod = lod[lod > thr]
      # Record whether the QTL was detected.
      if(length(lod) > 0) {
        locs = snps[match(names(lod), snps[,1]),3]
        detected[j] = any(abs(locs - sim.qtl[1,2]) <= 5.0, na.rm = T)
      } # if(length(lod) > 0)
    } # if(!is.na(sim.qtl[1,2])
  } # for(j)

  results2$pvar[i]  = mean(pvar, na.rm = T)
  results2$power[i] = mean(detected, na.rm = T)

} # for(i)

save(results2, file = "power.sim.loco.results.Rdata")


###
# Make a figure comparing the results.
setwd("/hpcdata/dgatti/PowerSim/")
load("power.sim.qtlrel.results.Rdata")
load("power.sim.loco.results.Rdata")

png("power.sim.png", width = 2000, height = 1000, res = 200)
layout(matrix(1:2, 1, 2))
plot(0, 0, col = 0, xlim = c(0, 0.35), ylim = c(0, 1), xlab = "% Variance",
     ylab = "Power", main = "QTLRel")
ss = sort(as.numeric(unique(results$sample)))
for(i in 1:length(ss)) {
  rng = which(results$sample == ss[i])
  points(results$pvar[rng], results$power[rng], pch = 16, col = i)
}
legend("bottomright", pch = 16, col = 1:length(ss), legend = ss)
plot(0,	0, col = 0, xlim = c(0,	0.35), ylim = c(0, 1), xlab = "% Variance",
     ylab = "Power", main = "LOCO")
ss = sort(as.numeric(unique(results2$sample)))
for(i in 1:length(ss)) {
  rng = which(results2$sample == ss[i])
  points(results2$pvar[rng], results2$power[rng], pch = 16, col = i)
}
legend("bottomright", pch = 16, col = 1:length(ss), legend = ss)
dev.off()
