library(DOQTL)
setwd("/hpcdata/dgatti/PowerSims/")

# Load in null phenotypes.
y = read.csv("null.phenotypes.for.typeI.error.csv")
rownames(y) = y[,1]
y = as.matrix(y[,-1])

# Load in all probs.
load("Petr/DO.gen8.1129.samples.founder.probs.Rdata")

# Load in the MUGA SNPs.
load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))
snps = muga_snps[muga_snps[,1] %in% dimnames(probs)[[3]],]

stopifnot(all(rownames(y) == rownames(probs)))
stopifnot(all(snps[,1] == dimnames(probs)[[3]]))

# Run permutations on 200 samples.
set.seed(99227471)
subset = sample(1:nrow(y), 200)
perms = scanone.eqtl(expr = y[subset,], probs = probs[subset,,], snps = snps)
save(perms, file = "qtlrel/perms.ss200.Rdata")


