################################################################################
# Gather the TB traits in Gillian Beamer's DO data.
# Daniel Gatti
# dan.gatti@jax.org
# Nov. 21, 2016
################################################################################
library(DOQTL)
setwd("/hpcdata/dgatti/Tufts/")

# Load phenotypes.
pheno = read.csv("data/GB_Tufts_phenotypes.csv")
rownames(pheno) = pheno[,1]

# Load genotypes and subset.
load("/hpcdata/dgatti/GigaMUGA/HMM/allele/gigamuga_allele_founder_probs.Rdata")

# Keep only the Tufts data.
probs = probs[grep("^GLB", rownames(probs)),,]
rownames(probs) = sub("^GLB", "", rownames(probs))
rownames(probs) = sub("Tufts_GB(\\.)?", "", rownames(probs))
rownames(probs) = sub("07.2J.DO|2_", "", rownames(probs))
rownames(probs) = sub("^_", "", rownames(probs))

# We need to remove the incorrect first run samples from 173 to 184.
probs = probs[-grep("\\.1$", rownames(probs)),,]

# Subset the phenotype and probs data to contain the same samples in the
# same order.
samples = intersect(rownames(pheno), rownames(probs))
pheno = pheno[samples,]
probs = probs[samples,,]
pheno[,7:27] = apply(pheno[,7:27], 2, as.numeric)
#save(pheno, probs, file = "data/GB_Tufts_pheno_probs.Rdata")

# Make a few diagnostic plots. The mice were given different amounts of TB.
for(i in 7:27) {
  boxplot(I(pheno[,i]+1) ~ pheno[,6], main = colnames(pheno)[i], log = "y")
  locator(1)
}

# We will rankZ transform each phenotype and use Mtb initial dose as
# a covariate. All mice were females.

# I will not use the unexposed mice.
pheno = pheno[pheno$Mtb.initial.dose > 0,]
probs = probs[rownames(pheno),,]

pheno.rz = apply(pheno[,7:27], 2, rankZ)

# Make the kinship matrices.
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
snps = GM_snps[,1:4]
K = kinship.probs(probs = probs, snps = snps, bychr = TRUE)

# Set up covariates.
covar = model.matrix(~Mtb.initial.dose, data = pheno)
colnames(covar)[1] = "sex"

# Save all of this data to use in mapping.
save(pheno, pheno.rz, covar, probs, K, snps, 
     file = "data/GB_Tufts_mapping_input.Rdata")


# Correlation plot between phenotypes.
library(GGally)
png("figures/pheno_cor.png", width = 3000, height = 3000,
    res = 128)
ggpairs(data.frame(pheno.rz[,-c(2:3, 8:9)]))
dev.off()
