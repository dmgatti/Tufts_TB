################################################################################
# WGCNA on Tufts DO TB samples. Using the data with probes containing SNPs
# removed, normalized by the BU bioinformatics core with no batch correction.
# Daniel Gatti
# July 6 2020
# dmgatti@coa.ede
################################################################################
options(stringsAsFactors = FALSE)
library(sva)
library(pcaMethods)
library(WGCNA)

base_dir    = '/media/dmgatti/hdb/projects/TB/'
expr_dir    = file.path(base_dir, 'data', 'expression')
results_dir = file.path(base_dir, 'results', 'WGCNA')
expr_file   = file.path(expr_dir, 'tb_expr.rds')
annot_file  = file.path(expr_dir, 'tb_expr_annot.rds')
covar_file  = file.path(expr_dir, 'tb_expr_covar.rds')
pheno_file  = file.path(base_dir, 'data', 'phenotypes', 'tufts_do_tb_pheno_cleaned.csv')

# Read in expression data.
expr  = readRDS(expr_file)
annot = readRDS(annot_file)
covar = readRDS(covar_file)
covar$batch = factor(covar$batch)

pheno = read.csv(pheno_file)
pheno$mouse = as.character(pheno$mouse)
pheno = pheno[pheno$mouse %in% covar$sample,]

# See if ComBat will help.
pca_before = pca(expr, nPcs = 50)
plot(loadings(pca_before),       pch = 16, col = covar$batch)
plot(loadings(pca_before)[,2:3], pch = 16, col = covar$batch)

# Combat batch normalization with only batch.
expr_cb = ComBat(dat = expr, batch = covar$batch)
pca_after = pca(expr_cb, nPcs = 50)
plot(loadings(pca_after),       pch = 16, col = covar$batch)
plot(loadings(pca_after)[,2:3], pch = 16, col = covar$batch)

# WGCNA wants samples in rows and genes in columns.
expr = t(expr)
expr_cb = t(expr_cb)

# The scale independence plot doesn't look great. But ComBat makes it look a bit better.
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(expr_cb, powerVector = powers, verbose = 5)

par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2",
     type = "n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2], labels = powers, 
     cex = cex1, col = "red")
abline(h = 0.90, col = "red")

plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")


net = blockwiseModules(expr_cb, power = 12, TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, saveTOMFileBase = "do_tb_TOM",verbose = 3)

mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05)

# This totally didn't work and I'm worried about the expression data now.


