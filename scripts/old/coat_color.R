################################################################################
# Check the coat color of the Tuft's mice.
# Daniel Gatti
# dan.gatti@jax.org
# June 22, 2017
################################################################################
options(stringsAsFactors = F)
library(DOQTL)

setwd("/hpcdata/dgatti/Tufts/")

pheno = read.csv("data/GB_Tufts_phenotypes.csv")
rownames(pheno) = pheno[,1]

library(rhdf5)
h5filename = "/hpcdata/gac/derived/CGD_DO_Genoprobs/GigaMUGA_hap_probs_v1.h5"
grp = h5read(file = h5filename, name = "/GLB")
dimnames(grp$probs) = list(grp$samples, grp$founders, grp$markers)
probs = grp$probs
rm(grp)

samples = intersect(rownames(pheno), rownames(probs))
pheno = pheno[samples,]
probs = probs[samples,,]
stopifnot(rownames(pheno) == rownames(probs))

coat = pheno$Coat.color
coat[grep("br", coat)] = "agouti"
coat[grep("bl", coat)] = "black"
coat[grep("wh", coat)] = "albino"
names(coat) = rownames(pheno)

load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
snps = GM_snps[,1:4]
snps = snps[dimnames(probs)[[3]],]

coat.chk = check.do.coat.color(probs = probs, markers = snps, coat = coat)
coat.chk = data.frame(coat.chk, pheno$Coat.color)
write.csv(coat.chk, file = "results/Tufts_coat_color_check.csv", quote = F,
          row.names = F)

# Map coat color.
library(qtl2)
library(qtl2convert)

probs = probs_doqtl_to_qtl2(probs = probs, map = snps, chr_column = "chr", 
        pos_column = "pos", marker_column = "marker")

coat = factor(coat, levels= c("albino", "agouti", "black"))
coat = matrix(as.numeric(coat) - 2, ncol = 1, dimnames = list(rownames(pheno), "coat"))

map = map_df_to_list(map = snps, pos_column = "pos")

qtl = scan1(genoprobs = probs, pheno = coat, cores = 5)
plot_scan1(qtl, map)

chr = 2
coef2 = scan1coef(genoprobs = probs[,chr], pheno = coat)
plot_coef(coef2, map = map[[chr]], columns = 1:8, col = CCcolors,
          scan1_output = qtl)

chr = 7
coef7 = scan1coef(genoprobs = probs[,chr], pheno = coat)
plot(coef7, map = map[[chr]], columns = 1:8, col = CCcolors, scan1_output = qtl)


