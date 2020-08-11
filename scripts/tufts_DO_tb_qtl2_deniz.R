##############################
# Tufts TB QTL Mapping.
# DMG
# Dec. 30, 2019
##############################
library(tidyverse)
library(qtl2)
library(qtl2convert)

base_dir = "/media/dmgatti/hdb/projects/TB/"

# Read in the U. Wisc. remapped markers.
markers = readr::read_csv('https://raw.githubusercontent.com/kbroman/MUGAarrays/master/UWisc/gm_uwisc_v1.csv',
                          col_types = 'ccnnncclllnnc') %>%
          filter(!is.na(chr)) %>%
          as.data.frame()
colnames(markers)[3:4] = c('pos', 'cM')
rownames(markers) = markers$marker

map  = qtl2convert::map_df_to_list(map = markers, pos_column = 'pos')
gmap = qtl2convert::map_df_to_list(map = markers, pos_column = 'cM')

# Read in cleaned phenotypes.
pheno    = read.csv(str_c(base_dir, "data/phenotypes/tb_pheno.csv"))
pheno_rz = read.csv(str_c(base_dir, "data/phenotypes/tb_pheno_rz.csv"))
rownames(pheno)    = pheno$mouse
rownames(pheno_rz) = pheno_rz$mouse

orig_pheno = pheno

# Read in probs.
probs = readRDS('/media/dmgatti/hdb/projects/TB/data/tufts_doqtl_probs_from_jax.rds')
b6_names = grep('^B6', dimnames(probs)[[3]])
dimnames(probs)[[3]][b6_names] = gsub('\\.', '-', dimnames(probs)[[3]][b6_names])
common_markers = intersect(rownames(markers), dimnames(probs)[[3]])
markers = markers[rownames(markers) %in% common_markers,]
probs   = probs[,,dimnames(probs)[[3]] %in% common_markers]
probs = probs[,,rownames(markers)]

print(paste(nrow(markers), 'Markers'))

# Synch samples.
samples2keep = intersect(rownames(pheno), rownames(probs))
pheno    = pheno[samples2keep,]
pheno_rz = pheno_rz[samples2keep,]
probs    = probs[samples2keep,,]
pheno$dose    = factor(pheno$dose,    levels = c(16, 28, 97, 127))
pheno_rz$dose = factor(pheno_rz$dose, levels = c(16, 28, 97, 127))

print(paste(nrow(pheno), 'Samples'))

# Convert to qtl2 format and resave.
probs = qtl2convert::probs_doqtl_to_qtl2(probs = probs, map = markers, pos_column = "pos")
map   = qtl2convert::map_df_to_list(map = markers, pos_column = "pos")
K = qtl2::calc_kinship(probs = probs, type = "loco", cores = 8, quiet = FALSE)
covar = model.matrix(~expt+dose, data = pheno)[,-1,drop = FALSE]

save(pheno, pheno_rz, covar, probs, map, K, markers, file = str_c(base_dir, "data/phenotypes/GB_Tufts_qtl2_input.Rdata"))

#############################
# Start here once the above has been run.
library(tidyverse)
library(qtl2)
library(qtl2convert)

base_dir = "/home/dmgatti/Documents/data/TB/"

load(file = str_c(base_dir, "data/phenotypes/GB_Tufts_qtl2_input.Rdata"))
peaks = NULL

# Read in sample subset.
training_samples = read.csv(paste0(base_dir, 'data/trainIDs/trainIDs_sel_0.csv'))
training_samples = training_samples$IDs

##########
# QTL Mapping.

result_dir = str_c(base_dir, "results/qtl2/dose_as_factor/")

for(i in 8:ncol(pheno_rz)) {
  
  pheno_name = colnames(pheno_rz)[i]
  
  samples2use = which(!is.na(pheno_rz[,pheno_name]))
  samples2use = intersect(samples2use, training_samples)
  
  print(str_c(pheno_name, " : ", length(samples2use)))
  
  tmpK = K
  for(j in 1:length(K)) {
    tmpK[[j]] = K[[j]][samples2use, samples2use]
  }
  tmp_covar = covar[samples2use,,drop = FALSE]
  tmp_covar = tmp_covar[,colSums(tmp_covar) > 0, drop = FALSE]

  print("    QTL")  
  lod = qtl2::scan1(genoprobs = probs[samples2use,], 
                    pheno     = pheno_rz[samples2use, pheno_name, drop = FALSE], 
                    kinship   = tmpK, 
                    addcovar  = tmp_covar, 
                    cores = 8)
  
  png(str_c(result_dir, pheno_name, "_qtl.png"), width = 1000, height = 800, res = 128)
  qtl2::plot_scan1(lod, map = map, main = pheno_name)
  dev.off()
  
  peaks = find_peaks(lod, map, threshold = 6, prob = 0.95)
  
  save(lod, file = str_c(result_dir, pheno_name, "_qtl2.Rdata"))
  
} # for(i)
