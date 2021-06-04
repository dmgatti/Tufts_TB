################################################################################
# Tufts TB QTL Mapping. Image data from Metin's group.
# DMG
# May 20, 2020
##############################
library(biomaRt)
library(survival)
library(tidyverse)
library(qtl2)

base_dir   = '/media/dmgatti/hdb/projects/TB'
input_file = file.path(base_dir, 'data', 'tufts_do_tb_qtl2_input.Rdata') 
result_dir = file.path(base_dir, 'results', 'qtl2', 'gen_factor2')

# Set up biomaRt.
ensembl = useEnsembl(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl')
#bmfilt = listFilters(ensembl)
#bmattr = listAttributes(ensembl)

# Read in the DOQTL formatted data.
load(input_file)

# CC SNP and gene database files (from Karl).
ccsnpdb = "/media/dmgatti/hda/data/MUGA/cc_variants.sqlite"
mgidb   = "/media/dmgatti/hda/data/MUGA/mouse_genes_mgi.sqlite"
snp_func  = create_variant_query_func(dbfile = ccsnpdb)
gene_func = create_gene_query_func(dbfile = mgidb)
ensembl = 93

peaks = NULL

##########
# Heritability.

addcovar = model.matrix(~gen, data = covar)[,-1,drop = FALSE]

herit = data.frame(pheno = colnames(pheno_rz),
                   herit = 0)

K_all = calc_kinship(probs, type = 'overall', cores = 4)

for(i in 1:ncol(pheno_rz)) {
  
  herit$herit[i] = qtl2::est_herit(pheno = pheno_rz[,i,drop = FALSE], 
                                   kinship = K_all, addcovar = addcovar)
  
} # for(i)

write_csv(herit, file.path(result_dir, 'tufts_do_tb_herit.csv'))

##########
# QTL Mapping with coef calculations.

addcovar = model.matrix(~gen, data = covar)[,-1,drop = FALSE]

for(i in 1:ncol(pheno_rz)) {
  
  pheno_name = colnames(pheno_rz)[i]
  
  samples2use = which(!is.na(pheno_rz[,pheno_name]))
  
  print(str_c(pheno_name, " : ", length(samples2use)))
  
  tmpK = K
  for(j in 1:length(K)) {
    tmpK[[j]] = K[[j]][samples2use, samples2use]
  }
  tmp_covar = addcovar[samples2use,,drop = FALSE]
  tmp_covar = tmp_covar[,colSums(tmp_covar) > 0, drop = FALSE]
  
  print("    QTL")  
  lod = qtl2::scan1(genoprobs = probs[samples2use,], 
                    pheno     = pheno_rz[samples2use, pheno_name, drop = FALSE], 
                    kinship   = tmpK, 
                    addcovar  = tmp_covar, 
                    cores     = 10)
  
  png(file.path(result_dir, str_c(pheno_name, "_qtl.png")), width = 1000, height = 800, res = 128)
  qtl2::plot_scan1(lod, map = map, main = pheno_name)
  dev.off()
  
  current_peaks = find_peaks(lod, map, threshold = 6, prob = 0.95)
  
  # If the confidence interval takes up more than 20 Mb, use the peak LOD +/- 5 Mb.
  coefs = NA
  assocs = NA
  if(nrow(current_peaks) > 0) {
    
    for(j in 1:nrow(current_peaks)) {
      if(current_peaks$ci_hi[j] - current_peaks$ci_lo[j] > 20) {
        current_peaks$ci_hi[j] = current_peaks$pos[j] + 5
        current_peaks$ci_lo[j] = max(current_peaks$pos[j] - 5, 3)  # This protects us from going off the beginning of the chr.
      } # if(current_peaks$ci_hi[j] - current_peaks$ci_lo[j])
    } # for(j)
    
    current_peaks = data.frame(current_peaks, 
                               as.data.frame(matrix(0, nrow = nrow(current_peaks), ncol = 8, 
                                                    dimnames = list(rownames(current_peaks), LETTERS[1:8]))))
    
    coefs = vector("list", nrow(current_peaks))
    names(coefs) = current_peaks$chr
    assocs = vector("list", nrow(current_peaks))
    names(assocs) = current_peaks$chr
    
    # For each peak, plot coefs and association mapping.
    for(j in 1:nrow(current_peaks)) {
      
      # Get chr, start and end.
      curr_chr = current_peaks$chr[j]
      print(str_c("    chr ", curr_chr))
      start = current_peaks$ci_lo[j]
      end   = current_peaks$ci_hi[j]
      
      # BLUP mapping on one chromosome.
      print(paste("    coef", j))
      coefs[[j]] = qtl2::scan1blup(genoprobs = probs[samples2use, curr_chr], 
                                   pheno     = pheno_rz[samples2use, pheno_name, drop = FALSE], 
                                   kinship   = K[[curr_chr]], 
                                   addcovar  = tmp_covar,
                                   se = TRUE,
                                   cores = 5)
      stopifnot(!is.nan(coefs[[j]]))
      
      png(file.path(result_dir, str_c(pheno_name, "_coef_chr", curr_chr, ".png")), width = 1000, height = 800, res = 128)
      plot_coefCC(coefs[[j]], map, scan1_output = lod, top_panel_prop = 0.6, main = pheno_name)
      dev.off()
      
      # Add the coefs to the peak table.
      mkr = find_marker(map, chr = curr_chr, pos = current_peaks$pos[j])
      current_peaks[j, LETTERS[1:8]] = coefs[[j]][mkr, LETTERS[1:8]]
      
      # Association mapping.
      print(paste("    assoc", j))
      assocs[[j]] = scan1snps(genoprobs = probs[,curr_chr], pheno = pheno_rz[,pheno_name, drop = FALSE], 
                              kinship = tmpK[[curr_chr]], addcovar = tmp_covar, map = map, 
                              chr = curr_chr, start = max(1, start - 1), end = end + 1, 
                              query_func = snp_func, keep_all_snps = TRUE, cores = 10)
      
      genes = gene_func(chr = curr_chr, start = start - 1, end = end + 1)
      png(file.path(result_dir, str_c(pheno_name, "_assoc_chr", curr_chr, ".png")), width = 2400, height = 1600, res = 300)
      plot_snpasso(scan1output = assocs[[j]]$lod, snpinfo = assocs[[j]]$snpinfo, genes = genes, drop_hilit = 1, 
                   top_panel_prop = 0.3, main = pheno_name, colors = 'black')
      dev.off()
      
      # Get the top SNPs.
      top = top_snps(assocs[[j]]$lod, assocs[[j]]$snpinfo)
      write_csv(top, path = file.path(result_dir, str_c(pheno_name, "_assoc_chr", curr_chr, "_top_snps.csv")))
      
    } # for(j)
    
    # Add current peaks to overall peaks.
    peaks = rbind(peaks, current_peaks)
    
  } # if(nrow(current_peaks > 0))
  
  save(lod, coefs, assocs, file = file.path(result_dir, str_c(pheno_name, "_qtl2.Rdata")))
  
  write_csv(peaks, file = file.path(result_dir, "tb_qtl_peaks.csv"))
  
} # for(i)

write_csv(peaks, file = file.path(result_dir, "tb_qtl_peaks.csv"))

################################################################################
# Survival mapping
# Set up survival model.
pheno_surv = data.frame(mouse = as.character(covar$mouse),
                        surv = Surv(time = covar$euth_day, event = (covar$euth == 'y') * 1))
rownames(pheno_surv) = pheno_surv$mouse
pheno_surv = pheno_surv[complete.cases(pheno_surv),]

# Also try survival model with longer living mice censored at 60 days.
time  = covar$euth_day
event = (covar$euth == 'y') * 1
wh = which(time > 60)
time[wh] = 60
event[wh] = 0
pheno_surv_60 = data.frame(mouse = as.character(covar$mouse),
                           surv = Surv(time = time, event = event))
rownames(pheno_surv_60) = pheno_surv_60$mouse
pheno_surv_60 = pheno_surv_60[complete.cases(pheno_surv_60),]

samples = intersect(rownames(pheno_surv), rownames(addcovar))
pheno_surv = pheno_surv[samples,]
probs = probs[samples,]
for(i in seq_along(K)) {
  K[[i]] = K[[i]][samples,samples]
}

# Null model log_likelihood.
null_ll = coxph(pheno_surv$surv ~ addcovar)$loglik[2]

lod = vector('list', length(probs))
names(lod) = names(probs)

coefs = vector('list', length(probs))
names(coefs) = names(probs)

# For each chromosome...
for(chr in names(probs)) {
  
  print(paste('chr', chr))
  
  # Get genoprobs for current samples.
  # Remove 'A' as basis.
  pr = probs[[chr]][,-1,]
  print(paste('   ', dim(pr)[3], 'markers'))
  
  lod_chr  = matrix(0, nrow = dim(pr)[3], ncol = 1, dimnames = list(dimnames(pr)[[3]], "survival"))
  coef_chr = matrix(0, nrow = dim(pr)[3], ncol = 8, dimnames = list(dimnames(pr)[[3]], LETTERS[1:8]))
  b_h = paste0('pr[, , j]', LETTERS[2:8])
  
  for(j in 1:dim(pr)[3]) {
    
    if(j %% 100 == 0) print(j)
    
    mod = coxph(pheno_surv$surv ~ addcovar + pr[,,j])
    lod_chr[j,1] = mod$loglik[2]
    coef_chr[j,] = c(0, coef(mod)[b_h])
    
  } # for(j)
  
  lod[[chr]]  = (lod_chr - null_ll) / log(10)
  coefs[[chr]] = coef_chr - rowMeans(coef_chr)
  class(coefs[[chr]]) = c('scan1coef', 'scan1', 'matrix')
  
  rm(pr, lod_chr, coef_chr)
  gc()
  
} # for(i)

# Recombine chromosomes.
tmp = lod[[1]]
for(i in 2:length(lod)) {
  tmp = rbind(tmp, lod[[i]])
}

lod = tmp
class(lod) = c('scan1', 'matrix')

save(lod, coefs, file = file.path(result_dir, 'survival_qtl2.Rdata'))

png(file.path(result_dir, "survival_qtl.png"), width = 1000, height = 800, res = 128)
plot(lod, map, main = 'Survival')
dev.off()

png(file.path(result_dir, "survival_coef_chr1.png"), width = 1000, height = 800, res = 128)
plot_coefCC(coefs[[1]], map, scan1_output = lod, main = 'Survival')
dev.off()

png(file.path(result_dir, "survival_coef_chr4.png"), width = 1000, height = 800, res = 128)
plot_coefCC(coefs[[4]], map, scan1_output = lod, main = 'Survival')
dev.off()

png(file.path(result_dir, "survival_coef_chr5.png"), width = 1000, height = 800, res = 128)
plot_coefCC(coefs[[5]], map, scan1_output = lod, main = 'Survival')
dev.off()

png(file.path(result_dir, "survival_coef_chr12.png"), width = 1000, height = 800, res = 128)
plot_coefCC(coefs[[12]], map, scan1_output = lod, main = 'Survival')
dev.off()

# Add survival peaks to survival table.
current_peaks = find_peaks(lod, map, threshold = 7, prob = 0.95)
current_peaks = data.frame(current_peaks, 
                           as.data.frame(matrix(0, nrow = nrow(current_peaks), ncol = 8, 
                                        dimnames = list(rownames(current_peaks), LETTERS[1:8]))))

for(j in 1:nrow(current_peaks)) {
  
  curr_chr = current_peaks$chr[j]
  mkr = find_marker(map, chr = curr_chr, pos = current_peaks$pos[j])
  current_peaks[j, LETTERS[1:8]] = coefs[[curr_chr]][mkr, LETTERS[1:8]]
  
} # for(j)


peaks = rbind(peaks, current_peaks)
readr::write_csv(peaks, path = file.path(result_dir, "tb_qtl_peaks.csv"))

################################################################################
# Survival association mapping.
# Read in linkage mapping from above. (lod & coefs)
load(file = file.path(result_dir, 'survival_qtl2.Rdata'))
rm(coefs)

# Set up survival model.
pheno_surv = data.frame(mouse = as.character(covar$mouse),
                        surv = Surv(time = covar$euth_day, event = (covar$euth == 'y') * 1))
rownames(pheno_surv) = pheno_surv$mouse
pheno_surv = pheno_surv[complete.cases(pheno_surv),]

samples = intersect(rownames(pheno_surv), rownames(addcovar))
pheno_surv = pheno_surv[samples,]
probs = probs[samples,]
for(i in seq_along(K)) {
  K[[i]] = K[[i]][samples,samples]
}

# Null model log_likelihood.
null_ll = coxph(pheno_surv$surv ~ addcovar)$loglik[2]

# Get peak locations for each chromosome.
surv_peaks = find_peaks(lod, map, threshold = 6, prob = 0.95)
surv_peaks$chr = as.character(surv_peaks$chr)

assocs = vector('list', nrow(surv_peaks))
names(assocs) = as.character(surv_peaks$chr)

# For each peak...
for(i in 1:nrow(surv_peaks)) {
  
  chr   = surv_peaks$chr[i]
  start = surv_peaks$ci_lo[i]
  end   = surv_peaks$ci_hi[i]
  print(paste('chr', chr))
  
  width = end - start
  if(width > 10) {
    start = max(3, surv_peaks$pos[i] - 5)
    end   = surv_peaks$pos[i] + 5
  } # if(width > 10)
  
  # Get SNP probs for current samples.
  snpinfo = snp_func(chr, start, end)
  snpinfo = subset(snpinfo, type == 'snp')
  snpinfo = qtl2:::index_snps(map, snpinfo)
  snppr   = genoprob_to_snpprob(probs, snpinfo)[[1]][,'B',]
  
  print(paste('   ', ncol(snppr), 'snps'))
  
  lod_chr  = matrix(0, nrow = ncol(snppr), ncol = 1, dimnames = list(colnames(snppr), "survival"))
  
  # fit Cox-PH model at each SNP.
  for(j in 1:ncol(snppr)) {
    
    if(j %% 100 == 0) print(j)
    
    mod = coxph(pheno_surv$surv ~ addcovar + snppr[,j])
    lod_chr[j,1] = mod$loglik[2]
    
  } # for(j)
  
  assocs[[chr]]  = (lod_chr - null_ll) / log(10)
  class(assocs[[chr]]) = c('scan1', 'matrix')
  
  # Get genes in interval.
  genes = gene_func(chr = chr, start = start, end = end)
  
  # Get GO categories from Biomart.
  bm = getBM(attributes = c('mgi_id', 'ensembl_gene_id', 'go_id', 'name_1006',  'namespace_1003'), 
             filters = 'mgi_id', 
             values  = genes$gene_id, 
             mart = ensembl)
  bm = subset(bm, namespace_1003 %in% c('biological_process', 'molecular_function'))
  keywords = paste(c('macrophage', 'cytokine', 'interferon', 'interleukin', 'immune', 'Gram', 
                     'antigen', 'T cell', 'B cell', 'virus', 'viral'),
                   collapse = '|')
  bm = subset(bm, base::grepl(keywords, bm$name_1006))
  genes$color = ifelse(genes$gene_id %in% bm$mgi_id, 'red', 'black')
  genes = genes[order(genes$start, genes$stop),]
  
  png(file.path(result_dir, paste0('survival_assoc_chr', chr, '.png')), 
      width = 2400, height = 1600, res = 300)
  plot_snpasso(assocs[[chr]], snpinfo, genes = genes, drop_hilit = 1, 
               top_panel_prop = 0.3, main = 'Survival', colors = genes$color)
  dev.off()
  
  rm(snppr, lod_chr)
  gc()
  
} # for(i)

save(lod, assocs, file = file.path(result_dir, 'survival_assoc.Rdata'))


##############
# Permutations.

perm_file = file.path(result_dir, "perms.rds")
perms = NULL

if(file.exists(file.path(result_dir, "perms.rds"))) {
  
  perms = readRDS(perm_file)
  
} else {
  
  pheno_name = "lung_cxcl5"
  
  samples2use = which(!is.na(pheno_rz[,pheno_name]))
  
  print(paste0(pheno_name, " : ", length(samples2use)))
  
  tmpK = K
  for(j in 1:length(K)) {
    tmpK[[j]] = K[[j]][samples2use, samples2use]
  }
  tmp_covar = addcovar[samples2use,,drop = FALSE]
  tmp_covar = tmp_covar[,colSums(tmp_covar) > 0,drop = FALSE]
  
  
  perms = qtl2::scan1perm(genoprobs = probs[samples2use,], 
                          pheno     = pheno_rz[samples2use, pheno_name, drop = FALSE], 
                          kinship   = tmpK, 
                          addcovar  = tmp_covar, 
                          cores     = 4,
                          n_perm    = 1000,
                          quiet     = FALSE) 
  saveRDS(perms, file = perm_file)
  
  rm(tmpK, tmp_covar, samples2use, pheno_name)
  
} # else

# Redraw QTL plots with threshold.
perms = readRDS(perm_file)
rdata_files = dir(result_dir, pattern = 'Rdata$')
rdata_files = rdata_files[-grep('survival', rdata_files)]

for(f in rdata_files) {
  # This loads in 'lod', 'coefs', & 'assocs'.
  load(file.path(result_dir, f))
  pheno_name = sub('_qtl2.Rdata$', '', f)
  png(file.path(result_dir, sub('2.Rdata$', '.png', f)), width = 1000, height = 800, res = 128)
  plot(lod, map, main = pheno_name)
  abline(h = quantile(perms, probs = c(0.95, 0.80)), col = c('red', 'orange'), lwd = 2)
  dev.off()
} # for(f)

##################################################################################
# Survival permutations.
nperm = 1000

surv_perms = rep(0, nperm)

for(p in 1:nperm) {
  
  print(paste('perm', p))
  
  # Permute the samples.
  new_order  = sample(1:nrow(pheno_surv))
  pheno_surv = pheno_surv[new_order,]
  addcovar   = addcovar[new_order,]
  
  # Null model log_likelihood.
  null_ll = coxph(pheno_surv$surv ~ addcovar)$loglik[2]
  
  # For each chromosome...
  for(chr in names(probs)) {
    
    print(paste('   chr', chr))
    
    # Get genoprobs for current samples.
    pr = probs[[chr]]
    
    lod_chr = rep(0, nrow = dim(pr)[3])
    
    for(j in 1:dim(pr)[3]) {
      
      lod_chr[j] = coxph(pheno_surv$surv ~ addcovar + pr[,,j])$loglik[2]
      
    } # for(j)
    
    max_lod = (max(lod_chr) - null_ll) / log(10)
    
    surv_perms[p] = max(surv_perms[p], max_lod)
    
  } # for(chr)
} # for(p)

saveRDS(surv_perms, file = file.path(result_dir, 'perms_survival.rds'))

# Plot the survival QTL with the threshold.
surv_perms = readRDS(file = file.path(result_dir, 'perms_survival.rds'))
# This load in 'lod' & 'coefs'
load(file = file.path(result_dir, 'survival_qtl2.Rdata'))

png(file.path(result_dir, "survival_qtl.png"), width = 1000, height = 800, res = 128)
plot(lod, map, main = 'Survival')
abline(h = quantile(surv_perms, probs = 0.95), col = 'red')
dev.off()


####################
# Remap s100a8 with the SNP at the peak QTL regressed out.

i = which(colnames(pheno_rz) == "lung_s100a8")
pheno_name = colnames(pheno_rz)[i]

samples2use = which(!is.na(pheno_rz[,pheno_name]))

print(str_c(pheno_name, " : ", length(samples2use)))

tmpK = K
for(j in 1:length(K)) {
  tmpK[[j]] = K[[j]][samples2use, samples2use]
}

# New covars.
load(file.path(result_dir, paste0(pheno_name, "_qtl2.Rdata")))
top = top_snps(assocs[['3']]$lod, assocs[['3']]$snpinfo)
top = top[top$lod == max(top$lod),]
s100a8_gt = genoprob_to_snpprob(genoprobs = probs, snpinfo = assocs[['3']]$snpinfo)[[1]]
s100a8_gt = as.matrix(s100a8_gt[,1,top$snp_id[1]])

tmp_covar = merge(addcovar, s100a8_gt, by = 'row.names')
rownames(tmp_covar) = tmp_covar$Row.names
tmp_covar = tmp_covar[,-1]

tmp_covar = tmp_covar[samples2use,,drop = FALSE]
tmp_covar = tmp_covar[,colSums(tmp_covar) > 0, drop = FALSE]

print("    QTL")  
lod = qtl2::scan1(genoprobs = probs[samples2use,], 
                  pheno     = pheno_rz[samples2use, pheno_name, drop = FALSE], 
                  kinship   = tmpK, 
                  addcovar  = tmp_covar, 
                  cores = 10)

png(file.path(result_dir, str_c(pheno_name, "_regress_gt_qtl.png")), width = 1000, height = 800, res = 128)
qtl2::plot_scan1(lod, map = map, main = str_c(pheno_name, " genotype regressed out"))
dev.off()

####################
# Remap serum Cxcl5 with the SNP at the peak QTL regressed out.

i = which(colnames(pheno_rz) == "serum_cxcl5")
pheno_name = colnames(pheno_rz)[i]
chr = '5'

samples2use = which(!is.na(pheno_rz[,pheno_name]))

print(str_c(pheno_name, " : ", length(samples2use)))

tmpK = K
for(j in 1:length(K)) {
  tmpK[[j]] = K[[j]][samples2use, samples2use]
}

# New covars.
load(file.path(result_dir, paste0(pheno_name, "_qtl2.Rdata")))
top = top_snps(assocs[[chr]]$lod, assocs[[chr]]$snpinfo)
top = top[top$lod == max(top$lod),]
cxcl5_gt = genoprob_to_snpprob(genoprobs = probs, snpinfo = assocs[[chr]]$snpinfo)[[1]]
cxcl5_gt = as.matrix(cxcl5_gt[,1,top$snp_id[1]])

tmp_covar = merge(addcovar, cxcl5_gt, by = 'row.names')
rownames(tmp_covar) = tmp_covar$Row.names
tmp_covar = tmp_covar[,-1]

tmp_covar = tmp_covar[samples2use,,drop = FALSE]
tmp_covar = tmp_covar[,colSums(tmp_covar) > 0, drop = FALSE]

# Get the baseline assoc plot.
start = 75
end   = 93
assoc0 = qtl2::scan1snps(genoprobs = probs[samples2use,], 
                         map       = map,
                         pheno     = pheno_rz[samples2use, pheno_name, drop = FALSE], 
                         kinship   = tmpK, 
                         addcovar  = addcovar, 
                         query_func = snp_func,
                         chr = chr, start = start, end = end,
                         cores = 4)

png(file.path(result_dir, str_c(pheno_name, "_regress_gt_qtl_baseline.png")), width = 1000, height = 800, res = 128)
qtl2::plot_snpasso(assoc0$lod, assoc0$snpinfo, main = str_c(pheno_name),
                   show_all_snps = TRUE, minlod = 0.1)
usr = par('usr')
dev.off()

# Regressing out SNP around 89.5 Mb (NOD & CAST)
assoc1 = qtl2::scan1snps(genoprobs = probs[samples2use,],
                         map       = map,
                         pheno     = pheno_rz[samples2use, pheno_name, drop = FALSE], 
                         kinship   = tmpK, 
                         addcovar  = tmp_covar, 
                         query_func = snp_func,
                         chr = chr, start = start, end = end,
                         cores = 4)

png(file.path(result_dir, str_c(pheno_name, "_regress_gt_qtl_89mb.png")), width = 1000, height = 800, res = 128)
qtl2::plot_snpasso(assoc1$lod, assoc1$snpinfo, main = str_c(pheno_name, " genotype regressed out"),
                   show_all_snps = TRUE, ylim = usr[3:4])
abline(v = 89.5, col = 'red', lwd = 2)
dev.off()

# Regressing out SNP around 8 Mb (NOD & CAST & 129)
top = top_snps(assoc0$lod, assoc0$snpinfo, drop = 2.5)
top = subset(top, pos <= 76)
top = top[top$lod == max(top$lod),]
cxcl5_gt = genoprob_to_snpprob(genoprobs = probs, snpinfo = assocs[[chr]]$snpinfo)[[1]]
cxcl5_gt = as.matrix(cxcl5_gt[,1,top$snp_id[1]])

tmp_covar = merge(addcovar, cxcl5_gt, by = 'row.names')
rownames(tmp_covar) = tmp_covar$Row.names
tmp_covar = tmp_covar[,-1]

tmp_covar = tmp_covar[samples2use,,drop = FALSE]
tmp_covar = tmp_covar[,colSums(tmp_covar) > 0, drop = FALSE]

assoc2 = qtl2::scan1snps(genoprobs = probs[samples2use,],
                         map       = map,
                         pheno     = pheno_rz[samples2use, pheno_name, drop = FALSE], 
                         kinship   = tmpK, 
                         addcovar  = tmp_covar, 
                         query_func = snp_func,
                         chr = chr, start = start, end = end,
                         cores = 4)

png(file.path(result_dir, str_c(pheno_name, "_regress_gt_qtl_77mb.png")), width = 1000, height = 800, res = 128)
qtl2::plot_snpasso(assoc2$lod, assoc2$snpinfo, main = str_c(pheno_name, " genotype regressed out"),
                   show_all_snps = TRUE, ylim = usr[3:4])
abline(v = 77.13, col = 'red', lwd = 2)
dev.off()



####################
# Create QTL heatmap
rdata_files = dir(result_dir, pattern = 'Rdata$')
qtl = NULL

for(f in rdata_files) {
  tmp = load(file.path(result_dir, f))
  if(is.null(qtl)) {
    qtl = lod
  } else {
    qtl = cbind(qtl, lod)
  }
} # for(f)

# Cluster QTL.
cl = hclust(as.dist(1.0 - cor(qtl)), method = "average")
qtl = qtl[,cl$order]

# Synch up marker names. This should have been done in setup scripts.
for(i in seq_along(map)) {
  map[[i]] = map[[i]][names(map[[i]]) %in% rownames(qtl)]
} # for(i)

qtl = cbind(marker = unlist(sapply(map, names)), 
            chr    = rep(names(map), sapply(map, length)), 
            pos    = unlist(map), 
            qtl)

qtl = qtl %>% 
  as.data.frame() %>% 
  mutate(chr = factor(chr, levels = c(1:19, 'X'))) %>%
  pivot_longer(cols = lung_il10:lung_cxcl2, names_to = 'pheno', values_to = 'lod') %>%
  mutate(lod   = as.numeric(lod),
         pos   = as.numeric(pos)) %>% 
  pivot_wider(names_from = pheno, values_from = lod)

# All QTL plot.
png(file.path(result_dir, "tb_all_qtl.png"), width = 4000, height = 3000, res = 256)
qtl %>%
  pivot_longer(cols = lung_il10:lung_cxcl2, names_to = 'pheno', values_to = 'lod') %>%
  mutate(pheno = factor(pheno, levels = cl$labels[cl$order])) %>%
  ggplot(aes(pos, lod)) +
  geom_line() +
  geom_abline(slope = 0, intercept = 7.3, color = 'red') +
  geom_abline(slope = 0, intercept = 6.0, color = 'orange') +
  facet_grid(pheno~chr, scales = 'free_x') +
  theme(panel.spacing.x = unit(0, 'lines'),
        panel.spacing.y = unit(0.1, 'lines'),
        axis.text.x = element_text(angle = 90, vjust = 0))
dev.off()

# QTL heatmap function.
qtl_heatmap = function(lod) {
  
  # Split markers and LOD.
  mkr = lod[,1:3]
  lod = as.matrix(lod[,-(1:3)])
  
  # Cluster phenotypes.
  cl = hclust(as.dist(1.0 - cor(lod)), method = 'average')
  
  # Get unique choromosomes.
  unique_chr = distinct(qtl, chr) %>% pull(chr)
  
  # Apportion 1000 markers proportionally to each chromosome.
  nm = 1000
  max_pos  = sapply(split(mkr$pos, mkr$chr), max)
  prop_chr = max_pos / sum(max_pos)
  num_mkr  = round(nm * prop_chr)
  # Add a few to get up to 1000.
  num_mkr[1]   = num_mkr[1] + 1
  num_mkr[2]   = num_mkr[2] + 1
  num_mkr['X'] = num_mkr['X'] + 1
  
  # Interpolate a new marker map.
  new_map = mapply(function(m, num) { approx(x = m, n = num)$y }, m = map, num = num_mkr)
  
  # Interpolated QTL map and LODs.
  new_qtl = NULL
  
  for(chr in names(new_map)) {
    
    # Get breakpoints and midpoints.
    brks = cut(map[[chr]], new_map[[chr]])
    mids = sapply(split(map[[chr]], brks), mean)
    
    # Get LOD for this chromosome.
    curr_lod = lod[mkr$chr == chr,]
    
    # Make new data frame for this chromosome.
    tmp = data.frame(chr = rep(chr, length(mids)),
                     pos = mids,
                     matrix(0, length(mids), ncol(curr_lod), dimnames = list(NULL, colnames(curr_lod))))
    
    for(j in colnames(curr_lod))  {
      tmp[[j]] = sapply(split(curr_lod[,j], brks), max)
    } # for(j)
    
    new_qtl = rbind(new_qtl, tmp)
    
  } # for(chr)
  
  colors = c('black', 'grey15', 'grey30', 'grey50', 'red', 'orange', 'yellow')
  
  new_qtl %>%
    mutate(chr = factor(chr, levels = c(1:19, 'X'))) %>%
    gather(pheno, lod, -(chr:pos)) %>%
    mutate(pheno = factor(pheno, levels = cl$labels[cl$order])) %>%
    ggplot() +
    geom_tile(aes(x = pos, y = 1, color = lod, fill = lod), width = 5) +
    scale_color_gradientn(colors = colors, values = scales::rescale(c(0, 5, 6, 7, 8, max(lod) + 0.1), to = c(0,1))) +
    scale_fill_gradientn(colors  = colors, values = scales::rescale(c(0, 5, 6, 7, 8, max(lod) + 0.1), to = c(0,1))) +
    facet_grid(pheno ~ chr, scales = 'free_x') +
    theme(panel.spacing = unit(0, 'lines'),
          axis.text.x = element_text(angle = 90, vjust = 0),
          axis.text.y = element_blank(),
          axis.ticks  = element_blank(),
          axis.title  = element_blank(),
          strip.text.x = element_text(size = 24),
          strip.text.y = element_text(size = 12, angle = 0))
  
} # qtl_heatmap()

png(file.path(result_dir, "tb_qtl_heatmap.png"), width = 4000, height = 2000, res = 256)
qtl_heatmap(qtl)
dev.off()

### QTL plot correlation.
qtl_files = dir(results_dir, pattern = '_qtl2.Rdata$')
qtl_files = qtl_files[-grep('PC', qtl_files)]

qtlmat = NULL

for(f in qtl_files) {
  load(file.path(results_dir, f))
  if(is.null(qtlmat)) {
    qtlmat = lod
  } else {
    qtlmat = cbind(qtlmat, lod)
  } # else
} # for(f)

colnames(qtlmat) = sub('_qtl2\\.Rdata$', '', qtl_files)
qtlmat = qtlmat[,colnames(qtlmat) != 'survival_60day']
#qtlmat = qtlmat / matrix(apply(qtlmat, 2, max), nrow(qtlmat), ncol(qtlmat), byrow = T)
#qtlmat[qtlmat < 0.4] = NA
qtlcor = cor(qtlmat, use = 'pairwise')
qltcl = hclust(as.dist(1.0 - qtlcor), method = 'average')
qtlcor = qtlcor[rev(qltcl$order), rev(qltcl$order)]

png(file.path(results_dir, 'qtl_cor_heatmap.png'), width = 2000, height = 2000, res = 300)
heatmap(qtlcor)
dev.off()

png(file.path(results_dir, 'qtl_corrplot.png'), width = 2000, height = 2000, res = 200)
corrplot.mixed(qtlcor, lower = 'number', upper = 'ellipse', order = 'original', tl.pos = 'lt')
dev.off()


phenocor = cor(pheno_rz, use = 'pairwise')
heatmap(phenocor, symm = T)

