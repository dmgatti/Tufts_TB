##############################
# Tufts TB QTL Mapping.
# DMG
# Feb. 23, 2020
##############################
options(stringsAsFactors = FALSE)
library(readxl)
library(tidyverse)
library(qtl2)

base_dir   = '/media/dmgatti/hdb/projects/TB'
input_file = file.path(base_dir, 'data', 'tufts_do_tb_qtl2_input.Rdata') 
result_dir = file.path(base_dir, 'results', 'qtl2', 'gen_factor')

# Read in the DOQTL formatted data.
load(input_file)

# CC SNP and gene database files (from Karl).
ccsnpdb = "/media/dmgatti/hda/data/MUGA/cc_variants.sqlite"
mgidb   = "/media/dmgatti/hda/data/MUGA/mouse_genes_mgi.sqlite"
snp_func  = create_variant_query_func(dbfile = ccsnpdb)
gene_func = create_gene_query_func(dbfile = mgidb)
ensembl = 93

peaks = NULL

# Read in image data.
img_data = read_xlsx(path = file.path(base_dir, 'data', 'image_measurements.xlsx'))
colnames(img_data) = c('mouse', 'prSS', 'necr')
plot(necr ~ prSS, data = img_data)

# Should probably be doing logistic regression this one...
qqnorm(log(img_data$prSS / (1 - img_data$prSS)))
qqline(log(img_data$prSS / (1 - img_data$prSS)))

qqnorm(img_data$necr)
qqline(img_data$necr)

# Transform necrosis ratio to rankZ.
rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

img_data$necr_rz = rankZ(img_data$necr)
img_data$necr_log = log(img_data$necr + 0.0009)

# Add necrosis ratio to the phenotypes.
pheno = data.frame(mouse = rownames(pheno), pheno)
pheno = merge(pheno, img_data, by = 'mouse', all.x = TRUE)

rownames(pheno) = pheno$mouse
pheno = pheno[,!colnames(pheno) %in% c('mouse', 'prSS', 'necr')]
pheno = as.matrix(pheno)


##########
# QTL Mapping with coef calculations.

addcovar = model.matrix(~gen, data = covar)[,-1,drop = FALSE]

for(i in 1:ncol(pheno)) {
  
  pheno_name = colnames(pheno)[i]
  
  samples2use = which(!is.na(pheno[,pheno_name]))
  
  print(paste0(pheno_name, " : ", length(samples2use)))
  
  tmpK = K
  for(j in 1:length(K)) {
    tmpK[[j]] = K[[j]][samples2use, samples2use]
  }
  tmp_covar = addcovar[samples2use,,drop = FALSE]
  tmp_covar = tmp_covar[,colSums(tmp_covar) > 0, drop = FALSE]

  print("    QTL")  
  lod = qtl2::scan1(genoprobs = probs[samples2use,], 
                    pheno     = pheno[samples2use, pheno_name, drop = FALSE], 
                    kinship   = tmpK, 
                    addcovar  = tmp_covar, 
                    cores = 10)
  
  png(file.path(result_dir, paste0(pheno_name, "_qtl.png")), width = 1000, height = 800, res = 128)
  qtl2::plot_scan1(lod, map = map, main = pheno_name)
  dev.off()
  
  current_peaks = find_peaks(lod, map, threshold = 6, prob = 0.95)
  
  # If the confidence interval takes up more than 20 Mb, use the peak LOD +/- 5 Mb.
  coefs  = NA
  assocs = NA
  if(nrow(current_peaks) > 0) {
    
    for(j in 1:nrow(current_peaks)) {
      if(current_peaks$ci_hi[j] - current_peaks$ci_lo[j] > 20) {
        current_peaks$ci_hi[j] = current_peaks$pos[j] + 5
        current_peaks$ci_lo[j] = max(current_peaks$pos[j] - 5, 3)  # This protects us from going off the beginning of the chr.
      } # if(current_peaks$ci_hi[j] - current_peaks$ci_lo[j])
    } # for(j)
    
    peaks = rbind(peaks, current_peaks)
  
    coefs = vector("list", nrow(current_peaks))
    names(coefs) = current_peaks$chr
    assocs = vector("list", nrow(current_peaks))
    names(assocs) = current_peaks$chr
  
    for(j in 1:nrow(current_peaks)) {
    
      curr_chr = current_peaks$chr[j]
      print(paste0("    chr ", curr_chr))
      start = current_peaks$ci_lo[j]
      end   = current_peaks$ci_hi[j]
  
      print(paste("    coef", j))
      coefs[[j]] = qtl2::scan1blup(genoprobs = probs[samples2use, curr_chr], 
                            pheno     = pheno[samples2use, pheno_name, drop = FALSE], 
                            kinship   = K[[curr_chr]], 
                            addcovar  = tmp_covar,
                            se = TRUE,
                            cores = 10)
      stopifnot(!is.nan(coefs[[j]]))
    
      png(file.path(result_dir, paste0(pheno_name, "_coef_chr", curr_chr, ".png")), width = 1000, height = 800, res = 128)
      plot_coefCC(coefs[[j]], map, scan1_output = lod, top_panel_prop = 0.6, main = pheno_name)
      dev.off()
    
      print(paste("    assoc", j))
      assocs[[j]] = scan1snps(genoprobs = probs[,curr_chr], pheno = pheno[,pheno_name, drop = FALSE], kinship = tmpK[[curr_chr]], addcovar = tmp_covar,
                        map = map, chr = curr_chr, start = max(1, start - 1), end = end + 1, query_func = snp_func, 
                        keep_all_snps = TRUE, cores = 10)
    
      genes = gene_func(chr = curr_chr, start = start - 1, end = end + 1)
      png(file.path(result_dir, paste0(pheno_name, "_assoc_chr", curr_chr, ".png")), width = 2400, height = 1600, res = 300)
      plot_snpasso(scan1output = assocs[[j]]$lod, snpinfo = assocs[[j]]$snpinfo, genes = genes, drop_hilit = 1, 
                   top_panel_prop = 0.3, main = pheno_name, colors = 'black')
      dev.off()
    
      # Get the top SNPs.
      top = top_snps(assocs[[j]]$lod, assocs[[j]]$snpinfo)
      write_csv(top, path = file.path(result_dir, paste0(pheno_name, "_assoc_chr", curr_chr, "_top_snps.csv")))

    } # for(j)
  
  } # if(nrow(current_peaks > 0))

  save(lod, coefs, assocs, file = file.path(result_dir, paste0(pheno_name, "_qtl2.Rdata")))
  
} # for(i)

write_csv(peaks, path = file.path(result_dir, "tb_qtl_peaks.csv"))

#############
# GWAS
for(i in 1:ncol(pheno)) {
  
  pheno_name = colnames(pheno)[i]
  
  samples2use = which(!is.na(pheno[,pheno_name]))
  
  print(paste0(pheno_name, " : ", length(samples2use)))
  
  tmpK = K
  for(j in 1:length(K)) {
    tmpK[[j]] = K[[j]][samples2use, samples2use]
  }
  tmp_covar = addcovar[samples2use,,drop = FALSE]
  tmp_covar = tmp_covar[,colSums(tmp_covar) > 0, drop = FALSE]
  
  ts = NULL
  
  pdf(file.path(result_dir, paste0(pheno_name, "_gwas.pdf")))
  
  for(chr in names(probs)) {
  
    print(chr)
    assoc = scan1snps(genoprobs = probs[,chr], 
                      map       = map,
                      pheno     = pheno[,pheno_name, drop = FALSE],
                      kinship   = tmpK[[chr]],
                      addcovar  = tmp_covar,
                      chr       = chr,
                      start     = 0,
                      end       = 200,
                      query_func = snp_func,
                      keep_all_snps = TRUE,
                      cores      = 4)
    
      plot_snpasso(assoc$lod, assoc$snpinfo, main = paste(pheno_name, ' : CHR', chr))
  
      ts = rbind(ts, top_snps(assoc$lod, assoc$snpinfo, chr = chr, drop = 0.5, show_all_snps = TRUE))
      
   } # for(chr)
  
  dev.off()
  
  write_csv(ts, file.path(result_dir, paste0(pheno_name, "_gwas.csv")))
  
  rm(tmpK, tmp_covar)
  gc()
  
} # for(i)



##############
# Permutations.

pheno_name = "cxcl5"

samples2use = which(!is.na(pheno_rz[,pheno_name]))

print(paste0(pheno_name, " : ", length(samples2use)))

tmpK = K
for(j in 1:length(K)) {
  tmpK[[j]] = K[[j]][samples2use, samples2use]
}
tmp_covar = covar[samples2use,,drop = FALSE]
tmp_covar = tmp_covar[,colSums(tmp_covar) > 0,drop = FALSE]

perms = qtl2::scan1perm(genoprobs = probs[samples2use,], 
                        pheno = pheno_rz[samples2use,pheno_name, drop = FALSE], 
                        kinship = tmpK, 
                        addcovar = tmp_covar, 
                        cores = 10,
                        n_perm = 1000)
saveRDS(perms, file = file.path(result_dir, "perms.rds"))


####################
# Remap s100a8 with the SNP at the peak QTL regressed out.

i = which(colnames(pheno_rz) == "s100a8")
pheno_name = colnames(pheno_rz)[i]

samples2use = which(!is.na(pheno_rz[,pheno_name]))

print(paste0(pheno_name, " : ", length(samples2use)))

tmpK = K
for(j in 1:length(K)) {
  tmpK[[j]] = K[[j]][samples2use, samples2use]
}

# New covars.
load(paste0(result_dir, pheno_name, "_qtl2.Rdata"))
top = top_snps(assocs[['3']]$lod, assocs[['3']]$snpinfo)
top = top[top$lod == max(top$lod),]
s100a8_gt = genoprob_to_snpprob(genoprobs = probs, snpinfo = assocs[['3']]$snpinfo)[[1]]
s100a8_gt = s100a8_gt[,1,top$snp_id[1]]

tmp_covar = cbind(covar, s100a8_gt)

tmp_covar = tmp_covar[samples2use,,drop = FALSE]
tmp_covar = tmp_covar[,colSums(tmp_covar) > 0, drop = FALSE]

print("    QTL")  
lod = qtl2::scan1(genoprobs = probs[samples2use,], 
                  pheno     = pheno_rz[samples2use, pheno_name, drop = FALSE], 
                  kinship   = tmpK, 
                  addcovar  = tmp_covar, 
                  cores = 10)

png(file.path(result_dir, paste0(pheno_name, "_regress_gt_qtl.png")), width = 1000, height = 800, res = 128)
qtl2::plot_scan1(lod, map = map, main = paste0(pheno_name, "genotype regressed out"))
dev.off()


####################
# Create QTL heatmap
rdata_files = dir(result_dir, pattern = 'Rdata$')
qtl = NULL

for(f in rdata_files) {
  tmp = load(paste0(result_dir, f))
  if(is.null(qtl)) {
    qtl = lod
  } else {
    qtl = cbind(qtl, lod)
  }
} # for(f)

# Cluster QTL.
cl = hclust(as.dist(1.0 - cor(qtl)), method = "average")
qtl = qtl[,cl$order]
qtl = cbind(markers[,1:3], qtl)

# All QTL plot.
png(file.path(result_dir, "tb_all_qtl.png"), width = 4000, height = 3000, res = 256)
qtl %>%
  mutate(chr = factor(chr, levels = c(1:19, 'X'))) %>%
  gather(pheno, lod, -(marker:pos)) %>%
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
  max_pos  = sapply(map, max)
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
      scale_color_gradientn(colors = colors, values = scales::rescale(c(0, 5.5, 6, 6.5, 7.0, max(lod) + 0.1), to = c(0,1))) +
      scale_fill_gradientn(colors  = colors, values = scales::rescale(c(0, 5.5, 6, 6.5, 7.0, max(lod) + 0.1), to = c(0,1))) +
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

### Lung Damage
