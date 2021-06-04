

# Add thresholds to plots.
# rds_files = dir(result_dir, pattern = 'rds$')
# rds_files = rds_files[-grep('perms', rds_files)]
# for(f in rds_files) {
#   
#   gwas = readRDS(file.path(result_dir, f))
#   
#   png(file.path(result_dir, sub('rds$', 'png', f)), width = 2000, height = 1000, res = 128)
#   plot_snpasso(gwas$lod, gwas$snpinfo, main = sub('_gwas\\.rds$', '', f))
#   abline(h = thr, col = 'red', lwd = 2)
#   dev.off()
#   
# } # for(f)



# Run detailed GWAS at peaks with LOD > 5.5
top_snp_files = dir(result_dir, pattern = 'chr.\\.csv$')

for(f in top_snp_files) {
  
  ts = readr::read_csv(file.path(result_dir, f))
  
  pheno_name = str_replace(f, '_gwas_chr.\\.csv$', '')
  
  samples2use = which(!is.na(pheno_rz[,pheno_name]))
  
  print(paste0(pheno_name, " : ", length(samples2use)))
  
  kinship = K
  for(j in 1:length(kinship)) {
    kinship[[j]] = kinship[[j]][samples2use, samples2use]
  }
  addcovar_tmp = addcovar[samples2use,,drop = FALSE]
  addcovar_tmp = addcovar_tmp[,colSums(addcovar_tmp) > 0, drop = FALSE]


  chr = str_split(f, pattern = '_')[[1]]
  chr = chr[str_detect(chr, 'chr')]
  chr = str_replace_all(chr, 'chr|\\.csv$', '')
  
  ctr = ts %>% 
          slice_max(lod, n = 1) %>% 
          summarize(pos = median(pos)) %>% 
          pull(pos)
  
  assoc = scan1snps(genoprobs = probs[samples2use, chr], 
                    map       = map,
                    pheno     = pheno_rz[samples2use, pheno_name, drop = FALSE],
                    kinship   = kinship[[chr]],
                    addcovar  = addcovar_tmp,
                    chr       = chr,
                    start     = 0,
                    end       = 200,
                    query_func = snp_func,
                    keep_all_snps = TRUE,
                    cores     = 1)

} # for(f)


  

