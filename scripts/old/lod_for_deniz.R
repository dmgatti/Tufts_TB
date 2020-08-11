base_dir  = '/media/dmgatti/hdb/projects/TB/results/qtl2/dose_as_factor/'

files = dir(base_dir, patter = '_qtl2.Rdata', full.names = TRUE)

snp_url = 'http://csbio.unc.edu/MUGA/snps.gigamuga.Rdata'
load(url(snp_url))

result = NULL

for(f in files) {
  load(f)
  if(is.null(result)) {
    result = lod
  } else {
    result = data.frame(result, lod)
  }
}
result = data.frame(marker = rownames(result), result)

result = merge(snps[,c(2,1,4)], result, by = 'marker', all.x = FALSE, all.y = TRUE)
result = result[order(result$chr, result$pos),]
result = result[,-grep('^PC', colnames(result))]

write.csv(result, file = paste0(base_dir, '../tufts_tb_qtl2_lod.csv'), 
          row.names = FALSE, quote = FALSE)

perms = readRDS(paste0(base_dir, 'perms.rds'))
write.csv(perms, file = paste0(base_dir, '../tufts_tb_qtl2_perms.csv'),
          row.names = FALSE, quote = FALSE)


