options(stringsAsFactors = FALSE)
library(AnnotationHub)
ens_ver = 93
hub = AnnotationHub()

# AnnotationHub timed out on me.
#ens_93 = names(hub)[hub$title == 'Mus_musculus.GRCm38.93.gtf']
#ensembl = hub[[ens_93]]

# Downloaded by hand. from http://ftp.ensembl.org/pub/current_gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz .
ensembl = read.table('/home/dmgatti/Documents/data/ensembl/Mus_musculus.GRCm38.93.gtf.gz',
                     sep = '\t', skip = 5)
ensembl = subset(ensembl, V3 == 'gene')

# Process column 9.
c9 = strsplit(ensembl[[9]], split = ';')
c9 = lapply(c9, function(z) {trimws(z)})

cn = sapply(strsplit(c9[[1]], ' '), '[', 1)
c9 = lapply(c9, strsplit, split = ' ')
c9 = lapply(c9, unlist)
c9 = lapply(c9, function(z) { z[c(2, 4, 6, 8, 10)] })
new_c9 = matrix(unlist(c9), nrow = nrow(ensembl), ncol = length(cn), dimnames = list(NULL, cn),
                byrow = TRUE)
rm(c9)
ensembl = data.frame(ensembl[,1:8], new_c9)

# Add gene symbols to association mapping output.
target_dir = '/home/dmgatti/Documents/data/TB/results/qtl2/dose_as_factor'
suffix     = '_top_snps.csv$'

files = dir(target_dir, pattern = suffix, full.names = TRUE)

for(f in files) {
  snps = read.csv(f)
  left = snps[,1:6]
  right = snps[,7:ncol(snps)]
  symbols = ensembl$gene_name[base::match(snps$ensembl_gene, ensembl$gene_id)]
  snps = data.frame(left, symbol = symbols, right)
  write.csv(snps, file = f, quote = FALSE, row.names = FALSE)
} # for(f)


