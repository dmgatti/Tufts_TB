################################################################################
# Tufts TB QTL Harvest.
# DMG
# May 20, 2020
##############################

# Hack to fix not being able to load org.Mm.eg.db.
# See https://github.com/rstudio/rstudio/issues/9219
options(connectionObserver = NULL)

library(biomaRt)
library(survival)
library(Mus.musculus)
library(qtl2)

base_dir   = '/media/dmgatti/hdb/projects/TB'
input_file = file.path(base_dir, 'data', 'tufts_do_tb_qtl2_input.Rdata') 
expr_file  = file.path(base_dir, 'data', 'expression', 'tb_expr.csv')
result_dir = file.path(base_dir, 'results', 'qtl2', 'gen_factor2')
peak_file  = file.path(result_dir, 'tb_qtl_peaks.csv')
perm_file  = file.path(result_dir, 'perms.rds')

# Set up biomaRt.
ensembl = useEnsembl(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl')
#bmfilt = listFilters(ensembl)
#bmattr = listAttributes(ensembl)

tx = TxDb.Mmusculus.UCSC.mm10.knownGene

# Read in the DOQTL formatted data.
load(input_file)

# CC SNP and gene database files (from Karl).
ccsnpdb = "/media/dmgatti/hda/data/MUGA/cc_variants.sqlite"
mgidb   = "/media/dmgatti/hda/data/MUGA/mouse_genes_mgi.sqlite"
snp_func  = create_variant_query_func(dbfile = ccsnpdb)
gene_func = create_gene_query_func(dbfile = mgidb)
ensembl = 93

addcovar = model.matrix(~gen, data = covar)[,-1,drop = FALSE]

# Read expression data.
expr = read.csv(expr_file)
colnames(expr)[1] = 'entrez'
expr$entrez = as.character(expr$entrez)

# Get the genes that we can find from the Bioconductor TxDB.
annot = select(tx, keys = expr$entrez, columns = c('TXNAME', 'TXCHROM', 'TXSTART', 'TXEND'), keytype = 'GENEID')
annot$TXCHROM = sub('^chr', '', annot$TXCHROM)

expr = expr[expr$entrez %in% annot$GENEID,]

# Read QTL peaks.
peaks = read.csv(peak_file)
peaks = GRanges(seqnames = peaks$chr, 
                ranges = IRanges(start = peaks$ci_lo * 1e6, end = peaks$ci_hi * 1e6),
                mcols  = peaks[,c(2, 5:15)])
colnames(mcols(peaks)) = sub('^mcols\\.', '', colnames(mcols(peaks)))
perms = readRDS(perm_file)
thr = quantile(perms, probs = 0.95)

# Keep peaks with lod > 7.76.
top_peaks = subset(peaks, lod >= thr)

# For each top peak, see if there are other QTL at lower thresholds within 
# the QTL interval.
for(i in 1:length(top_peaks)) {
  
  pheno_name = top_peaks$lodcolumn[i]
  chr   = as.character(seqnames(top_peaks)[i])
  start = start(top_peaks)[i]
  end   = end(top_peaks)[i]

  tmp = subsetByOverlaps(peaks, top_peaks[i])
  
} # for(i)



