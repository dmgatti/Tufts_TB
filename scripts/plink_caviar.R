# PLINK/CAVIAR analysis.

qtl_util_file = '/media/dmgatti/hdb/projects/TB/scripts/qtl2_utilities.R'
data_file     = '/media/dmgatti/hdb/projects/TB/data/tufts_do_tb_qtl2_input.Rdata'
assoc_file    = '/media/dmgatti/hdb/projects/TB/manuscripts/manuscript1/figures/PC1_assoc_chr1.rds'
out_prefix    = '/media/dmgatti/hdb/projects/TB/manuscripts/manuscript1/results/caviar/pc1'

source(qtl_util_file)

# Read in genoprobs, covar and pheno data.
load(data_file)

# Read in association mapping results.
assoc = readRDS(assoc_file)

pheno_sex = data.frame(sex   = rep('F', nrow(pheno)), 
                       pheno = pheno[,'lung_cxcl1'])

# PED file
ped = snps2plink(probs = probs, snpinfo = assoc$snpinfo, pheno_sex = pheno_sex)

write.table(ped, file = paste0(out_prefix, '.ped'), sep = ' ', quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# MAP file
plink_map = snps2plink(map)

write.table(plink_map, file = paste0(out_prefix, '.map'), sep = ' ', quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# COVAR file
plink_covar = data.frame(rownames(covar),    # Family ID
                         rownames(covar),    # Individual ID
                         model.matrix(~gen, data = covar)[,-1]) # Generation

write.table(plink_covar, file = paste0(out_prefix, '.covar'), sep = ' ', quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# Fake Z-score file
write.table(assoc$lod[,1], file = paste0(out_prefix, '.covar'), sep = ' ', quote = FALSE,
            row.names = FALSE, col.names = FALSE)
