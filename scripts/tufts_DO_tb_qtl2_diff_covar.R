##############################
# Tufts TB QTL Mapping.
# DMG
# Feb. 23, 2020
##############################
library(qtl2)

base_dir    = '/media/dmgatti/hdb/projects/TB'
input_file  = file.path(base_dir, 'data', 'tufts_do_tb_qtl2_input.Rdata') 
results_dir = file.path(base_dir, 'results', 'qtl2', 'try_different_covars')

# Read in the DOQTL formatted data.
load(input_file)

# CC SNP and gene database files (from Karl).
ccsnpdb = "/media/dmgatti/hda/data/MUGA/cc_variants.sqlite"
mgidb   = "/media/dmgatti/hda/data/MUGA/mouse_genes_mgi.sqlite"
snp_func  = create_variant_query_func(dbfile = ccsnpdb)
gene_func = create_gene_query_func(dbfile = mgidb)
ensembl = 93

peaks = NULL

########################################
# Create overall kinship matrix.
allK = calc_kinship(probs = probs, type = 'overall', cores = 4)

########################################
# QTL mapping with generation as factor.

# Crate additive covariates.
addcovar = model.matrix(~gen, data = covar)[,-1, drop = FALSE]

# Estimate heritability.
hsq = est_herit(pheno = pheno, kinship = allK, addcovar = addcovar, cores = 8)

# Perform genome scans.
qtl = scan1(genoprobs = probs, pheno = pheno, kinship = K, addcovar = addcovar, cores = 8)

# Save data.
save(qtl, hsq, file = file.path(results_dir, 'lod_gen_factor.Rdata'))

########################################
# QTL mapping with generation and run as factor.

# Crate additive covariates.
addcovar = model.matrix(~gen:aero_run, data = covar)[,-1, drop = FALSE]

# Estimate heritability.
hsq = est_herit(pheno = pheno, kinship = allK, addcovar = addcovar, cores = 8)

# Perform genome scans.
qtl = scan1(genoprobs = probs, pheno = pheno, kinship = K, addcovar = addcovar, cores = 8)

# Save data.
save(qtl, hsq, file = file.path(results_dir, 'lod_gen_run_factor.Rdata'))

########################################
# QTL mapping with mtb_dose as numeric.

# Create additive covariates.
addcovar = model.matrix(~mtb_dose, data = covar)[,-1, drop = FALSE]

# Estimate heritability.
hsq = est_herit(pheno = pheno, kinship = allK, addcovar = addcovar, cores = 8)

# Perform genome scans.
qtl = scan1(genoprobs = probs, pheno = pheno, kinship = K, addcovar = addcovar, cores = 8)

# Save data.
save(qtl, hsq, file = file.path(results_dir, 'lod_dose_numeric.Rdata'))

########################################
# QTL mapping with mtb_dose as factor.

# Create additive covariates.
covar$mtb_dose = factor(covar$mtb_dose)
addcovar = model.matrix(~mtb_dose, data = covar)[,-1, drop = FALSE]

# Estimate heritability.
hsq = est_herit(pheno = pheno, kinship = allK, addcovar = addcovar, cores = 8)

# Perform genome scans.
qtl = scan1(genoprobs = probs, pheno = pheno, kinship = K, addcovar = addcovar, cores = 8)

# Save data.
save(qtl, hsq, file = file.path(results_dir, 'lod_dose_factor.Rdata'))

########################################
# QTL mapping with euth as factor.

# Create additive covariates.
addcovar = model.matrix(~euth, data = covar)[,-1, drop = FALSE]

# Estimate heritability.
hsq = est_herit(pheno = pheno, kinship = allK, addcovar = addcovar, cores = 8)

# Perform genome scans.
qtl = scan1(genoprobs = probs, pheno = pheno, kinship = K, addcovar = addcovar, cores = 8)

# Save data.
save(qtl, hsq, file = file.path(results_dir, 'lod_euth_factor.Rdata'))

########################################
# QTL mapping with gen/cage as factor.

# Create additive covariates.
addcovar = model.matrix(~gen:cage, data = covar)[,-1, drop = FALSE]

# Estimate heritability.
hsq = est_herit(pheno = pheno, kinship = allK, addcovar = addcovar, cores = 8)

# Perform genome scans.
qtl = scan1(genoprobs = probs, pheno = pheno, kinship = K, addcovar = addcovar, cores = 8)

# Save data.
save(qtl, hsq, file = file.path(results_dir, 'lod_gen_cage_factor.Rdata'))


########################################
# Load in the data and evaluate the results.

load(file = file.path(results_dir, 'lod_gen_factor.Rdata'))
qtl_gen_f = qtl
hsq_gen_f = hsq

load(file.path(results_dir, 'lod_gen_run_factor.Rdata'))
qtl_gen_run_f = qtl
hsq_gen_run_f = hsq

load(file = file.path(results_dir, 'lod_dose_numeric.Rdata'))
qtl_dose_n = qtl
hsq_dose_n = hsq

load(file = file.path(results_dir, 'lod_dose_factor.Rdata'))
qtl_dose_f = qtl
hsq_dose_f = hsq

load(file = file.path(results_dir, 'lod_euth_factor.Rdata'))
qtl_euth_f = qtl
hsq_euth_f = hsq

load(file = file.path(results_dir, 'lod_gen_cage_factor.Rdata'))
qtl_gen_cage_f = qtl
hsq_gen_cage_f = hsq

hsq = data.frame(phenotype    = names(hsq_gen_f),
                 gen_factor   = hsq_gen_f,
                 gen_run_factor = hsq_gen_run_f,
                 dose_numeric = hsq_dose_n,
                 dose_factor  = hsq_dose_f,
                 euth_factor  = hsq_euth_f,
                 gen_cage_factor = hsq_gen_cage_f)

peaks = vector('list', 6)
names(peaks) = c('gen_f', 'gen_run', 'dose_n', 'dose_f', 'euth_f', 'gen_cage_f')
peaks[[1]] = find_peaks(qtl_gen_f, map, threshold = 7)
peaks[[2]] = find_peaks(qtl_gen_run_f, map, threshold = 7)
peaks[[3]] = find_peaks(qtl_dose_n, map, threshold = 7)
peaks[[4]] = find_peaks(qtl_dose_f, map, threshold = 7)
peaks[[5]] = find_peaks(qtl_euth_f, map, threshold = 7)
peaks[[6]] = find_peaks(qtl_gen_cage_f, map, threshold = 7)

peaks[[1]]$lodindex = 'gen_f'
peaks[[2]]$lodindex = 'gen_run_f'
peaks[[3]]$lodindex = 'dose_n'
peaks[[4]]$lodindex = 'dose_f'
peaks[[5]]$lodindex = 'euth_f'
peaks[[6]]$lodindex = 'gen_cage_f'

for(i in 2:length(peaks)) {
  peaks[[1]] = rbind(peaks[[1]], peaks[[i]])
}
peaks = peaks[[1]]

table(peaks$lodcolumn, peaks$chr)


