# Two-way association mapping.
# First prototpye code.
library(qtl2)

# CC SNP and gene database files (from Karl).
ccsnpdb = "/media/dmgatti/hda/data/MUGA/cc_variants.sqlite"
mgidb   = "/media/dmgatti/hda/data/MUGA/mouse_genes_mgi.sqlite"
snp_func  = create_variant_query_func(dbfile = ccsnpdb)
gene_func = create_gene_query_func(dbfile = mgidb)
ensembl = 93

scan2_snps = function(genoprobs, map, pheno, kinship = NULL, addcovar = NULL, Xcovar = NULL,
                      intcovar = NULL, weights = NULL, reml = TRUE, model = c("normal", "binary"),
                      query_func = NULL, chr = NULL, start = NULL, end = NULL, snpinfo = NULL,
                      batch_length = 20, keep_all_snps = FALSE, cores = 1, ...) {
  
  # We need to have chr, start and end.
  if(is.null(chr) | is.null(start) | is.null(end)) {
    stop("Please include chr, start and end (in Mb).")
  }

  # Get a common snpinfo for this interval. This speeds things up by skipping the disk access to get the SNPs every time.
  snpinfo = query_func(chr, start, end)

  # Get unique SNP probs.
  snpinfo = index_snps(map, snpinfo)
  snp_pr  = genoprob_to_snpprob(genoprobs, snpinfo)
  num_sdps = dim(snp_pr[[1]])[3]
  
  local_K = kinship[[chr]]
      
  # Fit the single locus model.
  lod1 = scan1(genoprobs = snp_pr, pheno = pheno, kinship = local_K, addcovar = addcovar,
               Xcovar = Xcovar, intcovar = intcovar, weights = weights, reml = reml, model = model,
               cores = cores)

  output = matrix(0, num_sdps, num_sdps)
  
  # For each SNP in snp_pr:
  #   add one SNP to the model
  #   compare it to the base model.
  #   decide what to keep.
  for(i in 1:num_sdps) {

    print(i)
    local_covar = cbind(addcovar, snp_pr[[1]][,1,i])
    lod2 = scan1(genoprobs = snp_pr, pheno = pheno, kinship = local_K, addcovar = local_covar,
                 Xcovar = Xcovar, intcovar = intcovar, weights = weights, reml = reml, model = model, cores = cores)

    output[i,] = lod2[,1]
        
  } # for(i)

  
} # scan2_snps()

# Testing
base_dir = "/media/dmgatti/hdb/projects/TB/"
load(file = paste0(base_dir, "data/phenotypes/GB_Tufts_qtl2_input.Rdata"))

# Get 1st PC of these 5 phenotypes.
library(pcaMethods)
chr15_pheno = pheno[,c("cxcl1", "cxcl2", "cxcl5", "s100a8", "tnf")]
chr15_pheno = apply(log1p(chr15_pheno), 2, scale)
chr15_pheno[chr15_pheno[,"tnf"] < -2.8,"tnf"] = NA
chr15_pca = pca(chr15_pheno, nPcs = 4, method = "bpca")
chr15_pheno = scores(chr15_pca)
rownames(chr15_pheno) = pheno$mouse

chr   = 15
start = 20
end   = 30

scan2_snps = function(genoprobs = probs, map = map, pheno = chr15_pheno[,"PC1", drop = FALSE], kinship = K, addcovar = covar, Xcovar = NULL,
                      intcovar = NULL, weights = NULL, reml = TRUE, model = c("normal", "binary"),
                      query_func = NULL, chr = NULL, start = NULL, end = NULL, snpinfo = NULL,
                      batch_length = 20, keep_all_snps = FALSE, cores = 8)

x = apply(output, 1, max)
y = apply(output, 2, max)
wh = which(x > 6)






  