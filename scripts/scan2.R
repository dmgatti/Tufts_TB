################################################################################
# Scan2 for haplotypes and SNPs.
# Daniel Gatti
# Dec 28, 2020

options(stringsAsFactors = FALSE)
library(qtl2)

# Scan2 for haplotypes.
# probs: qtl2-style genoprobs object.
# pheno: numeric matrix containing one phenotype.
# kinship: qtl2-style LOCO numeric matrices.
# addcovar: numeric matrix of covariates generated using model.matrix().
# map: qtl2-style marker map.
# pos: data.frame containing 2 rows and 3 columns: chr, start, end.
# cores: number of cores to use in scan1().
scan2 = function(probs, pheno, kinship = NULL, addcovar = NULL, map, pos, cores = 1) {
  
  # Get genoprobs on first chr.
  map_rng = map[[pos$chr[1]]][map[[pos$chr[1]]] >= pos$start[1] & map[[pos$chr[1]]] <= pos$end[1]]
  pr_chr1 = probs[,pos$chr[1]]
  pr_chr1[[1]] = pr_chr1[[1]][,,names(map_rng)]
  
  # Get genoprobs on second chr.
  map_rng = map[[pos$chr[2]]][map[[pos$chr[2]]] >= pos$start[2] & map[[pos$chr[2]]] <= pos$end[2]]
  pr_chr2 = probs[,pos$chr[2]]
  pr_chr2[[1]] = pr_chr2[[1]][,,names(map_rng)]
  
  n_markers = dim(pr_chr2[[1]])[3]
  
  # Results matrix. Put chr1 in rows and chr2 in columns. 
  result = matrix(0, n_markers, n_markers, dimnames = list(dimnames(pr_chr1)[[3]],
                                                           dimnames(pr_chr2)[[3]]))
  
  lod_chr1 = scan1(genoprobs = pr_chr1, pheno = pheno[,1,drop = FALSE], kinship = K[[pos$chr[1]]],
                   addcovar = addcovar, cores = 2)
  #lod_chr2 = scan1(genoprobs = pr_chr2, pheno = pheno[,1,drop = FALSE], kinship = K[[pos$chr[2]]],
  #                 addcovar = addcovar, cores = 2)
  
  # Subset pr_chr2 to just keep what we need: the probs without 'A'. (multicolinearity)
  pr_chr2 = pr_chr2[[1]][,LETTERS[2:8],]
  
  for(m in 1:n_markers) {
    
    if(m %% 10 == 0) print(paste(m, 'of', n_markers))
    addcovar_2 = cbind(addcovar, pr_chr2[,,m])
    intcovar  = 
    
    # TBD: Which kinship matrix do we use for kinship? Right now I use the first one.
    lod = scan1(genoprobs = pr_chr1, pheno = pheno, kinship = K[[pos$chr[1]]],
                 addcovar = addcovar_2, cores = cores)
  
    result[,m] = lod[,1]  

  } # for(m)

  # TBD: Not sure if this gives the additive LOD or not.
  result = lod_chr1[,1] - result

} # scan2()





### Test code
base_dir    = '/media/dmgatti/hdb/projects/TB'
input_file  = file.path(base_dir, 'data', 'tufts_do_tb_qtl2_input.Rdata') 
ccsnpdb = "/media/dmgatti/hda/data/MUGA/cc_variants.sqlite"
snp_func  = create_variant_query_func(dbfile = ccsnpdb)
load(input_file)
addcovar = model.matrix(~gen, data = covar)[,-1]
pheno = read.csv('/media/dmgatti/hdb/projects/TB/simulations/qtl_sims_1/simulated_phenotypes.csv')
rownames(pheno) = pheno[,1]
pheno = pheno[,-1]
pos = data.frame(chr   = c(5, 5),
                 start = c(111, 111),
                 end   = c(123, 123))

res = scan2(probs = probs, pheno = pheno[,1,drop = FALSE], kinship = K, addcovar = addcovar, 
            map = map, pos = pos, cores = 2)
