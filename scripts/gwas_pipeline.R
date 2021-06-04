################################################################################
# GWAS pipeline for qtl2
# Daniel Gatti
# May 20, 2020
# dmgatti@coa.edu
################################################################################
options(stringsAsFactors = FALSE)

library(tidyverse)
library(AnnotationHub)
require(BiocParallel)
library(rhdf5)
require(survival)
library(qtl2)

################################################################################
# genoprobs: qtl2 style genoprobs.
# map: qtl2 style physical marker map.
# pheno: phenotype matrix.
# kinship: qtl2 style kinship matrices.
# snp_file: full path to the CC variant file.
# thr: LOD threshold above which to write out top SNPs.
# cores: number of cores to use.
# out_dir: full path to output directory.
# verbose: boolean that is TRUE if progress should be shown.
run_gwas = function(genoprobs, map, pheno, kinship = NULL, addcovar = NULL, 
                    intcovar = NULL, snp_file, thr = 4, cores = 1, out_dir,
                    verbose = FALSE) {

  # Create SNP query function.
  snp_func  = create_variant_query_func(dbfile = snp_file)
  
  for(i in 1:ncol(pheno)) {
    
    pheno_name = colnames(pheno)[i]
    
    # Subset samples.
    samples2use = which(!is.na(pheno[,pheno_name]))
    samples2use = intersect(samples2use, rownames(addcovar))
    samples2use = intersect(samples2use, rownames(genoprobs[[1]]))
    
    if(verbose) {
      print(paste0(pheno_name, " : ", length(samples2use)))
    } # if(verbose)
    
    tmpK = K
    for(j in 1:length(K)) {
      tmpK[[j]] = K[[j]][samples2use, samples2use]
    }
    
    tmpcovar = NULL
    if(!is.null(addcovar)) {
      tmpcovar = addcovar[samples2use,,drop = FALSE]
      tmpcovar = tmpcovar[,colSums(tmpcovar) > 0, drop = FALSE]
    }
    
    tmpintcovar = NULL
    if(!is.null(intcovar)) {
      tmpintcovar = intcovar[samples2use,,drop = FALSE]
    }
    
    tmpprobs = genoprobs[samples2use,]
    tmppheno = pheno[samples2use, pheno_name, drop = FALSE]
    
    if(verbose) {  
      print('   mapping')
    }
    assoc = scan1snps(genoprobs = tmpprobs, 
                        map       = map,
                        pheno     = tmppheno[,pheno_name, drop = FALSE],
                        kinship   = tmpK,
                        addcovar  = tmpcovar,
                        intcovar  = tmpintcovar,
                        query_func = snp_func,
                        keep_all_snps = FALSE,
                        cores      = cores)
      
    if(verbose) {  
      print('   ploting')
    }
    png(file.path(out_dir, paste0(pheno_name, '_gwas.png')), width = 2000, height = 1000, res = 128)
    plot_snpasso(assoc$lod, assoc$snpinfo, main = pheno_name)
    abline(h = thr, col = 'red', lwd = 2)
    dev.off()
      
    if(verbose) {  
      print('   writing')
    }
    
    saveRDS(assoc, file = file.path(out_dir, paste0(pheno_name, '_gwas.rds')))
    
    topsnps = top_snps(assoc$lod, assoc$snpinfo, drop = 1.5)
    topsnps = subset(topsnps, lod >= thr)
    
    if(nrow(topsnps) > 0) {
       readr::write_csv(topsnps, file.path(out_dir, paste0(pheno_name, '_gwas_chr', chr,'.csv')))
    } # if(nrow(topsnps) > 0)
    
    rm(assoc, topsnps, tmpprobs, tmpK)  
    gc()
    
  } # for(i)
  
} # run_gwas()


################################################################################
# GWAS permutations.
# Pass in a phenotype with a single column.
run_gwas_perms = function(genoprobs, map, pheno, kinship = NULL, addcovar = NULL, 
                          intcovar = NULL, snp_file, 
                          cores = 1, out_dir, nperm = 1) {

  pheno_name = colnames(pheno)[1]
  
  samples2use = which(!is.na(pheno[,pheno_name]))
  
  print(paste0(pheno_name, " : ", length(samples2use)))
  
  pheno = pheno[samples2use,,drop = FALSE]
  for(j in 1:length(kinship)) {
    kinship[[j]] = kinship[[j]][samples2use, samples2use]
  }
  genoprobs = genoprobs[samples2use,]
  if(!is.null(addcovar)) {
    addcovar = addcovar[samples2use,,drop = FALSE]
    addcovar = addcovar[,colSums(addcovar) > 0, drop = FALSE]
  }
  if(!is.null(intcovar)) {
    intcovar = intcovar[samples2use,,drop = FALSE]
  }
  
  # Create nperm phenotypes.
  pheno_perm = matrix(sapply(1:nperm, function(z) { sample(pheno[,1]) }), nrow = nrow(pheno), ncol = nperm, 
                      dimnames = list(rownames(pheno), paste0('p', 1:nperm)))
  
  perms = matrix(0, nrow = length(genoprobs), ncol = nperm,
                 dimnames = list(names(genoprobs), paste0('p', 1:nperm)))
  
  for(chr in names(genoprobs)) {

      print(paste('CHR', chr))
      assoc = scan1snps(genoprobs = genoprobs[,chr], 
                        map       = map,
                        pheno     = pheno_perm,
                        kinship   = kinship[[chr]],
                        addcovar  = addcovar,
                        intcovar  = intcovar,
                        chr       = chr,
                        start     = 0,
                        end       = 200,
                        query_func = snp_func,
                        keep_all_snps = FALSE,
                        cores      = cores)
      perms[chr,] = apply(assoc$lod, 2, max)
      gc()

  } # for(chr)
    
  perms = apply(perms, 2, max)
  saveRDS(perms, file = file.path(out_dir, 'gwas_perms.rds'))

  return(perms)
  
} # run_gwas_perms()



################################################################################
# Plot gwas
# h5filename: full path to HDF5 file containing GWAS results.
# thr: numeric value of significance threshold.
# title: character with title.
# ensembl_ver: ensembl version.
# min_lod: floating point number below which points will not be plotted.
plot_gwas = function(h5filename, thr = NULL, title = '', ensembl_ver = 93, min_lod = 1) {
  
  #genes = get_ensembl(ensembl_ver)
  
  # Get the chromosomes from the file. (Called 'groups' in hdf5 lingo)
  grps    = h5ls(file = h5filename, recursive = FALSE)
  chr_num = as.numeric(sub('^chr', '', grps$name))
  grps    = grps[order(chr_num),]
  
  chr = NULL
  pos = NULL
  lod = NULL
  chrlen = rep(0, nrow(grps))
  names(chrlen) = grps$name
  
  for(i in 1:nrow(grps)) {
    
    print(i)
    
    lod_all  = h5read(h5filename, paste(grps$name[i], 'lod',     sep = '/'))
    info_all = h5read(h5filename, paste(grps$name[i], 'snpinfo', sep = '/'))
    snp_map  = qtl2:::snpinfo_to_map(info_all)
    snp_exp  = qtl2:::expand_snp_results(lod_all, snp_map, info_all)
    
    pos_lod = matrix(c(snp_exp$map[[1]], snp_exp$lod[,1]), ncol = 2, dimnames = list(NULL, c('pos', 'lod')))
    
    rm(lod_all, info_all, snp_map, snp_exp)
    
    chrlen[i] = max(pos_lod[,'pos'])
    # Create 100 breaks per Mb.
    brk   = seq(0, chrlen[i] + 1, 0.01)
    brk   = cut(pos_lod[,'pos'], breaks = brk)
    brk   = factor(as.numeric(brk))
    
    pos_i = sapply(split(pos_lod[,'pos'], brk), mean)
    lod_i = sapply(split(pos_lod[,'lod'], brk), max)
    
    chr = c(chr, rep(grps$name[i], length(pos_i)))
    pos = c(pos, pos_i)
    lod = c(lod, lod_i)
    
  } # for(i)
  
  data = tibble(chr, pos, lod) %>%
           mutate(chr = str_replace(chr, '^chr', ''),
                  chr = factor(chr, levels = c(1:19, 'X')))
  ggplot(data) +
    geom_point(aes(pos, lod)) +
    {
      if(!is.null(thr)) geom_hline(aes(yintercept = thr), color = 'red') 
    } +
    scale_x_continuous(breaks = c(0, 50, 100, 150)) +
    facet_wrap(~chr, nrow = 1, scales = 'free_x') +
    labs(title = title) +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
          panel.spacing.x = unit(0, 'lines'))
  
} # plot_gwas()


################################################################################
get_ensembl = function(ensembl_ver) {
  
  # Get the Ensembl GTF and just keep the gene locations.
  hub = AnnotationHub()
  ensembl = hub[[names(hub)[hub$title == paste0('Mus_musculus.GRCm38.', ensembl_ver, '.chr.gtf')]]]
  genes   = subset(ensembl, type == 'gene')
  
  # Reformat genes to satisfy qtl2.
  genes = data.frame(chr = seqnames(genes), start = start(genes) * 1e-6, stop = end(genes) * 1e-6, 
                     strand = strand(genes), Name = genes$gene_name, ensembl = genes$gene_id)
  
  return(genes)
  
} # get_ensembl()


################################################################################
# Harvest gwas peaks.
# h5filename: full path to HDF5 file containing GWAS results.
# thr: numeric value of significance threshold.
# ensembl_ver: ensembl version.
harvest_gwas = function(h5filename, thr = NULL, ensembl_ver = 93) {

  genes = get_ensembl(ensembl_ver)
  
  # Get the chromosomes from the file. (Called 'groups' in hdf5 lingo)
  grps    = h5ls(file = h5filename, recursive = FALSE)
  chr_num = as.numeric(sub('^chr', '', grps$name))
  grps    = grps[order(chr_num),]
  
  result = NULL
  chrlen = rep(0, nrow(grps))
  names(chrlen) = grps$name
  
  for(i in 1:nrow(grps)) {

    print(i)
    
    lod_all  = h5read(h5filename, paste(grps$name[i], 'lod',     sep = '/'))
    info_all = h5read(h5filename, paste(grps$name[i], 'snpinfo', sep = '/'))
    snp_map  = qtl2:::snpinfo_to_map(info_all)
    snp_exp  = qtl2:::expand_snp_results(lod_all, snp_map, info_all)
    snp_exp  = cbind(info_all, lod = snp_exp$lod)
    snp_exp  = subset(snp_exp, lod >= thr)

    snp_exp = left_join(snp_exp, select(genes, ensembl, Name), by = c('ensembl_gene' = 'ensembl'))
    result = bind_rows(result, snp_exp)
    
  } # for(i)
  
  result = select(result, snp_id:ensembl_gene, symbol = Name, consequence:lod)
  
  return(result)

} # harvest_gwas()


################################################################################
# Given a set of GWAS results, keep the results above a given threshold and
# summarize the width on each chromosome and the genes that intersect with
# SNPs.
# results: data.frame with phenotype in the first column, as proficed by harvest_gwas().
# thr: float that is the LOD threshold above which values will be kept.
summarize_gwas = function(results, thr = NULL) {

  if(!is.null(thr)) {
    results = filter(results, lod >= thr)
  } # if(!is.null(thr))
  
  res = results %>% 
          group_by(phenotype, chr) %>%
          summarize(proximal = min(pos, na.rm = TRUE),
                    distal   = max(pos, na.rm = TRUE),
                    ensembl  = str_c(base::unique(ensembl_gene), sep = ';', collapse = ';'),
                    symbol   = paste(base::unique(symbol),       sep = ';', collapse = ';')) %>%
          mutate(ensembl = str_replace(ensembl, 'NA', ''),
                 symbol  = str_replace(symbol,  'NA', ''),
                 ensembl = str_replace_all(ensembl, regex('^;'), ''),
                 symbol  = str_replace_all(symbol,  regex('^;'), ''),
                 ensembl = str_replace_all(ensembl, regex(';;'), ';'),
                 symbol  = str_replace_all(symbol,  regex(';;'), ';'))

   return(res)
  
} # summarize_gwas()


################################################################################
# Arguments:

# cores:     Number of cores to use in attempting parallel computation. 
#            Not implemented yet.


# genoprobs: qtl2-style genoprobs object with 8 founder allele probabilities.
#            List containing 20 elements, each of which is a 3 dimensional
#            numeric array with samples in rows, founders in columns and
#            markers in slices.
# map: qtl2 style physical marker map.
# pheno:     data.frame containing at least two columns:
#            'survival': Number of days that mouse lived.
#            'event':    Whether the death was observed or the mouse was euthanized.
#                        0 means that the mouse was euthanized at a time point.
#                        1 means that their death was observed.
# kinship: list of qtl2 style kinship matrices.
# addcovar:  Matrix containing additive covariates for the mapping model.
# snp_file: full path to the CC variant file.
# cores: number of cores to use.
# out_dir: full path to output directory.
# verbose: boolean that is TRUE if progress should be shown.
# Returns: Nothing. Writes results to HDF5 file in out_dir.
gwas_surv = function(genoprobs, map, pheno, kinship = NULL, addcovar = NULL, 
                     snp_file, cores = 1, out_dir, verbose = FALSE) {
  
  pheno_surv = Surv(time = pheno$survival, event = pheno$event)

  # Create SNP query function.
  snp_func  = create_variant_query_func(dbfile = snp_file)
  
  # Create HDF5 file.
  h5filename = file.path(out_dir, paste0('survival_gwas.h5'))
  
  if(file.exists(h5filename)) {
    file.remove(h5filename)
  }
  h5createFile(h5filename)
  
  # Null model log_likelihood.
  null_ll = coxph(pheno_surv ~ addcovar)$loglik[2]
  
  # For each chromosome...
  for(chr in names(probs)) {
    
    print(paste('chr', chr))
    
    # Get genoprobs for current samples.
    pr = probs[[chr]][rownames(pheno),-1,]
    
    snps   = snp_func(chr = chr, start = 0, end = 200)
    snps   = index_snps(map, snps)
    snp_pr = genoprob_to_snpprob(genoprobs, snps)[[1]]
    
    print(paste('   ', dim(snp_pr)[3], 'SDPs   ', nrow(snps), 'SNPs'))
    
    lod = matrix(0, nrow = dim(snp_pr)[3], ncol = 1, dimnames = list(dimnames(snp_pr)[[3]], "survival"))

    for(j in 1:dim(snp_pr)[3]) {

      if(j %% 100 == 0) print(j)

      lod[j,1] = coxph(pheno_surv ~ addcovar + snp_pr[,1,j])$loglik[2]

    } # for(j)
    
    lod = (lod - null_ll) / log(10)
    
    h5createGroup(h5filename, paste0('chr', chr))

    h5write(lod,  h5filename, paste0('chr', chr, '/lod'))
    h5write(snps, h5filename, paste0('chr', chr, '/snpinfo'))
    
    rm(pr, snps, snp_pr, lod)
    gc()

  } # for(i)
  
  H5close()
  
} # survscan()

### Test code.
#covar = covar[rownames(pheno),,drop = FALSE]

#qtl_surv = suppressWarnings(survscan(pheno_surv, probs, covar))

#qtl_surv = qtl_surv[rownames(qtl),,drop=FALSE]

#qtl = cbind(qtl, qtl_surv)

