################################################################################
# Survival QTL maping in the DO using qtl2.
# Daniel Gatti
# dmgatti@coa.edu
# Mar. 6, 2020
################################################################################

# Required libraries.
require(survival)
require(BiocParallel)

# Arguments:
# genoprobs: qtl2-style genoprobs object with 8 founder allele probabilities.
#            List containing 20 elements, each of which is a 3 dimensional
#            numeric array with samples in rows, founders in columns and
#            markers in slices.
# pheno:     data.frame containing at least two columns:
#            'survival': Number of days that mouse lived.
#            'event':    Whether the death was observed or the mouse was euthanized.
#                        0 means that the mouse was euthanized at a time point.
#                        1 means that their death was observed.
# addcovar:  Matrix containing additive covariates for the mapping model.
# cores:     Number of cores to use in attempting parallel computation. 
#            Not implemented yet.

scan1_surv = function(genoprobs, pheno, addcovar = NULL, cores = 1) {

  pheno_surv = Surv(time = pheno$survival, event = pheno$event)

  lod = vector("list", length(probs))
  
  # Null model log_likelihood.
  null_ll = coxph(pheno_surv ~ covar)$loglik[2]

  # For each chromosome...
  for(i in seq_along(probs)) {
    
    print(i)
    pr = probs[[i]][rownames(pheno),-1,]
    
    lod[[i]] = matrix(0, nrow = dim(pr)[3], ncol = 1, dimnames = list(dimnames(pr)[[3]], "survival"))
    
    for(j in 1:dim(pr)[[3]]) {
      
      lod[[i]][j,1] = coxph(pheno_surv ~ covar + pr[,,j])$loglik[2]
      
    } # for(j)
    
    lod[[i]] = (lod[[i]] - null_ll) / log(10)
    
  } # for(i)
  
  markers = unlist(sapply(probs, function(z) { dimnames(z)[[3]] }))
  
  retval = matrix(unlist(lod), ncol = 1, dimnames = list(markers, "survival"))
  class(retval) = c("scan1", "matrix")
  
  return(retval)
  
} # survscan()

### Test code.
#covar = covar[rownames(pheno),,drop = FALSE]

#qtl_surv = suppressWarnings(survscan(pheno_surv, probs, covar))

#qtl_surv = qtl_surv[rownames(qtl),,drop=FALSE]

#qtl = cbind(qtl, qtl_surv)

