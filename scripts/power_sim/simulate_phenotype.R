################################################################################
# Simulate a phenotype that has 2 QTL.  Allow the user to specify the founder strain 
# split, the location of the QTL and the effect size as a proportion of the total phenotype
# variance.
# Daniel Gatti
# dan.gatti@jax.org
# Oct. 19, 2016
################################################################################
# Arguments:
# probs: 3 dimensional numeric arrary of haplotype probabilities; samples in rows,
#        8 DO founders in columns and markers in dim 3.
# markers: data.frame containing the marker ID, Chr and Mb position in columns 1:3.
# chr: numeric vector of chromosome locations for QTL.
# pos: numeric vector of Mb locations for QTL.
# eff: numeric vector of effect sizes for each QTL; must be between 0 and 1 and sum of
#        of all effects must be < 1.
# spl: numeric vector of minor allele frequencies. We will split the QTL effects by choosing
#           founders at random. Each value must be between 1 and 4. Optionally, a character
#           vector containing the exact founders to split by. Ex: c("A", "ACDE")
simulate.phenotype = function(probs, markers, chr = c(1,1), pos = c(50, 55), eff = c(0.1, 0.1), 
                              spl = c(1, 4)) {

  # Get the probs at the requested positions.
  markers = markers[markers[,2] %in% chr,]
  rownames(markers) = markers[,1]
  keep = rep(0, length(chr))
  for(i in 1:length(chr)) {
    ss = markers[markers[,2] == chr[i],]
    wh = which.min(abs(pos[i] - ss[,3]))
    keep[i] = ss[wh,1]
  } # for(i)
  markers = markers[keep,]
  probs = probs[,,markers[,1]]

  # If the user hasn't given us founder splits, create them.
  splits = as.character(spl)
  if(is.numeric(spl)) {
    for(i in 1:length(chr)) {
      splits[i] = sample(LETTERS[1:8], spl[i])
    } # for(i)
  } # if(is.numeric(spl))
  split = strsplit(splits, split = "")
  
  # Simulate the phenotype.
  eff = matrix(0, nrow(probs), ncol(probs))
  for(i in 1:length(chr)) {

    samples = rowSums(probs[,splits[i],i] > 0)
    eff[samples,] = probs[samples,,i] * 2 * eff[i]

  } # for(i)
  
# TBD: Need to work on getting the effects scaled to be a proportion of variance.

  # Simulate the phenotype.
  pheno = eff + matrix(rnorm(nrow(eff) * ncol(eff)), nrow(eff), ncol(eff))

} # simulate.phenotype()