################################################################################
# Utility scripts to help with QTL analysis.
# Dec. 20, 2020
# Daniel Gatti
# dmgatti@coa.edu
################################################################################
# snps2plink: convert DO SNPs to PLINK ped format.
# map2plink: convert qtl2-style map to PLINK map format.
# snp_sdp_plot: Given a set of SNPs with LOD scores, get the founder SDP and
#               plot it as a matrix.

options(stringsAsFactors = FALSE)

library(qtl2)

################################################################################
# snps2plink: convert DO SNPs to PLINK ped format. This does not get the true
#             alleles from the Sanger sequence file. It writes out A & C for 
#             all alleles. This is to use for LD and mapping calculations.
# Arguments: probs: qtl2-style genoprobs object.
#            snpinfo: qtl2-style snpinfo object for the interval of interest
#            from scan1snps().
#            pheno_sex: data.frame containing two columns names 'sex' and 
#                       'pheno'. Sample IDs must be in rownames. Sex is a 
#                       character vector containing 'M', 'F' or NA. Pheno is
#                       numeric vector.
# Returns: data.frame in PLINK ped format, ready to be written out to a text
#          file.
snps2plink = function(probs, snpinfo, pheno_sex) {
  
  # Get SNP probs. We assume one chromosome here.
  pr = genoprob_to_snpprob(genoprobs = probs, snpinfo = snpinfo)[[1]][,1,]
  pr = round(2 * pr)
  
  # Convert to letters.
  tmp = matrix('AC', nrow(pr), ncol(pr), dimnames = dimnames(pr))
  tmp[pr == 0] = 'AA'
  tmp[pr == 2] = 'CC'
  rm(pr)
  
  # Create the 6 header columns: family, individual, fater, mother, sex, phenotype.
  sex_conv = ifelse(pheno_sex$sex == 'F', 2, 1)
  retval = data.frame(rownames(pheno_sex),     # family
                      rownames(pheno_sex),     # individual
                      rep(0, nrow(pheno_sex)), # father
                      rep(0, nrow(pheno_sex)), # mother)
                      sex_conv,                # sex
                      pheno_sex$pheno,         # phenotype
                      tmp)                     # genotypes
  
  return(retval)
  
} # snps2plink()


################################################################################
# map2plink: convert qtl2-style map to PLINK map format.
# Arguments: snpinfo: qtl2-style snpinfo object form scan1snps().
# Returns: data.frame containing four columns: chr, variant ID, posiiton in 
#          cM, position in base pairs.
map2plink = function(snpinfo) {
  
  retval = data.frame(snpinfo$chr,            # chr
                      snpinfo$snp_id,         # SNP ID
                      rep(0, nrow(snpinfo)),  # cM
                      snpinfo$pos * 1e6)      # pos
    
  return(retval)

} # map2plink()


################################################################################
# Marker correlation.
# From: https://stackoverflow.com/questions/1897704/angle-between-two-vectors-in-r
cos_vec = function(x, y) {
  dot.prod = x %*% y 
  norm.x = norm(x, type = "2")
  norm.y = norm(y, type = "2")
  dot.prod / (norm.x * norm.y)
} # cos_vec()

marker_corr = function(probs, map, chr) {
  
  pr = probs[[chr]]
  mp = map[[chr]]
  
  n_samples  = dim(pr)[1]
  n_founders = dim(pr)[2]
  n_markers  = dim(pr)[3]
  pr = matrix(pr, nrow = n_samples * n_founders, ncol = n_markers)
  
  # Remove markers at the same location.
  keep = which(!duplicated(mp))
  mp = mp[keep]
  pr = pr[,keep]
  
  colnames(pr) = mp
  
  return(cor(pr))

} # marker_corr()


################################################################################
# Founder SNP SDP plot.
# topsnps: data.frame from qtl2::top_snps() to be plotted.
# rng: numerice vector containing start and end in Mb for plot.
snp_sdp_plot = function(topsnps, rng) {
  
  # Get SDPs.
  founders = c("A_J", "C57BL_6J", "129S1_SvImJ", "NOD_ShiLtJ", "NZO_HlLtJ", "CAST_EiJ", "PWK_PhJ", "WSB_EiJ")
  # Using image() o na matrix generated a plot that was too dense.
  # I'm plotting a fake matrix of the same size as spds and then drawing the SNPs.
  fakemat = matrix(0, nrow = nrow(topsnps), ncol = 8)

  par(plt = c(0.15, 0.99, 0.15, 0.99))
  image(x = topsnps$pos, y = 1:8, fakemat, breaks = c(-0.5, 0.5), col = 'grey90', 
        xlim = rng, yaxt = 'n', ylab = '', bg = 'grey90')
  usr = par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90')
  mtext(founders, side = 2, line = 0.1, at = 1:8, las = 1)
  abline(h = 1:7 + 0.5, col = 'white')
  for(i in 1:length(founders)) {
    
    wh = which(topsnps[,founders[i]] > 1)
    if(length(wh) > 0) {
      pos = rep(topsnps$pos[wh], each = 3)
      y   = rep(c(i, i + 1, i + 2) - 0.5, length(pos))
      df = data.frame(x = pos, y = y)
      df[seq(3, nrow(df), 3),] = c(NA, NA)
      lines(df$x, df$y, col = qtl2::CCcolors[i])
    } # if(length(wh) > 0)
    
  } # for(i)
  
} # snp_sdp_plot()

