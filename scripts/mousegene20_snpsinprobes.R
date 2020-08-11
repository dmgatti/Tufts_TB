# Handle SNPs in probes on Mouse Gene 2.0 ST array.
options(stringsAsFactors = FALSE)
library(mogene20stprobeset.db)
library(GenomicRanges)
library(qtl2)
probe_file  = '/home/dmgatti/Documents/data/microarray/affymetrix/mogene20st/MoGene-2_0-st-v1-mm10-probe-tab/MoGene-2_0-st-v1.mm10.probe.tab'
probeset_file = '/home/dmgatti/Documents/data/microarray/affymetrix/mogene20st/MoGene-2_0-st-v1.na36.mm10.probeset.csv'
sanger_file = '/media/dmgatti/hda/data/MUGA/cc_variants.sqlite'
base_dir = '/media/dmgatti/hdb/projects/TB/'

snp_func = qtl2::create_variant_query_func(dbfile = sanger_file)

probes = read.delim(probe_file, comment.char = '#')
probes = subset(probes, startsWith(probes$seqname, 'chr'))
probes$seqname = sub('^chr._', '', probes$seqname)
#probes$seqname = sub('_random$', '.1', probes$seqname)
#probes$seqname = sub('^chrUn_', '', probes$seqname)
probes$seqname = sub('^chr', '', probes$seqname)
probes$seqname = sub('M', 'MT', probes$seqname)
probes = GRanges(seqnames = probes$seqname,
                 ranges = IRanges(start = probes$start,
                                  end   = probes$stop),
                 strand = probes$strand,
                 mcols = probes[,c(1:2,10)])

probeset = read.csv(probeset_file, comment.char = '#')

result = NULL

for(chr in seqlevels(probes)) {

  print(chr)  
  snps = snp_func(chr, 0, 200e6)
  snps = GRanges(seqnames = snps$chr,
                 ranges = IRanges(start = snps$pos * 1e6, width = 1),
                 mcols = snps[,c(1, 4, 6:ncol(snps))])
  
  ol = findOverlaps(probes, snps)
  result = rbind(result, data.frame(probes[queryHits(ol),], snps[subjectHits(ol),]))
  
} # for(chr)

saveRDS(result, file = paste0(base_dir, 'results/qtl2/eqtl/snps_in_probes.rds'))

# Connect the probes with SNPs to the probesets and Entrez IDs.
genes = select(mogene20stprobeset.db, keys = as.character(result$mcols.Transcript.Cluster.ID), keytype = 'PROBEID', columns = c('PROBEID', 'ENSEMBL', 'ENTREZID', 'GENENAME'))
genes$PROBEID = as.character(genes$PROBEID)

result = merge(result, genes, by.x = 'mcols.Transcript.Cluster.ID', by.y = 'PROBEID', all.x = TRUE)




