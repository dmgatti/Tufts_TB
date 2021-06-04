################################################################################
# Read in TB gene sets from GSEA and get their mouse genome locations and
# the LOD score for the M Tb burden phenotype.
# Will pass these to Anna and Matt at JAX.
#
# Daniel Gatti
# dmgtti@coa.edu
# Mar. 15, 2011
################################################################################

# Libraries.
library(biomaRt)
library(GenomicRanges)
library(tidyverse)

# Setup
base_dir = '/media/dmgatti/hdb/projects/TB'
working_dir = file.path(base_dir, 'results', 'tb_gene_list')
gwas_dir = file.path(base_dir, 'results', 'qtl2', 'gwas')
gs_file  = file.path(working_dir, 'gsea_tb_gene_lists.csv')

# Setup Biomart.
# Using GRCm38 and Ensembl 102
# We have human genes and need mouse orthologs and then mouse locations.
# Human mart
hmart = biomaRt::useEnsembl('genes', 'hsapiens_gene_ensembl', version = 102)

# Mouse mart
mmart = biomaRt::useEnsembl('genes', 'mmusculus_gene_ensembl', version = 102)

attribs = listAttributes(mmart)
filters = listFilters(mmart)

# Read in gene sets.
genes = read_csv(gs_file) %>% 
         distinct()
colnames(genes) = c('name', 'entrezid', 'symbol', 'desc')

# Get the mouse homologs. We have EntrezIDs, but biomart returns Ensembl IDs.
mgenes = biomaRt::getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'mmusculus_homolog_ensembl_gene'), 
                        filters = 'entrezgene_id', 
                        values = genes$entrezid, 
                        mart = hmart, verbose = TRUE)

# Look up location of genes by Entrez Gene ID.
res = biomaRt::getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'start_position', 'end_position'), 
                     filters = 'ensembl_gene_id', 
                     values = mgenes$mmusculus_homolog_ensembl_gene, 
                     mart = mmart, verbose = TRUE)

write_csv(res, file = file.path(working_dir, 'mouse_tb_genes.csv'))

# Read in the M Tb. burden GWAS and get the LOD at the SNP nearest the gene start.
genes = GRanges(seqnames = res$chromosome_name,
                ranges   = IRanges(start = res$start_position, end = res$end_position),
                ensembl  = res$ensembl_gene_id,
                symbol   = res$external_gene_name)

write.csv(as.data.frame(genes), file = file.path(working_dir, 'tb_genes.csv'), row.names = FALSE,
          quote = FALSE)

# Mtb burden
gwas = readRDS(file.path(gwas_dir, 'mtb_burden_gwas.rds'))

# Expand the SNPs to cover the whole genome.
snp_map = qtl2:::snpinfo_to_map(gwas$snpinfo)
tmp = qtl2:::expand_snp_results(gwas$lod, snp_map, gwas$snpinfo)
lod = GRanges(seqnames = rep(names(tmp$map), sapply(tmp$map, length)),
              ranges = IRanges(start = unlist(tmp$map) * 1e6, width = 1),
              snp = unlist(sapply(tmp$map, names)),
              lod = tmp$lod)

all_genes = biomaRt::getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'start_position', 'end_position'),
                     mart = mmart, verbose = TRUE)
all_genes = GRanges(seqnames = all_genes$chromosome_name,
                    ranges   = IRanges(start = all_genes$start_position, end = all_genes$end_position),
                    ensembl  = all_genes$ensembl_gene_id,
                    symbol   = all_genes$external_gene_name)

gene_snp_ol = findOverlaps(query = all_genes, subject = lod)
gene_snp_ol = split(gene_snp_ol, queryHits(gene_snp_ol))
names(gene_snp_ol) = all_genes$ensembl[as.numeric(names(gene_snp_ol))]

gene_ol = lapply(gene_snp_ol, subjectHits)
gene_ol = sapply(gene_ol, function(z) { max(lod$lod.mtb_burden[z], na.rm = TRUE) })
gene_ol = data.frame(ensembl = names(gene_ol), lod = gene_ol)

all_genes = subset(all_genes, ensembl %in% gene_ol$ensembl)
all_genes = as.data.frame(all_genes)
all_genes = merge(all_genes, gene_ol, by = 'ensembl')

write.csv(all_genes, file = file.path(working_dir, 'all_genes_lod.csv'), row.names = FALSE,
          quote = FALSE)









