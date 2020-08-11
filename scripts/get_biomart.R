# Get Go annotations for genes.
library(biomaRt)

listMarts()

mart = useMart("ENSEMBL_MART_ENSEMBL")

listDatasets(mart)

mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

afilt = listFilters(mart)

aattr = listAttributes(mart)

genes = select(x = mart,
               columns = c("chromosome_name", "start_position", "end_position", "strand", "ensembl_gene_id", "external_gene_name"), 
               keytype = "chromosome_name", 
               keys = "15")
genes = subset(genes, start_position > 20e6 & end_position <= 30e6)

attrib = 
filt = 
values = 


goids = select(x = mart,
               columns = c('ensembl_gene_id', 'external_gene_name', 'go_id', 'name_1006', 'definition_1006'),
               keytype = 'ensembl_gene_id', 
               keys = genes$ensembl_gene_id)





