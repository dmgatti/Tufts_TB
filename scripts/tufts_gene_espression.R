# Tufts Lung Gene Expression
options(stringsAfFactors = FALSE)
library(readxl)
library(tidyverse)
library(broom)

base_dir   = '/home/dmgatti/Documents/data/TB/'
pheno_dir  = str_c(base_dir, 'data/phenotypes/')
report_dir = str_c(base_dir, 'data/expression/reports/')

rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

# Read in phenotype data.
pheno = read_xlsx(str_c(pheno_dir, 'JDO lung phenotypes 081117.xlsx')) %>%
          select(`Mouse #`, `Mtb initial dose`, `Survival`:`Lung S100A8 (pg/ml)`) %>%
          rename(mouse = `Mouse #`,
                 dose  = `Mtb initial dose`,
                 survival = Survival,
                 euth     = `Euthanized due to morbidity`,
                 wt_loss  = `% wt loss`,
                 mtb_burden = `Lung Mtb burden`,
                 cxcl5      = `Lung CXCL5 pg/ml`,
                 cxcl2      = `Lung CXCL2 pg/ml`,
                 cxcl1      = `Lung CXCL1 pg/ml`,
                 perc_normal = `% Normal lung`,
                 ifng       = `Lung IFNg pg/ml`,
                 tnf        = `Lung TNF pg/ml`,
                 il12       = `Lung IL-12 pg/ml`,
                 il10       = `Lung IL-10 pg/ml`,
                 mmp8       = `Lung MMP8 (pg/ml)`,
                 vegf       = `Lung VEGF (pg/ml)`,
                 s100a8     = `Lung S100A8 (pg/ml)`) %>%
          filter(dose > 0) %>%
          gather(pheno, value, -(mouse:euth)) %>%
          mutate(value = parse_number(value),
                 value = log1p(value)) %>%
          spread(pheno, value) %>%
          mutate(mouse = as.character(mouse))

# Read in expression data.
expr1 = read_xlsx(str_c(report_dir, '2017-05-26_Beamer_analysis.xlsx'), sheet = '2017-05-26_Beamer_analysis', 
                  range = 'A2:CM25208') %>%
        select(`Brainarray probeset ID`:`KEGG Pathway(s)`, starts_with('Resis'), starts_with('Susc'), starts_with('Super')) %>%
        gather(sample, expr, -(`Brainarray probeset ID`:`KEGG Pathway(s)`)) %>%
        separate(sample, into = c('group', 'sample'))

expr2 = read_xlsx(str_c(report_dir, '2016-12-09_Beamer_analysis (2) - gb.xlsx'), sheet = 'analysis', 
                  range = 'A2:CG25208') %>%
        select(`Brainarray probeset ID`:`KEGG Pathway(s)`, starts_with('Resis'), starts_with('Susc'), starts_with('Super')) %>%
        gather(sample, expr, -(`Brainarray probeset ID`:`KEGG Pathway(s)`)) %>%
        separate(sample, into = c('group', 'sample'))

expr = bind_rows(expr1, expr2)
rm(expr1, expr2)

# Reshape data.:
annot = expr %>%
          select(`Brainarray probeset ID`:`KEGG Pathway(s)`)
expr = expr %>%
         select(sample, group, `Mouse Entrez Gene ID`, expr) %>%
         rename(entrez   = `Mouse Entrez Gene ID`) %>%
         spread(entrez, expr)

expr = full_join(select(pheno, mouse, dose), expr, by = c('mouse' = 'sample')) %>%
         select(-group) %>%
         gather(entrez, expr, -mouse, -dose) %>%
         group_by(dose) %>%
         nest()

#### STOPPED HERE ####



pheno = subset(pheno, pheno$mouse %in% expr$sample)

expr = expr %>%
         arrange(sample)
pheno = pheno %>%
          arrange(mouse)

stopifnot(all(pheno$mouse == expr$sample))

exprmat  = as.matrix(expr[,-(1:2)])
rownames(exprmat) = expr$sample
exprmat = apply(exprmat, 2, rankZ)
phenomat = as.matrix(pheno[,-(1:4)])
rownames(phenomat) = pheno$mouse
phenomat = apply(phenomat, 2, rankZ)

pecor = cor(phenomat, exprmat, use = 'pair')
df    = nrow(phenomat) - 2
pet   = pecor * sqrt((df) / (1.0 - pecor^2))
pep   = 2 * pt(abs(pet), df = df, lower.tail = FALSE)
peq   = matrix(p.adjust(pep), nrow(pep), ncol(pep), dimnames = dimnames(pep))

rowSums(peq < 0.05)
# Worried about large number of genes correlated with vegf.
# It seems to be due to strong correlation with dose.

g = apply(peq, 1, function(z) { colnames(peq)[z < 0.05] })

gl = g[[1]]
for(i in 2:length(g)) {
  gl = intersect(gl, g[[i]])
  print(length(gl))
}


