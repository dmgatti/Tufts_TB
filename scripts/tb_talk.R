options(stringsAsFactors = FALSE)
library(qtl2)
library(survival)
library(survminer)
library(readxl)
library(tidyverse)
base_dir = '/media/dmgatti/hdb/projects/TB'
pheno_file = 'Beamer_TB_Phenotype_data_20201124.xlsx'

# Read in phenotype data
pheno = readxl::read_xlsx(file.path(base_dir, 'data', 'phenotypes', pheno_file), 
                          sheet = 'Batch & phenotype data', na = 'na')

pheno = pheno %>% 
  select(mouse       = `Mouse #`,
         coat        = `Coat color`,
         strain      = `Strain`,
         gen         = `DOB or generation`,
         expt        = `Experiment (Batch)`,
         mach_type   = `Type of aerosol machine`,
         gigamuga    = `GigaMUGA`,
         cage        = `Cage number`,
         aero_run    = `Aerosol run no`,
         mtb_dose    = `Mtb initial dose`,
         susc_class  = `Susceptibility Class`,
         il10_assay    = `IL-10 Assay Source`,
         il12_assay    = `IL-12 Assay Source`,
         euth_day    = `Day of Euthanasia`,
         euth        = `Euthanized due to pulmonary TB`,
         pct_wt_loss = `% wt loss cf peak`,
         mtb_burden  = `Lung Mtb burden`,
         lung_cxcl5  = `Lung CXCL5 pg/ml`,
         lung_cxcl2  = `Lung CXCL2 pg/ml`,
         lung_cxcl1  = `Lung CXCL1 pg/ml`,
         lung_tnf    = `Lung TNF pg/ml`,
         lung_mmp8   = `Lung MMP8 pg/ml`,
         lung_s100a8 = `Lung S100A8 pg/mL`,
         lung_s100a9 = `Lung S100A9 pg/mL`,
         lung_resistin = `Lung resistin pg/ml`,
         lung_il1a     = `Lung IL-1a pg/mL`,
         lung_il1b     = `Lung IL-1B pg/mL`,
         lung_il6      = `Lung IL-6 pg/mL`,
         lung_ifng     = `Lung IFNg pg/ml`,
         lung_il12     = `Lung IL-12 pg/ml`,
         lung_il17     = `Lung IL-17 pg/ml`,
         lung_il10     = `Lung IL-10 pg/ml`,
         lung_il4      = `Lung IL-4 pg/mL`,
         lung_vegf     = `Lung VEGF pg/ml`,
         lung_ccl3     = `Lung CCL3 pg/ml`,
         lung_ccl4     = `Lung CCL4 pg/ml`,
         lung_ccl5     = `Lung CCL5 pg/ml`,
         lung_anti_esat6_cfp10 = `Lung anti-ESAT6:CFP10 pg/mL`,
         lung_anti_mtb_cw_igg  = `Lung anti-M.tb CW IgG \"ng/mL\"`,
         lung_anti_mtb_cfp_igg = `Lung anti-M.tb CFP IgG \"ng/mL\"`,
         lung_anti_crypto_igg  = `Lung anti-crypto IgG \"ng/mL\"`,
         blood_esat6_ifng      = `Blood ESAT-6 specific IFNg pg/ml`,
         blood_esat6_tnf       = `Blood ESAT-6 specific TNF pg/ml`,
         blood_esat6_il12p70   = `Blood ESAT-6 specific IL-12p40/p70 pg/ml`,
         blood_esat6_il10      = `Blood ESAT-6 specific IL-10 pg/ml`,
         blood_esat6_il2       = `Blood ESAT-6 specific IL-2 pg/ml`,
         blood_pct_neut        = `% blood neutrophils`,
         serum_cxcl5           = `Serum CXCL5 pg/ml`,
         serum_cxcl2           = `Serum CXCL2 pg/ml`,
         serum_cxcl1           = `Serum CXCL1 pg/ml`,
         serum_tnf             = `Serum TNF pg/ml`,
         seum_mmp8             = `Serum MMP8 pg/ml`,
         serum_s100a8          = `Serum S100A8 pg/ml`,
         serum_resistin        = `Serum resistin pg/ml`,
         serum_il1b            = `Serum IL-1B pg/ml`,
         serum_il6             = `Serum IL-6 pg/mL`,
         serum_ifng            = `Serum IFNg pg/mL`,
         serum_il12p70         = `Serum IL-12p70 pg/ml`,
         serum_il17            = `Serum IL-17 pg/ml`,
         serum_il10            = `Serum IL-10 pg/mL`,
         serum_il4             = `Serum IL-4 pg/mL`,
         serum_vegf            = `Serum VEGF pg/ml`,
         serum_ccl3            = `Serum CCL3 pg/ml`,
         serum_ccl4            = `Serum CCL4 pg/ml`,
         serum_ccl5            = `Serum CCL5 pg/ml`,
         serum_anti_mtb_cw_igg = `Serum anti-M.tb CW IgG \"ng/mL\"`,
         serum_anti_mtb_cfp_igg = `Serum anti-M.tb CFP IgG \"ng/mL\"`,
         serum_anti_crypto_igg = `Serum anti-crypto IgG \"ng/mL\"`)

pheno = pheno %>% 
          mutate(euth = if_else(str_detect(euth, '^n'), 'n', euth),
                 euth = if_else(str_detect(euth, '^y'), 'y', euth))

# Get BL6 and plot survival.
founders = pheno %>% 
             filter((strain == 'J:DO' | strain == 'C57BL/6J') & !is.na(strain))

surv = Surv(founders$euth_day, event = founders$euth == 'y')

surv = data.frame(surv, strain = founders$strain)

fit = survfit(surv ~ strain, data = surv)

ggsurvplot(fit, data = surv, risk.table = TRUE)

png(file.path(base_dir, 'figures', 'tb_survival.png'), width = 2000, height = 1200, res = 256)
par(bg = 'grey90', mgp = c(2, 0.5, 0))
cols = c('black', 'tomato2') 
plot(fit, col = NA, lwd = 1.5,  las = 1, main = 'TB survival',
     xlab = 'Days', ylab = '% Surviving')
abline(h = 0:5/5, col = 'white')
abline(v = 0:3 * 100, col = 'white')
lines(fit, col = cols, lwd = 1.5)
legend('bottomleft', legend = c('C57BL/6J', 'J:DO'), col = cols, lwd = 1.5)
dev.off()


