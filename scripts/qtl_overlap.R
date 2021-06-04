################################################################################
# Plot of QTL overlap per chr.
# Daniel Gatti
# 2021-02-14
# dmgatti@coa.edu
################################################################################
library(GenomicRanges)
library(tidyverse)

# Set up base directories.
base_dir = '/media/dmgatti/hdb/projects/TB/'
qtl_dir  = file.path(base_dir, 'results', 'qtl2', 'gen_factor2')
fig_dir  = file.path(base_dir, 'figures')
qtl_file = file.path(qtl_dir, 'tb_qtl_peaks.csv')

peaks = read_csv(qtl_file)

png(file.path(fig_dir, 'qtl_overlap.png'), width = 800, height = 800)
print(ggplot(peaks) +
        geom_segment(aes(x = ci_lo, xend = ci_hi, y = lodcolumn, yend = lodcolumn)) +
        geom_point(aes(pos, lodcolumn)) +
        facet_wrap(~chr)) +
        labs(title = 'DO TB QTL Overlap', x = '', y = '')
dev.off()

peaks_gr = GRanges(seqnames = peaks$chr, 
                   ranges   = IRanges(start = peaks$ci_lo, end = peaks$ci_hi))
ol = findOverlaps(peaks_gr, peaks_gr)

result = select(peaks[queryHits(ol),], lodcolumn:ci_hi) %>% 
            rename_with(.fn = ~str_c(., '_q'), .cols = everything()) %>% 
            bind_cols(select(peaks[subjectHits(ol),], lodcolumn:ci_hi) %>% 
                      rename_with(.fn = ~str_c(., '_s'), .cols = everything())) %>% 
            filter(lodcolumn_q != lodcolumn_s) %>% 
            arrange(chr_q, pos_q)




