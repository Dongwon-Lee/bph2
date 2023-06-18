# R script to classify genes based on their expression profiles
# across the four BP tissues 
#
# - Dongwon Lee
args<-commandArgs(trailingOnly=T)

adrenalfn<-args[1]
arteryfn<-args[2]
heartfn<-args[3]
kidneyfn<-args[4]

library('tidyverse')

d1<-read_tsv(adrenalfn, col_names=c('chrom', 'start', 'end', 'gene')) %>% 
    mutate(adrenal=TRUE)
d2<-read_tsv(arteryfn, col_names=c('chrom', 'start', 'end', 'gene')) %>%
    mutate(artery=TRUE)
d3<-read_tsv(heartfn, col_names=c('chrom', 'start', 'end', 'gene')) %>%
    mutate(heart=TRUE)
d4<-read_tsv(kidneyfn, col_names=c('chrom', 'start', 'end', 'gene')) %>%
    mutate(kidney=TRUE)

df<-d1 %>%
    full_join(d2) %>%
    full_join(d3) %>%
    full_join(d4)

df<-df %>%
    mutate(adrenal=ifelse(is.na(adrenal), FALSE, adrenal)) %>%
    mutate(artery=ifelse(is.na(artery), FALSE, artery)) %>%
    mutate(heart=ifelse(is.na(heart), FALSE, heart)) %>%
    mutate(kidney=ifelse(is.na(kidney), FALSE, kidney)) %>%
    mutate(ngene=adrenal+artery+heart+kidney)

df<-df %>% arrange(chrom, start)

write.table(df, 'bptissues_genes_combined_top10k_50kb.bed',
             col.names=F,
             row.names=F,
             quote=F,
             sep='\t')

df %>% filter(ngene==4) %>%
write.table('bptissues_genes_ubiq_top10k_50kb.bed',
             col.names=F,
             row.names=F,
             quote=F,
             sep='\t')

df %>% filter(ngene==1 & artery) %>%
write.table('bptissues_genes_artery_spec_top10k_50kb.bed',
             col.names=F,
             row.names=F,
             quote=F,
             sep='\t')

df %>% filter(ngene==1 & adrenal) %>%
write.table('bptissues_genes_adrenal_spec_top10k_50kb.bed',
             col.names=F,
             row.names=F,
             quote=F,
             sep='\t')

df %>% filter(ngene==1 & heart) %>%
write.table('bptissues_genes_heart_spec_top10k_50kb.bed',
             col.names=F,
             row.names=F,
             quote=F,
             sep='\t')

df %>% filter(ngene==1 & kidney) %>%
write.table('bptissues_genes_kidney_spec_top10k_50kb.bed',
             col.names=F,
             row.names=F,
             quote=F,
             sep='\t')

df %>% filter(ngene==2 & adrenal & heart) %>%
write.table('bptissues_genes_ad_he_top10k_50kb.bed',
             col.names=F,
             row.names=F,
             quote=F,
             sep='\t')

df %>% filter(ngene==2 & adrenal & artery) %>%
write.table('bptissues_genes_ad_ar_top10k_50kb.bed',
             col.names=F,
             row.names=F,
             quote=F,
             sep='\t')

df %>% filter(ngene==2 & adrenal & kidney) %>%
write.table('bptissues_genes_ad_ki_top10k_50kb.bed',
             col.names=F,
             row.names=F,
             quote=F,
             sep='\t')

df %>% filter(ngene==2 & heart & artery) %>%
write.table('bptissues_genes_he_ar_top10k_50kb.bed',
             col.names=F,
             row.names=F,
             quote=F,
             sep='\t')

df %>% filter(ngene==2 & heart & kidney) %>%
write.table('bptissues_genes_he_ki_top10k_50kb.bed',
             col.names=F,
             row.names=F,
             quote=F,
             sep='\t')

df %>% filter(ngene==2 & artery & kidney) %>%
write.table('bptissues_genes_ar_ki_top10k_50kb.bed',
             col.names=F,
             row.names=F,
             quote=F,
             sep='\t')

df %>% filter(ngene==3 & adrenal & heart & artery) %>%
write.table('bptissues_genes_ad_he_ar_top10k_50kb.bed',
             col.names=F,
             row.names=F,
             quote=F,
             sep='\t')

df %>% filter(ngene==3 & adrenal & heart & kidney) %>%
write.table('bptissues_genes_ad_he_ki_top10k_50kb.bed',
             col.names=F,
             row.names=F,
             quote=F,
             sep='\t')

df %>% filter(ngene==3 & adrenal & artery & kidney) %>%
write.table('bptissues_genes_ad_ar_ki_top10k_50kb.bed',
             col.names=F,
             row.names=F,
             quote=F,
             sep='\t')

df %>% filter(ngene==3 & heart & artery & kidney) %>%
write.table('bptissues_genes_he_ar_ki_top10k_50kb.bed',
             col.names=F,
             row.names=F,
             quote=F,
             sep='\t')
