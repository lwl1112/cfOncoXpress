# Weiling Li (wel4007@med.cornell.edu)
# ISMB2025 proceedings: Fig.6

setwd('/Users/wel4007/Documents/Github/cfOncoPath/GSEAdata/output/')
library(readr)
library(forcats)
library(dplyr)
library(ggplot2)
library(ggpubr)

### CRPC-AR
#import TSV file into data frame
log2TPMminmax_AR_pbmc_pos <- read_tsv('AR_hallmark_CRPC_pathways.Gsea.1737563629877/gsea_report_for_1_1737563629877.tsv')
significantAR <- log2TPMminmax_AR_pbmc_pos[log2TPMminmax_AR_pbmc_pos$'FDR q-val'<=0.1, ]
colorsAR <- c(replicate(dim(significantAR)[1], 'black'))
ARdata<- data.frame( name=c(significantAR$NAME),val=c(significantAR$NES))

ar=ARdata %>%
  mutate(name = fct_reorder(name, val)) %>%
  ggplot( aes(x=name, y=val)) +
  geom_bar(stat="identity", fill=colorsAR, alpha=.6, width=.4) +
  coord_flip() +
  ylab("Normalized Enrichment Score") + xlab("Pathway")+ylim(-3,3)+
  ggtitle("CRPC-AR vs. healthy controls")+
  theme_bw()+theme(plot.title=element_text(hjust=0.5))#,text = element_text(size = 15))
  
### CRPC-NE  
log2TPMminmax_NE_pbmc_pos <- read_tsv('NE_hallmark_CRPC_pathways.Gsea.1737397149906/gsea_report_for_2_1737397149906.tsv')
significantNE <- log2TPMminmax_NE_pbmc_pos[log2TPMminmax_NE_pbmc_pos$'FDR q-val'<=0.1, ]
colorsNE <- c(replicate(dim(significantNE)[1], 'black'))
NEdata<- data.frame( name=c(significantNE$NAME),val=c(significantNE$NES))

ne=NEdata %>%
  mutate(name = fct_reorder(name, val)) %>%
  ggplot( aes(x=name, y=val)) +
  geom_bar(stat="identity", fill=colorsNE, alpha=.6, width=.4) +
  coord_flip() +
  ylab("Normalized Enrichment Score") + xlab("Pathway")+ylim(-3,3)+
  ggtitle("CRPC-NE vs. healthy controls")+
  theme_bw()+theme(plot.title=element_text(hjust=0.5))#,text = element_text(size = 15))



ggarrange(ar, ne, ncol=1)


############


ggsave(
   filename='Figure6.jpg',
   plot = last_plot(),
   width = 20, height = 20
 )
