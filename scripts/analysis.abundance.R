library(grDevices)
library(ggplot2)
library(plyr)
library(reshape2)
library(gplots)
library(dplyr)
library(tidyr)
library(readr)
source("analyze.functions.R")



boxplots.prevalence_correlation = function()
{  
  set.seed(0)
  growth_matrix = read_tsv("../data/supplementary/growth_matrix.tab")
  colnames(growth_matrix) = gsub("^M", "", colnames(growth_matrix))
  growth_matrix = growth_matrix[!grepl("15|16", colnames(growth_matrix))]

  growth = growth_matrix %>% gather(medium, growth, 2:ncol(growth_matrix))
  abundant_species = read.table("../data/abundant_species.tab", sep="\t", quote="", header=F, stringsAsFactors=F, na.strings="")$V1
  species = read_tsv("../data/organisms.tab")
  species = species %>% filter(species %in% abundant_species) 
  abundances = read_tsv("../data/screenG_tax_info_specI_clusters.tab")
  
  specI <- species %>% 
    left_join(abundances %>% select(NT_code, SpecI_ID), by=c("nid"="NT_code")) %>% 
    select(species, SpecI_ID)
  specI_growth = growth_matrix %>% left_join(specI, by=c("Species"="species")) %>% filter(!is.na(SpecI_ID))
  
  specI_growth_matrix = as.matrix(subset(specI_growth,,-c(Species, SpecI_ID)))
  rownames(specI_growth_matrix) = specI_growth$SpecI_ID
  specI_growth_matrix[specI_growth_matrix==0] = NA
  specI_growth_matrix[apply(specI_growth_matrix, 2, rank) > 30] = NA
  
  all_abundances = read.delim("../data/screenG_rel_ab_specI_clusters.tab")
  all_abundances = all_abundances[ specI_growth$SpecI_ID, ]
  
  all_abundances_0na <- all_abundances
  all_abundances_0na[ all_abundances == 0 ] <- NA
  
  cor_matrix.resample.long = data.frame()
  for(sample.i in  1:20) {
    all_abundances_0na.resample = apply(all_abundances_0na, 2, sample)
    cor_matrix.resample = cor(all_abundances_0na.resample, specI_growth_matrix, use="p", method = "s")
    cor_matrix.resample = cor_matrix.resample[,1:(ncol(cor_matrix.resample)-1)]  # remove M16
    cor_matrix.resample.long = rbind(cor_matrix.resample.long, melt(cor_matrix.resample))
  }
  
  colnames(cor_matrix.resample.long) = c("subject", "media", "correlation")
  cor_matrix.resample.long$dataset = "background"
  
  cor_matrix = cor(all_abundances_0na, specI_growth_matrix, use="p", method = "s")
  cor_matrix = cor_matrix[,1:(ncol(cor_matrix)-1)]  # remove M16
  
  cor_matrix.long = melt(cor_matrix)
  colnames(cor_matrix.long) = c("subject", "media", "correlation")
  cor_matrix.long$dataset = "signal"
  cor_matrix.long = rbind(cor_matrix.long, cor_matrix.resample.long)
  
  cor_matrix.long.pvalues = data.frame()
  for(m in colnames(cor_matrix)) 
  {
    if(all(is.na(cor_matrix.long$correlation[cor_matrix.long$media==m]))) next
    m.test = t.test(cor_matrix.long$correlation[cor_matrix.long$media==m], cor_matrix.resample.long$correlation[cor_matrix.resample.long$media==m])
    cor_matrix.long.pvalues = rbind(cor_matrix.long.pvalues, data.frame(media=m, pvalue=m.test$p.value, pvalue.adj=m.test$p.value * ncol(cor_matrix)))
  }
  
  cor_matrix.long.pvalues.export = cor_matrix.long.pvalues
  cor_matrix.long.pvalues.export$media = gsub("^([0-9])", "M\\1", cor_matrix.long.pvalues.export$media)
  write.table(cor_matrix.long.pvalues.export, file="../data/supplementary/abundance_correlations.tab", sep="\t", quote=F, row.names=F, na="")
  
  
  media_order = c("mGAM", "WCA", "BHI++", "GMM", "9", "8", "11", "10", "7", "5", "4", "3", "2", "14", "13", "1")
  cor_matrix.long.pvalues$media = factor(cor_matrix.long.pvalues$media, media_order)
  cor_matrix.long$media = factor(cor_matrix.long$media, media_order)
  
  pdf("../report/media_correlation.pdf", paper="a4r", height=8.27, width=11.69)
  ggplot(cor_matrix.long) +
    geom_boxplot(aes(x=media, y=correlation, fill=dataset)) +
    geom_segment(aes(x=as.numeric(media)-0.3, xend=as.numeric(media)+0.3), y=0.78, yend=0.78, data=subset(cor_matrix.long.pvalues, pvalue < 0.05)) +
    geom_text(aes(x=media, label="**"), y=0.8, data=subset(cor_matrix.long.pvalues, pvalue < 0.05)) +
    geom_hline(yintercept=0, lty=2) +
    scale_fill_manual(values=c(background="#FFFFFF", signal="#D4D4D4")) +
    scale_y_continuous(breaks=seq(-1, 1, 0.25)) +
    theme_slim() +
    coord_flip() +
    labs(x="", y="correlation")
  dev.off()
}
  
