library(ggplot2)
library(plyr)
library(RColorBrewer)
library(gridExtra)
source("utils/ggplot.R")

boxplots_and_scater.pH_overview = function()
{
  media.annotation = read.delim("../data/media_names.tab", sep="\t", na="", quote="", stringsAsFactors=F, header=T)
  media.annotation = media.annotation[order(media.annotation$Order),]
  
  organisms = read.table("../data/organisms.tab", sep="\t", quote="", header=T, stringsAsFactors=T, na.strings="")
  
  pH_breaks = seq(3.75, 8.25, 0.5)
  pH_breaks2 = c(3.5, 5, 6, 7, 8.5)
  pH_labels2 = c("(3.5,5]"="<=5", "(5,6]"="5-6", "(6,7]"="6-7", "(7,8.5]"=">7")
  curves.a = read.table("../data/curves_annotation.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="")
  curves.a = subset(curves.a, !is.na(Species) & is.na(ConditionSpecies))
  curves.a$MaxOD = curves.a$MaxOD #-curves.a$BlankOD
  curves.pH = ddply(subset(curves.a, !is.na(pH_48.mean), c(Species, Media, pH_48.mean)), .(Species, Media), summarize, pH=median(pH_48.mean, na.rm=T))
  curves.a = ddply(curves.a, .(Species, Media), summarize, MaxOD=median(MaxOD, na.rm=T))
  curves.a = merge(curves.a, curves.pH, by=c("Species", "Media"))
  curves.a = merge(curves.a, media.annotation[,c("ShortName", "Type")], by.x="Media", by.y="ShortName")
  curves.a = merge(curves.a, organisms[,c("species", "genus", "phylum")], by.x="Species", by.y="species")
  curves.a = ddply(curves.a, .(phylum), mutate, phylum=paste0(phylum, " (", length(phylum), ")"))
  curves.a = ddply(curves.a, .(genus), mutate, genus=paste0(genus, " (", length(genus), ")"))
  curves.a$genus = factor(curves.a$genus, unique(as.character(curves.a$genus)))
  curves.a$phylum = factor(curves.a$phylum, unique(as.character(curves.a$phylum)))
  colnames(curves.a)[colnames(curves.a)=="Type"] = "MediaType"
  curves.a$pH_group = cut(curves.a$pH,breaks=pH_breaks)
  curves.a$pH_group = factor(pH_breaks[as.numeric(curves.a$pH_group)] + 0.25)
  curves.a$pH_group2 = as.character(cut(curves.a$pH,breaks=pH_breaks2))
  curves.a$pH_group2 = factor(pH_labels2[curves.a$pH_group2], levels=pH_labels2)
  curves.a$color = colorRampPalette(RColorBrewer::brewer.pal(9, "Greens"))(30)[as.numeric(factor(curves.a$Species))]
  
  pdf("../report/pH_growth_effect.pdf", height=5, width=6.941)
  curves.a.phylum_c = ddply(curves.a, .(phylum, pH_group2), summarize, Count=length(phylum), MaxOD=median(MaxOD, na.rm=T))
  
  curves.a_fit = with(curves.a, expand.grid(MediaType=unique(MediaType), pH=c(T,F)))
  curves.a_fit$text = curves.a_fit$p.value = curves.a_fit$rho = curves.a_fit$slope = curves.a_fit$x = curves.a_fit$y = NA
  for(i in 1:nrow(curves.a_fit)) {
    curves.a.d = subset(curves.a, MediaType==curves.a_fit$MediaType[i] & (pH>=7) == curves.a_fit$pH[i])
    curves.a.d_cor = cor.test(curves.a.d$pH, curves.a.d$MaxOD)
    curves.a.d_lm =  lm(MaxOD ~ pH, curves.a.d)
    curves.a_fit$slope[i] = curves.a.d_lm$coefficients["pH"]
    curves.a_fit$rho[i] = curves.a.d_cor$estimate
    curves.a_fit$p.value[i] = curves.a.d_cor$p.value
    curves.a_fit$x[i] = 4
    curves.a_fit$y[i] = ifelse(curves.a_fit$pH[i], 1.3, 1.2)
    curves.a_fit$n[i] = nrow(curves.a.d)
    curves.a_fit$text[i] = with(curves.a_fit[i,], paste0(
      "R2 (", ifelse(curves.a_fit$pH[i], ">=7", "<7"), ")=",
      round(rho,2), 
      "(p<=", sprintf("%1.0e", p.value), ", ", round(p.value,3), "; n=", n,")"
    ))
  }
  gridExtra::grid.arrange(
    ggplot(curves.a, aes(x=pH, y=MaxOD)) +
      geom_point(alpha=0.3, size=0.75) +
      geom_smooth(aes(x=pH), method="lm", formula=y~splines::bs(x, knots=c(3, 7, 8), degree=1)) +
      geom_text(aes(label=text, x=x, y=y), data=curves.a_fit, hjust=0, size=2) +
      geom_vline(xintercept=7, lty=2, color="#FF0000") +
      facet_wrap(~MediaType, ncol=4) +
      coord_cartesian(ylim=c(0,1.5)) +
      theme_bw(base_size=16)
  ,  
    ggplot(curves.a) +
      geom_boxplot(aes(x=phylum, y=MaxOD, fill=pH_group2), alpha=0.3, outlier.size = NA) +
      geom_point(aes(x=phylum, y=MaxOD, fill=pH_group2), position=position_jitterdodge(dodge.width=0.75, jitter.width=0.8), shape=21, data=curves.a, size=0.75) +
      geom_text(aes(x=phylum, y=MaxOD, fill=pH_group2, label=Count), position=position_dodge(width=0.75), data=curves.a.phylum_c) +
      scale_fill_discrete(drop=FALSE) +
      scale_color_discrete(drop=FALSE) +
      scale_x_discrete(drop=FALSE) +
      labs(x="") +
      theme_bw(base_size=16) +
      guides(fill=FALSE)
  )
  dev.off()
}
