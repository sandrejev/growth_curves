library(gridExtra)
library(plyr)
library(ggplot2)
library(reshape2)
library(grid)

passages_reproducibility = function() {
  media_names = read.delim("../data/media_names.tab")
  media_names$MediaType = media_names$Type
  
  res.everything = read.delim("../data/curves_annotation.tab", stringsAsFactors=F, na="")
  media.order = rev(c("1", "13", "14", "15 A", "15 B", "16", "2", "3", "4", "5", "7", "10", "11", "8", "9", "GMM", "BHI++", "WCA", "mGAM"))
  

  th = 0.15
  res.everything.f = subset(res.everything, Volume==9 & !is.na(Species) & !is.na(Media))
  res.everything.f$MaxODRaw = res.everything.f$MaxOD
  res.everything.f$MaxOD = res.everything.f$MaxOD - res.everything.f$BlankOD
  res.everything.f$MaxODNegative = res.everything.f$MaxOD - res.everything.f$BlankOD
  res.everything.f$MaxOD[res.everything.f$MaxOD < 0] = 0
  res.everything.f = ddply(res.everything.f, .(File, Species, Media, Passage), summarize, MaxOD=mean(MaxOD), MaxODNegative=mean(MaxODNegative), MaxODRaw=mean(MaxODRaw), MaxTime=mean(MaxTime), BlankOD=mean(BlankOD), Class=Class[1])
  res.combinations = merge(subset(res.everything.f, Passage==1), subset(res.everything.f, Passage==3), by=c("Species", "Media"), suffixes=c(".1", ".3"))
  res.combinations = merge(res.combinations, media_names[,c("ShortName", "MediaType")], by.x="Media", by.y="ShortName")
  res.combinations$unreproducible.3 = with(res.combinations, ifelse(MaxOD.1>=th & MaxOD.3<=th & MaxOD.3 < MaxOD.1-th , "unreproducible", "reproducible"))
  res.combinations$unreproducible.1 = with(res.combinations, ifelse(MaxOD.3>=th & MaxOD.1<=th & MaxOD.1 < MaxOD.3-th , "unreproducible", "reproducible"))
  res.combinations$unreproducible = with(res.combinations, ifelse(unreproducible.1=="unreproducible" | unreproducible.3=="unreproducible", "unreproducible", "reproducible"))
  res.combinations$MaxOD.diff = with(res.combinations, MaxODNegative.3-MaxODNegative.1)
  res.combinations$Media = factor(res.combinations$Media, media.order)
  res.combinations = res.combinations[order(res.combinations$Media),]
  res.combinations.unrep = subset(res.combinations, unreproducible.3=="unreproducible")
  res.combinations.unrep = res.combinations.unrep[order(res.combinations.unrep$MediaType),]
  
  x = res.combinations
  x$effect = "none"
  x$effect[x$unreproducible.3=="unreproducible"] = "reduction"
  x$effect[x$unreproducible.1=="unreproducible"] = "increase"
  write.table(x[,c("Media", "Species", "effect")], file="../data/generated/passage_effect.tab", sep="\t", quote=F, row.names=F, na="")
  
  cor.res = with(subset(res.combinations, unreproducible!="anything"), cor.test(MaxOD.1, MaxOD.3))
  cor.pval = paste0("1e", round(log10(cor.res$p.value)+1, 0))
  cor.n = sum(res.combinations$unreproducible!="anything")
  cor.R2 = round(cor.res$estimate^2, 2)
  
  set.seed(2)
  p.correlation = ggplot(res.combinations) +
    geom_abline(intercept=0, slope=1, color="#999999", size=2) +
    geom_point(aes(MaxOD.1, MaxOD.3, color=factor(MediaType), shape=unreproducible), size=3, alpha=0.7) +
    geom_text(x=0.25, y=1.75, label=paste0("R2=", cor.R2, " (p<", cor.pval, ", n=", cor.n, ")"), hjust=0, vjust=0, size=5, data=data.frame(x=c(1))) +
    geom_smooth(aes(MaxOD.1, MaxOD.3), method="lm") +
    geom_path(aes(x=x, y=y), data=data.frame(x=c(-2, th, 2*th, Inf), y=c(-2-th, 0, th, th)), color="#FF0000") +
    geom_path(aes(x=x, y=y), data=data.frame(x=c(-2-th, 0, th, th), y=c(-2, th, 2*th, Inf)), color="#FF0000") +
    labs(x=expression(MaxOD[Passage1]), y=expression(MaxOD[Passage3])) +
    coord_cartesian(xlim=c(0,2), ylim=c(0,2)) +
    theme_bw(base_size=16)  +
    theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
    scale_shape_manual(values=c(reproducible=21, unreproducible=22))
  p.correlation
  
  res.combinations.unrep_text = ddply(res.combinations.unrep, .(Species, MediaType), function(z) {
    d = ddply(z, .(Media), summarize, 
              count=length(Media),
              Media=gsub(" ", "", Media[1]),
              details=paste0(ifelse(grepl("^[0-9]", Media), "M", ""), Media, " [", paste(round(MaxOD.1, 2), collapse=" "), "]"),
              MaxOD.1=mean(MaxOD.1))
    
    data.frame(count=sum(d$count), MedianMaxOD.1=median(d$MaxOD.1), details=paste0(d$details, collapse=", "))
  })
  res.combinations.unrep_text = ddply(res.combinations.unrep_text, .(Species), mutate, cumcount=cumsum(count), cumcount2=cumcount-count+0.05)
  p.unrep_cases = ggplot(res.combinations.unrep) +
    geom_bar(aes(x=Species, y=1, fill=MediaType), stat="identity", position="stack") +
    geom_text(aes(x=Species, y=cumcount2, color=MediaType, label=details), stat="identity", hjust=0, size=4, data=res.combinations.unrep_text) +
    scale_fill_manual(values=c(Defined="#e6ee9c", Minimal="#c5e1a5", Mucin="#a5d6a7", Rich="#80cbc4")) +
    scale_color_manual(values=c(Defined="#827717", Minimal="#33691e", Mucin="#1b5e20", Rich="#004d40")) +
    scale_y_continuous(breaks=0:max(ddply(res.combinations.unrep, .(Species), nrow)$V1)) +
    labs(y="Not reproduced in passage 3 (count)", x="") +
    theme_bw(base_size = 30) +
    theme(panel.grid.major.y=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor.y=element_blank(), panel.grid.minor.x=element_blank(), legend.position = c(1, 1), legend.justification = c(1, 1)) +
    coord_flip()
  
  p4 = tableGrob(res.combinations.unrep[,c("Species", "Media", "MaxOD.1", "MaxOD.3")], rows=NULL, theme=ttheme_minimal(core=list(fg_params=list(cex=0.5)), colhead=list(fg_params=list(cex=0.5))))
  p4 = gtable::gtable_add_grob(p4, grobs = segmentsGrob(name='segment', y1 = unit( 0, 'npc' ), gp=gpar(lty=1, lwd=1)), t=1, l=1, r=ncol(p4))
  p4 = gtable::gtable_add_grob(p4, grobs = segmentsGrob(name='segment', x1 = unit( 0, 'npc' ), gp=gpar(lty=1, lwd=1)), t=1, l=2, r=2, b=nrow(p4))
  
  
  res.combinations.f = subset(res.combinations, unreproducible!="X")
  res.combinations.pval = ddply(res.combinations.f, .(Media), function(z) {
    data.frame(
      pvalue.paired=wilcox.test(z$MaxOD.3, z$MaxOD.1, paired=T)$p.value,
      n=nrow(z)
    )
  })
  subset(res.combinations.pval, pvalue.paired<0.05)
  res.combinations.pval$text = with(res.combinations.pval, ifelse(pvalue.paired<0.05, ifelse(pvalue.paired<0.01, paste0("p<",paste0("10", ceiling(log10(pvalue.paired)))), paste0("p=",round(pvalue.paired, 2))), ""))
  res.combinations.pval$pvalue.paired.adj = p.adjust(res.combinations.pval$pvalue.paired)
  write.table(res.combinations.pval, file="../data/generated/media_effect_on_passages.tab", sep="\t", quote=F, row.names=F, na="")
  
  p.media_effect = ggplot(res.combinations.f) +
    geom_boxplot(aes(x=Media, y=MaxOD.diff), fill="#FEEC66") + 
    #geom_text(aes(x=Media, label="*"), y=1.38, data=subset(res.combinations.pval, pvalue.paired < 0.05 & pvalue.paired.adj > 0.05)) +
    #geom_text(aes(x=Media, label="**"), y=1.38, data=subset(res.combinations.pval, pvalue.paired < 0.05 & pvalue.paired.adj < 0.05)) +
    geom_text(aes(x=Media, label=text), y=1.38, data=subset(res.combinations.pval, pvalue.paired < 0.05), hjust=0) +
    geom_segment(aes(x=as.numeric(Media)-0.3, xend=as.numeric(Media)+0.3), y=1.35, yend=1.35, data=subset(res.combinations.pval, pvalue.paired < 0.05)) +
    geom_rect(aes(xmin=as.numeric(Media)-0.45, xmax=as.numeric(Media)+0.45), ymin=1, ymax=1.3, fill="#FFFFFFBB", color="#CCCCCC", data=res.combinations.pval) +
    geom_segment(aes(x=as.numeric(Media)-0.45, xend=as.numeric(Media)+0.45), y=1, yend=1.3, color="#666666AA", data=res.combinations.pval) +
    geom_point(aes(x=as.numeric(Media)-0.35+MaxOD.3/2, y=MaxOD.1/5+1.025), size=0.1) +
    theme_light(base_size=30) + 
    theme(panel.grid.major.y=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor.y=element_blank(), panel.grid.minor.x=element_blank(), legend.position = c(1, 1), legend.justification = c(1, 1)) +
    scale_y_continuous(breaks=seq(-1, 1, by=0.5)) +
    labs(y=expression(MaxOD[Passage3]-MaxOD[Passage1]), x="") +
    coord_flip()
  p.media_effect
  
  
  growth.table = with(res.combinations, table(
    Passage1=factor(ifelse(MaxOD.1>=th | unreproducible.1=="reproducible", "Growth", "No Growth")),
    Passage3=factor(ifelse(MaxOD.3>=th | unreproducible.3=="reproducible", "Growth", "No Growth"))
  ))
  growth.table.df = dcast(data.frame(growth.table), Passage3 ~ Passage1, value.var="Freq")
  colnames(growth.table.df)[1] = ""
  p.table = tableGrob(growth.table.df, rows=NULL, theme=ttheme_minimal(core = list(fg_params=list(cex=2.5)), colhead=list(fg_params=list(cex=2.5))))
  p.table = gtable::gtable_add_grob(p.table, grobs = segmentsGrob(name='segment', y1 = unit( 0, 'npc' ), gp=gpar(lty=1, lwd=1)), t=1, l=1, r=ncol(p.table))
  p.table = gtable::gtable_add_grob(p.table, grobs = segmentsGrob(name='segment', x1 = unit( 0, 'npc' ), gp=gpar(lty=1, lwd=1)), t=1, l=2, r=2, b=nrow(p.table))
  
  pdf("../report/passages_reproducibility.pdf", height=20, width=24)
  grid.arrange(
    arrangeGrob(p.correlation, p.unrep_cases, heights=c(0.7, 0.3)),
    arrangeGrob(p.media_effect), 
    widths=c(0.5, 0.4)
  )
  grid.arrange(arrangeGrob(p.table))
  dev.off()
  
  
  
  
  ###################################
  # Unreproducible (3) curves
  ###################################
  
  curves = read.delim("../data/curves.tab", header=T, sep="\t", stringsAsFactors=F, na="")
  curves = subset(curves, !is.na(Species) & Volume==9 & !is.na(Species))
  curves = merge(curves, media_names[, c("ShortName", "MediaType")], by.x="Media", by.y="ShortName")
  
  curves.unrep = subset(curves, paste(Species, Media, File) %in% paste(res.combinations.unrep$Species, res.combinations.unrep$Media, res.combinations.unrep$File.1) | paste(Species, Media, File) %in% paste(res.combinations.unrep$Species, res.combinations.unrep$Media, res.combinations.unrep$File.3))
  ggplot(curves.unrep) +
    geom_line(aes(x=Time/3600, y=OD, group=paste(File, Well), color=factor(Passage))) +
    geom_point(aes(y=MaxODRaw.3, x=MaxTime.3/3600, color="3"), data=res.combinations.unrep, alpha=0.3) +
    geom_hline(aes(yintercept=BlankOD.3, color="3"), data=res.combinations.unrep, alpha=0.8) +
    geom_point(aes(y=MaxODRaw.1, x=MaxTime.1/3600, color="1"), data=res.combinations.unrep, alpha=0.3) +
    geom_hline(aes(yintercept=BlankOD.1, color="1"), data=res.combinations.unrep, alpha=0.8) +
    facet_wrap( ~ Species + Media) 
  
  ###################################
  # Unreproducible (1) curves
  ###################################
  res.combinations.unrep1 = subset(res.combinations, unreproducible.1=="unreproducible")
  curves.unrep1 = merge(subset(curves, Passage==3), res.combinations.unrep1, by.x=c("File", "MediaType"), by.y=c("File.3", "MediaType"), suffixes=c("", ".unreproducible"))
  curves.unrep1 = subset(curves.unrep1, !grepl("NoGrowth", Class.1))
  curves.unrep1 = ddply(curves.unrep1, .(File, TechnicalReplicates, Media, Species, Species.unreproducible, Media.unreproducible, Time), summarize, OD=mean(OD))
  curves.unrep1$Source = ""
  curves.unrep1$Source[with(curves.unrep1, Species != Species.unreproducible)] = "Contaminant species"
  curves.unrep1$Source[with(curves.unrep1, Species == Species.unreproducible)] = "Unreproducible species"
  curves.unrep1$Source[with(curves.unrep1, Species == Species.unreproducible & Media == Media.unreproducible)] = "Unreproducible experiment"
  curves.unrep1$Source[with(curves.unrep1, Species != Species.unreproducible & Media == Media.unreproducible)] = "Contaminant experiment"
  ggplot(curves.unrep1) +
    geom_line(aes(x=Time/3600, y=OD, group=paste(File, TechnicalReplicates), color=Source), alpha=0.7) +
    scale_color_manual(values=c("Contaminant species"="#42a5f5", "Unreproducible species"="#8bc34a", "Unreproducible experiment"="#ff3d00", "Contaminant experiment"="#d500f9")) +
    facet_wrap( ~ Species.unreproducible + Media.unreproducible) 
  
}
