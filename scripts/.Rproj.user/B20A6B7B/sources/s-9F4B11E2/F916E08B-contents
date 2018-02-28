library(ggplot2)
library(reshape2)
library(plyr)
library(gplots)
library(RColorBrewer)
source("utils/ggplot.R")
source("analyze.functions.R")
library(gridExtra)
library(igraph)

cluster = function()
{
  media.annotation = read.delim("../data/media_names.tab", sep="\t", na="", quote="", stringsAsFactors=F, header=T)
  media.annotation = media.annotation[order(media.annotation$Order),]
  media.names = sapply(media.annotation$ShortName, function(x) media.annotation$FullName[match(x, media.annotation$ShortName)])
  media.names2 = sapply(media.annotation$ShortName, function(x) media.annotation$FullName2[match(x, media.annotation$ShortName)])
  media.general = media.annotation$ShortName[media.annotation$IsGeneral==1]
  media.rich = media.annotation$ShortName[media.annotation$IsRich==1]
  
  curves = read.table("../data/curves.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="")
  curves.f = subset(curves, !is.na(Species) & is.na(ConditionSpecies) & !grepl("_|15|16", Media))
  #curves.f = subset(curves.f, File=="20170609_Monoculture_Revision_Plate1-12_Species31-53_plus_empty_plate11")
  curves.f = ddply(curves.f, .(File, Species, Time, Media), summarize, OD=mean(OD, na.rm=T))
  
  # Is growing
  curves.a = read.table("../data/curves_annotation.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="")
  curves.a = subset(curves.a, !is.na(Species) & is.na(ConditionSpecies) & !grepl("_|15|16", Media))
  curves.a = ddply(curves.a, .(File, Species, Media), summarize, IsGrowing=any(!grepl("NoGrowth|Unrep|Undef", Class)))
  curves.f = merge(curves.f, curves.a[,c("File", "Species", "Media", "IsGrowing")], by=c("File", "Species", "Media"))
  
  curves.f_comb = merge(curves.f, curves.f, by=c("File", "Media", "Time"), suffixes=c(".1", ".2"))
  curves.f_comb$Same = substr(curves.f_comb$Species.1, 1, 7) == substr(curves.f_comb$Species.2, 1, 7)
  order.cols = c("Species.1", "Species.2", "OD.1", "OD.2", "IsGrowing.1", "IsGrowing.2")
  curves.f_comb[,order.cols] = t(apply(curves.f_comb[,order.cols], 1, function(z) {
    z.n = 3
    z[rep(0:(z.n-1), each=2)*2 + rep(order(z[1:2]),z.n)]
  }))
  curves.f_comb$IsGrowing.1 = grepl("TRUE", curves.f_comb$IsGrowing.1)
  curves.f_comb$IsGrowing.2 = grepl("TRUE", curves.f_comb$IsGrowing.2)
  curves.f_comb$OD.1 = as.numeric(curves.f_comb$OD.1)
  curves.f_comb$OD.2 = as.numeric(curves.f_comb$OD.2)
  curves.f_comb = curves.f_comb[with(curves.f_comb, order(File, Species.1, Species.2, Time)),]
  curves.f_comb$Time = curves.f_comb$Time/3600
  curves.f_dist = ddply(curves.f_comb, .(File, Same, Species.1, Species.2), function(z) {
    #z = subset(curves.f_comb, File=="150122_Big_growth_curves_1_plate_1" & Species.1=="E. coli ED1a" & Species.2=="L. plantarum")
    z = subset(z, (IsGrowing.1 | IsGrowing.2) & !(Media%in%media.rich) & !is.na(Time))
    if(!nrow(z)) return(data.frame(distance=NA, n=0))
    
    #ggplot(z)+
    #  geom_line(aes(Time, OD.1, color=Species.1, group=paste(File, Species.1, Media))) +
    #  geom_line(aes(Time, OD.2, color=Species.2, group=paste(File, Species.2, Media))) +
    #  facet_wrap(~Media)
    
    #ggplot(z)+
    #  geom_line(aes(Time, OD.diff, color=paste(Species.1, Species.2), group=paste(File, Species.1, Species.2, Media))) +
    #  facet_wrap(~Media)
    
    
    z$OD.diff = abs(z$OD.1-z$OD.2)^2
    z.integral = integrate(approxfun(z$Time,z$OD.diff), 1, max(z$Time), subdivisions=1000L, stop.on.error=F)
    r = data.frame(distance=z.integral$value/max(z$Time), n=length(unique(z$Media)))
    r$distance.norm = r$distance/r$n
    r
  })
  
  curves.f_dist$SameGenus = ifelse(curves.f_dist$Same, "Same genus", "Different genus")
  ggplot(curves.f_dist) +
    geom_density(aes(x=log2(distance.norm), fill=SameGenus), alpha=0.7) +
    geom_vline(xintercept=-9, color="#FF0000") +
    scale_x_continuous(breaks=-15:0) +
    scale_fill_brewer("Greys") +
    theme_bw()

  
  curves.f_dist.f = subset(curves.f_dist, !is.na(distance) & n>1 & !Same & log2(distance.norm)< -9)
  curves.f_dist.f$File = as.factor(as.character(curves.f_dist.f$File))
  curves.f_dist.f$FileN = as.numeric(curves.f_dist.f$File)
  curves.f_dist.f$Var.1 = paste(curves.f_dist.f$FileN, "-", curves.f_dist.f$Species.1)
  curves.f_dist.f$Var.2 = paste(curves.f_dist.f$FileN, "-", curves.f_dist.f$Species.2)
  curves.f_dist.f = curves.f_dist.f[,c("Var.1", "Var.2", "Species.1", "Species.2", "File", "FileN", "distance", "n", "distance.norm")]
  #write.table(curves.f_dist.f, file="../data/suspected_contamination.tab", sep="\t", quote=F, row.names=F, na="")
  
  curves.f.contaminated = melt(curves.f_dist.f, id.vars="File", measure.vars=c("Species.1", "Species.2"), value.name="Species")
  curves.f.contaminated = merge(curves.f, curves.f.contaminated, by=c("File", "Species"))
  curves.f.contaminated$Media = factor(curves.f.contaminated$Media, sort(unique(as.numeric(curves.f.contaminated$Media))))
  pp = list()
  for(file in unique(curves.f_dist.f$File)) {
    pp[[length(pp)+1]] = 
      ggplot(subset(curves.f.contaminated, File %in% file))+
      geom_line(aes(Time/3600, OD, color=Species, group=paste(File, Species, Media))) +
      labs(x="", y="", title=file) + 
      scale_x_continuous(breaks=c(), labels=c()) +
      facet_wrap(~Media, nrow=1) +
      theme_classic() +
      theme(strip.text.x=element_text(angle=0, size=rel(0.5)))
  }
  pp[["ncol"]] = 1
  
  pdf("../report/suspected_contamination.pdf", paper="a4r", height=10, width=10)
  do.call(grid.arrange, pp)
  dev.off()
  
  g = graph_from_data_frame(curves.f_dist.f, directed=F, vertices=NULL)
  E(g)$width = (1 - curves.f_dist.f$distance.norm)*7
  l <- layout.fruchterman.reingold(g)*5
  plot(g, vertex.size=1, layout=l)

  ddply(curves.f_dist.f, .(FileN), summarize, File=File[1])

}

x= function()
{
  curves.a = read.table("../data/curves_annotation.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="")
  curves.a = subset(curves.a, !is.na(Species) & is.na(ConditionSpecies))
  curves.merge.final = cuves.merge_annotations2(subset(curves.a, Volume<10))
  curves.merge.review = cuves.merge_annotations2(subset(curves.a, Volume<9))
  
  
  pdf("../report/replicates_improvement_histogram.pdf", paper="a4r", height=5, width=5)
  ggplot() +
    geom_histogram(aes(x=Replicates, fill="Final"), bins=14, data=curves.merge.final, alpha=0.5)  +
    geom_histogram(aes(x=Replicates, fill="Before review"), bins=14, data=curves.merge.review, alpha=0.5) +
    geom_vline(aes(color="Final"), xintercept=median(curves.merge.final$Replicates)) +
    geom_vline(aes(color="Before review"), xintercept=median(curves.merge.review$Replicates)) +
    labs(title=paste0("More than 2: ", 
                      round(mean(curves.merge.review$Replicates>=3)*100, 0), "% / ", round(mean(curves.merge.final$Replicates>=3)*100, 0), "%")) +
    theme_bw()
  
  dev.off()
  
}

reproducibility_emptywells = function()
{
  curves = read.table("../data/curves.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="")
  curves.empty = ddply(subset(curves, Passage==1 & !is.na(Media)), .(File, Row, Col, Well, TechnicalReplicates, Species, Media, ConditionSpecies), summarize, MaxOD=max(OD, na.rm=T), BlankOD=min(OD, na.rm=T))
  curves.empty.f = subset(curves.empty, is.na(Species) & is.na(ConditionSpecies))
  with(curves.empty.f, table(MaxOD - BlankOD <= 0.15))
  curves.empty.ff = subset(curves.empty.f, MaxOD - BlankOD > 0.15)
  
  pdf("../report/contaminated_empty_wells.pdf", paper="a4r", height=8.27, width=11.69)
  curves.f = subset(curves, paste(File, Well) %in% paste(curves.empty.ff$File, curves.empty.ff$Well))
  ggplot(curves.f) +
    geom_line(aes(x=Time, y=OD, group=Well)) +
    facet_wrap(Well~File) +
    theme_bw(base_size=18)
  dev.off()
}

reproducibility_quantitative = function()
{
  media.annotation = read.delim("../data/media_names.tab", sep="\t", na="", quote="", stringsAsFactors=F, header=T)
  media.annotation = media.annotation[order(media.annotation$Order),]
  media.names = sapply(media.annotation$ShortName, function(x) media.annotation$FullName[match(x, media.annotation$ShortName)])
  media.files = sapply(media.annotation$ShortName, function(x) media.annotation$Filename[match(x, media.annotation$ShortName)])
  
  variables = c("MaxOD", "StatOD", "Rate") #, "OvergrowthOD")
  curves.a = read.table("../data/curves_annotation.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="")
  curves.a = subset(curves.a, !is.na(Species) & is.na(ConditionSpecies))
  curves.a = curves.a[!duplicated(curves.a[,c("Species", "Volume", "Media")]),]
  curves.a = merge(curves.a, cuves.merge_annotations2(curves.a)[,c("Species", "Media", "Growing", "Replicates")], by=c("Species", "Media"), all.x=T)
  curves.a = subset(curves.a, 
                    grepl("NoGrowth", Class) & (is.na(curves.a$Growing) | Growing/Replicates < 0.5) | 
                      !grepl("NoGrowth", Class) & !is.na(Growing) & !is.na(Replicates) & Growing/Replicates > 0.5)
  
  curves.a$MediaName = media.names[curves.a$Media]
  curves.a$MediaFile = media.files[curves.a$Media]
  curves.a$MaxOD = curves.a$MaxOD - curves.a$BlankOD
  curves.a$MaxOD[curves.a$MaxOD < 0] = 0
  curves.a$StatOD = curves.a$StatOD - curves.a$BlankOD
  curves.a$StatOD[curves.a$StatOD < 0] = 0
  curves.a$Rate[curves.a$Rate < 0] = 0
  curves.a.f = subset(curves.a, !(Media %in% c("15 A", "15 B", "16")))
  
  pdf("../report/rep_quantitative_summary2.pdf", paper="a4r", height=8.27, width=11.69)
  for(var in variables) {
    curves.a.rep = ddply(curves.a.f, .(Media, MediaName, ConditionSpecies, Species), compatible.replicates, corrected=F, column=var)
    curves.a.rep_cor = ggplot.cor_data(curves.a.rep, mapping=aes(x, y), method="pearson")
    curves.a.rep_cor$max = apply(curves.a.rep_cor, 1, function(x) max(as.numeric(x[c("right", "top_ci")])))
    print(ggplot(curves.a.rep, aes(x, y)) +
            geom_abline(intercept=0, slope=1,  size=0.5, alpha=0.8, color="#000000") +
            geom_point() +
            geom_smooth(method="lm") +
            geom_point(aes(max, max), data=curves.a.rep_cor, alpha=0) +
            geom_text(aes(x=max, y=0, label=short_str), curves.a.rep_cor, hjust=1, vjust=0, size=12, color="#D73027") +
            labs(title=paste0("Reproducibility(", var, ")")) +
            theme_classic(base_size=14))
  }
  dev.off()
  
}

