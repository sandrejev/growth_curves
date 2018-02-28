library(plyr)
library(reshape2)
library(ggplot2)
library(ape)
library(pheatmap)
source("analyze.functions.R")


phylumcol = c("Bacteroidetes"="#8DD3C7", "Firmicutes"="#FFFFB3", "Proteobacteria"="#BEBADA", "Actinobacteria"="#FB8072", "Verrucomicrobia"="#80B1D3", "Fusobacteria"="#B3DE69")

species.general = c("A. odontolyticus", "B. adolescentis", "B. animalis subsp. lactis BI-07", "B. animalis subsp. lactis BL-04", 
                    "B. caccae", "B. clarus", "B. dorei", "B. fragilis", "B. fragilis enterotoxigenic", "B. fragilis HM-20", "B. fragilis HM-709", 
                    "B. fragilis HM-710", "B. fragilis HM-711", "B. fragilis HM-712", "B. fragilis HM-713", "B. fragilis HM-714", "B. hansenii", 
                    "B. longum subsp. infantis", "B. obeum", "B. ovatus", "B. stercoris", "B. thetaiotaomicron", "B. uniformis", "B. uniformis HM-715", 
                    "B. uniformis HM-716", "B. vulgatus", "B. vulgatus HM-720", "B. xylanisolvens", "C. aerofaciens", "C. comes", "C. perfringens C36", 
                    "C. perfringens S107", "C. ramosum", "C. saccharolyticum", "D. formicigenerans", "D. piger", "E. coli CFT073", "E. coli E2348/69", 
                    "E. coli ED1a", "E. coli H10407", "E. coli HM605", "E. coli IAI1", "E. coli UTI89", "E. eligens", "E. lenta", "E. limosum", "E. rectale", 
                    "E. siraeum", "F. nucleatum subsp. animalis", "F. nucleatum subsp. nucleatum", "L. fermentum", "L. gasseri", "L. lactis", "L. paracasei", 
                    "L. plantarum", "L. salivarius", "L. vaginalis", "P. difficile", "P. distasonis", "P. merdae", "R. gnavus", "R. intestinalis", 
                    "S. flexneri", "S. parasanguinis", "S. salivarius", "S. sonnei", "S. typhimurium ATCC14028", "S. typhimurium LT2", "V. cholerae A1552", 
                    "V. cholerae N16961", "V. parvula", "Y. pseudotuberculosis")

media.annotation = read.delim("../data/media_names.tab", sep="\t", na="", quote="", stringsAsFactors=F, header=T)
media.annotation = media.annotation[order(media.annotation$Order),]
media.names = sapply(media.annotation$ShortName, function(x) media.annotation$FullName[match(x, media.annotation$ShortName)])
media.names2 = sapply(media.annotation$ShortName, function(x) media.annotation$FullName2[match(x, media.annotation$ShortName)])
media.general = media.annotation$ShortName[media.annotation$IsGeneral==1]
media.rich = media.annotation$ShortName[media.annotation$IsRich==1]
media.defined = as.character(c(2:5, 7, 10:11))
media.files = sapply(media.annotation$ShortName, function(x) media.annotation$Filename[match(x, media.annotation$ShortName)])


organisms.ann = read.delim("../data/organisms.tab", sep="\t", na="", quote="", stringsAsFactors=F, header=T)
rownames(organisms.ann) = organisms.ann$species


table.replicates_number = function() {
  curves.a = read.table("../data/curves_annotation.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="")
  curves.a = subset(curves.a, !is.na(Species) & is.na(ConditionSpecies) & Passage==1 & !grepl("_", Media) & !(grepl("Unrep|Undef", Class)))
  
  curves.a$TechnicalReplicatesCount = sapply(curves.a$TechnicalReplicates, function(x) {
    x.chars = unlist(strsplit(gsub("([0-9])", "", x), "-"))
    x.nums = unlist(strsplit(gsub("([A-Z])", "", x), "-"))
    s.size1 = strtoi(charToRaw(x.chars[2]),16L) - strtoi(charToRaw(x.chars[1]),16L) + 1
    s.size2 = as.numeric(x.nums[2]) - as.numeric(x.nums[1]) + 1
    pmax(s.size1, s.size2)
  })
  curves.tech = ddply(curves.a, .(Species, Media), summarize, TechnicalReplicatesCount=sum(TechnicalReplicatesCount, na.rm=T))
  
  growth_matrix.q1 = cuves.merge_annotations2(curves.a)
  growth_matrix.q1 = merge(growth_matrix.q1, curves.tech, by=c("Species", "Media"))
  growth_matrix.q1$ReplicatesText = with(growth_matrix.q1, paste0(Growing, "/", Replicates, " (", TechnicalReplicatesCount, ")"))
  growth_matrix.q.replicates = dcast(growth_matrix.q1, Species ~ Media, value.var="ReplicatesText")
  growth_matrix.q.replicates = growth_matrix.q.replicates[,c("Species", "GMM", "BHI++", "WCA", "mGAM", as.character(1:5), as.character(7:11), "13", "14", "15 A", "15 B", "16")]
  growth_matrix.q.replicates[is.na(growth_matrix.q.replicates)] = "0 (0)"
  colnames(growth_matrix.q.replicates)[6:20] = paste0("M", colnames(growth_matrix.q.replicates)[6:20])
  
  write.table(growth_matrix.q.replicates, file="../data/supplementary/replicates_number.tab", sep="\t", quote=F, row.names=F, na="")
}

table.growth_matrix = function()
{
  curves.a.all = read.table("../data/curves_annotation.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="")
  curves.a = subset(curves.a.all, !is.na(Species) & is.na(ConditionSpecies) & Passage==1 & !grepl("_", Media))
  
  growth_matrix.q1 = cuves.merge_annotations2(curves.a)
  growth_matrix.q1$Media = gsub("^([0-9])", "M\\1", growth_matrix.q1$Media)
  growth_matrix.q1$MaxOD_with_sd = paste0(round(growth_matrix.q1$MaxOD,3), ifelse(growth_matrix.q1$MaxOD>0, paste0(" (", round(growth_matrix.q1$MaxOD_sd,3), ")"), ""))
  
  growth_matrix.q.raw = dcast(growth_matrix.q1, Species ~ Media, value.var="MaxOD")
  growth_matrix.q.export = growth_matrix.q.raw[,c("Species", "GMM", "BHI++", "WCA", "mGAM", paste0("M", c(1:5,7:11, 13:14, "15 A", "15 B", "16")))]
  write.table(growth_matrix.q.export, file="../data/supplementary/growth_matrix.tab", sep="\t", quote=F, row.names=F, na="")
  
  growth_matrix.q.raw = dcast(growth_matrix.q1, Species ~ Media, value.var="MaxOD_with_sd")
  growth_matrix.q.export = growth_matrix.q.raw[,c("Species", "GMM", "BHI++", "WCA", "mGAM", paste0("M", c(1:5,7:11, 13:14, "15 A", "15 B", "16")))]
  write.table(growth_matrix.q.export, file="../data/supplementary/growth_matrix_maxod-sd.tab", sep="\t", quote=F, row.names=F, na="")
}

histogram.overview = function()
{
  growth_matrix.q.raw = read.table("../data/supplementary/growth_matrix.tab", sep="\t", quote="", header=T, na="", check.names=F)
  rownames(growth_matrix.q.raw) = growth_matrix.q.raw$Species
  growth_matrix.q.raw$Species = NULL
  growth_matrix.q.raw = t(growth_matrix.q.raw)
  
  media.hist = data.frame(media=rownames(growth_matrix.q.raw), species.n=rowSums(!is.na(growth_matrix.q.raw) &  growth_matrix.q.raw > 0))
  media.hist$media_type = "Other"
  media.hist$media_type[media.hist$media %in% c("GMM", "BHI++", "WCA", "mGAM")] = "Rich"
  media.hist$media_type[media.hist$media %in% c("M8", "M9")] = "Mucin"
  media.hist$media_type[media.hist$media %in% paste0("M", c(2:7, 10:11))] = "Defined"
  media.hist$media_type[media.hist$media %in% c("M1", "M13", "M14", "M15 A", "M15 B", "M16")] = "Minimal"
  media.hist$media_type = factor(media.hist$media_type, c("Minimal", "Defined", "Mucin", "Rich"))
  species.hist = data.frame(species=colnames(growth_matrix.q.raw), media.n=colSums(!is.na(growth_matrix.q.raw) &  growth_matrix.q.raw > 0))
  
  pdf("../report/growth_capacity_histogram.pdf", height=4.5, width=7.5)
  gridExtra::grid.arrange(  
    ggplot(media.hist) +
      geom_histogram(aes(x=species.n, fill=media_type), binwidth=5, position="stack", color="#000000") +
      theme_slim() +
      scale_fill_manual(values=RColorBrewer::brewer.pal(5, "Greys")[c(1,2,4,5)]) +
      coord_cartesian(xlim=c(-2.5, 100)) +
      scale_x_continuous(breaks=seq(0, 100, 25)) +
      theme(legend.position = c(0.1, 0.8)) +
      labs(x = "Number of species growing in media", y="")
    ,
    ggplot(species.hist) +
      geom_histogram(aes(x=media.n), binwidth=1, position="stack", color="#000000", fill="#D4D4D4") +
      theme_slim() +
      labs(x = "Number of media supporting species", y="")
    ,
    ncol=1)
  dev.off()
}

heatmap.obsulute = function()
{
  growth_matrix.q.raw = read.table("../data/supplementary/growth_matrix.tab", sep="\t", quote="", header=T, na="", check.names=F)
  rownames(growth_matrix.q.raw) = growth_matrix.q.raw$Species
  colnames(growth_matrix.q.raw) = gsub("^M", "", colnames(growth_matrix.q.raw))
  growth_matrix.q.raw$Species = NULL
  growth_matrix.q.raw = t(growth_matrix.q.raw)
  
  
  growth_matrix.q.species_growth = colMeans(growth_matrix.q.raw[media.general,], na.rm=T)
  growth_matrix.q.media_growth = rowMeans(growth_matrix.q.raw[,species.general], na.rm=T)
  growth_matrix.q.raw = growth_matrix.q.raw[,order(growth_matrix.q.species_growth)]
  growth_matrix.q.species_growth = colMeans(growth_matrix.q.raw[media.general,], na.rm=T)
  growth_matrix.q.media_growth = rowMeans(growth_matrix.q.raw[,species.general], na.rm=T)
  
  growth_matrix.q = growth_matrix.q.raw
  growth_matrix.q.rownames_raw = rownames(growth_matrix.q.raw)
  growth_matrix.q.rownames = media.names2[growth_matrix.q.rownames_raw]
  growth_matrix.q.pal = colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu")[2:9])
  growth_matrix.q.breaks50 = c(0, 0.01, seq(0.02, max(growth_matrix.q.raw, na.rm=T), length.out = 48))
  growth_matrix.q.colnames = colnames(growth_matrix.q)
  
  rownames(growth_matrix.q) = growth_matrix.q.rownames
  colnames(growth_matrix.q) = growth_matrix.q.colnames

  growth_matrix.q.count_raw = dcast(growth_matrix.q1, Media ~ Species, value.var="Count")
  growth_matrix.q.count = growth_matrix.q.count_raw
  growth_matrix.q.count[is.na(growth_matrix.q.count)] = ""
  rownames(growth_matrix.q.count) = growth_matrix.q.count$Media
  growth_matrix.q.count$Media = NULL
  growth_matrix.q.count = growth_matrix.q.count[rownames(growth_matrix.q.raw),colnames(growth_matrix.q.raw)]
  
  #
  # Annotations
  #
  growth_matrix.q.col_growth.color = growth_matrix.q.pal(100)[sapply(growth_matrix.q.species_growth, function(x) { which.min(abs(x - seq(0, max(growth_matrix.q, na.rm=T), length.out=100))) })]
  growth_matrix.q.row_growth.color = growth_matrix.q.pal(100)[sapply(growth_matrix.q.media_growth, function(x) { which.min(abs(x - seq(0, max(growth_matrix.q, na.rm=T), length.out=100))) })]
  
  growth_matrix.q.col.ann = cbind(data.frame(avg.OD=growth_matrix.q.species_growth, row.names=colnames(growth_matrix.q)), organisms.ann[growth_matrix.q.colnames, ])
  growth_matrix.q.row.ann = data.frame(
    avg.OD=growth_matrix.q.media_growth, 
    row.names=rownames(growth_matrix.q)[sapply(names(growth_matrix.q.media_growth), function(x) grep(paste0("^", gsub("[^A-Za-z0-9]", "", x), "\\b"), rownames(growth_matrix.q)))], 
    stringsAsFactors=F)
  
  anncol = c("avg.OD", "phylum")
  bincol = c(Yes="#1F78B4", No="#A6CEE3")
  growth_matrix.q.ann.colors = list(
    avg.OD=growth_matrix.q.pal(100), 
    ARRANGEMENT_SINGLES=bincol, 
    ARRANGEMENT_CLUSTER=bincol,
    phylum=phylumcol,
    pathogen=bincol, 
    motility=bincol, 
    media_utilisation=colorRampPalette(RColorBrewer::brewer.pal(9, "Greens"))(100),
    sort=colorRampPalette(RColorBrewer::brewer.pal(9, "Greens"))(100), 
    gut_prevalence=colorRampPalette(RColorBrewer::brewer.pal(9, "Greens"))(100)
  )[anncol]
  
  
  growth_matrix.q.order = growth_matrix.q[,order(organisms.ann$sort[match(growth_matrix.q.colnames, organisms.ann$species)])]
  growth_matrix.q.count.order = growth_matrix.q.count[,order(organisms.ann$sort[match(growth_matrix.q.colnames, organisms.ann$species)])]
  pheatmap(t(growth_matrix.q.order), cluster_rows=F,
           cluster_cols=F, gaps_row=rep(0,4), gaps_col=c(0, 0, rep(which(!duplicated(media.annotation$Type))[-1]-2, each=2)),
           color=c("#CCCCCC", growth_matrix.q.pal(48)), breaks=growth_matrix.q.breaks50, border_color="#FFFFFF",
           main="Growth heatmap (OD)", key.title="", key.ylab="", key.xlab="OD", 
           fontsize=18, cellwidth=18, cellheight=12, height=21, width=14, fontsize_row=11, fontsize_col=16,
           annotation_col=growth_matrix.q.row.ann,  
           annotation_row=growth_matrix.q.col.ann[,anncol], annotation_color=growth_matrix.q.ann.colors,
           filename="../report/growth_matrix_sorted.pdf")
  
  
  growth_matrix.q.row_dist = dist.species(growth_matrix.q)
  pheatmap(mat=t(growth_matrix.q), #display_numbers=t(growth_matrix.q.count), 
           clustering_distance_rows=growth_matrix.q.row_dist, clustering_method="average",
           cluster_cols=F, gaps_row=rep(0,4), gaps_col=c(0, 0, rep(which(!duplicated(media.annotation$Type))[-1]-2, each=2)),
           color=c("#CCCCCC", growth_matrix.q.pal(48)), breaks=growth_matrix.q.breaks50, border_color="#FFFFFF", 
           main="Growth heatmap (OD)", key.title="", key.ylab="", key.xlab="OD",
           fontsize=18, cellwidth=18, cellheight=12, height=21, width=14, fontsize_row=11, fontsize_col=16,
           fontsize_number=8,
           annotation_col=growth_matrix.q.row.ann,  
           annotation_row=growth_matrix.q.col.ann[,anncol], annotation_color=growth_matrix.q.ann.colors,
           filename="../report/growth_matrix_clustered.pdf")
  # 
  # species.defined_preference = c("A. odontolyticus","B. adolescentis","B. dorei","B. longum subsp. infantis","C. comes","D. formicigenerans","E. eligens","E. limosum","P. distasonis","S. parasanguinis","S. typhimurium LT2")
  # growth_matrix.q.pref = growth_matrix.q[,species.defined_preference]
  # growth_matrix.q.pref.row_dist = dist.species(growth_matrix.q.pref)
  # pheatmap(mat=t(growth_matrix.q.pref), #display_numbers=t(growth_matrix.q.count), 
  #          clustering_distance_rows=growth_matrix.q.pref.row_dist, clustering_method="average",
  #          cluster_cols=F, gaps_row=rep(0,4), gaps_col=c(0, 0, rep(which(!duplicated(media.annotation$Type))[-1]-2, each=2)),
  #          color=c("#CCCCCC", growth_matrix.q.pal(48)), breaks=growth_matrix.q.breaks50, border_color="#FFFFFF", 
  #          main="Growth heatmap (OD)", key.title="", key.ylab="", key.xlab="OD",
  #          fontsize=18, cellwidth=18, cellheight=12, height=21, width=14, fontsize_row=11, fontsize_col=16,
  #          fontsize_number=8,
  #          annotation_col=growth_matrix.q.row.ann,  
  #          annotation_row=growth_matrix.q.col.ann[,anncol], annotation_color=growth_matrix.q.ann.colors,
  #         filename="../report/growth_matrix_defined_preference.pdf")
  # 
  # species.bacteroides = organisms.ann$species[organisms.ann$genus=="Bacteroides"]
  # growth_matrix.q.bact = growth_matrix.q[,species.bacteroides]
  # growth_matrix.q.bact.row_dist = dist.species(growth_matrix.q.bact)
  # pheatmap(mat=t(growth_matrix.q.bact), #display_numbers=t(growth_matrix.q.count), 
  #          clustering_distance_rows=growth_matrix.q.bact.row_dist, clustering_method="average",
  #          cluster_cols=F, gaps_row=rep(0,4), gaps_col=c(0, 0, rep(which(!duplicated(media.annotation$Type))[-1]-2, each=2)),
  #          color=c("#CCCCCC", growth_matrix.q.pal(48)), breaks=growth_matrix.q.breaks50, border_color="#FFFFFF", 
  #          main="Growth heatmap (OD)", key.title="", key.ylab="", key.xlab="OD",
  #          fontsize=18, cellwidth=18, cellheight=12, height=21, width=14, fontsize_row=11, fontsize_col=16,
  #          fontsize_number=8,
  #          annotation_col=growth_matrix.q.row.ann,  
  #          annotation_row=growth_matrix.q.col.ann[,anncol], annotation_color=growth_matrix.q.ann.colors,
  #          filename="../report/growth_matrix_bacteroides.pdf")
  # 
  # species.mucin = c("A. muciniphila", "R. torques")
  # growth_matrix.q.mucin = growth_matrix.q[,species.mucin]
  # growth_matrix.q.mucin.row_dist = dist.species(growth_matrix.q.mucin)
  # pheatmap(mat=t(growth_matrix.q.mucin), #display_numbers=t(growth_matrix.q.count), 
  #          clustering_distance_rows=growth_matrix.q.mucin.row_dist, clustering_method="average",
  #          cluster_cols=F, gaps_row=rep(0,4), gaps_col=c(0, 0, rep(which(!duplicated(media.annotation$Type))[-1]-2, each=2)),
  #          color=c("#CCCCCC", growth_matrix.q.pal(48)), breaks=growth_matrix.q.breaks50, border_color="#FFFFFF", 
  #          main="Growth heatmap (OD)", key.title="", key.ylab="", key.xlab="OD",
  #          fontsize=18, cellwidth=18, cellheight=12, height=21, width=14, fontsize_row=11, fontsize_col=16,
  #          fontsize_number=8,
  #          annotation_col=growth_matrix.q.row.ann,  
  #          annotation_row=growth_matrix.q.col.ann[,anncol], annotation_color=growth_matrix.q.ann.colors,
  #          filename="../report/growth_matrix_mucin.pdf")
}

heatmap.curves = function()
{
  curves = read.table("../data/curves.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="")
  curves.f.medias = c("1", "13", "14", "15 A", "15 B", "16", "2", "3", "4", "5", "7", "10", "11", "8", "9", "GMM", "BHI++", "WCA", "mGAM")
  curves.f = subset(curves, !is.na(Species) & is.na(ConditionSpecies) & !is.na(Media) & Passage==1 & Media %in% curves.f.medias)
  curves.f$MediaType = "Other"
  curves.f$MediaType[curves.f$Media %in% media.rich] = "Rich"
  curves.f$MediaType[curves.f$Media %in% c("8", "9")] = "Mucin"
  curves.f$MediaType[curves.f$Media %in% as.character(c(2:7, 10:11))] = "Defined"
  curves.f$MediaType[curves.f$Media %in% c("1", "13", "14", "15 A", "15 B", "16")] = "Minimal"
  curves.f$MediaType = factor(curves.f$MediaType, c("Minimal", "Defined", "Mucin", "Rich"))
  curves.f$Media = factor(curves.f$Media, )
  curves.f = ddply(curves.f, .(File), mutate, AllSpecies=paste(unique(Species), collapse=","))
  curves.f = curves.f[with(curves.f, order(Species, Media, Time)),]
  curves.f$File2 = paste0(curves.f$File, "\n", substring(curves.f$AllSpecies, 1, 50))
  curves.f$Species = gsub("\n$", "", gsub("^([^ ]+ [^ ]+ [^ ]+ [^ ]+)", "\\1\n", gsub("^([^ ]+ [^ ]+)", "\\1\n", curves.f$Species)))
  
  pdf("../report/curves_matrix.pdf", height=60, width=15)
  ggplot(curves.f) +
    geom_rect(aes(fill=MediaType), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.1, data=unique(curves.f[,c("Species", "Media", "MediaType")])) +
    geom_line(aes(x=Time/3600, y=OD, group=paste(File, TechnicalReplicates, Well), color=factor(Volume)), size=0.5) +
    facet_grid(Species ~ Media) +
    theme_classic(base_size=7) + 
    theme(strip.text.y = element_text(angle = 0))
  dev.off()
  
  pdf("../report/curves_matrix_file.pdf", height=60, width=15)
  contaminants = read.table("../data/suspected_contamination.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="")
  contaminants = merge(contaminants, unique(curves.f[,c("File", "File2", "Media")]))
  ggplot(curves.f[!duplicated(curves.f$File, curves.f$Media),]) +
    geom_rect(aes(fill=MediaType), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.1, data=unique(curves.f[,c("File2", "Media", "MediaType")])) +
    geom_rect(aes(fill="Contamination"), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0.2, alpha=0.1, data=contaminants) +
    geom_line(aes(x=Time/3600, y=OD, group=paste(File, TechnicalReplicates, Well), color=Species), size=0.1) +
    facet_grid(File2 ~ Media) +
    theme_classic(base_size=7) + 
    theme(strip.text.y = element_text(angle = 0, family="Times", size=rel(1.4))) +
    guides(color=element_blank())
  dev.off()
  
}

heatmap.relative = function()
{
  set.seed(0)
  
  curves.a = read.table("../data/curves_annotation.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="")
  curves.a = subset(curves.a, !is.na(Species) & is.na(ConditionSpecies) & !is.na(Media) & Passage==1 & !grepl("_", Media))
  #curves.a = subset(curves.a, !is.na(Species) & is.na(ConditionSpecies) & Passage==1 & !(grepl("Unrep|Undef", Class)))
  
  curves.a.wide = ddply(curves.a, .(File, Species), function(z) {
    z.3 = subset(z, Media==3)
    if(!nrow(z.3)) return(NULL)
    z.other = subset(z, Media != 3)
    
    maxod.other = z.other$MaxOD-z.other$BlankOD
    if(any(maxod.other<0)) maxod.other[maxod.other<0] = jitter(mean(maxod.other[maxod.other>0], na.rm=T))
    maxod.3 = z.3$MaxOD-z.3$BlankOD
    if(any(maxod.3<0)) maxod.3[maxod.3<0] = jitter(mean(maxod.other[maxod.other>0], na.rm=T))
    maxod.fold = log2(maxod.other/maxod.3)

    ret = data.frame(Media=z.other$Media, MaxOD.fold=maxod.fold, MaxOD.3=maxod.3, MaxOD.other=maxod.other, Class.other=z.other$Class, Class.3=z.3$Class)
    ret
  })


  curves.a.fold = ddply(curves.a.wide, .(Species, Media), function(z) {
    z.f = !grepl("NoGrowth", z$Class.other) | !grepl("NoGrowth", z$Class.3)
    if(mean(z.f) > 0.5) {
      z.ff = z.f & is.finite(z$MaxOD.fold) & !grepl("Unfin", z$Class.other)
      maxod.pval = c()
      maxod.fold = c()
      if(sum(z.ff) > 1) {
        maxod.pval = t.test(z$MaxOD.other[z.ff], z$MaxOD.3[z.ff], paired=T)$p.value
        maxod.fold = mean(z$MaxOD.fold[z.ff])
      }
      
      z.ff = z.f & is.finite(z$MaxOD.fold)
      if(sum(z.ff) > 1) {
        maxod.pval = t.test(z$MaxOD.other[z.ff], z$MaxOD.3[z.ff], paired=T)$p.value
        maxod.fold = mean(z$MaxOD.fold[z.ff])
      }
      
      if(length(maxod.pval) == 0) {
        return(data.frame(MaxOD.fold=mean(z$MaxOD.fold[z.f]), MaxOD.pval=NA))
      }
      
      return(data.frame(MaxOD.fold=maxod.fold[which.min(maxod.pval)], MaxOD.pval=min(maxod.pval)))
    }
    z.f = grepl("NoGrowth", z$Class.other) & grepl("NoGrowth", z$Class.3)
    if(mean(z.f) > 0.5) return(data.frame(MaxOD.fold=0, MaxOD.pval=1))
    
    return(NULL)
  })

  curves.a.fold$MaxOD.pval.adj = p.adjust(curves.a.fold$MaxOD.pval)
  curves.a.rep = ddply(curves.a, .(Species, Media), summarize, replicates=paste0(sum(!grepl("NoGrowth", Class)), "/",length(File)))
  curves.a.fold = merge(curves.a.fold, curves.a.rep, by=c("Species", "Media"))
  write.table(curves.a.fold, file="../data/generated/curves_annotation_rel.tab", sep="\t", quote=F, row.names=F, na="")
  
  curves.a.fold$MaxOD.text = with(curves.a.fold, ifelse(MaxOD.fold>0, "+", ifelse(MaxOD.fold<0, "-", "")))
  curves.a.fold$MaxOD.pval.text = with(curves.a.fold, ifelse(!is.na(MaxOD.pval) & MaxOD.pval <= 0.05, "*", ""))
  curves.a.fold$text = with(curves.a.fold, paste0(MaxOD.text, replicates, MaxOD.pval.text))
    
  growth_matrix = dcast(curves.a.fold, Species ~ Media, value.var="MaxOD.fold")
  growth_matrix = merge(growth_matrix, organisms.ann[,c("species", "sort")], by.x="Species", by.y="species")
  growth_matrix = growth_matrix[order(growth_matrix$sort),]
  rownames(growth_matrix) = growth_matrix$Species
  growth_matrix$Species = NULL
  growth_matrix$sort = NULL
  growth_matrix = growth_matrix[,c("1", "13", "14", "15 A", "15 B", "16", "2", "4", "5", "7", "10", "11", "8", "9", "GMM", "BHI++", "WCA", "mGAM")]
  growth_matrix = data.matrix(growth_matrix)
  growth_matrix.vals = growth_matrix[is.numeric(growth_matrix) & is.finite(growth_matrix)]

  growth_matrix.text = dcast(curves.a.fold, Species ~ Media, value.var="text")
  rownames(growth_matrix.text) = growth_matrix.text$Species
  growth_matrix.text$Species = NULL
  growth_matrix.text = as.matrix(growth_matrix.text)
  
  growth_matrix.ann = subset(cuves.merge_annotations2(curves.a, nogrowth=T), Media==9, c(Species, MaxOD))
  colnames(growth_matrix.ann) = c("species", "M9")
  growth_matrix.ann = merge(growth_matrix.ann, organisms.ann, by="species")[,c("species", "M9", "phylum")]
  rownames(growth_matrix.ann) = growth_matrix.ann$species
  growth_matrix.ann$species = NULL
  
  growth_matrix.pal_max = 5
  growth_matrix.pal_breaks = seq(-max(abs(growth_matrix.vals)), max(abs(growth_matrix.vals)), length.out=109)
  growth_matrix.pal = colorRampPalette(c(rev(RColorBrewer::brewer.pal(8, "PuBu")[1:7]), rep("white", 1), RColorBrewer::brewer.pal(8, "Reds")[1:7]))(sum(growth_matrix.pal_breaks>=-growth_matrix.pal_max & growth_matrix.pal_breaks<=growth_matrix.pal_max))
  growth_matrix.pal = c(rep(growth_matrix.pal[1], sum(growth_matrix.pal_breaks<(-growth_matrix.pal_max))), growth_matrix.pal)
  growth_matrix.pal = c(growth_matrix.pal, rep(growth_matrix.pal[length(growth_matrix.pal)], sum(growth_matrix.pal_breaks>growth_matrix.pal_max)))
  growth_matrix.ann.colors = list(
    M9=colorRampPalette(RColorBrewer::brewer.pal(9, "Greys")[2:9])(100),
    phylum = c("Bacteroidetes"="#8DD3C7", "Firmicutes"="#FFFFB3", "Proteobacteria"="#BEBADA", "Actinobacteria"="#FB8072", "Verrucomicrobia"="#80B1D3", "Fusobacteria"="#B3DE69")
  )

  fig4a.species = c("B. obeum", "B. uniformis HM−715", "B. uniformis HM−716", "B. fragilis HM−20", "B. ovatus", "B. eggerthii", 
                    "B. clarus", "F. nucleatum subsp. nucleatum", "B. thetaiotaomicron", "P. merdae", "B. xylanisolvens",
                    "R. intestinalis", "C. saccharolyticum", "L. sakei subsp. sakei", "L. paracasei", "L. salivarius", "L. vaginalis", 
                    "L. ruminis", "C. perfringens C36")
  fig4a.media = c("5", "10", "11")
  fig4a.matrix = growth_matrix[fig4a.species,fig4a.media]
  fig4a.matrix_text = growth_matrix.text[fig4a.species,fig4a.media]
  growth_matrix.ann.colors.f.subset = 1:(max(growth_matrix.ann[fig4a.species, "M9"])/max(growth_matrix.ann$M9)*100)
  growth_matrix.ann.colors.f = growth_matrix.ann.colors
  growth_matrix.ann.colors.f$M9 = growth_matrix.ann.colors.f$M9[growth_matrix.ann.colors.f.subset]
  pheatmap(fig4a.matrix, main="Relative effect on growth (OD)", 
           cluster_rows=F,cluster_cols=F,
           annotation_row=growth_matrix.ann, 
           annotation_colors=growth_matrix.ann.colors.f,
           color=growth_matrix.pal, breaks=growth_matrix.pal_breaks, border_color="#CCCCCC",
           fontsize=18, cellwidth=18, cellheight=12, height=21, width=14, fontsize_row=11, fontsize_col=16,
           #display_numbers=fig4a.matrix_text, fontsize_number=8,
           filename="../report/fig4a_media_influence_scfa.pdf")

  figS4b.media = c("8")
  fig4b.species = list(
    top_left = c("E. coli CFT073", "E. coli E2348/69", "E. coli ED1a", "E. coli H10407", "E. coli HM605", "E. coli IAI1", "E. coli UTI89", "C. comes", "L. paracasei", "P. difficile"),
    bottom_left = c("B. stercoris", "R. hominis", "B. longum subsp. longum", "C. aerofaciens", "C. bolteae"),
    top_right = c("A. muciniphila", "R. torques", "C. perfringens S107", "B. hansenii"),
    bottom_right = c("B. uniformis HM-716", "P. capillosus", "B. ovatus", "B. clarus", "B. thetaiotaomicron", "B. xylanisolvens", "P. merdae", "R. intestinalis", "B. caccae", "B. fragilis HM-20", "B. fragilis HM-711", "B. fragilis enterotoxigenic")
  )
  for(pos in names(fig4b.species))
  {
    species = fig4b.species[[pos]]
    figS4b.matrix = growth_matrix[species,figS4b.media, drop=F]
    figS4b.matrix_text = growth_matrix.text[figS4b.species,figS4b.media, drop=F]
    
    growth_matrix.ann.colors.f.subset = 1:(max(growth_matrix.ann[species, "M9"])/max(growth_matrix.ann$M9)*100)
    growth_matrix.ann.colors.f = growth_matrix.ann.colors
    growth_matrix.ann.colors.f$M9 = growth_matrix.ann.colors.f$M9[growth_matrix.ann.colors.f.subset]
    
    pheatmap(figS4b.matrix, main="Relative effect on growth (OD)", 
             cluster_rows=F,cluster_cols=F,
             annotation_row=growth_matrix.ann, 
             annotation_colors=growth_matrix.ann.colors.f,
             color=growth_matrix.pal, breaks=growth_matrix.pal_breaks, border_color="#CCCCCC",
             fontsize=18, cellwidth=18, cellheight=12, height=21, width=14, fontsize_row=11, fontsize_col=16,
             #display_numbers=figS4a.matrix_text, fontsize_number=8,
             filename=paste0("../report/fig4b_", pos, ".pdf"))
  }  
}

itol.phylogenetic_tree = function()
{
  organisms.ann = read.delim("../data/organisms.tab", sep="\t", na="", quote="", stringsAsFactors=F, header=T)
  
  growth_matrix = read.table("../data/supplementary/growth_matrix.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="", check.names=F)
  colnames(growth_matrix) = gsub("^M", "", colnames(growth_matrix))
  rownames(growth_matrix) = gsub("[^A-Za-z0-9]+", "_", growth_matrix$Species)
  growth_matrix$species = NULL
  
  growth_matrix.dist = dist.species(t(growth_matrix))
  growth_matrix.hclust = hclust(growth_matrix.dist, method="average")
  growth_matrix_phylo = as.phylo(growth_matrix.hclust)
  dir.create("../report/itol", showWarnings = FALSE)
  write.tree(growth_matrix_phylo, "../report/itol/growth_tree.newick")

  organisms.ann$phylum.color = phylumcol[organisms.ann$phylum]
  organisms.ann$phylum.color_type = "range"
  rownames(organisms.ann) = gsub("[^A-Za-z0-9]+", "_", organisms.ann$species)
  
  writeLines(paste0(
    "TREE_COLORS
    SEPARATOR TAB
    DATASET_LABEL\tSPECIES_ANNOTATION
    DATA"), "../report/itol/phylogenetic_color.txt")
  write.table(subset(organisms.ann, , c(phylum.color_type, phylum.color, phylum)), file="../report/itol/phylogenetic_color.txt", append=T, sep="\t", quote=F, col.names=F, row.names=T, na="")
}

plots.media_preferrence = function()
{
  media.defined = as.character(c(2:5, 7, 10:11))
  abundant_species = read.table("../data/abundant_species.tab", sep="\t", quote="", header=F, stringsAsFactors=F, na.strings="")$V1
  organisms.ann = read.delim("../data/organisms.tab", sep="\t", na="", quote="", stringsAsFactors=F, header=T)

  # Load abundance average data (only species passing 1% rel. abundance and 10% prevalence threshold)
  # Aggregate abundance data by species level (merge strains)
  organisms.gut = read.table("../data/screenG_tax_info_specI_clusters.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="")
  organisms.gut = merge(organisms.gut, organisms.ann[!is.na(organisms.ann$nid), c("species", "lineage_species", "nid")], by.x="NT_code", by.y="nid")
  organisms.gut = ddply(organisms.gut, .(lineage_species, SpecI_ID), summarize, Prevalence=mean(Prevalence))
  organisms.gut = ddply(organisms.gut, .(lineage_species), summarize, Prevalence=mean(Prevalence))
  rownames(organisms.gut) = organisms.gut$lineage_species
#  dim(organisms.gut)

  # Load species monoculture growth data
  # Aggregate growth data by species level (merge strains)
  growth_matrix = read.table("../data/supplementary/growth_matrix.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="", check.names=F)
  colnames(growth_matrix)[1] = "species"
  growth_matrix = merge(growth_matrix, organisms.ann[, c("species", "lineage_species")], by="species")
  growth_matrix = ddply(growth_matrix, .(lineage_species), function(x) {
    x = x[,setdiff(colnames(x), c("species", "lineage_species"))]
    x.mean = matrix(rep(colMeans(x>0), each=nrow(x)), ncol=ncol(x))
    x[x==0 & x.mean>0.5] = NA
    x = apply(x, 2, median, na.rm=T)
    x[is.nan(x)] = NA
    x
  })
  colnames(growth_matrix) = gsub("^M", "", colnames(growth_matrix))
  rownames(growth_matrix) = growth_matrix$lineage_species
  growth_matrix$lineage_species = NULL
  growth_matrix = growth_matrix[,c(media.defined, media.rich)]
  
  # Create ranked and scaled growth matrix
  growth_matrix.rank = apply(growth_matrix, 2, function(x) { 
    x.res = length(x) - rank(x, ties.method="min") + 1
    x.res
  })
  
  growth_matrix.rank_table = data.frame(lineage_species=rownames(growth_matrix.rank))
  rownames(growth_matrix.rank_table) = growth_matrix.rank_table$lineage_species
  growth_matrix.rank_table$rich_medrank = apply(growth_matrix.rank[, media.rich], 1, median, na.rm=T)
  growth_matrix.rank_table$defined_medrank = apply(growth_matrix.rank[, media.defined], 1, function(z) mean(sort(z)[3:5], na.rm=T))
  growth_matrix.rank_table$preference = with(growth_matrix.rank_table, log2(rich_medrank / defined_medrank))
  growth_matrix.rank_table$preferred_media = "no preferrence"
#  growth_matrix.rank_table$preferred_media[growth_matrix.rank_table$preference <= -1] = "rich"
  growth_matrix.rank_table$preferred_media[growth_matrix.rank_table$preference >= 1] = "defined"
  growth_matrix.rank_table$preferred_media = factor(growth_matrix.rank_table$preferred_media, c("rich", "no preferrence", "defined"))
  growth_matrix.rank_table$Prevalence = organisms.gut[rownames(growth_matrix.rank_table), "Prevalence"]
  growth_matrix.rank_table$AverageDefinedMaxOD = rowMeans(growth_matrix[, as.character(intersect(colnames(growth_matrix), c(1:8, 10:11)))])[rownames(growth_matrix.rank_table)]
  
  growth_matrix.rank_table.export = subset(growth_matrix.rank_table, preferred_media=="defined", c(lineage_species, Prevalence))
  colnames(growth_matrix.rank_table.export)[1] = "Species"
  write.table(growth_matrix.rank_table.export, file="../data/supplementary/defined_media_preferrence.tab", sep="\t", quote=F, row.names=F, na="")

  dens = density(growth_matrix.rank_table$preference, na.rm=T)
  df = data.frame(x=dens$x, y=dens$y)
  df$preferred_media = "no preferrence"
  #df$preferred_media[df$x <= -1] = "rich"
  df$preferred_media[df$x >= 1] = "defined"
 
  pdf("../report/media_preferrence.pdf", height=4.5, width=7.5)
  defined_growth_pval = with(growth_matrix.rank_table, t.test(AverageDefinedMaxOD[preferred_media=='defined'], AverageDefinedMaxOD[preferred_media!='defined']))$p.value
  ggplot(growth_matrix.rank_table) +
    geom_density(aes(fill=preferred_media, x=AverageDefinedMaxOD), alpha=0.3) + 
    theme_slim() +
    labs(fill="Preferred media", x="Average MaxOD on defined media") +
    scale_fill_manual(values=c("rich"="#0000004d", "no preferrence"="#bdbdbd4d", "defined"="#0000004d")) +
    theme(legend.position = c(0.2, 0.8))
  
  ggplot(df, aes(x, y)) +
    geom_line() +
    geom_ribbon(aes(ymin=0, ymax=y, fill=preferred_media)) +
    geom_text(aes(x=-2.5, y=0.4, label="Rich media preferred"), size=5) +
    geom_text(aes(x=2.5, y=0.4, label="Defined media preferred"), size=5) +
    labs(title="Media preferrence", x="rich(%)/defined(%)", y="") + 
    scale_fill_manual(values=c("rich"="#0000004d", "no preferrence"="#FFFFFF00", "defined"="#d4d4d4ff")) +
    theme_slim()
  
  set.seed(0)
  growth_matrix.rank_table.f = subset(growth_matrix.rank_table, !is.na(Prevalence))
  growth_matrix.rank_table.metahitsum = ddply(growth_matrix.rank_table.f, .(preferred_media), summarize, N=length(preferred_media), MeanPrevalence=mean(Prevalence))
  growth_matrix.rank_table.ttest = with(growth_matrix.rank_table.f, t.test(Prevalence[preferred_media=="defined"], Prevalence[preferred_media!="defined"]))
  ggplot(growth_matrix.rank_table.f) +
    geom_boxplot(aes(preferred_media, Prevalence), fill="#d4d4d4ff") +
    geom_jitter(aes(preferred_media, Prevalence), fill="#d4d4d4ff") +
    geom_text(aes(preferred_media, MeanPrevalence, label=N), data=growth_matrix.rank_table.metahitsum) +
    theme_slim() +
    labs(title="Species prevalence in MetaHIT", y="Prevalence in MetaHIT, %", x="") +
    coord_cartesian(ylim=c(0, 1))
  dev.off()
}
  
  