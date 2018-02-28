library(ggplot2)
library(plyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(beeswarm)
library(VennDiagram)
source("analyze.functions.R")


plot.mucin_validation_data = function() {
  medias = "^(9|8|8_salic075|9_salic075)$"
  species = rev(c("A. muciniphila", "B. thetaiotaomicron", "B. stercoris", "B. clarus", "C. comes", "R. torques"))
  
  growth_matrix = read.table("../data/supplementary/growth_matrix.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="", check.names=F)
  growth_matrix = melt(growth_matrix, variable.name="Media", value.name="IsGrowing")
  growth_matrix$IsGrowing = growth_matrix$IsGrowing>0
  growth_matrix$Media = gsub("^M", "", growth_matrix$Media)
  
  curves.a = read.table("../data/curves_annotation.tab", sep="\t", quote="", header=T, na.strings="")
  curves.a$MaxOD = curves.a$MaxOD-curves.a$BlankOD
  curves.a = merge(curves.a, growth_matrix, by=c("Species", "Media"), all.x=T) 
  curves.a = subset(curves.a, grepl(medias, Media) & Species %in% species)
  curves.a = subset(curves.a, grepl("salic", Media) | IsGrowing == !grepl("NoGrowth", Class) & !grepl("Unrep|Undef", Class))
  curves.a.f.raw = curves.a

  curves.a.f.raw$MediaOrigin = factor(gsub("_.*", "", curves.a.f.raw$Media))
  curves.a.f.raw$MediaOriginDirection = as.numeric(curves.a.f.raw$MediaOrigin)*2 - 3
  curves.a.f.raw$Mucin = ifelse(grepl("salic", curves.a.f.raw$Media), "Salic Acid", "Unpurified Mucin")
  
  curves.a.f.raw = ddply(curves.a.f.raw, .(MediaOrigin, Species), mutate, 
                     MaxTotal=max(MaxOD, na.rm=T), 
                     Pvalue=t.test(MaxOD[!grepl("salic", Media)], mu=MaxOD[grepl("salic", Media)])$p.value, 
                     PvalueText=ifelse(Pvalue<0.05, "*", ""), 
                     n1=sum(grepl("salic", Media)), 
                     n2=sum(!grepl("salic", Media)),
                     PvalueText=ifelse(Pvalue>0.01, as.character(round(Pvalue,2)), paste0("<10", ceiling(log10(Pvalue))))
                     )
  curves.a.f.raw$Species = as.character(curves.a.f.raw$Species)
  curves.a.f.raw = curves.a.f.raw[order(match(curves.a.f.raw$Species, species)),]
  curves.a.f.raw$Species = factor(curves.a.f.raw$Species, unique(as.character(curves.a.f.raw$Species)))
  curves.a.f = ddply(curves.a.f.raw, .(Species, Media, MediaOrigin, MediaOriginDirection, Mucin, Pvalue, PvalueText, MaxTotal), summarize, MaxOD.se=sd(MaxOD, na.rm=T)/sqrt(length(MaxOD)), MaxOD=mean(MaxOD, na.rm=T), n1=mean(n1), n2=mean(n2))
    
  unique(subset(curves.a.f, , c(Species, MediaOrigin, Pvalue, PvalueText)))
  
  pdf("../report/salic_acid_with_points.pdf", paper="a4r", height=8.27, width=11.69)
  write.table(unique(subset(curves.a.f, , c(Species, MediaOrigin, Pvalue))), file="../data/generated/sialic_acid.tab", sep="\t", quote=F, row.names=F, na="")
  ggplot(curves.a.f) +
    geom_bar(aes(Species, MaxOD*MediaOriginDirection, fill=Mucin), stat="identity", position="dodge") +
    geom_point(aes(Species, MaxOD*MediaOriginDirection, fill=Mucin, color=factor(MediaOriginDirection)), stat="identity", position=position_jitterdodge(dodge.width=1, jitter.width=0.4), data=curves.a.f.raw, size=0.7) +
    geom_text(aes(Species, y=MediaOriginDirection*1.5, label=paste0("p=",PvalueText)), stat="identity", position=position_dodge(width=1)) +
    geom_hline(yintercept=0, size=2, color="#FFFFFF") +
    geom_hline(yintercept=0, size=0.5, color="#666666") +
    geom_hline(yintercept=c(-.15, .15), size=0.5, color="#FF0000", lty=2) +
    geom_errorbar(aes(Species, ymin=MediaOriginDirection*(MaxOD), ymax=MediaOriginDirection*(MaxOD+MaxOD.se), fill=Mucin), stat="identity", position=position_dodge(width=1), width=0.3) +
    labs(x="", y="MaxOD") +
    scale_fill_manual(values=c("Salic Acid"="#666666", "Unpurified Mucin"="#CCCCCC")) +
    coord_flip(ylim=c(-1.5, 1.5)) +
    scale_y_continuous(breaks=seq(-1.5, 1.5, 0.5)) +
    theme_bw(base_size=16) 
  dev.off()
}


plot.mucin_figure = function() {
  mucin.known = c("A. muciniphila", "R. torques", "B. thetaiotaomicron", "B. fragilis", "B. vulgatus", "B. caccae", "B. uniformis", "R. gnavus", "B. longum subsp. infantis", "B. longum subsp. longum")
  
  
  all_species_star  = c("A. muciniphila*", "R. torques*", "C. perfringens S107", "C. ramosum", "C. comes", "A. odontolyticus", 
                        "C. aerofaciens", "B. clarus", "F. nucleatum subsp. nucleatum", "L. paracasei", "B. thetaiotaomicron*", 
                        "B. fragilis*", "B. vulgatus*", "B. caccae* different strain", "B. eggerthii", "B. uniformis*", "B. xylanisolvens", "P. capillosus", 
                        "R. gnavus*", "B. longum subsp. infantis*", "B. longum subsp. longum* different strain", 
                        
                        "C. perfringens C36", 
                        "B. fragilis HM-20", "B. fragilis HM-709",  "B. fragilis HM-710", "B. fragilis HM-711", "B. fragilis HM-712", "B. fragilis HM-713", "B. fragilis HM-714", "B. fragilis enterotoxigenic", 
                        "B. vulgatus HM-720", 
                        "B. uniformis HM-716", "B. uniformis HM-715", 
                        "B. ovatus", "B. dorei", "B. stercoris", "B. coprocola",
                        "P. merdae", "P. distasonis",
                        "B. hansenii", "B. obeum",                  
                        "F. nucleatum subsp. animalis", "F. nucleatum subsp. vincentii", 
                        "E. coli CFT073", "E. coli UTI89", "E. coli ED1a", "E. coli E2348/69",  "E. coli IAI1", "E. coli H10407", "E. coli HM605",  "S. flexneri", "S. typhimurium LT2", "S. sonnei", "Y. pseudotuberculosis", 
                        "C. bolteae", "C. leptum", "C. saccharolyticum", 
                        "P. difficile", 
                        "E. eligens", "E. limosum", "E. rectale", "E. siraeum", 
                        "S. parasanguinis", "S. salivarius",                   
                        "L. delbrueckii subsp. delbrueckii", "L. ruminis", "L. acidophilus", "L. fermentum", "L. gasseri", "L. plantarum", "L. sakei subsp. sakei", "L. salivarius", "L. vaginalis", "L. lactis",   
                        
                        "B. animalis subsp. lactis BL-04", "B. animalis subsp. lactis BI-07", "B. adolescentis",
                        
                        "E. lenta", "V. parvula", "B. crossotus", 
                        "V. cholerae N16961", 
                        "P. melaninogenica", 
                        "A. putredinis", 
                        "D. piger", "R. intestinalis", "R. hominis",
                        "P. copri", 
                        "D. formicigenerans",  "O. splanchnicus", 
                        "H. parainfluenzae", "A. shahii", 
                        "R. bromii")
  
  mucin_species_star = c("A. muciniphila*", "R. torques*", "C. perfringens S107", "C. ramosum", "C. comes", "A. odontolyticus", 
                         "C. aerofaciens", "B. clarus", "F. nucleatum subsp. nucleatum", "L. paracasei", "B. thetaiotaomicron*",
                         "B. fragilis*", "B. vulgatus*", "B. caccae* different strain", "B. eggerthii", "B. uniformis*", "B. xylanisolvens", "P. capillosus", 
                         "R. gnavus*", "B. longum subsp. infantis*", "B. longum subsp. longum* different strain")
  
  mucin_species = gsub("\\*.*", "", mucin_species_star)
  all_species = gsub("\\*.*", "", all_species_star)
  
  re_mucin_species = c(mucin_species, "perfringens")
  re_cazy_mucin = "^(GH2_20_42|GH33|CBM32|CBM40|CBM51|GH84_85|GH89|GH29_95|GH98|GH101_129)$"
  
  organisms = read.table("../data/organisms.tab", na.strings="", header=T, sep="\t", stringsAsFactors=F)
  organisms = subset(organisms, !is.na(assembly))

  ##################
  # Read data
  ##################
  organisms.cazy = data.frame()
  for(f in list.files("../data/generated/cazydb_results/")) {
    f.cazy = read.table(paste0("../data/generated/cazydb_results/", f), na.strings="", header=T, sep="\t", stringsAsFactors=F, quote="")
    f.cazy$assembly = gsub("_protein.cazy", "", f)
    organisms.cazy = rbind(organisms.cazy, f.cazy)
  }
  organisms.cazy = merge(organisms.cazy, organisms[,c("assembly", "species")], by="assembly")
  organisms.cazy$hit = gsub(".hmm", "", organisms.cazy$hit)
  organisms.cazy$hit = gsub("^(GH2|GH20|GH42)$", "GH2_20_42", organisms.cazy$hit)
  organisms.cazy$hit = gsub("^(GH84|GH85)$", "GH84_85", organisms.cazy$hit)
  organisms.cazy$hit = gsub("^(GH29|GH95)$", "GH29_95", organisms.cazy$hit)
  organisms.cazy$hit = gsub("^(GH101|GH129)$", "GH101_129", organisms.cazy$hit)
  organisms.cazy = organisms.cazy[,c("species", "assembly", "hit", "evalue")]
  #organisms.cazy  = rbind(organisms.cazy, organisms.nan)
  organisms.cazy = ddply(organisms.cazy, .(species, assembly, hit), function(z) z[which.min(z$evalue),])
  organisms.cazy$is_mucin_degrader = sapply(organisms.cazy$species, function(x) any(sapply(re_mucin_species, grepl, x)))
  organisms.cazy$in_mucin_table = organisms.cazy$species %in% mucin_species
  organisms.cazy$is_mucin_enzyme = grepl(re_cazy_mucin, organisms.cazy$hit)
  organisms.cazy$species_text = organisms.cazy$species
  organisms.cazy$species_text = all_species_star[match(organisms.cazy$species, all_species)]
  organisms.cazy$species_text = factor(organisms.cazy$species_text, rev(all_species_star))
  organisms.cazy$species = factor(organisms.cazy$species, rev(all_species))
  organisms.cazy = subset(organisms.cazy, grepl("GH", hit))
  #write.table(organisms.cazy, file="../data/supplementary/cazy.tab", sep="\t", quote=F, row.names=F, na="")
  organisms.cazy = subset(organisms.cazy, evalue < 1e-5)
  
  
  #
  # CAZy genes boxplot 
  #
  genes_cutoff = 3
  MaxOD8_effect_cutoff = log2(2)
  MaxOD9_cutoff = 0.05
  reps_cutoff = 0.5
  
  organisms.cazy_count = ddply(organisms.cazy, .(species, species_text, is_mucin_degrader, in_mucin_table), summarize, mucin_genes=sum(is_mucin_enzyme))
  organisms.cazy_count$type = ifelse(organisms.cazy_count$is_mucin_degrader, "Mucin degrader", "Other")
  
  curves_annotation = read.table("../data/curves_annotation.tab", na.strings="", header=T, sep="\t", stringsAsFactors=F)
  curves_annotation = subset(curves_annotation, !is.na(Species) & is.na(ConditionSpecies) & Passage==1)
  
  curves_annotation.9 = subset(curves_annotation, !is.na(Species) & is.na(ConditionSpecies) & !is.na(Media) & Media==9)
  curves_annotation.9 = ddply(curves_annotation.9, .(Species), summarize, MaxOD9=median(MaxOD-BlankOD), growth9=sum(!grepl("NoGrowth", Class) & MaxOD9>MaxOD9_cutoff), replicates9=length(Class))
  curves_annotation.9$hit9 = with(curves_annotation.9, (growth9 / replicates9) > reps_cutoff & MaxOD9>=MaxOD9_cutoff)
  curves_annotation.8 = subset(read.table("../data/generated/curves_annotation_rel.tab", na.strings="", header=T, sep="\t", stringsAsFactors=F), Media==8, c(Species, Media, MaxOD.fold, MaxOD.pval, replicates))
  curves_annotation.8$growth8 = as.numeric(gsub("/.*", "", curves_annotation.8$replicates))
  curves_annotation.8$replicates8 = as.numeric(gsub("[^/]+/", "", curves_annotation.8$replicates))
  curves_annotation.8$replicates = NULL
  colnames(curves_annotation.8)[3:4] = c("MaxOD8_effect", "MaxOD8_effect_pval")
  curves_annotation.8$hit8 = with(curves_annotation.8, MaxOD8_effect >= MaxOD8_effect_cutoff & (growth8/replicates8) > reps_cutoff)

  mucin.influence.f = merge(organisms.cazy_count, curves_annotation.8, by.x="species", by.y="Species")
  mucin.influence.f = merge(mucin.influence.f, curves_annotation.9, by.x="species", by.y="Species")
  mucin.influence.f$genes_threshold = mucin.influence.f$mucin_genes >= genes_cutoff
  #mucin.influence.f$hit9[mucin.influence.f$species %in% c("E. coli UTI89", "E. coli HM605", "S. salivarius", "C. aerofaciens")] = F
  mucin.influence.f$is_known = with(mucin.influence.f, species %in% mucin.known)
  mucin.influence.f$pwpch = 16
  mucin.influence.f$pwpch[mucin.influence.f$is_known] = 47+1:sum(mucin.influence.f$is_known)
  mucin.influence.f$species = as.character(mucin.influence.f$species)
  mucin.influence.f$species_text = as.character(mucin.influence.f$species_text)
  
  
  data = with(mucin.influence.f, data.frame(
    species=c(species, species[!is_known], species[is_known], species[!genes_threshold], species[genes_threshold], species[hit8 & hit9], species[!hit8 & !hit9]),
    species_text=c(species_text, species_text[!is_known], species_text[is_known], species_text[!genes_threshold], species_text[genes_threshold], species_text[hit8 & hit9], species_text[!hit8 & !hit9]),
    MaxOD8_effect=c(MaxOD8_effect, MaxOD8_effect[!is_known], MaxOD8_effect[is_known], MaxOD8_effect[!genes_threshold], MaxOD8_effect[genes_threshold], MaxOD8_effect[hit8 & hit9], MaxOD8_effect[!hit8 & !hit9]),
    mucin_genes=c(mucin_genes, mucin_genes[!is_known], mucin_genes[is_known], mucin_genes[!genes_threshold], mucin_genes[genes_threshold], mucin_genes[hit8 & hit9], mucin_genes[!hit8 & !hit9]),
    dataset=c(rep("All", length(is_known)), rep("Not known", sum(!is_known)), rep("Known", sum(is_known)), rep(paste("<", genes_cutoff, "genes"), sum(!genes_threshold)), rep(paste(">=", genes_cutoff, "genes"), sum(genes_threshold)), rep("M8 and M9", sum(hit8 & hit9)), rep("not (M8 or M9)", sum(!hit8 & !hit9)))
  ), stringsAsFactors)
  data = merge(data, mucin.influence.f[,c("species", "MaxOD8_effect_pval")], by="species")
  data$is_known = data$species %in% mucin.known
  data$dataset = as.character(data$dataset)
  data = data[nrow(data):1,]
  
  data.pvalues = do.call(rbind, apply(combn(unique(data$dataset), 2), 2, function(x) {
    d.1 = subset(data, dataset==x[1])
    d.2 = subset(data, dataset==x[2])
    
    m8.res = t.test(d.1$MaxOD8_effect, d.2$MaxOD8_effect)
    genes.res = wilcox.test(d.1$mucin_genes, d.2$mucin_genes)
    data.frame(dataset1=x[1], dataset2=x[2], MaxOD8_effect_pval=m8.res$p.value, mucin_genes_pval=genes.res$p.value)
  }))
  
  pdf("../report/mucin_barplots.pdf", paper="a4r", height=8.27, width=11.69)
  pars=list(par(mar=c(5,10,3,2)))
  datasets.f = c("All", "Known", "not (M8 or M9)", "M8 and M9")
  data.f = subset(data, (dataset %in% datasets.f))
  data.f$dataset = factor(data.f$dataset, rev(datasets.f))
  boxplot(mucin_genes ~ dataset, data.f, horizontal=T, las=2, outline=F, ylab="", xlab="# of mucin genes", pars=pars)
  beeswarm(jitter(mucin_genes, factor=0.7) ~ dataset, data.f, add=T, pwcol=ifelse(data.f$is_known, "red", "black"), horizontal=T, pars=pars)
  pval = subset(data.pvalues, dataset2=="Known" & dataset1=="not (M8 or M9)")
  text(6,2.5, paste0("p < 10", ceiling(log10(pval$mucin_genes_pval))))
  dev.off()
  
  
  #
  # Hits Venn diagram 
  #
  mucin.influence.f$species_text_rep = with(mucin.influence.f, paste0(species_text, " [", 
    ifelse(hit8, paste0(growth8, "(", round(MaxOD8_effect, 1), ")"), ""),
    ifelse(hit8 & hit9, "-", ""), 
    ifelse(hit9, paste0(growth9, "(", round(MaxOD9, 1), ")"), ""),  
    "/", 
    replicates9, "]"))
  
  mucin.influence.f$species_text_rep = gsub(" different strain", "", mucin.influence.f$species_text_rep)
  hits.cazy = with(subset(mucin.influence.f, genes_threshold), species_text_rep)
  hits.8 = with(subset(mucin.influence.f, hit8), species_text_rep)
  hits.9 = with(subset(mucin.influence.f, hit9), species_text_rep)
  hits.known = paste0(mucin.known, "*")
  
  venn.species = list(hits.cazy, hits.8, hits.9, hits.known)
  names(venn.species) = c(paste(">=", genes_cutoff, "mucin degradation\nassociated genes (CAZy)"), "M8 hits: dGMM+LAB plus mucin", "M9 hits: dGMM+LAB only mucin", "Known species")
  grey.palette = colorRampPalette(c("#F0F0F0", "#969696", "#252525"))
  
  venn.polygon = venn.diagram(
    venn.species[1:3],
    fill=grey.palette(3),
    filename=NULL,
    cex=1.6,
    main.cex=1,
    cat.cex=1,
    col="transparent", main="")
  
  pdf("../report/mucin_venn.pdf",  height=8, width=8)
  plot.new()
  grid.draw(venn.polygon)
  intersect.all = intersect(intersect(hits.8, hits.9), hits.cazy)
  intersect.89 = setdiff(intersect(hits.8, hits.9), hits.cazy)
  intersect.9 = setdiff(setdiff(hits.9, hits.8), hits.cazy)
  intersect.8 = setdiff(setdiff(hits.8, hits.9), hits.cazy)
  intersect.8cazy = intersect(setdiff(hits.8, hits.9), hits.cazy)
  intersect.9cazy = intersect(setdiff(hits.9, hits.8), hits.cazy)
  text(0.05,0.6, paste(c("Degraders with known pathways:", intersect.all), collapse="\n"), adj=c(0,1))
  text(0.6,0.6, paste(c("Degraders with unknown pathways:", intersect.89), collapse="\n"), adj=c(0,1))
  text(-0.01,0.2, paste(c("8 only:", intersect.8), collapse="\n"), adj=c(0,1))
  text(0.3,0.37, paste(c("8+cazy:", intersect.8cazy), collapse="\n"), adj=c(0,1))
  text(0.55,0.37, paste(c("9+cazy:", intersect.9cazy), collapse="\n"), adj=c(0,1))
  text(0.8,0.2, paste(c("9 only:", intersect.9), collapse="\n"), adj=c(0,1))
  dev.off()
}
