library(ggplot2)
library(plyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(shiny)
library(VennDiagram)
library(R.utils)
library(beeswarm)
source("analyze.functions.R")

tigrfam_cutoff = 1e-5


#
# Read data
#
universe = read.table("../data/agora_reactions.tab", na.strings="", header=T, sep="\t", stringsAsFactors=F)
universe = ddply(universe, .(ReactionID, ReactionName, Subsystem, Equation), summarize, EC=unique(unlist(strsplit(EC, ","))))
universe.subsystems = ddply(universe, .(Subsystem), summarize, SubsystemExpectedCount=length(Subsystem))
universe.enzymes = universe$EC[!is.na(universe$EC)]
universe.enzymes = unique(unlist(sapply(universe.enzymes, strsplit, ",")))
universe.enzymes = universe.enzymes[grepl("\\d+\\.\\d+\\.\\d+\\.\\d+", universe.enzymes)]

organisms = read.table("../data/organisms.tab", na.strings="", header=T, sep="\t", stringsAsFactors=F)
organisms.f = organisms[complete.cases(organisms[, c("assembly", "agora", "species")]), c("assembly", "agora", "species", "family")]

medias = c("M3_M4", "M5", "M7", "M11", "M2", "M1", "M13", "M14")
growth = read.table("../report/tables/S4_growth_matrix.tab", na.strings="", header=T, sep="\t", stringsAsFactors=F)
growth$M3_M4 = rowMeans(growth[,c("M3", "M4")])
growth = growth[,c("Species", medias)]
growth$is_growing = rowMeans(growth[,medias]) > 0
growth = merge(growth, organisms, by.x="Species", by.y="species")

tigrfam_info = read.table("../data/TIGRFAMs_15.0_INFO.tsv", na.strings="", header=T, sep="\t", stringsAsFactors=F)
tigrfam_info = tigrfam_info[!duplicated(tigrfam_info$accession),]
tigrfam_info = ddply(tigrfam_info, .(accession), function(z) data.frame(enzyme=unique(unlist(strsplit(z$enzyme, " ")))))
tigrfam_info.enzymes = unique(with(tigrfam_info, enzyme[!is.na(enzyme) & grepl("\\d+\\.\\d+\\.\\d+\\.\\d+", enzyme)]))
compared.enzymes = setdiff(intersect(tigrfam_info.enzymes, universe.enzymes), "3.1.3.5")


# Load annotations TIGRFAM annotations for organisms
if(!file.exists("../tmp/organisms.tigrfam.raw.rda")) {
  organisms.tigrfam.raw = data.frame()
  for(f in list.files("../data/tigrfam15_results/")) {
    f.tigrfam = read.table(paste0("../data/tigrfam15_results/", f), na.strings="", header=T, sep="\t", stringsAsFactors=F, quote="")
    f.tigrfam$assembly = gsub("_protein.tigrfam", "", f)
    organisms.tigrfam.raw = rbind(organisms.tigrfam.raw, f.tigrfam)
  }

  organisms.tigrfam.raw = merge(organisms.tigrfam.raw, tigrfam_info, by="accession")
  organisms.tigrfam.raw = subset(organisms.tigrfam.raw, evalue < tigrfam_cutoff)
  organisms.tigrfam.raw$evalue_log10 = abs(log10(organisms.tigrfam.raw$evalue))
  organisms.tigrfam.raw_sum = ddply(organisms.tigrfam.raw, .(assembly, query), summarize, evalue_log10.max=max(evalue_log10))
  organisms.tigrfam.raw = merge(organisms.tigrfam.raw, organisms.tigrfam.raw_sum, by=c("assembly", "query"))
  save(organisms.tigrfam.raw, file="../tmp/organisms.tigrfam.raw.rda")
} else {
  load("../tmp/organisms.tigrfam.raw.rda")
}

organisms.tigrfam = subset(organisms.tigrfam.raw, evalue_log10 > 0.5*evalue_log10.max & !is.na(enzyme) & !grepl("-", enzyme))
organisms.tigrfam = merge(organisms.tigrfam, organisms.f, by="assembly")
organisms.tigrfam = organisms.tigrfam[with(organisms.tigrfam, order(agora, enzyme, evalue)),]
organisms.tigrfam = organisms.tigrfam[!duplicated(organisms.tigrfam[,c("agora", "enzyme")]),]
organisms.tigrfam = merge(organisms.tigrfam, unique(subset(universe, !is.na(EC), c(EC, Subsystem))), by.x="enzyme", by.y="EC", all.x=T)
organisms.tigrfam = ddply(organisms.tigrfam, .(agora, enzyme), function(z) {
  z.ret = z[1,]
  z.ret$Subsystem = paste(na.omit(z$Subsystem), collapse="; ")
  z.ret
})
organisms.tigrfam$Subsystem[organisms.tigrfam$Subsystem==""] = NA

#
# Read gapfills
#
gapfills = read.table("../report/tables/S8_gapfills.tsv", na.strings="", header=T, sep="\t", stringsAsFactors=F)
gapfills = merge(organisms[,c("species", "agora", "assembly")], gapfills, by.x="agora", by.y="ModelID")
gapfills$SubsystemGroup = gapfills$Subsystem
gapfills$SubsystemGroup[grepl("B2|B6|B12|Biotin|Thiamine|Folate", gapfills$Subsystem)] = "Vitamin metabolism*"
gapfills$SubsystemGroup[grepl("Arginine|Glycine|Valin|Lysine|Tryptophan|Alanine|Methionine|Glutamate|Phenylalanine|Histidine", gapfills$Subsystem)] = "Amino acid metabolism*"
gapfills$SubsystemGroup[grepl("Fatty acid", gapfills$Subsystem)] = "Fatty acid metabolism*"
gapfills$SubsystemGroup[grepl("Galactose|Glycolisis|Starch", gapfills$Subsystem)] = "Carbon metabolism*"
gapfills.gaps_export = subset(gapfills, !is.na(MediumID) & !is.na(Growing) & Growing, c(species, MediumID, agora, EC, ReactionID, ReactionName, Subsystem))

write.table(gapfills.gaps_export, file="../report/tables/S8_gapfills_export.tab", sep="\t", quote=F, row.names=F, na="")


gapfills.x = ddply(gapfills, .(species, MediumID, agora, assembly, ReactionID, ReactionName, Subsystem, Growing), function(z) data.frame(EC=unique(unlist(strsplit(z$EC, " |,")))))
gapfills.y = ddply(gapfills, .(species, MediumID, agora, assembly, ReactionID, ReactionName, Subsystem, SubsystemGroup, Growing), function(z) data.frame(EC=unique(unlist(strsplit(z$EC, " |,")))))
dim(gapfills.x)
dim(gapfills.y)

#
# Split EC numbers
#
gapfills$EC[grepl("2.4.2.82.4.2.7", gapfills$EC)] = "2.4.2.8,2.4.2.7"
gapfills = ddply(gapfills, .(species, MediumID, agora, assembly, ReactionID, ReactionName, Subsystem, SubsystemGroup, Growing), function(z) data.frame(EC=unique(unlist(strsplit(z$EC, " |,")))))
gapfills.gaps = subset(gapfills, !is.na(MediumID) & !is.na(Growing) & Growing)
gapfills.orig = subset(gapfills, is.na(MediumID))

species.summary = ddply(gapfills, .(species), summarize, any(Growing[is.na(MediumID)]))
species.summary = ddply(gapfills, .(species), summarize, 
                        is_agora        = any(!is.na(agora)), # Has AGORA model
                        is_gapfilled    = is_agora && any(!is.na(MediumID)),  # This model was gapfilled because +InVivo/-InSilico
                        is_successfull  = !is_gapfilled || is_gapfilled && any(!is.na(Growing) & Growing)) # This model was successfully gapfilled +InVivo/-InSilico/+Growing
species.summary = merge(species.summary, growth[,c("Species", "is_growing")], by.x="species", by.y="Species")

species.summary_pie = data.frame()
species.summary_pie = rbind(species.summary_pie, data.frame(variable="Not growing on tested medias",  value=sum(!species.summary$is_growing)))
species.summary_pie = rbind(species.summary_pie, data.frame(variable="No AGORA model",                value=sum(!species.summary$is_agora & species.summary$is_growing)))
species.summary_pie = rbind(species.summary_pie, data.frame(variable="AGORA model can't be imrpoved", value=sum(!species.summary$is_gapfilled & species.summary$is_growing & species.summary$is_agora)))
species.summary_pie = rbind(species.summary_pie, data.frame(variable="Gapfilling failed",             value=sum(!species.summary$is_successfull & species.summary$is_gapfilled & species.summary$is_growing & species.summary$is_agora)))
species.summary_pie = rbind(species.summary_pie, data.frame(variable="Successfully gapfilled",        value=sum(species.summary$is_successfull & species.summary$is_gapfilled & species.summary$is_growing & species.summary$is_agora)))


gapfill.overview_before = c(
  "Valid models*"=nrow(subset(species.summary, !is_gapfilled & is_growing & is_agora)),
  "Incomplete models"=nrow(subset(species.summary, is_gapfilled & is_growing & is_agora)))
gapfill.overview_after = c(
  "Valid models*"=nrow(subset(species.summary, is_successfull & is_growing & is_agora)),
  "Incomplete models"=nrow(subset(species.summary, !is_successfull & is_growing & is_agora)))
gapfill.overview_ggplot = data.frame(
  when=factor(rep(c("before", "after"), each=2), c("before", "after")),
  subset=factor(rep(c("Valid models*", "Incomplete models"), 2)),
  count=c(gapfill.overview_before, gapfill.overview_after)
)
gapfill.overview_ggplot = ddply(gapfill.overview_ggplot, .(when), mutate, cumsum=cumsum(count) - count/2, ytext=cumsum - count/2)
gapfill.overview_ggplot_valid = subset(gapfill.overview_ggplot, subset=="Valid models*")

pdf("../report/gapfill_summary.pdf", paper="a4", height=11.69, width=8.27)
ggplot(gapfill.overview_ggplot) +
  geom_bar(aes(x=when, y=count, fill=subset), stat="identity", position="stack", width=0.5) +
  geom_text(aes(x=when, y=cumsum, label=count, color=subset), stat="identity", position="stack") +
  geom_polygon(fill="#666666",
               x=with(gapfill.overview_ggplot_valid, c(when, rev(when))+0.25*c(1,-1,-1,1)),
               y=with(gapfill.overview_ggplot_valid, c(count, 0, 0))) +
  scale_fill_manual(values=c("#BBBBBB", "#666666"))+
  scale_color_manual(values=c("#666666", "#BBBBBB"))+
  labs(x="", fill="", y="", title="Improved AGORA models") +
  guides(color=F) +
  theme_classic(base_size=14)  +
  theme(axis.text=element_blank(), axis.ticks=element_blank())
dev.off()


#
# HEATMAP with all gapfilled SUBSYSTEMs
#
gapfills.orig.f = subset(gapfills.orig, EC %in% compared.enzymes & !is.na(agora) & agora %in% organisms.tigrfam$agora)
gapfills.orig.f = gapfills.orig.f[!duplicated(gapfills.orig.f[,c("agora", "assembly", "species","Subsystem", "EC")]),]
gapfills.orig.f$enzyme_agora = gapfills.orig.f$EC
gapfills.gaps.f = subset(gapfills.gaps, EC %in% compared.enzymes & !is.na(agora) & agora %in% organisms.tigrfam$agora)
gapfills.gaps.f = gapfills.gaps.f[!duplicated(gapfills.gaps.f[,c("agora", "assembly", "species","Subsystem", "EC")]),]
gapfills.gaps.f$enzyme_gf = gapfills.gaps.f$EC
organisms.tigrfam.f = subset(organisms.tigrfam, enzyme %in% compared.enzymes & agora %in% gapfills.gaps.f$agora)
organisms.tigrfam.f = organisms.tigrfam.f[!duplicated(organisms.tigrfam.f[,c("agora", "assembly", "species","enzyme")]),]
organisms.tigrfam.f$enzyme_tigrfam = organisms.tigrfam.f$enzyme
organisms.tigrfam.f$Subsystem_tigrfam = organisms.tigrfam.f$Subsystem

gapfills.gaps_tigrfam = merge(gapfills.orig.f, gapfills.gaps.f, by.x=c("agora", "assembly", "species","EC"), by.y=c("agora", "assembly", "species","EC"), all=T, suffixes=c("_agora", "_gf"))
gapfills.gaps_tigrfam = merge(gapfills.gaps_tigrfam, organisms.tigrfam.f, by.x=c("assembly", "species", "EC"), by.y=c("assembly","species", "enzyme"), all=T, suffixes=c("", "_tigrfam"))
gapfills.gaps_tigrfam$ReactionID = gapfills.gaps_tigrfam$ReactionID_agora
gapfills.gaps_tigrfam$ReactionID[is.na(gapfills.gaps_tigrfam$ReactionID)] = gapfills.gaps_tigrfam$ReactionID_gf[is.na(gapfills.gaps_tigrfam$ReactionID)]
gapfills.gaps_tigrfam$Subsystem[is.na(gapfills.gaps_tigrfam$Subsystem)] = gapfills.gaps_tigrfam$Subsystem_gf[is.na(gapfills.gaps_tigrfam$Subsystem)]
gapfills.gaps_tigrfam$Subsystem = gapfills.gaps_tigrfam$Subsystem_agora
gapfills.gaps_tigrfam$Subsystem[is.na(gapfills.gaps_tigrfam$Subsystem)] = gapfills.gaps_tigrfam$Subsystem_gf[is.na(gapfills.gaps_tigrfam$Subsystem)]
gapfills.gaps_tigrfam$Subsystem[is.na(gapfills.gaps_tigrfam$Subsystem)] = gapfills.gaps_tigrfam$Subsystem_tigrfam[is.na(gapfills.gaps_tigrfam$Subsystem)]
gapfills.gaps_tigrfam$SubsystemGroup = gapfills.gaps_tigrfam$Subsystem
gapfills.gaps_tigrfam$SubsystemGroup[grepl("B2|B6|B12|Biotin|Thiamine|Folate", gapfills.gaps_tigrfam$Subsystem)] = "Vitamin metabolism*"
gapfills.gaps_tigrfam$SubsystemGroup[grepl("Arginine|Glycine|Valin|Lysine|Tryptophan|Alanine|Methionine|Glutamate|Phenylalanine|Histidine", gapfills.gaps_tigrfam$Subsystem)] = "Amino acid metabolism*"
gapfills.gaps_tigrfam$SubsystemGroup[grepl("Fatty acid", gapfills.gaps_tigrfam$Subsystem)] = "Fatty acid metabolism*"
gapfills.gaps_tigrfam$SubsystemGroup[grepl("Galactose|Glycolisis|Starch", gapfills.gaps_tigrfam$Subsystem)] = "Carbon metabolism*"

#gapfills.gaps_tigrfam_subsystem = merge(gapfills.gaps_tigrfam, gapfills.subsystem_prop, by="Subsystem", all.x=T)
gapfills.gaps_tigrfam_subsystem = ddply(gapfills.gaps_tigrfam, .(species, SubsystemGroup), summarize, 
                                        ConfirmedProportion=mean(!is.na(enzyme_tigrfam[!is.na(enzyme_gf)])), 
                                        GapfilledCount=sum(!is.na(enzyme_gf)),
                                        MediaList=list(na.omit(unique(MediumID_gf))))
gapfills.gaps_tigrfam_subsystem$ConfirmedProportion[is.nan(gapfills.gaps_tigrfam_subsystem$ConfirmedProportion)] = 0

gapfills.gaps_tigrfam_subsystem.gp_enzymes = pheatmap_matrix(gapfills.gaps_tigrfam_subsystem, species ~ SubsystemGroup, column="GapfilledCount", na.value=0)
gapfills.gaps_tigrfam_subsystem.gp_enzymes = gapfills.gaps_tigrfam_subsystem.gp_enzymes[rowSums(gapfills.gaps_tigrfam_subsystem.gp_enzymes)>0, colSums(gapfills.gaps_tigrfam_subsystem.gp_enzymes)>0]
gapfills.gaps_tigrfam_subsystem.confirmed = pheatmap_matrix(gapfills.gaps_tigrfam_subsystem, species ~ SubsystemGroup, column="ConfirmedProportion", na.value=0)
gapfills.gaps_tigrfam_subsystem.confirmed = gapfills.gaps_tigrfam_subsystem.confirmed[rownames(gapfills.gaps_tigrfam_subsystem.gp_enzymes), colnames(gapfills.gaps_tigrfam_subsystem.gp_enzymes)]


pdf("../report/gapfill.pdf", paper="a4", height=11.69, width=8.27)
pheatmap(gapfills.gaps_tigrfam_subsystem.gp_enzymes, display_numbers=pheatmap_prop_text(gapfills.gaps_tigrfam_subsystem.confirmed, percent=T),
         cluster_rows=F, fontsize_row=6, fontsize_number=6, clustering_method="average",
         clustering_distance_cols=dist(t(gapfills.gaps_tigrfam_subsystem.gp_enzymes)),
         color=c("#777777", colorRampPalette(brewer.pal(n=9, name="Blues"))(20)), 
         breaks=c(-0.1, seq(0, max(gapfills.gaps_tigrfam_subsystem.gp_enzymes), length.out=20)),
         border_color="#CACACA33", number_color="#88419D", drop_levels=F,
         main="Gapfilled SUBSYSTEMS reactions (confirmed with TIGRFAM)")

selected.subsystem = "Ubiquinone and other terpenoid-quinone biosynthesis"
gapfills.gaps_tigrfam_subsystem.summary = ddply(gapfills.gaps_tigrfam_subsystem, .(SubsystemGroup), function(z) {
  ret = data.frame(SpeciesGapfilledCount = length(unique(z$species[z$GapfilledCount>0])))
  ret$SpeciesGapfilledProportion=ret$SpeciesGapfilledCount/length(unique(gapfills.gaps$species)) 
  ret$SpeciesConfirmedCount=length(unique(z$species[z$ConfirmedProportion*z$GapfilledCount>0]))
  ret$SpeciesConfirmedProportion=ret$SpeciesConfirmedCount/ret$SpeciesGapfilledCount
  ret$ConfirmedProportion=ifelse(sum(z$GapfilledCount)>0, sum(z$ConfirmedProportion*z$GapfilledCount)/sum(z$GapfilledCount), 0) 
  ret$GapfilledCount=sum(z$GapfilledCount)
  ret$SubsystemGroupText=paste0(z$SubsystemGroup[1], " (species: ", ret$SpeciesConfirmedCount, "/", ret$SpeciesGapfilledCount, ", reactions:", ret$GapfilledCount, ")")
  ret$SelectedSubsystem=grepl(selected.subsystem, z$SubsystemGroup[1])
  
  z.medialist = table(unlist(z$MediaList))
  z.medialist = paste0(names(z.medialist), ":", z.medialist, collapse=", ")
  if(z.medialist==":") z.medialist = ""
  z.medialist = gsub("M3_M4", "M3", z.medialist)
  ret$MediaList = z.medialist

  ret
})

#    SpeciesGapfilledCount=length(unique(species[GapfilledCount>0])),
#    SpeciesGapfilledProportion=SpeciesGapfilledCount/length(unique(gapfills.gaps$species)), 
#    SpeciesConfirmedCount=length(unique(species[ConfirmedProportion*GapfilledCount>0])), 
#    SpeciesConfirmedProportion=SpeciesConfirmedCount/SpeciesGapfilledCount, 
#    ConfirmedProportion=ifelse(sum(GapfilledCount)>0, sum(ConfirmedProportion*GapfilledCount)/sum(GapfilledCount), 0), 
#    GapfilledCount=sum(GapfilledCount),
#    SubsystemGroupText=paste0(SubsystemGroup[1], " (species: ", SpeciesConfirmedCount, "/", SpeciesGapfilledCount, ", reactions:", GapfilledCount, ")"),
#    SelectedSubsystem=grepl(selected.subsystem, SubsystemGroup[1]),
#    MediaList=MediaList ), collapse=";"))
gapfills.gaps_tigrfam_subsystem.summary = gapfills.gaps_tigrfam_subsystem.summary[order(gapfills.gaps_tigrfam_subsystem.summary$ConfirmedProportion),]
gapfills.gaps_tigrfam_subsystem.summary$SubsystemGroupText = factor(gapfills.gaps_tigrfam_subsystem.summary$SubsystemGroupText, gapfills.gaps_tigrfam_subsystem.summary$SubsystemGroupText)
gapfills.gaps_tigrfam_subsystem.summary = subset(gapfills.gaps_tigrfam_subsystem.summary, GapfilledCount > 1 & ConfirmedProportion>0)
gapfills.gaps_tigrfam_subsystem.summary = melt(gapfills.gaps_tigrfam_subsystem.summary, measure.vars=c("ConfirmedProportion", "SpeciesGapfilledProportion"))
gapfills.gaps_tigrfam_subsystem.summary = gapfills.gaps_tigrfam_subsystem.summary[order(rev(gapfills.gaps_tigrfam_subsystem.summary$variable), gapfills.gaps_tigrfam_subsystem.summary$value),]

ggplot(gapfills.gaps_tigrfam_subsystem.summary) +
  geom_bar(aes(x=paste0(SubsystemGroupText, ": ", MediaList), y=value, fill=variable), stat="identity", position=position_dodge(width=0.6), width=0.6) +
  coord_flip() + theme_bw(base_size=14) +
  scale_fill_manual(values=c("#BBBBBB", "#666666"))+
  labs(x="", y="%, percentage", fill="") +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"),
        legend.position = c(0.8, 0.2))

gapfills.gaps_tigrfam_subsystem.summary.long = ddply(gapfills.gaps_tigrfam_subsystem.summary, .(SubsystemGroup) , summarize, MediaStr=unlist(strsplit(MediaList, ", ")))
gapfills.gaps_tigrfam_subsystem.summary.long$Media = gsub(":.*", "", gapfills.gaps_tigrfam_subsystem.summary.long$MediaStr)
gapfills.gaps_tigrfam_subsystem.summary.long$Count = as.numeric(gsub(".*:", "", gapfills.gaps_tigrfam_subsystem.summary.long$MediaStr))
gapfills.gaps_tigrfam_subsystem.summary.long = gapfills.gaps_tigrfam_subsystem.summary.long[!duplicated(gapfills.gaps_tigrfam_subsystem.summary.long$MediaStr),]
ggplot(gapfills.gaps_tigrfam_subsystem.summary.long) +
  geom_bar(aes(x=SubsystemGroup, color=Media, y=Count), stat="identity", position="stack") +
  theme_classic() + 
  scale_color_discrete(drop=FALSE)


gapfills.gaps_tigrfam.f = subset(gapfills.gaps_tigrfam, grepl(selected.subsystem, SubsystemGroup))
list.orig = unique(with(subset(gapfills.gaps_tigrfam.f, !is.na(enzyme_agora)), paste(species, enzyme_agora)))
list.gapfill = unique(with(subset(gapfills.gaps_tigrfam.f, !is.na(enzyme_gf)), paste(species, enzyme_gf)))
list.tigrfam = unique(with(subset(gapfills.gaps_tigrfam.f, !is.na(enzyme_tigrfam)), paste(species, enzyme_tigrfam)))
grey.palette = colorRampPalette(c("#F0F0F0", "#969696", "#252525"))

venn.polygon = venn.diagram(
  list("TIGRFAM"=list.tigrfam, "ORIGINAL"=list.orig, "GAPFILL"=list.gapfill),
  fill=grey.palette(3),
  filename=NULL,
  col="#00000000", main=selected.subsystem)
plot.new()
grid.draw(venn.polygon)
dev.off()


#
# Ubiquinone gapfill
#
triptophan.pathway = function() {
  #curves.f = subset(curves,
                    
  curves.plantarum = subset(curves, grepl("plantarum", Species) & !(Media %in% c("1", "15 A", "15 B", "16")))
  curves.plantarum = ddply(curves.plantarum, .(File, Volume, Media, Species, Time, TechnicalReplicates), summarize, OD=mean(OD))
  ggplot(curves.plantarum) +
    geom_line(aes(x=Time/3600, y=OD, group=paste(File, TechnicalReplicates), color=Volume)) + 
    facet_wrap(~Media) +
    theme_bw()
  
  x = gapfills
  x$Subsystem = gsub(" synthesis| metabolism|beta-", "", x$Subsystem)
  x.subsystem = with(ddply(x, .(Subsystem), summarize, Count=sum(!is.na(MediumID))), Subsystem[Count>0])
  x.species = with(ddply(x, .(species), summarize, Count=sum(!is.na(MediumID))), species[Count>0])
  x = x[order(x$ReactionID),]
  x = subset(x, species %in% x.species & Subsystem %in% x.subsystem)
  x = dcast(x, species ~ Subsystem, value.var="MediumID", fun.aggregate=function(z) ifelse(length(z)==0, "", ifelse(any(!is.na(z)), na.omit(z)[1], "+")))
  rownames(x) = x$species  
  x$species = NULL
  x.values = x
  map_dist = c("+"=0, "M3_M4"=1,"M5"=2, "M7"=3, "M11"=4, "M2"=5, "M1"=6, "M14"=7)
  x.values = apply(x.values, 1:2, function(z) map_dist[z])
  x.values[is.na(x.values)] = -1
  
  pheatmap(x.values, display_numbers=x)

  x = subset(gapfills, grepl("Glycine|Valin|Lysine|Tryptophan|Alanine|Methionine|Glutamate|Phenylalanine|Histidine|Ubiquinone|Shikimate", Subsystem))
  x$ReactionID = paste0(gsub(" metabolism|beta-", "", x$Subsystem), ": ", x$ReactionID, " (", x$EC, ")")
  x = x[order(x$ReactionID),]
  x = dcast(x, species ~ ReactionID, value.var="MediumID", fun.aggregate=function(z) ifelse(length(z)==0, "", ifelse(any(!is.na(z)), na.omit(z)[1], "+")))
  x = x[rowMeans(apply(x, 1:2, function(z) grepl("^M.*", z)))>0,]
  x.sum = x[1,,drop=F]
  x.sum[,-1] = colSums(apply(x[,-1], 1:2, function(z) grepl("^M.*", z)))
  x.sum[,1] = "Summary"
  x = rbind(x, x.sum)
  write.table(x, file="c:/Users/Sergej/Desktop/paula_AA.txt", sep="\t", quote=F, row.names=F, na="")

  x = gapfills
  x$EC = as.character(x$EC)
  x.long = ddply(x, .(species, ReactionID, Subsystem, MediumID), summarize, EC=unique(unlist(strsplit(EC, ",|;| "))))
  organisms.tigrfam$enzyme_tigrfam = as.character(organisms.tigrfam$enzyme)
  organisms.tigrfam$enzyme_tigrfam[is.na(organisms.tigrfam$enzyme_tigrfam)] = ""
  x = merge(x.long, subset(organisms.tigrfam,,-Subsystem), by.x=c("species", "EC"), by.y=c("species", "enzyme"), all.x=T)
  x$Subsystem = gsub(" synthesis| metabolism|beta-", "", x$Subsystem)
  x.subsystem = with(ddply(x, .(Subsystem), summarize, Gapfilled=any(!is.na(MediumID))), Subsystem[Gapfilled])
  x.species = with(ddply(x, .(species), summarize, Gapfilled=any(!is.na(MediumID))), species[Gapfilled])
  x = subset(x, species %in% x.species & Subsystem %in% x.subsystem)
  x$ReactionID = paste0(gsub(" metabolism|beta-", "", x$Subsystem), ": ", x$ReactionID, " (", x$EC, ")")
  x = x[order(x$ReactionID),]
  x.medium = dcast(x, ReactionID ~ species, value.var="MediumID", fun.aggregate=function(z) ifelse(length(z)==0, "", ifelse(any(!is.na(z)), na.omit(z)[1], "+")))
  rownames(x.medium) = x.medium$ReactionID
  x.medium$ReactionID = NULL
  x.medium = as.matrix(x.medium)
  x.enzyme = dcast(x, ReactionID ~ species, value.var="enzyme_tigrfam", fun.aggregate=function(z) {
    ifelse(length(z)==0, "", na.omit(z)[1])
  })
  rownames(x.enzyme) = x.enzyme$ReactionID
  x.enzyme$ReactionID = NULL
  x.enzyme = as.matrix(x.enzyme)
  x.enzyme[is.na(x.enzyme)] = ""
  x.final = matrix(paste(x.medium, x.enzyme), ncol=ncol(x.enzyme))
  colnames(x.final) = colnames(x.enzyme)
  rownames(x.final) = rownames(x.enzyme)
  View(x.final)
  write.table(x.final, file="c:/Users/Sergej/Desktop/gapfilled_reactions.txt", sep="\t", quote=F, row.names=T, col.names=T, na="")
}


#
# Ubiquinone gapfill
#
ubiquinone.pathway = function() {
  ubiquinone_reactions = c("R_5HLTDL", "R_5HXKYNDCL", "R_ACACT1r", "R_ANPRT", "R_ANS", "R_ANS2", "R_IGPS", "R_LTDCL", "R_PRAI", "R_TRPAS2","R_TRPS1","R_TRPS2r","R_TRPS3r")
  ubiquinone_reactions = intersect(ubiquinone_reactions, gapfills$ReactionID)
  x = ddply(subset(gapfills, ReactionID %in% ubiquinone_reactions), .(species), function(z) {
    z.ret = z[1,"species", drop=F]
    for(r_id in ubiquinone_reactions) {
      z.ret[[r_id]] = ifelse(r_id %in% z$ReactionID, ifelse(is.na(z$MediumID[which(r_id==z$ReactionID)[1]]), "Original", z$MediumID[which(r_id==z$ReactionID)[1]]), NA)
    }
    z.ret
  })
  
  x = x[rowMeans(is.na(x[,ubiquinone_reactions]))<2/3,]
  branch_reactions = c("R_DHNAOT", "R_AMMQT8r", "R_DHNAOPT", "R_PHYQS")
  x$R_Branch = ifelse(!is.na(x$R_DHNAOT) & !is.na(x$R_AMMQT8r) | !is.na(x$R_DHNAOPT) & !is.na(x$R_PHYQS), "Original", NA)
  
  shortcut_reactions = c("R_DHNAS", "R_NCOAH", "R_NPHS")
  x$R_Shortcut = ifelse(!is.na(x$R_DHNAS) & !is.na(x$R_NCOAH) | !is.na(x$R_NPHS), "Original", NA)
  
  ubiquinone_reactions_f = c(setdiff(setdiff(ubiquinone_reactions, branch_reactions), shortcut_reactions), c("R_Branch", "R_Shortcut"))
  x.gapfilled = x[rowSums(!is.na(x[,ubiquinone_reactions]) & x[,ubiquinone_reactions]!="Original")>0, c("species", ubiquinone_reactions)]
  x.gapfilled
  
  ubiquinone_enzymes = c("5.4.99.6", "4.2.99.20", "2.2.1.9", "4.2.1.113", "6.2.1.26", "4.1.3.36", "3.1.2.28", "2.5.1.74", "1.3.99.1")
  x.tigrfam = subset(tigrfam.df, species %in% x.gapfilled$species & enzyme %in% ubiquinone_enzymes)
  x.tigrfam = ddply(x.tigrfam, .(species), function(z) {
    z.ret = z[1,"species", drop=F]
    for(ec in ubiquinone_enzymes) {
      z.ret[[ec]] = ec %in% z$enzyme
    }
    z.ret
  })
  x.tigrfam
}


