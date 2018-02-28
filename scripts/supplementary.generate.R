source("analyze.functions.R")
library("plyr")
library("reshape2")
#
# Table S3. Annotated data
#
table.S3 = function() {
  curves.a = read.table("../data/curves_annotation.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="")
  curves.a = subset(curves.a, !is.na(Species) & is.na(ConditionSpecies))
  
  curves.a = subset(curves.a, is.na(ConditionSpecies))
  curves.a$TechnicalError = ifelse(grepl("Unrep|Undef", curves.a$Class) & curves.a$Passage!=2, 1, 0)
  curves.a$Unfinished = ifelse(grepl("Unfinished", curves.a$Class), 1, 0)
  curves.a$Unfinished[curves.a$TechnicalError==1] = NA
  curves.a$Growth = ifelse(!grepl("NoGrowth", curves.a$Class), 1, 0)
  curves.a$Growth[curves.a$TechnicalError==1] = NA
  curves.a$MaxOD[curves.a$TechnicalError==1] = NA
  curves.a$BlankOD[curves.a$TechnicalError==1] = NA
  curves.a$Rate[curves.a$TechnicalError==1] = NA
  curves.a$AUC8[curves.a$TechnicalError==1] = NA
  curves.a$AUC12[curves.a$TechnicalError==1] = NA
  curves.a$AUC18[curves.a$TechnicalError==1] = NA
  curves.a$AUC24[curves.a$TechnicalError==1] = NA
  curves.a = curves.a[,c("File", "TechnicalReplicates", "Media", "Passage", "Species", "TechnicalError", "Growth",  "Unfinished", "pH_48.mean", "BlankOD", "MaxOD", "Rate", "AUC8", "AUC12", "AUC18", "AUC24")]

  write.table(curves.a, file="../report/tables/S3_Annotated_data.tab", sep="\t", quote=F, row.names=F, na="")
}

#
# Table S4. Growth matrix
#
table.s4 = function()
{
  curves.a.all = read.table("../data/curves_annotation.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="")
  curves.a = subset(curves.a.all, !is.na(Species) & is.na(ConditionSpecies) & Passage==1 & !grepl("_", Media))

  growth_matrix.q1 = cuves.merge_annotations2(curves.a)
  growth_matrix.q1$Media = gsub("^([0-9])", "M\\1", growth_matrix.q1$Media)
  growth_matrix.q1$MaxOD_with_sd = paste0(round(growth_matrix.q1$MaxOD,3), ifelse(growth_matrix.q1$MaxOD>0, paste0(" (", round(growth_matrix.q1$MaxOD_sd,3), ")"), ""))
  growth_matrix.q.raw = dcast(growth_matrix.q1, Species ~ Media, value.var="MaxOD_with_sd")
  growth_matrix.q.export = growth_matrix.q.raw[,c("Species", "GMM", "BHI++", "WCA", "mGAM", paste0("M", c(1:5,7:11, 13:14, "15 A", "15 B", "16")))]
  write.table(growth_matrix.q.export, file="../report/tables/S4_growth_matrix2.tab", sep="\t", quote=F, row.names=F, na="")
  
  
  ##
  curves.new.raw = read.table("../data/S3_curves_annotation_new.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="")
  curves.new.raw = subset(curves.new.raw, !measurement.error & !is.na(species) & passage==1 & !grepl("_", media))
  curves.new = ddply(curves.new.raw, .(species, media), summarize, non_viable=sum(viable==0, na.rm=T), viable=sum(viable!=0, na.rm=T), replicates=viable+non_viable)
  
  curves.old.raw = read.table("../data/S3_curves_annotation_old.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="")
  curves.old.raw = subset(curves.old.raw, !is.na(Species) & is.na(ConditionSpecies) & !grepl("Unrep|Undefined", Shape))
  curves.old = ddply(curves.old.raw, .(Species, Media), summarize, non_viable=sum(grepl("NoGrowth", Shape), na.rm=T), viable=sum(!grepl("NoGrowth", Shape), na.rm=T), replicates=viable+non_viable)
  
  curves.ambiguous = merge(curves.old, curves.new, by.y=c("species", "media"), by.x=c("Species", "Media"), suffixes=c(".old", ".new"))
  curves.ambiguous$ambiguous.old = with(curves.ambiguous, replicates.old==1 | viable.old/replicates.old==0.5)
  curves.ambiguous$ambiguous.new = with(curves.ambiguous, replicates.new==1 | viable.new/replicates.new==0.5)
  curves.ambiguous$minority.grow = ifelse(with(curves.ambiguous, replicates.new==1 | viable.new/replicates.new<=0.5), "yes", "no")
  curves.ambiguous = subset(curves.ambiguous, , c(Species, Media, ambiguous.old, ambiguous.new, minority.grow))
  curves.ambiguous$Media = gsub(" ", "", gsub("^([0-9])", "M\\1", curves.ambiguous$Media))
  #####
  
  curves.additional = unique(subset(curves.a.all, grepl("evisio|mucin", File), c(Species, Media)))
  curves.additional$Media = gsub(" ", "", gsub("^([0-9])", "M\\1", curves.additional$Media))
  curves.additional$passaging = "yes"
  new = melt(growth_matrix.q.export, variable.name="Media", id.vars="Species", value.name="MaxOD")
  new$Media = gsub(" ", "", new$Media)
  newb = melt(growth_matrix.q.export, variable.name="Media", id.vars="Species", value.name="MaxOD.new_withblank")
  newb$Media = gsub(" ", "", newb$Media)
  old = read.table("../report/tables/S4_old.txt", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="")
  old = melt(old, variable.name="Media", id.vars="Species", value.name="MaxOD")
  old2new = merge(old, new, by=c("Species", "Media"), suffixes=c(".old", ".new"))
  old2new = merge(old2new, newb, by=c("Species", "Media"))
  old2new = merge(old2new, curves.additional, by=c("Species", "Media"), all.x=T)
  old2new = merge(old2new, curves.ambiguous, by=c("Species", "Media"))
  old2new$passaging[is.na(old2new$passaging)] = "no"
  old2new = subset(old2new, !grepl("15|16", Media))

  interesting = subset(old2new, passaging=="no" & minority.grow=="no" & abs(MaxOD.old - MaxOD.new) > 0.1)
  curves.a.old = curves.old.raw
  curves.a.old$Media = gsub("^([0-9])", "M\\1", curves.a.old$Media)
  curves.a.old.interesting = merge(curves.a.old, interesting[,c("Species", "Media")]) 
  curves.a.new = curves.a.all
  curves.a.new$Media = gsub("^([0-9])", "M\\1", curves.a.new$Media)
  curves.a.new.interesting = merge(curves.a.new, interesting[,c("Species", "Media")]) 
  
  ggplot(old2new) + 
    geom_abline(slope=1, color="#999999") +
    geom_point(aes(MaxOD.old, MaxOD.new, shape=minority.grow, color=passaging)) +
    facet_wrap(~Media) +
    guides(shape=guide_legend(title="Majority do not grow"), color=guide_legend(title="Additional experiments\nin responce to review #1")) + 
    labs(title="MaxOD before and after review #1", x="MaxOD before review #1", y="MaxOD after review #1") +
    scale_shape_manual(values=c("yes"=4, "no"=16)) +
    coord_cartesian(xlim=c(0, 1.5), ylim=c(0, 1.5)) +
    theme_bw(base_size=16)
  
  ggplot(subset(old2new, minority.grow == "no" & passaging == "no")) + 
    geom_density(aes(MaxOD.old - MaxOD.new))
    
}

#
# Table S5 Exact sample size
#
table.S5 = function() {
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
  View(growth_matrix.q.replicates)
  write.table(growth_matrix.q.replicates, file="../report/tables/S5_Exact_sample_size.tab", sep="\t", quote=F, row.names=F, na="")
}

table.S6_fig3a_pvalues = function() {
  
}
