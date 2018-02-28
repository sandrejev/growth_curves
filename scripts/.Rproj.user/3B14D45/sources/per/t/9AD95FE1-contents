library(XLConnect)
library(reshape2)
library(plyr)
library(ggplot2)
library(RSQLite)


media.annotation = read.delim("../data/media_names.tab", sep="\t", na="", quote="", stringsAsFactors=F, header=T)
media.annotation = media.annotation[order(media.annotation$Order),]
media.names = sapply(media.annotation$ShortName, function(x) media.annotation$FullName[match(x, media.annotation$ShortName)])
media.names2 = sapply(media.annotation$ShortName, function(x) media.annotation$FullName2[match(x, media.annotation$ShortName)])
media.general = media.annotation$ShortName[media.annotation$IsGeneral==1]
media.rich = media.annotation$ShortName[media.annotation$IsRich==1]
media.files = sapply(media.annotation$ShortName, function(x) media.annotation$Filename[match(x, media.annotation$ShortName)])

grey.palette = colorRampPalette(c("#F0F0F0", "#969696", "#252525"))

piechart.species_summary = function()
{
  species.summary = data.frame(name=">=1% rel.\nabundance\n>=10%\nprevalence", number=58)
  species.summary = rbind(species.summary, data.frame(name="pathogens", number=13))
  species.summary = rbind(species.summary, data.frame(name="probiotics", number=13))
  species.summary = rbind(species.summary, data.frame(name="linked to\ncolorectal cancer", number=3))
  species.summary = rbind(species.summary, data.frame(name="linked to\nfoodborne disease", number=2))
  species.summary = rbind(species.summary, data.frame(name="causing bacteremia\nor sepsis", number=2))
  species.summary = rbind(species.summary, data.frame(name="additional\nrepresentative\nof genus", number=3))
  species.summary = rbind(species.summary, data.frame(name="representating separate\nmetabolic clades", number=2))
  
  pdf("../report/species_selection_overview.pdf", height=6, width=8)
  pie(species.summary$number, labels=species.summary$name, col=grey.palette(nrow(species.summary)), clockwise=T, init.angle=90)
  dev.off()
}

piechart.media_summary = function()
{
  wb = XLConnect::loadWorkbook("../data/compounds.xlsx")
  df = XLConnect::readWorksheet(wb, sheet=1, header=T)
  
  result = list(cols=list(annotations=c("Class", "Compound")))
  result$cols$medias = setdiff(colnames(df), result$cols$annotations)
  df$Class[df$Class=="Buffer compounds"] = "Buffer"
  
  data.regexp = "([0-9]+(\\.[0-9]+)?) (.*)"
  data.unitdefs = list(
    "mg" = list(canonic="g", factor=1.0e-3),
    "g" = list(canonic="g", factor=1),
    "µg" = list(canonic="g", factor=1.0e-6),
    "mL" = list(canonic="L", factor=1.0e-3)
  )
  
  medias = intersect(colnames(df), paste0("M", 1:11))
  df = df[rowSums(!is.na(df[,medias]))>0,]
  data = df[,medias]
  data.numbers = apply(data, 1:2, function(x) as.numeric(gsub(data.regexp, "\\1", x)))
  data.units = apply(data, 1:2, function(x) gsub(data.regexp, "\\3", x))
  
  df.norm = apply(data.units, 1:2, function(x) ifelse(is.na(x), NA, data.unitdefs[[x]]$factor)) * data.numbers
  df.norm[is.na(df.norm)] = 0
  df.norm = cbind(df[,c("Class", "Compound")], df.norm)
  
  media_matrix.long = melt(df.norm, id.vars=c("Class", "Compound"), variable.name="Media", value.name="Amount")
  compound_sum = ddply(subset(media_matrix.long, Media=="M3"), .(Compound), summarize, Amount=sum(Amount, na.rm=T))
  rownames(compound_sum) = compound_sum$Compound
  
  media_matrix.long_s = ddply(media_matrix.long, .(Class, Media), function(z) {
    amount = mean(ifelse(compound_sum[z$Compound, "Amount"] == 0, as.numeric(z$Amount>0), z$Amount / compound_sum[z$Compound, "Amount"]))
    data.frame(Media=z$Media[1], Amount=amount)
  })
  
  media_colors = c("M1"="#F46D43", "M2"="#D53E4F", "M3"="#9E0142", "M4"="#FFDB00FF", "M5"="#49FF00FF", "M7"="#00FF92FF", "M8"="#0092FFFF", "M9"="#4900FFFF", "M10"="#FF00DBFF", "M11"="#FF00BFFF")
  media_colors = c("M1"="#F28260", "M2"="#D35665", "M3"="#9B2E5C", "M4"="#FFEB7F", "M5"="#A5FF7F", "M7"="#7FFFC7", "M8"="#7FC9FF", "M9"="#A37FFF", "M10"="#FF7FED", "M11"="#FF7FDF")
  media_colors = rev(media_colors)
  category_colors = c("Salts & Minerals"="#E6AB02", "Vitamins & Antioxidants"="#E7298A", "SCFA"="#D95F02", "Sugar"="#1B9E77", "Mucin"="#386CB0", 
                      "Amino acids"="#7570B3", "Nucleotids"="#A6761D", "Others"="#666666", "Buffer"="#66A61E")
  
  media_matrix.long_s$Amount[media_matrix.long_s$Amount==0] = NA
  media_matrix.long_s = media_matrix.long_s[order(as.numeric(gsub("M", "", media_matrix.long_s$Media))),]
  media_matrix.long_s$MediaName = media.names2[gsub("M", "", media_matrix.long_s$Media)]
  media_matrix.long_s$MediaName = factor(media_matrix.long_s$MediaName, rev(unique(media_matrix.long_s$MediaName)))
  media_matrix.long_s$Class = factor(media_matrix.long_s$Class, names(category_colors))
  media_matrix.long_s$MediaComposition = with(media_matrix.long_s, ifelse(Amount>=1 | (Media %in% c("M1", "M2", "M3")), "GMM+LAB", "Modified"))
  media_matrix.long_s$MediaComposition[media_matrix.long_s$Class=="Salts & Minerals" & media_matrix.long_s$Media %in% paste0("M", 5:11)] = "GMM+LAB"
  names(media_colors) = media.names2[gsub("M", "", names(media_colors))]
  
  pdf("../report/media_selection_overview.pdf", height=6, width=8)
  media_matrix.long_s_vline = unique(media_matrix.long_s[,"Class", drop=F])
  ggplot(media_matrix.long_s) +
    geom_hline(aes(yintercept=y), data=data.frame(y=seq(1, nlevels(media_matrix.long_s$MediaName), 2)), size=13, color="#F0F0F0") + 
    geom_hline(aes(yintercept=y), data=data.frame(y=0:nlevels(media_matrix.long_s$MediaName)+0.5), size=0.25, color="#D9D9D9", lty=2) + 
    geom_hline(yintercept=nlevels(media_matrix.long_s$MediaName)-2.5, size=2, color="#D9D9D9") + 
    geom_vline(aes(xintercept=x), data=data.frame(x=0:nlevels(media_matrix.long_s$Class)+0.5), size=0.25, color="#D9D9D9") + 
    geom_point(aes(x=as.numeric(Class), y=as.numeric(MediaName), size=Amount, color=Class, alpha=Amount), shape=20) + 
    geom_point(aes(x=as.numeric(Class), y=as.numeric(MediaName)), size=1, color="#FFFFFF", fill="#FFFFFF", shape=21, data=subset(media_matrix.long_s, MediaComposition=="Modified")) + 
    theme_classic(base_size=16) +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), axis.ticks=element_blank()) +
    labs(x="", y="") + 
    scale_fill_manual(values=category_colors, guide=F) + 
    scale_color_manual(values=category_colors, guide=F) + 
    scale_alpha_continuous(range=c(0.4, 0.8), guide=F) +
    scale_size_area(max_size=16, guide=F) +
    scale_shape_manual(values=c("GMM+LAB"=16, "Modified"=1), guide=F) +
    scale_x_continuous(breaks=1:nlevels(media_matrix.long_s$Class), labels=levels(media_matrix.long_s$Class), limits=c(0.95, nlevels(media_matrix.long_s$Class)+0.05)) +
    scale_y_continuous(breaks=1:nlevels(media_matrix.long_s$MediaName), labels=levels(media_matrix.long_s$MediaName), limits=c(0.95, nlevels(media_matrix.long_s$MediaName)+0.25))
  
    media.summary = data.frame(name=c("newly compounded defined media (10)", "described minimal and defined media (5)", "rich media (4)"))
    media.summary$number = as.numeric(gsub(".*\\((\\d+)\\)", "\\1", media.summary$name))
    pie(media.summary$number, labels=media.summary$name, col=grey.palette(3), clockwise=T, init.angle=90)
  dev.off()
}

cummulative_plot.enzymatic_coverage = function()
{
  organisms = read.table("../data/organisms.tab", na.strings="", header=T, sep="\t", stringsAsFactors=F)
  top95clusters = read.table("../data/screenG_tax_input_specI_clusters.tab", na.strings="", header=T, sep="\t", stringsAsFactors=F)
  top95genus = read.table("../data/screenG_genus_order.tab", na.strings="", header=F, sep="\t", stringsAsFactors=F)$V1
  db = dbConnect(SQLite(), dbname="../data/kegg.db")
  all.ec = dbGetQuery(db, paste0("SELECT DISTINCT ec4 AS ec4 FROM org2ec INNER JOIN tax2org ON (tax2org.org=org2ec.org) WHERE taxid IN (", paste(top95clusters$taxid, collapse=","), ")"))
  
  top95ggplot = data.frame()
  top95genus.ec = c()
  for(genus in top95genus) {
    unique(organisms$species[organisms$genus=="Bilophila"])
    taxid.g = unique(organisms$taxid[!is.na(organisms$taxid) & organisms$genus==genus])
    orgs.g = unique(organisms$nearest_kegg_org[!is.na(organisms$nearest_kegg_org) & organisms$genus==genus])
    
    #  query = paste0("SELECT DISTINCT ec4 FROM org2ec INNER JOIN tax2org ON (tax2org.org=org2ec.org) WHERE tax2org.taxid IN ('", paste(taxid.g, collapse="','"), "') OR org2ec.org IN ('", paste(orgs.g, collapse="','"), "')")
    query = paste0("SELECT DISTINCT ec4 FROM org2ec WHERE org2ec.org IN ('", paste(orgs.g, collapse="','"), "')")
    genus.ec = dbGetQuery(db, query)
    
    top95genus.ec = unique(c(top95genus.ec, genus.ec$ec4))
    top95ggplot = rbind(top95ggplot, data.frame(genus=genus, coverage=mean(all.ec$ec4 %in% genus.ec$ec4), cummulative=mean(all.ec$ec4 %in% top95genus.ec)))
  }
  top95ggplot$genus_i = as.numeric(top95ggplot$genus) - 0.5
  
  pdf("../report/species_selection_enzyme_coverage.pdf", height=6, width=8)
  par(las=2, mar=c(10,5,5,5))
  plot(cummulative ~ genus, top95ggplot, type="s", ylim=c(0, 1), xlab="", col="#A4A4A4")
  lines(cummulative ~ genus_i, top95ggplot, type="s")
  dev.off()
}