source("utils/gen5.R")
library(gtools)
library(reshape2)
library(plyr)
library(gridExtra)
library(grid)

curves2.read = function(dir, section="578")
{
  layout.map = read.delim(file.path(dir, "layout.map"), stringsAsFactors=F)  
  layout.map$Species[layout.map$Species==""] = NA
  layout.map$ConditionSpecies[layout.map$ConditionSpecies==""] = NA
  
  res = data.frame()
  files = list.files(dir, pattern="^[^$~]+\\.xlsx?$", full.names=T)
  for(path in files)
  {
    writeLines(paste0(which(path == files), "/", length(files), ": ", basename(path)))
    filename = gsub("(\\.xpt|\\.xlsx)", "", basename(path))
    
    layout.map.i = subset(layout.map, File==filename, -File)
    res.i = gen5.parse(path)$long
    res.i$File = filename
    if(!all(layout.map.i$Annotation %in% res.i$Annotation)) {
      stop(paste0("Some cells in file '", path, "' don't have annotations in layout.map"))
    }
    res.i = merge(layout.map.i, res.i, by=c("Annotation"))
    
    if(any(is.na(res.i$File))) print(path)
    
    res.i.prune = res.i[,c("File", "Batch", "Section", "Well", "Row", "Col", "TechnicalReplicates", "Annotation", "Media", "ConditionSpecies", "Species", "Time", "OD", "Excluded")]
    res = rbind(res, res.i.prune)
  }
  
  if(is.character(section)) res = subset(res, Section==section)
  
  res = res[order(res$File, res$Media, res$Species, res$TechnicalReplicates, res$Time),]
  res
}

curves3.read = function(dir, section="578")
{
  layout.map = read.delim(file.path(dir, "layout.map"), stringsAsFactors=F)  
  layout.map$Species[layout.map$Species==""] = NA
  layout.map$ConditionSpecies[layout.map$ConditionSpecies==""] = NA
  
  
  files = list.files(dir, pattern="^[^$~]+\\.xlsx?$", full.names=T)
  
  files.unknown = setdiff(gsub("(\\.xpt|\\.xlsx)", "", basename(files)), unique(layout.map$File))
  if(length(files.unknown)) {
    stop(paste0("Some files are missing from layout map: ", paste(files.unknown, collapse=", ")))
  }
  
  
  res = data.frame()
  for(path in files)
  {
    writeLines(paste0(which(path == files), "/", length(files), ": ", basename(path)))
    filename = gsub("(\\.xpt|\\.xlsx)", "", basename(path))
    
    layout.map.i = subset(layout.map, File==filename, -File)
    res.i = gen5.parse(path)$long
    res.i$File = filename
    if(!all(layout.map.i$Annotation %in% res.i$Annotation)) {
      stop(paste0("Some cells in file '", path, "' don't have annotations in layout.map"))
    }
    res.i = merge(layout.map.i, res.i, by=c("Annotation"))
    
    res.i.prune = res.i[,c("File", "Batch", "Section", "Well", "Row", "Col", "TechnicalReplicates", "Annotation", "Media", "ConditionSpecies", "Species", "Time", "OD", "Excluded")]
    res = rbind(res, res.i.prune)
  }
  
  if(is.character(section)) res = subset(res, Section==section)
  
  res = res[order(res$File, res$Media, res$Species, res$TechnicalReplicates, res$Time),]
  res
}

curves46.read = function(dir, section="578")
{
  layout.map = read.delim(file.path(dir, "layout.map"), stringsAsFactors=F, na.strings="")  
  organisms.map1 = read.delim(file.path(dir, "organisms.map"), stringsAsFactors=F) 
  organisms.map = organisms.map1$species
  names(organisms.map) = organisms.map1$filename
  
  files = list.files(dir, pattern="^[^$~]+\\.xlsx?$", full.names=T)
  species_fnames = gsub("CM\\d_\\w_(.*)\\.xlsx", "\\1", unname(sapply(files, basename)))
  
  if(!all(species_fnames %in% names(organisms.map))) {
    stop(paste0("Not all species are in organism map ", paste(setdiff(unique(species_fnames), names(organisms.map)), collapse=",")))
  }
  
  res = data.frame()
  for(path in files)
  {
    writeLines(paste0(which(path == files), "/", length(files), ": ", basename(path)))
    filename = gsub("(\\.xpt|\\.xlsx)", "", basename(path))
    species = organisms.map[gsub("CM\\d_\\w_", "", filename)]
    batch = as.numeric(gsub("CM(\\d).*", "\\1", filename))
    layout = gsub("CM\\d_(\\w)_.*", "\\1", filename)
    
    layout.map.i = subset(layout.map, Batch==batch & Layout==layout, )
    res.i = gen5.parse(path)$long
    if(is.character(section)) res.i = subset(res.i, Section==section)
    
    res.i$File = filename
    if(!all(layout.map.i$Annotation %in% res.i$Annotation)) {
      stop(paste0("Some cells in file '", path, "' don't have annotations in layout.map"))
    }
    res.i = merge(layout.map.i, res.i, by=c("Annotation"))
    res.i$Species[!is.na(res.i$Media)] = species
    
    coi = c("File", "Batch", "Section", "Well", "Row", "Col", "TechnicalReplicates", "Annotation", "Media", "ConditionSpecies", "Species", "Time", "OD", "Excluded")
    res.i.prune = res.i[,coi]
    res = rbind(res, res.i.prune)
  }
  
  res = res[order(res$File, res$Media, res$Species, res$TechnicalReplicates, res$Time),]
  res
}

curves78.read = function(dir, section="578")
{
  section = "578"
  layout.map = read.delim(file.path(dir, "layout.map"), stringsAsFactors=F, na.strings="")
  species = na.omit(unique(layout.map$Species))
  known_species = c("A. muciniphila", "A. odontolyticus", "A. putredinis", "B. adolescentis", "B. animalis subsp. lactis BI-07", "B. animalis subsp. lactis BL-04", 
                    "B. crossotus", "B. fragilis", "B. fragilis enterotoxigenic", "B. hansenii", "B. longum subsp. infantis", "B. longum subsp. longum", "B. thetaiotaomicron", 
                    "B. uniformis", "B. vulgatus", "C. bolteae", "C. comes", "C. perfringens C36", "C. perfringens S107", "C. ramosum", "C. saccharolyticum", "D. piger", 
                    "E. coli CFT073", "E. coli E2348/69", "E. coli ED1a", "E. coli H10407", "E. coli HM605", "E. coli IAI1", "E. coli UTI89", "E. lenta", "E. rectale", 
                    "E. siraeum", "F. nucleatum subsp. animalis", "F. nucleatum subsp. nucleatum", "F. nucleatum subsp. vincentii", "H. parainfluenzae", "L. gasseri", 
                    "L. lactis", "L. paracasei", "L. plantarum", "P. capillosus", "P. copri", "P. difficile", "P. melaninogenica", "R. bromii", "R. gnavus", "R. torques", 
                    "S. flexneri", "S. salivarius", "S. sonnei", "S. typhimurium ATCC14028", "S. typhimurium LT2", "V. cholerae A1552", "V. cholerae N16961", "V. parvula", 
                    "Y. pseudotuberculosis", "B. wadsworthia", "B. dorei", "B. caccae", "B. eggerthii", "B. clarus", "B. coprocola", "B. ovatus", "B. stercoris", 
                    "B. vulgatus HM-720", "B. fragilis HM-20", "B. fragilis HM-709", "B. fragilis HM-710", "B. fragilis HM-711", "B. fragilis HM-714", "B. fragilis HM-712", 
                    "B. fragilis HM-713", "B. xylanisolvens", "B. uniformis HM-715", "B. uniformis HM-716", "P. capillosus", "R. intestinalis", "E. limosum", "B. obeum", 
                    "A. shahii", "P. merdae", "S. parasanguinis", "C. aerofaciens", "S. parasanguinis", "C. aerofaciens", "P. distasonis", "E. eligens", "D. formicigenerans", 
                    "R. hominis", "O. splanchnicus", "C. leptum", "L. salivarius", "L. ruminis", "L. sakei subsp. sakei", "L. delbrueckii subsp. delbrueckii", "C. catus", 
                    "L. fermentum", "L. vaginalis", "L. acidophilus")
  files = list.files(dir, pattern="^[^$~]+\\.xlsx?$", full.names=T)
  
  if(!all(species %in% known_species)) {
    stop(paste0("Not all species are in organism map ", paste(setdiff(unique(species), known_species), collapse=",")))
  }
  
  res = data.frame()
  for(path in files)
  {
    writeLines(paste0(which(path == files), "/", length(files), ": ", basename(path)))
    filename = gsub("(\\.xpt|\\.xlsx)", "", basename(path))
    
    res.i_bck = gen5.parse(path)$long
    res.i_bck$File = gsub("(\\.xpt|\\.xlsx)", "", basename(res.i_bck$File))
    res.i = merge(res.i_bck, layout.map, by=c("File", "TechnicalReplicates"))
    if(nrow(res.i) != nrow(res.i_bck)) {
      stop(paste("Something was wrong with file", filename))
    }
    if(is.character(section)) res.i = subset(res.i, Section==section)
    
    res.i$ConditionSpecies = NA
    res.i$Volume = 8
    res.i$Batch = gsub("(\\d\\d\\d\\d)(\\d\\d)(\\d\\d)", "\\1-\\2-\\3", gsub("_.*", "", filename))
    
    coi = c("File", "Volume", "Batch", "Section", "Well", "Row", "Col", "TechnicalReplicates", "Annotation", "Media", "ConditionSpecies", "Species", "Time", "OD", "Excluded")
    if(!all(coi %in% colnames(res.i))) {
      stop(paste("Columns ", setdiff(coi, colnames(res.i)), "are missing"))
    }
    res.i.prune = res.i[,coi]
    res = rbind(res, res.i.prune)
  }

  res = res[order(res$File, res$Media, res$Species, res$TechnicalReplicates, res$Time),]
}

curves9.read = function(dir, section="578")
{
  layout.map = read.delim(file.path(dir, "layout.map"), stringsAsFactors=F, na.strings="")
  species = na.omit(unique(layout.map$Species))
  known_species = c("A. muciniphila", "A. odontolyticus", "A. putredinis", "B. adolescentis", "B. animalis subsp. lactis BI-07", "B. animalis subsp. lactis BL-04", 
                    "B. crossotus", "B. fragilis", "B. fragilis enterotoxigenic", "B. hansenii", "B. longum subsp. infantis", "B. longum subsp. longum", "B. thetaiotaomicron", 
                    "B. uniformis", "B. vulgatus", "C. bolteae", "C. comes", "C. perfringens C36", "C. perfringens S107", "C. ramosum", "C. saccharolyticum", "D. piger", 
                    "E. coli CFT073", "E. coli E2348/69", "E. coli ED1a", "E. coli H10407", "E. coli HM605", "E. coli IAI1", "E. coli UTI89", "E. lenta", "E. rectale", 
                    "E. siraeum", "F. nucleatum subsp. animalis", "F. nucleatum subsp. nucleatum", "F. nucleatum subsp. vincentii", "H. parainfluenzae", "L. gasseri", 
                    "L. lactis", "L. paracasei", "L. plantarum", "P. capillosus", "P. copri", "P. difficile", "P. melaninogenica", "R. bromii", "R. gnavus", "R. torques", 
                    "S. flexneri", "S. salivarius", "S. sonnei", "S. typhimurium ATCC14028", "S. typhimurium LT2", "V. cholerae A1552", "V. cholerae N16961", "V. parvula", 
                    "Y. pseudotuberculosis", "B. wadsworthia", "B. dorei", "B. caccae", "B. eggerthii", "B. clarus", "B. coprocola", "B. ovatus", "B. stercoris", 
                    "B. vulgatus HM-720", "B. fragilis HM-20", "B. fragilis HM-709", "B. fragilis HM-710", "B. fragilis HM-711", "B. fragilis HM-714", "B. fragilis HM-712", 
                    "B. fragilis HM-713", "B. xylanisolvens", "B. uniformis HM-715", "B. uniformis HM-716", "P. capillosus", "R. intestinalis", "E. limosum", "B. obeum", 
                    "A. shahii", "P. merdae", "S. parasanguinis", "C. aerofaciens", "S. parasanguinis", "C. aerofaciens", "P. distasonis", "E. eligens", "D. formicigenerans", 
                    "R. hominis", "O. splanchnicus", "C. leptum", "L. salivarius", "L. ruminis", "L. sakei subsp. sakei", "L. delbrueckii subsp. delbrueckii", "C. catus", 
                    "L. fermentum", "L. vaginalis", "L. acidophilus")
  files = list.files(dir, pattern="^[^$~]+\\.xlsx?$", full.names=T)
  
  if(!all(species %in% known_species)) {
    stop(paste0("Not all species are in organism map ", paste(setdiff(unique(species), known_species), collapse=",")))
  }
  
  res = data.frame()
  for(path in files)
  {
    writeLines(paste0(which(path == files), "/", length(files), ": ", basename(path)))
    filename = gsub("(\\.xpt|\\.xlsx)", "", basename(path))
    
    res.i_bck = gen5.parse(path)$long
    res.i_bck$File = gsub("(\\.xpt|\\.xlsx)", "", basename(res.i_bck$File))
    res.i = merge(res.i_bck, layout.map, by=c("File", "TechnicalReplicates"))
    
    if(nrow(res.i) != nrow(res.i_bck)) {
      stop(paste("Something was wrong with file", filename))
    }
    if(is.character(section)) res.i = subset(res.i, Section==section)
    
    res.i$ConditionSpecies = NA
    res.i$Volume = 9
    res.i$Batch = gsub("(\\d\\d\\d\\d)(\\d\\d)(\\d\\d)", "\\1-\\2-\\3", gsub("_.*", "", filename))
    
    coi = c("File", "Passage", "Volume", "Batch", "Section", "Well", "Row", "Col", "TechnicalReplicates", "Annotation", "Media", "ConditionSpecies", "Species", "Time", "OD", "Excluded")
    if(!all(coi %in% colnames(res.i))) {
      stop(paste("Columns ", setdiff(coi, colnames(res.i)), "are missing"))
    }
    res.i.prune = res.i[,coi]
    res = rbind(res, res.i.prune)
  }
  
  res = res[order(res$File, res$Media, res$Species, res$TechnicalReplicates, res$Time),]
}




od = function() {
  ann.cols = c("File", "Passage", "Volume", "Batch", "Media", "TechnicalReplicates", "ConditionSpecies", "Species")
  curves9 = curves9.read("../data/curves_raw/10")

  curves9.a_auto = ddply(curves9, ann.cols, summarize,
                             BlankOD=min(OD),
                             BlankTime=Time[which.min(OD)],
                             Class="NoGrowthAuto",
                         
                             TopOD=max(OD),
                             MaxOD=max(OD)-min(OD), MaxTime=max(Time)/3600,
                             StatOD=max(OD)-min(OD), StatTime=max(Time)/3600,
                             Rate=0, RateODIntercept=min(OD),
                             OvergrowthOD=max(OD)-min(OD), OvergrowthTime=max(Time)/3600,
                             AUC8=MaxOD*8, AUC12=MaxOD*8, AUC18=MaxOD*18, AUC24=MaxOD*24)
  curves9.a_empty = subset(curves9.a_auto, is.na(Species) | is.na(Media))
  curves9.a_nogrowth = subset(curves9.a_auto, MaxOD < 0.15)
  curves9.a_candidates = subset(curves9.a_auto, MaxOD >= 0.15)

  assist.a = read.delim("../data/curves_annotation.tab", quote="", na.strings="", stringsAsFactors=F)
  assist.af = subset(assist.a, !grepl("NoGrowth", Class) & !is.na(Species) & !is.na(Media) & is.na(ConditionSpecies))
  assist = read.delim("../data/curves.tab", quote="", na.strings="", stringsAsFactors=F)
  assist.f = merge(assist, assist.af[,c("File", "TechnicalReplicates")])
  assist.af$OD = assist.af$MaxOD
  assist.af$Time = assist.af$MaxTime
  assist.af = merge(assist.af, curves9.a_candidates[,c("File", "Media", "Species")], by=c("Media", "Species"))
  assist.af$File = assist.af$File.y

  curves9.f = merge(subset(curves9, !is.na(Species)), curves9.a_candidates)
  curves9.f = curves9.f[with(curves9.f, order(File, Species, Col, Row, Time)),]

  res.mu = gen5.annotate_coordinates(d=subset(curves9, !is.na(Species)), assist = res.all.assist, other=assist.f, aggregate=T, start=1, mu=F)
  colnames(res.mu)[6:7] = c("Rate", "RateODIntercept") 

  
  write.table(assist.a, file="../data/curves_annotation_10.tab", sep="\t", quote=F, row.names=F, na="")
  assist.a = read.delim("../data/curves_annotation.tab", quote="", na.strings="", stringsAsFactors=F)
  assist.a = rbind(assist.a, res.all)
  
  curves9$TT = gsub("[0-9]", "", curves9$TechnicalReplicates)
  curves9 = curves9[order(curves9$File, curves9$Well, curves9$Time),]
  curves9$MS = paste(curves9$Media, "\n", curves9$Spe)
  ggplot(subset(curves9, Media %in% c(3,8,9, "8_salic025", "8_salic075", "9_salic025", "9_salic075"))) +
    geom_line(aes(x=Time/3600, y=OD, group=paste(File, Well))) +
    facet_grid(Media ~ Species + TT) + 
    theme_bw()

  ggplot(subset(curves9, Media == "salic025" & grepl("theta", Species))) +
    geom_line(aes(x=Time/3600, y=OD, group=paste(File, Well))) +
    facet_wrap(File~TechnicalReplicates)
  
    #
 
  # Annotate shape
  question="(x)dead, (unrep)roducible, (u)nfinished, (undef)ined, (d)iauxic, (n)oisy, (o)vergrowth, (g)ood, (r)aising, (f)alling, (z)ero"
  classes=c(g="Good", x="NoGrowth", o="Overgrowth", d="Diauxic", n="Noisy", r="Raising", z="Zero", f="Falling", u="Unfinished", unrep="Unreproducible", undef="Undefined")
  classes.perm = function(cl, n) {
    unlist(apply(permutations(length(cl), n, names(cl)), 1, function(x, ccl) {
      ret = list()
      ret[[paste(x, collapse="")]] = paste(sort(ccl[x]), collapse=";")
      ret
    }, ccl=cl))
  }
  classes.multi = c(classes, classes.perm(classes, 2), classes.perm(classes, 3), classes.perm(classes, 4), classes.perm(classes, 5))
  
  res.shape = gen5.annotate_boolean(subset(curves9, !is.na(Species)), aggregate=T, classes=classes.multi, show.options=F, question="", size=nrow(unique(subset(curves9, !is.na(Species))[,c("File", "TechnicalReplicates")])), start=1)
  colnames(res.shape)[5] = "Class"
  save(res.shape_ng, file="res.shape_ng.rda")
  #save(curves9.f, file="curves9.f.rda")
}

correlation = function() {
  data7.ann = read.delim("data/curves7_annotation.tab", quote="", na.strings="", stringsAsFactors=F)
  
  data8.blank = read.delim("data/curves8_res.blank.tab", quote="", na.strings="", stringsAsFactors=F)
  data8.blank$BlankTime = data8.blank$Time
  data8.blank$BlankOD = data8.blank$OD
  data8.max = read.delim("data/curves8_res.max.tab", quote="", na.strings="", stringsAsFactors=F)
  data8.max$MaxTime = data8.max$Time
  data8.max$MaxOD = data8.max$OD
  data8.shape = read.delim("data/curves8_res.shape.tab", quote="", na.strings="", stringsAsFactors=F)
  data8.shape$Class = data8.shape$cl
  data8.stat = read.delim("data/curves8_res.stat.tab", quote="", na.strings="", stringsAsFactors=F)
  data8.stat$StatTime = data8.stat$Time
  data8.stat$StatOD = data8.stat$OD

  data8.ann = merge(data8.blank, data8.max[,c("File", "TechnicalReplicates", "MaxTime", "MaxOD")])
  data8.ann = merge(data8.ann, data8.stat[,c("File", "TechnicalReplicates", "StatTime", "StatOD")])
  data8.ann = merge(data8.ann, data8.shape[,c("File", "TechnicalReplicates", "Class")])

  data78.ann = merge(data7.ann, data8.ann, by=c("Media", "Species"), suffixes=c(".7", ".8"))
  data78.ann$Label = paste0(data78.ann$Media, ": ", data78.ann$Species)
  data78.ann = data78.ann[,c("Media", "Species", "Label", setdiff(colnames(data78.ann), c("Media", "Species", "Label")))]
  
  data78.ann.changed = subset(data78.ann, grepl("NoGrowth", Class.7) != grepl("NoGrowth", Class.8))
  View(data78.ann.changed[order(data78.ann.changed$Species, data78.ann.changed$Media), c("Species", "Media", "Class.7", "Class.8")])
  table(grepl("NoGrowth", data78.ann.changed$Class.7), grepl("NoGrowth", data78.ann.changed$Class.8))
  
  #f = which(with(data8.ann, grepl("20160613", File) & Species %in% "E. siraeum" & !(Media %in% rich.media)))
  #data8.ann[f, c("Media", "Class")]


  par(mfrow=c(1,3))
  data78.ann_lm = lm(MaxOD.7 ~ MaxOD.8, data78.ann)
  data78.ann_cor = cor(data78.ann$MaxOD.8, data78.ann$MaxOD.7)
  plot(MaxOD.7 ~ MaxOD.8, data78.ann, main="all MaxOD")
  abline(data78.ann_lm)
  legend("topleft", legend=round(data78.ann_cor^2, 2))
  
  
  data78.ann.gr = subset(data78.ann, !grepl("NoGrowth", Class.7) & !grepl("NoGrowth", Class.8))
  data78.ann.gr_lm = lm(MaxOD.7 ~ MaxOD.8, data78.ann.gr)
  data78.ann.gr_cor = cor(data78.ann.gr$MaxOD.8, data78.ann.gr$MaxOD.7)
  plot(MaxOD.7 ~ MaxOD.8, data78.ann.gr, main="growing MaxOD")
  abline(data78.ann.gr_lm)
  legend("topleft", legend=round(data78.ann.gr_cor^2, 2))
  
  data78.ann.fin = subset(data78.ann, !grepl("NoGrowth|Unfin", Class.7) & !grepl("NoGrowth|Unfin", Class.8))
  data78.ann.fin_lm = lm(MaxOD.7 ~ MaxOD.8, data78.ann.fin)
  data78.ann.fin_cor = cor(data78.ann.fin$MaxOD.8, data78.ann.fin$MaxOD.7)
  plot(MaxOD.7 ~ MaxOD.8, data78.ann.fin, main="finished MaxOD")
  abline(data78.ann.fin_lm)
  legend("topleft", legend=round(data78.ann.fin_cor^2, 2))
  par(mfrow=c(1,1))
}

annotate = function()
{
  question="(x)dead, (unrep)roducible, (u)nfinished, (undef)ined, (d)iauxic, (n)oisy, (o)vergrowth, (g)ood, (r)aising, (f)alling, (z)ero"
  classes=c(g="Good", x="NoGrowth", o="Overgrowth", d="Diauxic", n="Noisy", r="Raising", z="Zero", f="Falling", u="Unfinished", unrep="Unreproducible", undef="Undefined")
  classes.perm = function(cl, n) {
    unlist(apply(permutations(length(cl), n, names(cl)), 1, function(x, ccl) {
      ret = list()
      ret[[paste(x, collapse="")]] = paste(sort(ccl[x]), collapse=";")
      ret
    }, ccl=cl))
  }
  classes.multi = c(classes, classes.perm(classes, 2), classes.perm(classes, 3), classes.perm(classes, 4), classes.perm(classes, 5))
  
  # Annotate shape
  curves7 = read.delim("data/curves.tab", header=T, sep="\t", stringsAsFactors=F, na="")
  curves7_f = subset(curves7, !is.na(Species) & Volume==8)
  curves7_f = merge(curves7_f, ddply(res.shape_1, .(Species, Media), summarize, cl=paste(unique(cl), collapse=";")), all.x=T)
  curves7_f.max = ddply(curves7_f, .(File, TechnicalReplicates), summarize, max=max(OD, na.rm=T)-min(OD, na.rm=T))
  curves7_f = merge(curves7_f, curves7_f.max)
  curves7_f = subset(curves7_f, max >= 0.15)
  curves7_f = curves7_f[order(curves7_f$File, curves7_f$Species, curves7_f$Col, curves7_f$Row, curves7_f$Time),]
  res.shape = gen5.annotate_boolean(curves7_f, aggregate=T, classes=classes.multi, show.options=F, question="", size=nrow(unique(curves7_f[,c("File", "TechnicalReplicates")])), start=1)
  write.table(wells[,c("File", "TechnicalReplicates", "Media", "ConditionSpecies", "Species")], file="data/aaa.txt", sep="\t", quote=F, row.names=F, na="")
  
  
  
  pdf("report/curves8_growth_summary.pdf", paper="a4r", height=8.27, width=11.69)
  res.shape_all = read.delim("data/curves8_res.shape.tab", quote="", na.strings="", stringsAsFactors=F)
  curves.a = read.delim("data/curves_annotation.tab", quote="", na.strings="", stringsAsFactors=F)
  curves.a_species = ddply(curves.a, .(Species), summarize, batch="old", good=sum(!grepl("NoGrowth|Unfinished", Class)), grow=sum(!grepl("NoGrowth", Class)), all=length(Class))
  summary_species = ddply(res.shape_all, .(Species), summarize, batch="new", good=sum(!grepl("NoGrowth|Unfinished", cl)), grow=sum(!grepl("NoGrowth", cl)), all=length(cl))
  summary_species = rbind(summary_species, curves.a_species[curves.a_species$Species %in% summary_species$Species,])
  curves.a_media = ddply(curves.a, .(Media), summarize, batch="old", good=sum(!grepl("NoGrowth|Unfinished", Class)), grow=sum(!grepl("NoGrowth", Class)), all=length(Class))
  summary_media = ddply(res.shape_all, .(Media), summarize, batch="new", good=sum(!grepl("NoGrowth|Unfinished", cl)), grow=sum(!grepl("NoGrowth", cl)), all=length(cl))
  summary_media = rbind(summary_media, curves.a_media[curves.a_media$Media %in% summary_media$Media,])
  summary_files = ddply(res.shape_all, .(File), summarize, batch="new", good=sum(!grepl("NoGrowth|Unfinished", cl)), grow=sum(!grepl("NoGrowth", cl)), all=length(cl))
  
  
  grid.arrange(
    ggplot(subset(summary_files, batch=="new")) + 
      geom_bar(aes(x=reorder(File, good), y=grow/all, fill="new grow", group=batch), stat="identity", position="dodge", width=.7) + 
      coord_flip() + theme_minimal(base_size=8) + xlab("") + guides(fill="none")
    ,
    ggplot(subset(summary_species, batch=="new")) + 
      geom_bar(aes(x=reorder(Species, good), y=grow/all, fill=paste(batch, "grow"), group=batch), stat="identity", position="dodge", width=.7) + 
      geom_bar(aes(x=reorder(Species, good), y=good/all, fill=paste(batch, "good"), group=batch), stat="identity", position="dodge", width=.7) + 
      coord_flip() + theme_minimal(base_size=8) + xlab("") + guides(fill="none")
    , 
    ggplot(subset(summary_media, batch=="new")) + 
      geom_bar(aes(x=reorder(Media, good), y=grow/all, fill=paste(batch, "grow"), group=batch), stat="identity", position="dodge", width=.7) + 
      geom_bar(aes(x=reorder(Media, good), y=good/all, fill=paste(batch, "good"), group=batch), stat="identity", position="dodge", width=.7) + 
      coord_flip() + theme_minimal(base_size=8) + xlab("")
    ,ncol=3, top=textGrob("8: Percentage of species grow and percentage of finished measurements (only new data)")
  )
  
  
  grid.arrange(
    ggplot(summary_files) + 
      geom_bar(aes(x=reorder(File, good), y=grow/all, fill="new grow", group=batch), stat="identity", position="dodge", width=.7) + 
      geom_bar(aes(x=reorder(File, good), y=good/all, fill="new good", group=batch), stat="identity", position="dodge", width=.7) + 
      coord_flip() + theme_minimal(base_size=8) + xlab("") + guides(fill="none")
    ,
    ggplot(summary_species) + 
      geom_bar(aes(x=reorder(Species, good), y=grow/all, fill=paste(batch, "grow"), group=batch), stat="identity", position="dodge", width=.7) + 
      geom_bar(aes(x=reorder(Species, good), y=good/all, fill=paste(batch, "good"), group=batch), stat="identity", position="dodge", width=.7) + 
      coord_flip() + theme_minimal(base_size=8) + xlab("") + guides(fill="none")
    , 
    ggplot(summary_media) + 
      geom_bar(aes(x=reorder(Media, good), y=grow/all, fill=paste(batch, "grow"), group=batch), stat="identity", position="dodge", width=.7) + 
      geom_bar(aes(x=reorder(Media, good), y=good/all, fill=paste(batch, "good"), group=batch), stat="identity", position="dodge", width=.7) + 
      coord_flip() + theme_minimal(base_size=8) + xlab("")
    ,ncol=3, top=textGrob("8: Percentage of species grow and percentage of finished measurements (new vs old data)")
  )
  
  
  growth_ggplot = subset(res.shape_all, , c(Media, Species, max, cl))
  colnames(growth_ggplot) = c("Media", "Species", "MaxOD", "Class")
  growth_ggplot$Batch = "New"
  growth_ggplot2 = subset(curves.a, Media %in% c(1:15, "BHI++", "GMM", "mGAM", "WCA"), c("Media", "Species", "MaxOD", "BlankOD", "Class"))
  growth_ggplot2$MaxOD = growth_ggplot2$MaxOD - growth_ggplot2$BlankOD
  growth_ggplot2$BlankOD = NULL
  growth_ggplot2$Batch = "Old"  
  growth_ggplot_sp = intersect(growth_ggplot2$Species, growth_ggplot$Species)
  growth_ggplot2 = subset(growth_ggplot2, Species %in% growth_ggplot_sp)
  growth_ggplot = subset(growth_ggplot, Species %in% growth_ggplot_sp)
  growth_ggplot_all = rbind(growth_ggplot, growth_ggplot2)
  growth_ggplot_all = subset(growth_ggplot_all, !grepl("NoGrowth", Class) & !grepl("Unfinished", Class))
  growth_ggplot_all = subset(growth_ggplot, Species %in% growth_ggplot_sp)

  growth_ggplot_all.wide = merge(growth_ggplot, growth_ggplot2, by=c("Media", "Species"), suffixes=c("_new", "_old"))
  growth_ggplot_all.wide_qualitative = ddply(growth_ggplot_all.wide, .(Media), summarize, 
    P=sum(grepl("NoGrowth", Class_old) == F), 
    N=sum(grepl("NoGrowth", Class_old) == T), 
    TP=sum((grepl("NoGrowth", Class_new) == F) & (grepl("NoGrowth", Class_old) == F)), 
    FP=sum((grepl("NoGrowth", Class_new) == F) & (grepl("NoGrowth", Class_old) == T)), 
    TN=sum((grepl("NoGrowth", Class_new) == T) & (grepl("NoGrowth", Class_old) == T)), 
    FN=sum((grepl("NoGrowth", Class_new) == T) & (grepl("NoGrowth", Class_old) == F)),
    sensitivity = TP/P, specificity=TN/N)
  growth_ggplot_all.wide_qualitative$Media = factor(growth_ggplot_all.wide_qualitative$Media, c(15:1, "WCA", "mGAM", "GMM", "BHI++"))
  growth_ggplot_all.wide_qualitative.specificity = growth_ggplot_all.wide_qualitative
  growth_ggplot_all.wide_qualitative.specificity$type="specificity"
  growth_ggplot_all.wide_qualitative.specificity$value = growth_ggplot_all.wide_qualitative.specificity$specificity
  growth_ggplot_all.wide_qualitative.sensitivity = growth_ggplot_all.wide_qualitative
  growth_ggplot_all.wide_qualitative.sensitivity$type="sensitivity"
  growth_ggplot_all.wide_qualitative.sensitivity$value = growth_ggplot_all.wide_qualitative.sensitivity$sensitivity

  growth_ggplot_all.wide_grow = subset(growth_ggplot_all.wide, !grepl("NoGrowth", Class_new) & !grepl("Unfinished", Class_new) & !grepl("NoGrowth", Class_old) & !grepl("Unfinished", Class_old))
  growth_ggplot_all.cor_grow = ggplot.cor_data(growth_ggplot_all.wide_grow, mapping=aes(MaxOD_new, MaxOD_old), facets=~Media, method="spearman", facet_fun="facet_wrap", scales="free")
  growth_ggplot_all.cor_grow.sp = ggplot.cor_data(growth_ggplot_all.wide_grow, mapping=aes(MaxOD_new, MaxOD_old), facets=~Species, method="spearman", facet_fun="facet_wrap", scales="free")
  
  grid.arrange(
    ggplot(rbind(growth_ggplot_all.wide_qualitative.specificity, growth_ggplot_all.wide_qualitative.sensitivity)) +
      geom_bar(aes(Media, value, fill=type), stat="identity", position="dodge") + 
      theme_bw() + 
      ggtitle("Qualitative growth replicate of old data") +
      theme(axis.title=element_blank()) + 
      coord_flip()
    ,
    ggplot(growth_ggplot_all.wide_grow, aes(MaxOD_new, MaxOD_old)) + 
      geom_point() + 
      geom_abline(intercept=0, slope=1,  size=0.5, alpha=0.8, color="#000000") + 
      geom_smooth(aes(), method="lm") +
      geom_text(aes(x=left, y=top-0.1, label=short_str), growth_ggplot_all.cor_grow, hjust=0, vjust=0, size=3, color="#D73027") +
      facet_wrap(~Media, scales="free") + 
      theme_minimal(base_size=16) + 
      theme(axis.ticks=element_blank(), axis.text=element_blank()) + 
      ggtitle("MaxOD correlation (old vs new)"),
    widths=c(1, 2.5)
  )
  
  ggplot(growth_ggplot_all.wide_grow, aes(MaxOD_new, MaxOD_old)) + 
    geom_point() + 
    geom_abline(intercept=0, slope=1,  size=0.5, alpha=0.8, color="#000000") + 
    geom_smooth(aes(), method="lm") +
    geom_text(aes(x=left, y=top-0.1, label=short_str), growth_ggplot_all.cor_grow.sp, hjust=0, vjust=0, size=3, color="#D73027") +
    facet_wrap(~Species, scales="free") + 
    theme_minimal(base_size=16) + 
    theme(axis.ticks=element_blank(), axis.text=element_blank()) + 
    ggtitle("MaxOD correlation (old vs new)")
  dev.off()
  
  
  curves.a = read.delim("data/curves_annotation.tab", quote="", na.strings="", stringsAsFactors=F)
  curves.a7 = read.delim("data/curves7_annotation.tab", quote="", na.strings="", stringsAsFactors=F)
  curves.a = rbind(curves.a, curves.a7)
  res.shape = read.delim("data/curves8_res.shape.tab", quote="", na.strings="", stringsAsFactors=F)
  
  curves = read.delim("data/curves.tab", header=T, sep="\t", stringsAsFactors=F, na="")
  curves8 = subset(curves, !is.na(Species) & Volume==8)
  curves8 = merge(curves8, res.shape[,c("File", "TechnicalReplicates", "cl")])
  curves8_f = subset(curves8, !grepl("NoGrowth", cl))
  
  # Annotate max OD
  res.max.a = read.delim("data/curves7_annotation.tab", quote="", na.strings="", stringsAsFactors=F)
  res.max.a$OD = res.max.a$MaxOD
  res.max.a$Time = res.max.a$MaxTime
  res.max = gen5.annotate_coordinates(d=curves8_f, other=curves, assist=res.max.a, aggregate=T, start=1)
  res.max.ng = ddply(subset(curves8, grepl("NoGrowth", cl)), .(File, Media, ConditionSpecies, Species, TechnicalReplicates), summarize, OD=min(OD), Time=max(Time))
  write.table(rbind(res.max, res.max.ng), file="_data/curves8_res.blank.tab", sep="\t", quote=F, row.names=F, na="")
  
  # Annotate stationary phase OD
  res.stat2 = gen5.annotate_coordinates(curves7_f, aggregate=T, start=1, assist=res.stat)
  save(res.stat2, file="tmp/curves7/res.stat2.rda")
  write.table(res.stat2, file="tmp/curves2_stat_Sergej.tab", sep="\t", quote=F, row.names=F, na="")
    
  # Annotate growth rate
  res.rate = gen5.annotate_coordinates(d=curves8_f, aggregate=T, start=1, mu=T)
  res.rate.ng = ddply(subset(curves8, grepl("NoGrowth", cl)), .(File, Media, ConditionSpecies, Species, TechnicalReplicates), summarize, Rate=0, RateODIntercept=max(diff(sort(OD))))
  write.table(rbind(res.rate, res.rate.ng), file="data/curves8_res.rate.tab", sep="\t", quote=F, row.names=F, na="")
  
  # Annotate blank
  res.blank = gen5.annotate_coordinates(curves7_f, aggregate=T, start=1)
  save(res.blank, file="tmp/curves7/res.blank1.rda")
  write.table(res.blank, file="tmp/curves2_blank.tab", sep="\t", quote=F, row.names=F, na="")
  
  # Annotate max
  res.max = gen5.annotate_coordinates(curves7_f, aggregate=T, start=1)
  save(res.max, file="tmp/curves7/res.max.rda")
  write.table(res.max, file="tmp/curves2_max.tab", sep="\t", quote=F, row.names=F, na="")
  
  # Annotate overgrowth
  res.over = gen5.annotate_coordinates(curves7_f, aggregate=T, start=1)
  save(res.over, file="tmp/curves7/res.over.rda")
  write.table(res.over.rda, file="tmp/res.over.rda.tab", sep="\t", quote=F, row.names=F, na="")
  
  load("X:/Sergej/ConditionedMedia/tmp/final2/res.fit3_a.rda")
  table(res.fit$File == "GMM_1_C_ramosum")
  
  # Annotate lag
  res.lag = gen5.annotate_coordinates(curves7_f, aggregate=T, start=1)
  save(res.lag, file="tmp/curves7/res.lag.rda")
  
  res.max = read.delim("data/curves8_res.max.tab", quote="", na.strings="", stringsAsFactors=F)
  res.blank = read.delim("data/curves8_res.blank.tab", quote="", na.strings="", stringsAsFactors=F)
  res.over = read.delim("data/curves8_res.max.tab", quote="", na.strings="", stringsAsFactors=F)
  res.rate = read.delim("data/curves8_res.rate.tab", quote="", na.strings="", stringsAsFactors=F)
  res.shape = read.delim("data/curves8_res.shape.tab", quote="", na.strings="", stringsAsFactors=F)
  res.stat = read.delim("data/curves8_res.stat.tab", quote="", na.strings="", stringsAsFactors=F)
  
  res.control = ddply(subset(curves, !is.na(Species)), .(File, Media, ConditionSpecies, Species, TechnicalReplicates), summarize,
                      Volume=Batch[1], 
                      Class=ifelse(all(is.na(Species)), "NoGrowthAuto;Control", ifelse(mean(max(OD, na.rm=T) - min(OD, na.rm=T)) < 0.15, "NoGrowthAuto", "")),
                      Rate=0, RateODIntercept=0, BlankOD=min(OD, na.rm=T),
                      StatOD=max(OD, na.rm=T), StatTime=Time[which.max(OD)],
                      MaxOD=max(OD, na.rm=T), MaxTime=Time[which.max(OD)],
                      OvergrowthOD=max(OD, na.rm=T), OvergrowthTime=Time[which.max(OD)])
  res.control = subset(res.control, grepl("NoGrowthAuto", Class))
  
 
  res.batch_c = unique(curves8[,c("File", "TechnicalReplicates", "Batch")])
  colnames(res.batch_c)[colnames(res.batch_c) == "Batch"] = "Volume"
  
  res.shape_c = res.shape
  colnames(res.shape_c)[colnames(res.shape_c) == "cl"] = "Class"
  
  res.rate_c = res.rate
  #colnames(res.rate_c)[colnames(res.rate_c) == "Mu"] = "Rate"
  #colnames(res.rate_c)[colnames(res.rate_c) == "ODIntercept"] = "RateODIntercept"
  
  res.blank_c = res.blank
  colnames(res.blank_c)[colnames(res.blank_c) == "OD"] = "BlankOD"
  
  res.stat_c = res.stat
  colnames(res.stat_c)[colnames(res.stat_c) == "OD"] = "StatOD"
  colnames(res.stat_c)[colnames(res.stat_c) == "Time"] = "StatTime"
  
  res.max_c = res.max
  colnames(res.max_c)[colnames(res.max_c) == "OD"] = "MaxOD"
  colnames(res.max_c)[colnames(res.max_c) == "Time"] = "MaxTime"
  
  res.over_c = res.over
  colnames(res.over_c)[colnames(res.over_c) == "OD"] = "OvergrowthOD"
  colnames(res.over_c)[colnames(res.over_c) == "Time"] = "OvergrowthTime"
  
  dim(res.final)
  res.final = res.rate_c
  res.final = merge(res.final, res.batch_c[,c("File", "TechnicalReplicates", "Volume")])
  res.final = merge(res.final, res.shape_c[,c("File", "TechnicalReplicates", "Class")])
  res.final = merge(res.final, res.blank_c[,c("File", "TechnicalReplicates", "BlankOD")])
  res.final = merge(res.final, res.stat_c[,c("File", "TechnicalReplicates", "StatOD", "StatTime")])
  res.final = merge(res.final, res.max_c[,c("File", "TechnicalReplicates", "MaxOD", "MaxTime")])
  res.final = merge(res.final, res.over_c[,c("File", "TechnicalReplicates", "OvergrowthOD", "OvergrowthTime")])
  # res.final = rbind(res.final, res.control)
  write.table(res.final, file="data/curves8_annotation.tab", sep="\t", quote=F, row.names=F, na="")
  
  
  res.1 = read.delim("data/curves_annotation_2-6.tab", quote="", na.strings="", stringsAsFactors=F)
  res.7 = read.delim("data/old/curves7_annotation.tab", quote="", na.strings="", stringsAsFactors=F)
  res.8 = read.delim("data/old/curves8_annotation.tab", quote="", na.strings="", stringsAsFactors=F)
  
  res = rbind(res.1, res.7[,colnames(res.1)], res.8[,colnames(res.1)])
  write.table(res, file="data/curves_annotation.tab", sep="\t", quote=F, row.names=F, na="")
  
}

currate = function()
{
  load("tmp/curves.rda")
  curves.a = read.table("data/curves_annotation.tab", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="")
  
  var = "Stat"
  if(var == "Rate") 
  {
    var.od = "Rate"
    var.time = "RateODIntercept"
  } else {
    var.od = paste0(var, "OD")
    var.time = paste0(var, "Time")
    var = var.od
  }
  
  assist = read.table("C:\\Users\\Sergej\\Desktop\\new  1.txt", sep="\t", quote="", header=T, stringsAsFactors=F, na.strings="")
  assist = unique(subset(assist, Var==var, -Var))
  assist = merge(curves.a, assist)
  assist$OD.max = assist$MaxOD
  assist$Time.max = assist$MaxTime
  assist$OD.stat = assist$StatOD
  assist$Time.stat = assist$StatTime
  assist$Mu.rate = assist$Rate
  assist$ODIntercept.rate = assist$RateODIntercept
  
  curves.f = merge(curves, assist[,c("File", "TechnicalReplicates")])
  annotation = gen5.annotate_coordinates(curves.f, aggregate=T, start=1, mu=(var=="Rate"), assist=assist)
  
  curves.a.f = match(paste0(annotation$File, annotation$TechnicalReplicates), paste0(curves.a$File, curves.a$TechnicalReplicates))
  if(var == "Rate") 
  {
    curves.a[[var.od]][curves.a.f] = annotation$Mu
    curves.a[[var.time]][curves.a.f] = annotation$ODIntercept
  } else {
    curves.a[[var.od]][curves.a.f] = annotation$OD
    curves.a[[var.time]][curves.a.f] = annotation$Time
  }
  
  write.table(curves.a, file="data/curves_annotation.tab", sep="\t", quote=F, row.names=F, na="")  
}


curves4_6.curate_shape_annotation = function()
{
  #curves46 = curves46.read(dir="X:/Sergej/Curves/data/curves4-6")
  #save(curves46, file="tmp/curves46/curves46.rda")
  load("tmp/curves46/curves46.rda")

  
  
  
  curves_shape = read.table("data/curves3_shape_all.tab", stringsAsFactors=F, na.strings="", sep="\t", quote="", header=T)
  curves_shape = merge(curves_shape, unique(curves3[,c("File", "TechnicalReplicates", "Media")]), by=c("File", "TechnicalReplicates"))
  
  curves_shape$Class = unlist(lapply(strsplit(curves_shape$Class, ";"), function(x) paste(sort(x), collapse=";")))
  curves_shape$Class = factor(curves_shape$Class, names(table(curves_shape$Class))[order(table(curves_shape$Class), decreasing=T)])
  curves_shape.sum = ddply(curves_shape, .(File, Media, ConditionSpecies, Species, TechnicalReplicates), function(z) {
    z.res = data.frame(
      cl.prop = max(table(z$Class))/length(na.omit(z$Class)),
      cl.s1   = ifelse(any(z$Annotator=="Sergej" & z$Batch==1), as.character(z$Class[z$Annotator=="Sergej" & z$Batch==1]), NA),
      cl.s2   = ifelse(any(z$Annotator=="Sergej" & z$Batch==2), as.character(z$Class[z$Annotator=="Sergej" & z$Batch==2]), NA),
      cl.m1   = ifelse(any(z$Annotator=="Melanie" & z$Batch==1), as.character(z$Class[z$Annotator=="Melanie" & z$Batch==1]), NA),
      cl.m2   = ifelse(any(z$Annotator=="Melanie" & z$Batch==2), as.character(z$Class[z$Annotator=="Melanie" & z$Batch==2]), NA),
      cl.s3   = ifelse(any(z$Annotator=="Sergej" & z$Batch==3), as.character(z$Class[z$Annotator=="Sergej" & z$Batch==3]), NA),
      stringsAsFactors=F
    )
    
    if(!is.na(z.res$Class.s3)) { z.res$Class = z.res$Class.s3
    } else {
      if(z.res$Class.m1==z.res$Class.s1) { z.res$Class = z.res$Class.m1
      } else {
        if(!is.na(z.res$Class.m2) & !is.na(z.res$Class.s2) & z.res$Class.m2==z.res$Class.s2) { z.res$Class = z.res$Class.m2
        } else {
          if(!is.na(z.res$Class.m2) & z.res$Class.m2==z.res$Class.s1) { z.res$Class = z.res$Class.m2
          } else {
            if(!is.na(z.res$Class.s2) & z.res$Class.m1==z.res$Class.s2) { z.res$Class = z.res$Class.m1
            } else { z.res$Class = NA }
          }}}}
    
    z.res
  })
  
  #write.table(curves_shape.sum[,c("File", "Media", "ConditionSpecies", "Species", "TechnicalReplicates", "cl")], file="data/curves3_shape.tab", sep="\t", quote=F, na="", row.names=F, col.names=T)
  
  pdf("report/curves3/shape_reproducibility_detailed.pdf", paper="a4r", height=11.69, width=8.27)
  for(m in c("All", unique(curves_shape.sum$Media)))
  {
    print(m)
    plot.new()
    text(.5,.5, paste("Media", m), cex=4)
    
    curves.a.f = subset(curves_shape.sum, (m == "All" | Media == m) & is.na(ConditionSpecies))
    
    #####################################################################
    # Plot annotation reproducibility
    #####################################################################
    if(annotation.reproducibility)
    {
      d1.cols = c("File", "ConditionSpecies", "Species", "TechnicalReplicates")
      curves_shape.s = merge(subset(curves_shape, Annotator=="Sergej" & Batch==1), subset(curves_shape, Annotator=="Sergej" & Batch==2), by=d1.cols, suffixes=c(".s1", ".s2"))
      curves_shape.s.rep = table(Melanie=curves_shape.s$Class.s1, Sergej=curves_shape.s$Class.s2)
      heatmap.2(curves_shape.s.rep, cellnote=ifelse(curves_shape.s.rep>1,curves_shape.s.rep,""),
                Rowv=F, Colv=F, dendrogram="none", 
                key=F, trace="none", 
                notecol="#FFFFFF", col=greenred(11), breaks=c(-1, 0, 2^(0:9)), 
                margins=c(10,10), cexRow=.7, cexCol=.7, notecex=0.5,
                xlab="Sergej 1", ylab="Sergej 2", main=paste0("Raw growth type annotation reproducibility (Sergej B1/2 - ", m, ")"))
      
      curves_shape.m = merge(subset(curves_shape, Annotator=="Melanie" & Batch==1), subset(curves_shape, Annotator=="Melanie" & Batch==2), by=d1.cols, suffixes=c(".m1", ".m2"))
      curves_shape.m.rep = table(Melanie=curves_shape.m$Class.m1, Sergej=curves_shape.m$Class.m2)
      heatmap.2(curves_shape.m.rep, cellnote=ifelse(curves_shape.m.rep>1,curves_shape.m.rep,""),
                Rowv=F, Colv=F, dendrogram="none", 
                key=F, trace="none", 
                notecol="#FFFFFF", col=greenred(11), breaks=c(-1, 0, 2^(0:9)), 
                margins=c(10,10), cexRow=.7, cexCol=.7, notecex=0.5,
                xlab="Melanie 1", ylab="Melanie 2", main=paste0("Raw growth type annotation reproducibility (Melanie B1/2 - ", m, ")"))
      
      curves_shape.sm1 = merge(subset(curves_shape, Annotator=="Sergej" & Batch==1), subset(curves_shape, Annotator=="Melanie" & Batch==1), by=d1.cols, suffixes=c(".s1", ".m1"))
      curves_shape.sm1.rep = table(Melanie=curves_shape.sm1$Class.s1, Sergej=curves_shape.sm1$Class.m1)
      heatmap.2(curves_shape.sm1.rep, cellnote=ifelse(curves_shape.sm1.rep>1,curves_shape.sm1.rep,""),
                Rowv=F, Colv=F, dendrogram="none", 
                key=F, trace="none", 
                notecol="#FFFFFF", col=greenred(11), breaks=c(-1, 0, 2^(0:9)), 
                margins=c(10,10), cexRow=.7, cexCol=.7, notecex=0.5,
                xlab="Sergej 1", ylab="Melanie 1", main=paste0("Raw growth type annotation reproducibility (Sergej/Melanie - ", m, ")"))
      
      curves_shape.sm2 = merge(subset(curves_shape, Annotator=="Sergej" & Batch==2), subset(curves_shape, Annotator=="Melanie" & Batch==2), by=d1.cols, suffixes=c(".s2", ".m2"))
      curves_shape.sm2.rep = table(Melanie=curves_shape.sm2$Class.s2, Sergej=curves_shape.sm2$Class.m2)
      heatmap.2(curves_shape.sm2.rep, cellnote=ifelse(curves_shape.sm2.rep>1,curves_shape.sm2.rep,""),
                Rowv=F, Colv=F, dendrogram="none", 
                key=F, trace="none", 
                notecol="#FFFFFF", col=greenred(11), breaks=c(-1, 0, 2^(0:9)), 
                margins=c(10,10), cexRow=.7, cexCol=.7, notecex=0.5,
                xlab="Sergej 2", ylab="Melanie 2", main=paste0("Raw growth type annotation reproducibility (Sergej/Melanie - ", m, ")"))
    }
  }
  dev.off()
}


curves3.curate_shape_annotation = function()
{
  #curves3 = curves3.read(dir="data/curves3")
  #save(curves3, file="tmp/curves3/curves3.rda")
  load("tmp/curves3/curves3.rda")
  curves_shape = read.table("data/curves3_shape_all.tab", stringsAsFactors=F, na.strings="", sep="\t", quote="", header=T)
  curves_shape = merge(curves_shape, unique(curves3[,c("File", "TechnicalReplicates", "Media")]), by=c("File", "TechnicalReplicates"))
  
  curves_shape$Class = unlist(lapply(strsplit(curves_shape$Class, ";"), function(x) paste(sort(x), collapse=";")))
  curves_shape$Class = factor(curves_shape$Class, names(table(curves_shape$Class))[order(table(curves_shape$Class), decreasing=T)])
  curves_shape.sum = ddply(curves_shape, .(File, Media, ConditionSpecies, Species, TechnicalReplicates), function(z) {
    z.res = data.frame(
      cl.prop = max(table(z$Class))/length(na.omit(z$Class)),
      cl.s1   = ifelse(any(z$Annotator=="Sergej" & z$Batch==1), as.character(z$Class[z$Annotator=="Sergej" & z$Batch==1]), NA),
      cl.s2   = ifelse(any(z$Annotator=="Sergej" & z$Batch==2), as.character(z$Class[z$Annotator=="Sergej" & z$Batch==2]), NA),
      cl.m1   = ifelse(any(z$Annotator=="Melanie" & z$Batch==1), as.character(z$Class[z$Annotator=="Melanie" & z$Batch==1]), NA),
      cl.m2   = ifelse(any(z$Annotator=="Melanie" & z$Batch==2), as.character(z$Class[z$Annotator=="Melanie" & z$Batch==2]), NA),
      cl.s3   = ifelse(any(z$Annotator=="Sergej" & z$Batch==3), as.character(z$Class[z$Annotator=="Sergej" & z$Batch==3]), NA),
      stringsAsFactors=F
    )
    
    if(!is.na(z.res$Class.s3)) { z.res$Class = z.res$Class.s3
    } else {
      if(z.res$Class.m1==z.res$Class.s1) { z.res$Class = z.res$Class.m1
      } else {
        if(!is.na(z.res$Class.m2) & !is.na(z.res$Class.s2) & z.res$Class.m2==z.res$Class.s2) { z.res$Class = z.res$Class.m2
        } else {
          if(!is.na(z.res$Class.m2) & z.res$Class.m2==z.res$Class.s1) { z.res$Class = z.res$Class.m2
          } else {
            if(!is.na(z.res$Class.s2) & z.res$Class.m1==z.res$Class.s2) { z.res$Class = z.res$Class.m1
            } else { z.res$Class = NA }
          }}}}
    
    z.res
  })
  
  #write.table(curves_shape.sum[,c("File", "Media", "ConditionSpecies", "Species", "TechnicalReplicates", "cl")], file="data/curves3_shape.tab", sep="\t", quote=F, na="", row.names=F, col.names=T)
  
  pdf("report/curves3/shape_reproducibility_detailed.pdf", paper="a4r", height=11.69, width=8.27)
  for(m in c("All", unique(curves_shape.sum$Media)))
  {
    print(m)
    plot.new()
    text(.5,.5, paste("Media", m), cex=4)
    
    curves.a.f = subset(curves_shape.sum, (m == "All" | Media == m) & is.na(ConditionSpecies))
    
    #####################################################################
    # Plot annotation reproducibility
    #####################################################################
    if(annotation.reproducibility)
    {
      d1.cols = c("File", "ConditionSpecies", "Species", "TechnicalReplicates")
      curves_shape.s = merge(subset(curves_shape, Annotator=="Sergej" & Batch==1), subset(curves_shape, Annotator=="Sergej" & Batch==2), by=d1.cols, suffixes=c(".s1", ".s2"))
      curves_shape.s.rep = table(Melanie=curves_shape.s$Class.s1, Sergej=curves_shape.s$Class.s2)
      heatmap.2(curves_shape.s.rep, cellnote=ifelse(curves_shape.s.rep>1,curves_shape.s.rep,""),
                Rowv=F, Colv=F, dendrogram="none", 
                key=F, trace="none", 
                notecol="#FFFFFF", col=greenred(11), breaks=c(-1, 0, 2^(0:9)), 
                margins=c(10,10), cexRow=.7, cexCol=.7, notecex=0.5,
                xlab="Sergej 1", ylab="Sergej 2", main=paste0("Raw growth type annotation reproducibility (Sergej B1/2 - ", m, ")"))
      
      curves_shape.m = merge(subset(curves_shape, Annotator=="Melanie" & Batch==1), subset(curves_shape, Annotator=="Melanie" & Batch==2), by=d1.cols, suffixes=c(".m1", ".m2"))
      curves_shape.m.rep = table(Melanie=curves_shape.m$Class.m1, Sergej=curves_shape.m$Class.m2)
      heatmap.2(curves_shape.m.rep, cellnote=ifelse(curves_shape.m.rep>1,curves_shape.m.rep,""),
                Rowv=F, Colv=F, dendrogram="none", 
                key=F, trace="none", 
                notecol="#FFFFFF", col=greenred(11), breaks=c(-1, 0, 2^(0:9)), 
                margins=c(10,10), cexRow=.7, cexCol=.7, notecex=0.5,
                xlab="Melanie 1", ylab="Melanie 2", main=paste0("Raw growth type annotation reproducibility (Melanie B1/2 - ", m, ")"))
      
      curves_shape.sm1 = merge(subset(curves_shape, Annotator=="Sergej" & Batch==1), subset(curves_shape, Annotator=="Melanie" & Batch==1), by=d1.cols, suffixes=c(".s1", ".m1"))
      curves_shape.sm1.rep = table(Melanie=curves_shape.sm1$Class.s1, Sergej=curves_shape.sm1$Class.m1)
      heatmap.2(curves_shape.sm1.rep, cellnote=ifelse(curves_shape.sm1.rep>1,curves_shape.sm1.rep,""),
                Rowv=F, Colv=F, dendrogram="none", 
                key=F, trace="none", 
                notecol="#FFFFFF", col=greenred(11), breaks=c(-1, 0, 2^(0:9)), 
                margins=c(10,10), cexRow=.7, cexCol=.7, notecex=0.5,
                xlab="Sergej 1", ylab="Melanie 1", main=paste0("Raw growth type annotation reproducibility (Sergej/Melanie - ", m, ")"))
      
      curves_shape.sm2 = merge(subset(curves_shape, Annotator=="Sergej" & Batch==2), subset(curves_shape, Annotator=="Melanie" & Batch==2), by=d1.cols, suffixes=c(".s2", ".m2"))
      curves_shape.sm2.rep = table(Melanie=curves_shape.sm2$Class.s2, Sergej=curves_shape.sm2$Class.m2)
      heatmap.2(curves_shape.sm2.rep, cellnote=ifelse(curves_shape.sm2.rep>1,curves_shape.sm2.rep,""),
                Rowv=F, Colv=F, dendrogram="none", 
                key=F, trace="none", 
                notecol="#FFFFFF", col=greenred(11), breaks=c(-1, 0, 2^(0:9)), 
                margins=c(10,10), cexRow=.7, cexCol=.7, notecex=0.5,
                xlab="Sergej 2", ylab="Melanie 2", main=paste0("Raw growth type annotation reproducibility (Sergej/Melanie - ", m, ")"))
    }
  }
  dev.off()
}

currate.media_composition = function()
{
  ##################################################################
  #
  # Correlate Media composition
  #
  ##################################################################
  medias = read.delim("tmp/media.txt", quote="", na.strings="", stringsAsFactors=F)
  medias[is.na(medias)] = 0
  medias.w = data.frame(t(dcast2(medias, compounds ~ media, value.var="maxflux", fun.aggregate=function(z) z[1])))
  #plot(medias.w, col=ifelse(medias.w < medias.w$m13, "red", "black"))
  #abline(a=0, b=1)
  
  m2 = "m3"
  par(mfrow=c(2,3))
  for(m1 in setdiff(unique(medias$media), m2))
  {
    medias.lhs = subset(medias, media==m1, c(compounds, maxflux))
    medias.rhs = subset(medias, media==m2, c(compounds, maxflux))
    medias.p = merge(medias.lhs, medias.rhs, by="compounds", suffixes=paste0(".", c(m1, m2)), all=T)
    medias.p[is.na(medias.p)] = 0
    medias.p = medias.p[which(medias.p[[paste0("maxflux.", m1)]] < 8 & medias.p[[paste0("maxflux.", m2)]] < 8),]
    medias.f = medias.p[which(medias.p[[paste0("maxflux.", m1)]] < medias.p[[paste0("maxflux.", m2)]]),]
    plot(x=sqrt(medias.p[[paste0("maxflux.", m1)]]), y=sqrt(medias.p[[paste0("maxflux.", m2)]]), xlab=m1, ylab=m2, col=ifelse(medias.p$compounds %in% medias.f$compounds, "red", "black"), main=paste(m1, "/", m2))
    text(x=sqrt(medias.f[[paste0("maxflux.", m1)]]), y=sqrt(medias.f[[paste0("maxflux.", m2)]]), label=medias.f$compounds, pos=4)
    abline(a=0, b=1)
  }
  par(mfrow=c(1,1))
}