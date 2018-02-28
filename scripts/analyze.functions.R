library(plyr)
library(ggplot2)
library(grid)

dist.species = function(x) {
  x.rownames = gsub("([^: ]+).*", "\\1", rownames(x))
  x_rich.dist = dist(t(x[x.rownames %in% media.rich,]))
  x_def.dist = dist(t(x[(x.rownames %in% media.general) & !(x.rownames %in% media.rich),]))
  x.dist = (x_rich.dist + x_def.dist)/2
  return(x.dist)
}

cuves.merge_annotations = function(curves.a.f) 
{
  growth_matrix.q1 = ddply(curves.a.f, .(Media, Species), function(z) {
    z.f = subset(z, !grepl("NoGrowth", Class) & !grepl("Unreproducible", Class) & !grepl("Unfinished", Class))
    z.blank = mean(z$BlankOD, na.rm=T)
    n.finished = sum(with(z, !grepl("NoGrowth", Class) & !grepl("Unreproducible", Class) & !grepl("Unfinished", Class)))
    n.growing = sum(with(z, !grepl("NoGrowth", Class) & !grepl("Unreproducible", Class)))
    n.all = nrow(z)
    
    ascii = function(x) { strtoi(charToRaw(x),16L) }
    sym.l1 = sapply(gsub("^([A-Z]).*", "\\1", z$TechnicalReplicates), ascii)
    sym.l2 = sapply(gsub(".*([A-Z])[0-9]+$", "\\1", z$TechnicalReplicates), ascii)
    sym.n1 = as.numeric(gsub("^[A-Z]([0-9]+).*", "\\1", z$TechnicalReplicates))
    sym.n2 = as.numeric(gsub(".*[A-Z]([0-9]+)$", "\\1", z$TechnicalReplicates))
    n.technical = pmax(sym.l2-sym.l1, sym.n2 - sym.n1) + 1
    
    if(nrow(z.f)) {
      r = data.frame(MaxOD=mean(z.f$MaxOD, na.rm=T), MaxOD_sd=sd(z.f$MaxOD, na.rm=T), Rate=median(z.f$Rate, na.rm=T), Rate_sd=sd(z.f$Rate, na.rm=T))
      if(any(grepl("AUC18", colnames(z.f)))) {
        r$AUC18 = mean(z.f$AUC18, na.rm=T)
        r$AUC18_sd = sd(z.f$AUC18, na.rm=T)
      }
    } else {
      z.f = subset(z, !grepl("NoGrowth", Class) & !grepl("Unreproducible", Class))
      if(nrow(z.f)) {
        r = data.frame(MaxOD=mean(z.f$MaxOD, na.rm=T), MaxOD_sd=sd(z.f$MaxOD, na.rm=T), Rate=mean(z.f$Rate, na.rm=T), Rate_sd=sd(z.f$Rate, na.rm=T))
        if(any(grepl("AUC18", colnames(z.f)))) {
          r$AUC18 = mean(z.f$AUC18, na.rm=T)
          r$AUC18_sd = sd(z.f$AUC18, na.rm=T)
        }
      } else {
        r = data.frame(MaxOD=0, MaxOD_sd=NA, Rate=0, Rate_sd=NA)
        if(any(grepl("AUC18", colnames(z.f)))) {
          r$AUC18 = 0
          r$AUC18_sd = NA
        }
      }
    }
    
    r$Count = ""
    if(n.finished > 0 & n.finished < n.all) { 
      r$Count = paste0(n.finished, "/", n.all) 
    } else {
      if(n.growing > 0 & n.finished < n.all) { 
        r$Count = paste0(":", ifelse(n.growing == n.all, n.all, paste0(n.growing, "/", n.all)))
      }
    }
    
    r$Growing = n.growing
    r$Finished = n.finished
    r$Replicates = n.all
    r$TechnicalReplicates = sum(n.technical)
    r$MaxOD = r$MaxOD - z.blank
    if(any(r$MaxOD < 0)) r$MaxOD[r$MaxOD < 0] = 0
    if(any(grepl("AUC18", colnames(z.f))) & any(r$AUC18 < 0)) {
      r$AUC18[r$AUC18 < 0] = 0
    }
    
    return(r)
  })
  
  growth_matrix.q1
}

cuves.merge_annotations2 = function(curves.a.f, nogrowth=F) 
{
  growth_matrix.q1 = ddply(curves.a.f, .(Media, Species), function(z) {
    #z = subset(curves.a.f, Species=="R. bromii"  & Media %in% c(2))
    n.all = nrow(z)
    z = subset(z, !grepl("Unreproducible", Class))
    z.blank = mean(z$BlankOD, na.rm=T)
    growth.prop = mean(!grepl("NoGrowth", z$Class))
    n.finished = sum(with(z, !grepl("NoGrowth", Class) & !grepl("Unfinished", Class)))
    n.growing = sum(with(z, !grepl("NoGrowth", Class)))
    
    z.f_pre = subset(z, !grepl("NoGrowth", Class) & (growth.prop>0.5)) 
    z.f = subset(z.f_pre, !grepl("Unfinish", Class))
    if(!nrow(z.f)) { z.f = z.f_pre }
    if(nrow(z.f)) {
      r = data.frame(
        MaxOD=mean(z.f$MaxOD-z.f$BlankOD, na.rm=T), MaxOD_sd=sd(z.f$MaxOD-z.f$BlankOD, na.rm=T),
        MaxOD_blank=mean(z.f$MaxOD, na.rm=T), MaxOD_blank_sd=sd(z.f$MaxOD, na.rm=T),
        Rate=mean(z.f$Rate, na.rm=T), Rate_sd=sd(z.f$Rate, na.rm=T))
      if(any(grepl("AUC18", colnames(z.f)))) {
        r$AUC18 = mean(z.f$AUC18, na.rm=T)
        r$AUC18_sd = sd(z.f$AUC18, na.rm=T)
      }
    } else {        
      if(!nogrowth) {
        r = data.frame(MaxOD=0, MaxOD_sd=NA, Rate=0, Rate_sd=NA)
      } else {
        z.f = subset(z, grepl("NoGrowth", Class))
        r = data.frame(
          MaxOD=mean(z.f$MaxOD-z.f$BlankOD, na.rm=T), MaxOD_sd=sd(z.f$MaxOD-z.f$BlankOD, na.rm=T), 
          MaxOD_blank=mean(z.f$MaxOD, na.rm=T), MaxOD_blank_sd=sd(z.f$MaxOD, na.rm=T), 
          Rate=mean(z.f$Rate, na.rm=T), Rate_sd=sd(z.f$Rate, na.rm=T))
      }
      if(any(grepl("AUC18", colnames(z.f)))) {
        r$AUC18 = 0
        r$AUC18_sd = NA
      }
    }
    
    r$MaxOD[is.nan(r$MaxOD)] = NA
    r$AUC18[is.nan(r$AUC18)] = NA
    r$Count = paste0(n.finished, ifelse(n.growing - n.finished > 0, paste0("(", n.growing - n.finished, ")"), ""), "/", n.all)
    r$Growing = n.growing
    r$Finished = n.finished
    r$Replicates = n.all
    if(any(r$MaxOD < 0)) r$MaxOD[r$MaxOD < 0] = 0
    if(any(grepl("AUC18", colnames(z.f))) && (any(r$AUC18 < 0 | is.na(r$AUC18)))) {
      r$AUC18[r$AUC18 < 0] = 0
    }
    
    return(r)
  })
  
  growth_matrix.q1
}


t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  pval = NA
  if(df) pval = 2*pt(-abs(t),df)
  list(se=se, t=t, p.value=pval)
}

cuves.filter_annotations = function(curves.a.f, variable="MaxOD", FUN=mean) 
{
  res = ddply(curves.a.f, .(Media, Species), function(z) {
    #z = subset(curves.a, Species == "A. muciniphila" & Media == 3)
    z.f = subset(z, !grepl("NoGrowth", Class) & !grepl("Unreproducible", Class) & !grepl("Unfinished", Class))
    if(nrow(z.f)) {
      return(z.f)
    }
    
    z.f = subset(z, !grepl("NoGrowth", Class) & !grepl("Unreproducible", Class))
    if(nrow(z.f)) {
      return(z.f)
    }
    
    z.f = subset(z, !grepl("Unreproducible", Class))
    if(nrow(z.f)) {
      return(z.f)
    }
  })
  
  res
}

medium_8_effect = function(curves_annotation) {
  curves_annotation.8 = subset(curves_annotation, !is.na(Species) & is.na(ConditionSpecies) & !is.na(Media) & Media %in% c(3,8))
  curves_annotation.8 = ddply(curves_annotation.8, .(File, Species), function(x) {
    x.0 = subset(x, !grepl("Unrep", Class))
    if(!sum(x.0$Media==3) | !sum(x.0$Media==8)) {
      res = data.frame(MaxOD8_effect=NA, Class="Unreproducible")
      return(res)
    }
    
    
    x.1 = subset(x, !grepl("Unfin|Unrep|NoGrowth", Class))
    if(sum(x.1$Media==3) & sum(x.1$Media==8)) {
      res = data.frame(MaxOD8_effect=with(x.1, max(MaxOD[Media==8] / MaxOD[Media==3])), Class="Good")
      return(res)
    }
    
    x.2 = subset(x, !grepl("NoGrowth|Unrep|Unfin", Class))
    if(length(unique(x.2$Media)) > 0) {
      x.2 = subset(x, !grepl("Unrep", Class))
      res = with(x.2, max(MaxOD[Media==8] / MaxOD[Media==3]))
      res = res[which.max(abs(res))]
      if(res > 8) res = 8
      if(res < 1/8) res = 1/8
      res = data.frame(MaxOD8_effect=res, Class=ifelse(res>1, "Recovery", "Lethal"))
      return(res)
    }
    
    x.3 = subset(x, grepl("Unfin", Class))
    if(sum(x.3$Media==3) & sum(x.3$Media==8)) {
      res = data.frame(MaxOD8_effect=with(x.3, max(MaxOD[Media==8] / MaxOD[Media==3])), Class="Unfinished")
      return(res)
    }
    
    x.4 = subset(x, grepl("Unfin|NoGrowth", Class))
    if(sum(x.4$Media==3) & sum(x.4$Media==8)) {
      res = data.frame(MaxOD8_effect=with(x.4, max(MaxOD[Media==8] / MaxOD[Media==3])), Class="NoGrowth")
      return(res)
    }
    
    res = data.frame(MaxOD8_effect=NA, Class=NA)
  })
  
  curves_annotation.8 = ddply(curves_annotation.8, .(Species), function(x) {
    x = ddply(x, .(Species, Class), summarize, MaxOD8_effect=mean(MaxOD8_effect)) # MaxOD8_effect=MaxOD8_effect[which.max(abs(log2(MaxOD8_effect)))]
    # x = ddply(x, .(Species, Class), summarize, MaxOD8_effect=MaxOD8_effect[which.max(abs(log2(MaxOD8_effect)))])
    x[order(match(x$Class, c("Good", "Recovery", "Lethal", "Unfinished", "NoGrowth", "Unreproducible"))),][1,]
  })
  
  curves_annotation.8
}


pheatmap_matrix = function(x, formula, column="SubsystemProp", na.value=NA) {
  x = dcast(x, formula, value.var=column)
  col = as.character(formula)[2]
  rownames(x) = x[[col]]
  x[[col]] = NULL
  x = data.matrix(x)
  x[is.na(x)] = na.value
  x
}


pheatmap_prop_text = function(m, percent=F) {
  m[m==0 | is.nan(m)] = NA
  if(!percent) {
    m[!is.na(m)] = sprintf("%.2f", m[!is.na(m)])
    m[is.na(m)] = ""
    m[m == "1.00"] = "1"
    m[grepl("^0", m)] = gsub("^0", "", m[grepl("^0", m)])
    m[grepl("0+$", m)] = gsub("0+$", "", m[grepl("0+$", m)])
  } else {
    m[!is.na(m)] = paste0(round(m[!is.na(m)]*100, 0), "%")
    m[is.na(m)] = ""
  }
  
  m
}

align_matrices = function(src, template) {
  template[, setdiff(colnames(template), colnames(src))] = NA
  template[setdiff(rownames(template), rownames(src)), ] = NA
  i.col = intersect(colnames(template), colnames(src))
  i.row = intersect(rownames(template), rownames(src))
  template[i.row, i.col] = src[i.row, i.col]
  template
}

element_grob.element_custom <- function(element, ...)  {
  segmentsGrob(c(1,0), c(0,0), c(0,0), c(0,1), gp=gpar(lwd=1))
}

border_custom <- function(...){
  structure(list(...), # this ... information is not used, btw
            class = c("element_custom","element_blank", "element") # inheritance test workaround
  )
}

compatible.replicates = function(z, corrected=F, column="OD") 
{
  if(nrow(z) < 2) return(NULL)
  if(!("Class" %in% colnames(z))) stop("No 'Class' column")
  if(!(column %in% colnames(z))) stop(paste0("No '", column, "' column"))
  
  lhs = combn(1:nrow(z), 2)[1,]
  rhs = combn(1:nrow(z), 2)[2,]
  different.volumes = z$Volume[lhs] != z$Volume[rhs]
  lhs = lhs[different.volumes]
  rhs = rhs[different.volumes]
  z.f = grepl("NoGrowth", z$Class[lhs]) & grepl("NoGrowth", z$Class[rhs]) |
    (corrected | (!grepl("Overgrowth", z$Class[lhs]) & !grepl("Overgrowth", z$Class[rhs]))) & grepl("Good|Falling", z$Class[lhs]) & grepl("Good|Falling", z$Class[rhs]) | 
    grepl("Overgrowth", z$Class[lhs]) & grepl("Overgrowth", z$Class[rhs]) | 
    grepl("Raising", z$Class[lhs]) & grepl("Raising", z$Class[rhs]) |
    (corrected & grepl("Raising|Good|Falling|Overgrowth", z$Class[lhs]) & grepl("Raising|Good|Falling|Overgrowth", z$Class[rhs]))
  
  if(!sum(z.f)) return(NULL)
  
  data.frame(
    x=z[lhs, column], y=z[rhs, column], 
    rhs=gsub("[^A-Z;]", "", z$Class[rhs]), lhs=gsub("[^A-Z;]", "", z$Class[lhs]), 
    rhs.loc=paste0(z$File[rhs], ":", z$TechnicalReplicates[rhs]),
    lhs.loc=paste0(z$File[lhs], ":", z$TechnicalReplicates[lhs]),
    compatible=z.f)
}

theme_slim = function(base_size=NA, ...)
{
  if(is.na(base_size)) base_size = 16
  theme_bw(base_size=base_size) + #eliminates baground, gridlines, and chart border  
    theme(
      panel.border=border_custom(), 
      panel.grid.major=element_blank(), 
      panel.grid.minor=element_blank(), 
      legend.background=element_blank(),
      legend.key=element_blank(),
      legend.key.size = unit(0.3, "cm"))
}
