library(ggplot2)

ggplot.cor_data = function(d, mapping, facets=NULL, scales="fixed", facet_fun="none", method="spearman")
{
    x = as.character(mapping$x)
    y = as.character(mapping$y)
    
    facets.char = c()
    if(!is.null(facets)) {
        facets.char = rownames(attr(terms(facets), "factors"))
    }
    if(length(facets.char) == 0) fx = "," 
    if(length(facets.char) == 1) fx = facets.char[1]
    if(length(facets.char) > 1) fx = facets.char[2]
    if(facet_fun == "none") fx = "."

    if(length(facets.char) == 0) fy = "." 
    if(length(facets.char) == 1) fy = "." 
    if(length(facets.char) > 1) fy = facets.char[1]
    if(facet_fun == "none") fy = "."
    
    gr = c()
    if(fx != ".") { gr = c(gr, fx) 
    } else { fx = c() }
    if(fy != ".") { gr = c(gr, fy) 
    } else { fy = c() }
    
    if(facet_fun == "grid" & scales == "free") scales = "free_colrow"
    if(facet_fun == "grid" & scales == "free_x") scales = "free_col"
    if(facet_fun == "grid" & scales == "free_y") scales = "free_row"
    if(facet_fun == "none") scales = "free"
    
    horizontal = NULL
    if(scales %in% c("free_colrow", "free_col"))
    {
        horizontal = ddply(d, fx, function(z) { 
            data.frame(left = min(z[[x]], na.rm=T), right = max(z[[x]], na.rm=T))
        })
    }

    if(scales %in% c("fixed", "free_y", "free_row"))
    {
        left = min(d[[x]], na.rm=T)
        right = max(d[[x]], na.rm=T)
        horizontal = ddply(d, fx, function(z) { data.frame(left=left, right=right) })
    }
        
    
    vertical = NULL
    if(scales %in% c("free_colrow", "free_row"))
    {
        vertical = ddply(d, fy, function(z) { 
            data.frame(
                bottom = min(z[[y]], na.rm=T),
                top = max(z[[y]], na.rm=T))
        })
        
        vertical_ci = ddply(d, gr, function(z) {
            z.lm = lm(z[[y]] ~ z[[x]])
            z.pred = predict(z.lm, interval="confidence")
            
            if(sum(complete.cases(z.pred)) < 3) {
                bottom_ci = min(z[[y]], na.rm=T)
                top_ci = max(z[[y]], na.rm=T)
            } else {
                bottom_ci = min(z.pred[,"lwr"], na.rm=T)
                top_ci = max(z.pred[,"upr"], na.rm=T)
            }
    
            data.frame(bottom_ci=bottom_ci, top_ci=top_ci)
        })
        
        vertical_ci = ddply(vertical_ci, fy, summarize, bottom_ci=min(bottom_ci, na.rm=T), top_ci=max(top_ci, na.rm=T))
        vertical = ddply(merge(vertical, vertical_ci), fy, mutate, top_ci=pmax(top, top_ci), bottom_ci=pmin(bottom, bottom_ci))
    }

    if(scales %in% c("fixed", "free_y", "free_row"))
    {
        bottom = min(d[[y]], na.rm=T)
        top = max(d[[y]], na.rm=T)
        
        vertical_ci = ddply(d, gr, function(z) {
            z.lm = lm(z[[y]] ~ z[[x]])
            z.pred = predict(z.lm, interval="confidence")
            
            bottom_ci = NULL
            if(sum(complete.cases(z.pred)) < 3) {
                bottom_ci = min(z[[y]], na.rm=T)
                top_ci = max(z[[y]], na.rm=T)
            } else {
                bottom_ci = min(z.pred[,"lwr"], na.rm=T)
                top_ci = max(z.pred[,"upr"], na.rm=T)
            }
            
            data.frame(bottom_ci=bottom_ci, top_ci=top_ci)
        })
        
        vertical_ci = ddply(vertical_ci, c(), summarize, bottom_ci=min(bottom_ci, na.rm=T), top_ci=max(top_ci, na.rm=T))
        vertical = ddply(d, fy, function(z) { 
            data.frame(
                top=top, top_ci=pmax(top, vertical_ci$top_ci, na.rm=T), 
                bottom=bottom, bottom_ci=pmin(bottom, vertical_ci$bottom_ci, na.rm=T)) 
        })
    }

    d.group = ddply(d, gr, function(z) {
        z.res = data.frame()
        
        x.na = is.na(z[[x]])
        x.inf = !x.na & sapply(z[[x]], is.infinite)
        y.na = is.na(z[[y]])
        y.inf = !y.na & sapply(z[[y]], is.infinite)
        
        xy.na = x.na | y.na
        xy.inf = !(x.na | y.na) & (x.inf | y.inf)
        xy.infsum = sum(!x.na & !y.na & xor(x.inf, y.inf))
        
        z.x = z[[x]][!xy.na & !xy.inf]
        z.y = z[[y]][!xy.na & !xy.inf]
        
        if(length(z.x) < 3 || length(z.y) < 3 || sd(z.x, na.rm=T) == 0 || sd(z.y, na.rm=T) == 0) {
            z.res = data.frame(rho=NA, R2=NA, pval=NA, sq_dist=NA)
        } else {
            zz <<- z
            z.cor = cor.test(z[[y]], z[[x]], method=method, exact=F)        
            z.sq_dist = sum(1-pmin(z[[y]],z[[x]])/pmax(z[[y]],z[[x]]))/nrow(z)
            z.res = data.frame(rho=z.cor$estimate, R2=z.cor$estimate^2, pval=z.cor$p.value, sq_dist=z.sq_dist)
        }
        
        if(scales %in% c("free", "free_x")) {
            if(length(z.x) == 0) {
                z.res$left = 0
                z.res$right = 0
            } else {
                z.res$left = min(z.x, na.rm=T)
                z.res$right = max(z.x, na.rm=T)
            }
            
        } else {
            h = horizontal[1,]
            if(nrow(horizontal) > 1) { h = horizontal[horizontal[,1] == z[1,fx], ] }
            
            z.res$left = h$left
            z.res$right = h$right
        }
        
        if(scales %in% c("free", "free_y")) {
            if(length(z.y) == 0) {
                z.res$bottom = 0
                z.res$top = 0
            } else {
                z.res$bottom = min(z.y, na.rm=T)
                z.res$top = max(z.y, na.rm=T)
            }
            
            if(length(z.x) < 3 || length(z.y) < 3)
            {
                z.res$bottom_ci = z.res$bottom
                z.res$top_ci = z.res$top
            } else {
                z.lm = lm(z.y ~ z.x)
                z.pred = predict(z.lm, interval="confidence")
                if(sum(complete.cases(z.pred)) > 2) {
                    z.res$bottom_ci = min(z.pred[,"lwr"], na.rm=T)
                    z.res$top_ci = max(z.pred[,"upr"], na.rm=T)
                } else {
                    z.res$bottom_ci = z.res$bottom
                    z.res$top_ci = z.res$top
                }
            }
        } else {
            v = vertical[1,]
            if(nrow(vertical) > 1) { v = vertical[vertical[,1] == z[1,fy], ] }
            
            z.res$bottom = v$bottom
            z.res$top = v$top
            z.res$bottom_ci = v$bottom_ci
            z.res$top_ci = v$top_ci
        }
        
        z.res$infsum = xy.infsum
        z.res
    })
    
    d.group$top_ci = pmax(d.group$top_ci, d.group$top, na.rm=T)
    d.group$bottom_ci = pmin(d.group$bottom_ci, d.group$bottom_ci, na.rm=T)
    
    d.group$short_str = apply(d.group, 1, function(z) {
        pval = as.numeric(z["pval"])
        rho = as.numeric(z["R2"])
        infsum = as.numeric(z["infsum"])
        if(is.na(rho) || is.na(pval)) return("")
        
        p = ""
        if(pval <= 0.001) p = "***"
        if(pval <= 0.01) p = "**"
        if(pval <= 0.05) p = "*"
        if(pval > 0.05) p = paste0(" (p=", round(pval, 2), ")")
        rho_str = round(rho, 2)   
        
        infsum_str = ""
        if(infsum > 0) infsum_str = paste0(" ! ", infsum)
        paste0("R2=", rho_str, p, infsum_str)
    })
    
    rho_str = round(d.group$rho, 2)
    pval_str = format(round(d.group$pval, 3), scientific=T)
    infsum_str = ifelse(d.group$infsum > 0, paste0(" ! -", d.group$infsum), "")
    d.group$long_str = paste0("rho=", rho_str, " (P-value: ", pval_str, ")", infsum_str)
    
    d.group
}

ggplot.piechart = function(d, factor="factor", count="count", facet=NA, title="", guide="")
{
    if(is.na(facet))
    {
        d$pos = cumsum(d[[count]]) - d[[count]]/2
    } else {
        d = d[order(d[[facet]], d[[factor]]),]
        d$facet_tmp = d[[facet]]
        d = ddply(d, facet,  function(d1) {
            d1$pos = cumsum(d1[[count]]) - d1[[count]]/2
            d1
        })
    }
    
    
    p = ggplot(d, aes(x=1)) + 
        geom_bar(aes_string(y=count, fill=factor), stat="identity") +
        geom_text(aes_string(y="pos", label=count), x=1.3) + 
        coord_polar(theta="y") +
        theme_bw() +
        theme(axis.line=element_blank(), 
              axis.title=element_blank(), 
              axis.text=element_blank(),
              axis.ticks=element_blank(),
              panel.border=element_blank(), 
              panel.grid=element_blank()) +
        guides(fill=guide_legend(title=guide)) + 
        ggtitle(title)
    
    if(is.na(facet))
    {
        return(p)
    } else {
        ggfacet = facet_wrap(~facet_tmp)
        if(length(unique(d[[facet]])) == 2) {
            ggfacet$ncol = 1
        }
            
        return(p + ggfacet)
    }
}

ggplot.legend<-function(g){
    tmp <- ggplot_gtable(ggplot_build(g))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}