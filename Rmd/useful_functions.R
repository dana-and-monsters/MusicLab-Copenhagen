adj.cor = function(df, p.adjust = FALSE, p.adjust.method = "none", threshold = 1, cor.method = "kendall"){
  cor_test = rcor.test(df, p.adjust = p.adjust, p.adjust.method = p.adjust.method, method = cor.method)
  r.mat = cor_test$cor.mat     # matrix of coefficient values
  p.list = cor_test$p.values   # p.list will be 3 columns
  
  # initiate empty matrix for the p values only
  MAT = matrix(, nrow = ncol(df), ncol = ncol(df))
  diag(MAT) = 0
  
  # starting to convert p.list to p.values matrix
  for(ind in 1:(length(p.list[,1]))){
    var1 = p.list[ind, 1] # var1
    var2 = p.list[ind, 2] # var2
    p.value = p.list[ind, 3] # p value
    MAT[var1, var2] = p.value
    MAT[var2, var1] = p.value
  }
  rownames(MAT) = names(df)
  colnames(MAT) = names(df)
  # At this point, MAT has the p.values in a matrix
  
  # subset only the coefficients with p values < 0.05 (or threshold)
  subset = ifelse(MAT < threshold, r.mat, NA)
  rownames(subset) = names(df)
  colnames(subset) = names(df)
  
  output = list(adj.p.values = MAT, threshold.r = subset)
  return(output)
}


lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 2)))
  as.character(as.expression(eq));
}

#Define gppr_theme() function

# library(extrafont)
# font_import()
theme_DSQ <- function(){ 
  font <- "GeosansLight"   #assign font family up front
  
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      #panel.grid.major = element_blank(),    #strip major gridlines
      #panel.grid.minor = element_blank(),    #strip minor gridlines
      #axis.ticks = element_blank(),          #strip axis ticks
      
      #since theme_minimal() already strips axis lines, 
      #we don't need to do that again
      
      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 20,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = 14),               #font size
      
      plot.caption = element_text(           #caption
        family = font,            #font family
        size = 9,                 #font size
        hjust = 1),               #right align
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 10),               #font size
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 9),                #font size
      
      axis.text.x = element_text(            #margin for axis text
        margin=margin(5, b = 10))
      
      #since the legend often requires manual tweaking 
      #based on plot content, don't define it here
    )
}

cross.correlation = function(corData, colours_list, Title, Subtitle){
  # create corrected p-values
  corellation = adj.cor(corData, p.adjust = TRUE, p.adjust.method = "BH", threshold = 0.05, cor.method = "kendall")
  pval = corellation$adj.p.values  # shows the p.values in a matrix format
  thresR = corellation$threshold.r  # shows only the p value <0.05
  
  # convert to triangle
  thresR[lower.tri(thresR)]<- NA
  pval[lower.tri(pval)]<- NA
  #create melted values
  melted_cormat<-melt(thresR,value.name = "r")
  melted_pmat<-melt(pval, value.name = "p")
  corPlotData<-full_join(melted_cormat, melted_pmat, by = c("Var1", "Var2"))
  corPlotData$r[corPlotData$r==1]<-NA
  corPlotData.c<-na.omit(corPlotData)
  corPlotData.c$stars<-cut(corPlotData.c$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance
  corPlotData.c$r<-round(corPlotData.c$r,2)
  
  # Create a ggheatmap
  ggheatmap <- ggplot(corPlotData.c, aes(Var2, Var1, fill = r))+
    geom_tile(color = "white")+
    scale_fill_gradientn(
      colours = colours_list,
      limits = c(-1,1)
    )+
    theme_minimal()+ 
    theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45, size = 4))+
    theme(axis.text.y = element_text(size = 4))+
    labs(title = Title,
         subtitle = Subtitle)+
    coord_fixed()
  
  print(ggheatmap)
  
  ggheatmap2<-ggheatmap+
    geom_text(aes(Var2, Var1, label = paste(stars,r, sep ="\n")), color = "black", size = 1.3) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.8),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 5, barheight = 0.4,
                                 title.position = "top", title.hjust = 0.5, 
                                 title.theme = element_text(size = 5),
                                 label.theme = element_text(size = 5)))
  
  print(ggheatmap2)
}
