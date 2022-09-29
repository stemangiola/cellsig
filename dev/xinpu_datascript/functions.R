
plot_markers = function(my_ct1, my_ct2){
  ARMET::ARMET_ref %>% 
    filter(ct1==my_ct1 & ct2==my_ct2) %>% 
    nanny::nest_subset(data = -symbol) %>% 
    arrange(rank) %>% slice(1:10) %>% 
    mutate(symbol = factor(symbol, levels = .$symbol)) %>% 
    unnest(data) %>% 
    filter(`Cell type category` %in% c(my_ct1, my_ct2)) %>% 
    ggplot(aes(`Cell type category`, `count scaled bayes`+1)) +
    geom_boxplot() +
    geom_jitter() +
    facet_wrap(~symbol) +
    scale_y_log10()
}

ComputeMarkers <- function(sobj){
  
  markers.all <- list()
  
  
  Idents(sobj) <- sobj@meta.data[["cluster"]]
  
  sobj.markers <- FindAllMarkers(object= sobj, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
  
  # sobj.markers$AUROC <- NA
  # for(j in 1:nrow(sobj.markers)){
  #   sobj.markers$AUROC[j] <- roc.curve(scores.class0 = GetAssayData(sobj)[sobj.markers$gene[j],], weights.class0 = sobj@active.ident == sobj.markers$cluster[j])$auc
  # }
  
  top.10 <- sobj.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
  DoHeatmap(object = sobj, features = top.10$gene) + NoLegend()
  
}