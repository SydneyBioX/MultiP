## functions to calculate summaries of MultiP

library(table1)
library(ggnewscale)
library(data.tree)
library(DiagrammeR)
library(DiagrammeRsvg)
library(png)
library(grid)
library(gridExtra)
library(rsvg)
library(DT)
library(plyr)


classification_metrics <- function(fit, test){
  classes_levels <- unique(c(fit,test))
  if(length(classes_levels) == 1){
    return(1)
  } else{
    fit_factor <- factor(fit, levels = classes_levels)
    test_factor <- factor(test, levels = classes_levels)  
    cont_table <- table(fit_factor,test_factor)
    TN <- cont_table[1,1]
    TP <- cont_table[2,2]
    FN <- cont_table[1,2]
    FP <- cont_table[2,1]
    P <- TP + FN
    N <- FP + TN
    PP <- TP + FP
    PN <- FN + TN
    
    TPR <- TP/(TP+FN)
    TNR <- TN/(TN+FP)
    ACC <- (TP + TN)/(P+N)
    BA <- (TPR + TNR)/2
    Precision <- TP/PP
    Recall = TP/P
    Specificity <- TN/N
    F1 = 2*TP/(2*TP + FP + FN)
    
    metrics = data.frame(metric = c("ACC", "BA", "Precision", "Recall", "Specificity", "F1"),
                         value = c(ACC, BA, Precision, Recall, Specificity, F1))
    return(metrics)
  }
}

# Calculates weighted score for table
calculateScore <- function(data, properties, weights){
  originalLength = length(data)

  #for each property and weight, create a new column corresponding to its contribution to the total weight
  for(i in 1:length(properties)){
    colname = paste0(properties[i],".weight")
    if(properties_direction[i] == "high"){
      data[colname] <- data[properties[i]]/max(data[properties[i]]) 
    } else{
      data[colname] <- min(data[properties[i]])/data[properties[i]]
    }
  }
  
  weighted_scores <- colSums(apply(data[(originalLength+1):length(data)], 1, function(x) x * weights))
  data = data %>% mutate(Score = weighted_scores) # sum up all new cols
  
  return(data)
}

calculateScore_ranking <- function(data, properties, properties_direction, weights){
  n_vars = length(data)
  n_perms = nrow(data)
  
  #for each property and weight, create a new column corresponding to its rank
  for(i in 1:length(properties)){
    colname = paste0(properties[i],".weight")
    data = data %>% arrange(get(properties[i]))
    if(properties_direction[i] == "high"){
      data[colname] <- 1:n_perms
    } else{
      data[colname] <- n_perms:1
    }
  }
  
  weighted_scores <- colSums(apply(data[(n_vars+1):length(data)], 1, function(x) x * weights))
  data = data %>% mutate(Score = weighted_scores) # sum up all new cols
  
  return(data)
}

abbreviate <- function(sequences){
  sapply(sequences, function(x) strsplit(x, "-")[[1]] %>% substr(1,1) %>% toupper %>% paste(collapse = "-"))
}

# Summary table
summaryTable <- function(MultiPobject){

  tab = MultiPobject$TSER.summary %>% 
    dplyr::distinct() %>%
    calculateScore_ranking(MultiPobject$properties,MultiPobject$properties_direction,MultiPobject$weights) %>%
    select(Platform.Sequence, balancedAcc, TSA.retained, Cost.byPlatform,Cost.Total_str, Cost.Total2, Score) %>%
    arrange(-Score) 
  
  if(length(unique(tab$Cost.Total2)) == 1){
    tab <- select(tab, -Cost.Total2) %>% 
      DT::datatable(option=list(columnDefs=list(list(className="dt-center", targets=list(1,6))),
                                                                   pageLength=10),
                                                       colnames = c("Sequence","Balanced Acc", "Overall Acc", "Platform Costs", "Total Cost", "Score")) %>%  
      DT::formatRound(columns=c(2, 3, 6), digits=3) 
  } else{
    tab <- tab %>% DT::datatable(option=list(columnDefs=list(list(className="dt-center", targets=list(1,6))),
                                             pageLength=10),
                                 colnames = c("Sequence","Balanced Acc", "Overall Acc", "Platform Costs", "Total Cost", "Total Cost2", "Score")) %>%  
      DT::formatRound(columns=c(2, 3, 7), digits=3) 
    }
  
  return(tab)
}

topSequence <- function(MultiPobject){
  TSUM = MultiPobject$TSER.summary
  TSUM = calculateScore_ranking(TSUM,properties,weights)
  tab = TSUM %>% arrange(-Score)
  return(tab$Platform.Sequence[1])
}

topBAccSequence <- function(MultiPobject){
  TSUM = MultiPobject$TSER.summary
  tab = TSUM %>% arrange(-balancedAcc)
  return(tab$Platform.Sequence[1])
}

bubblePlot <- function(MultiPobject, nVis = 8){

  tser_sum = MultiPobject$TSER.summary %>%
    dplyr::distinct() %>% 
    calculateScore_ranking(MultiPobject$properties,MultiPobject$properties_direction,MultiPobject$weights) %>%
    arrange(-Score) %>%
    head(nVis) %>%
    arrange(Platform.Sequence)
  platformSequence = tser_sum$Platform.Sequence
  nperms = nrow(tser_sum)
  colours = RColorBrewer::brewer.pal(nperms, "Set2")[1:nperms]
  
  theme_set(theme_bw())
  g <- ggplot(tser_sum, aes(x = balancedAcc, y = Cost.Total, fill = factor(abbreviate(Platform.Sequence)), alpha = TSA.retained.weight)) 
  n_platforms_used <- 0
  #for each permutation we want to add concentric circles representing proportions of retention
  for(nperm in 1:nperms){
    platform_order = unlist(strsplit(platformSequence[nperm], "-"))
    df = MultiPobject$Results.stratification[[nperm]] %>%
      dplyr::mutate(platform = factor(platform, levels = platform_order)) %>%
      dplyr::filter(strata == "Retained") %>%
      dplyr::count(platform) %>%
      dplyr::mutate(size = cumsum(n)/sum(n)) %>%
      dplyr::mutate(alpha = seq(from = 0.9, to = 0.4, length.out = nrow(.)))
    
    n_platforms_used = max(n_platforms_used, nrow(df))
    
    #for each level of the permutation make a circle representing the relative proportion of the level (e.g. retained for each experiment)
    for(i in length(df$platform):1){
      g = g + annotate("point", 
                       x=tser_sum$balancedAcc[nperm], 
                       y= tser_sum$Cost.Total[nperm], 
                       size = df$size[i]*10,
                       colour = colours[nperm],
                       alpha = df$alpha[i])
    }
  }

  xmin = min(tser_sum$balancedAcc)
  xmax = max(tser_sum$balancedAcc)
  xdiff = xmax - xmin
  ymin = min(tser_sum$Cost.Total)
  ymax = max(tser_sum$Cost.Total)
  ydiff = ymax - ymin
  
  g = g + geom_point(x= -1, y = -1, size = 0) +
    coord_cartesian(xlim= c(xmin - xdiff/10,min(xmax + xdiff/10,1)),
                    ylim=c(ymin - ydiff/10,ymax + ydiff/10)) +
    scale_fill_manual(values = colours, name = "Sequence") + 
    scale_alpha(name='Platform', breaks=n_platforms_used:1, labels=paste("Platform", 1:n_platforms_used)) +
    guides(alpha=guide_legend(override.aes=list(size = 5, fill = "black")),
           fill = guide_legend(override.aes = list(size = 5, alpha = 1, colour = colours))) +
    labs(title="Pathway Comparison Plot") + xlab("Sequence Balanced Accuracy") + ylab("Total Cost") +
    theme(legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          axis.text=element_text(size=15),
          axis.title=element_text(size=15),
          plot.title=element_text(size=20))
  return(g)
}

strataPlot <- function(MultiPobject, sequence){
  nperm = which(MultiPobject$TSER.summary$Platform.Sequence == sequence)[1]
  stra = MultiPobject$Results.stratification[[nperm]]
  stra$platform = factor(stra$platform, levels = c(unique(stra$platform), "True Class"))
  stra = dplyr::arrange(stra, platform, true_class, sser) %>% filter(strata == "Retained")
  stra$ID = 1:nrow(stra)
  stra$colour = mapvalues(stra$true_class, from = sort(unique(stra$true_class)), to = c("#4DAF4A", "#984EA3"))
    
  n_platforms = length(unique(stra$platform))

  p = ggplot(mapping = aes(x = ID, y = platform), data = stra) +
    geom_tile(aes(fill= true_class)) +
    scale_fill_manual(values = c("#4DAF4A", "#984EA3"))  +
    labs(title=paste("Sequence:", abbreviate(sequence)),
         fill="True class",
         y="",
         x="")+
    guides(fill = guide_legend(title.position = "top")) +
    new_scale_fill() +
    geom_tile(aes(fill= 1-sser)) +
    scale_fill_gradient(low="#377EB8", high="#E41A1C") +
    labs(fill = "Accuracy") +
    guides(fill = guide_colorbar(title.position = "top")) +
    theme(panel.background = element_blank(),
          axis.text.x = element_blank(),
          aspect.ratio=1/4, 
          plot.title=element_text(face ="bold", size = 20),
          legend.title = element_text(face ="bold", size = 12),
          legend.text = element_text(size = 10),
          legend.position="bottom",
          axis.text=element_text(size=15)) +
    annotate("tile",
               x=stra$ID,
               y= n_platforms + 0.8,
               height = 0.6,
               fill = stra$colour)  +
    coord_cartesian(expand = FALSE)  +
    geom_hline(yintercept = n_platforms + 0.5, size = 1.5) 
    
  return(p)
}

GetEdgeLabel <- function(node){
  if(node$isRoot || node$nodeType == "Platform"){
    label = NULL
  } else {
    val = round((node$counter/node$sum)*100)
    label = paste0( (val),  "% (", node$counter, ")")
  }
  return(label)
}

GetFillColour <- function(node){
  if(node$isRoot || node$nodeType == "Platform"){
    colour = "#87C37D"
  } else if(node$nodeType == "Class1"){
    colour = "#ADCEE0"
  } else if(node$nodeType == "Class2"){
    colour = "salmon"
  } else {
    colour = "snow3"
  }
  return(colour)
}

GetNodeShape <- function(node){
  if(node$isRoot || node$nodeType == "Platform"){
    shape = "oval"
  } else {
    shape = "box"
  }
  return(shape)
}


flowChart <- function(MultiPobject, sequence){
  nperm = which(MultiPobject$TSER.summary$Platform.Sequence == sequence)[1]
  nplatforms=ncol(MultiPobject$myperms)
  nperms=nrow(MultiPobject$myperms)
  stra = MultiPobject$Results.stratification[[nperm]]
  n_total = sum(MultiPobject$Results.stratification[[nperm]]$strata == "Retained")
  
  platforms = unlist(strsplit(sequence, "-"))
  outcomes <- unique(stra$predict_class)
  class1 <- outcomes[1]
  class2 <- outcomes[2]
  
  treeRoot = Node$new(platforms[1])
  currNode = treeRoot
  
  for(nplatform in 1:nplatforms){
    
    retain_level = levels(stra$platformstrata)[((nplatform-1)*3)+1]
    progress_level = levels(stra$platformstrata)[((nplatform-1)*3)+2]
    noprocess_level = levels(stra$platformstrata)[((nplatform-1)*3)+3]
    
    n_retained = stra %>% filter(platformstrata == retain_level) %>% nrow()
    n_class1 = stra %>% filter(platformstrata == retain_level, predict_class == class1) %>% nrow()
    n_class2 = stra %>% filter(platformstrata == retain_level, predict_class == class2) %>% nrow()
    n_progress = stra %>% filter(platformstrata == progress_level) %>% nrow()
    n_noprocess = stra %>% filter(platformstrata == noprocess_level) %>% nrow() 
    
    platform = platforms[nplatform]
    
    if(n_retained>0){ # Class 1 node
      class1_preds = currNode$AddChild(class1, counter = n_class1, sum = n_total, nodeType = "Class1")
    }
    
    if(n_progress>0){ # To Progress node
      uncertain = currNode$AddChild("Uncertain", counter = n_progress, sum = n_total, nodeType = "Uncertain")
    }
    
    if(n_retained>0){ # Class 2 node
      class2_preds = currNode$AddChild(class2, counter = n_class2, sum = n_total, nodeType = "Class2")
    }
    
    if(n_noprocess>0){ # Not processed node
      notprog = currNode$AddChild(noprocess_level, counter = n_noprocess, sum = n_total, nodeType = "NotProcessed")
    }
    
    if(n_progress==0){
      break;
    } else{ # Platform node
      currNode = uncertain
      progress = currNode$AddChild(platforms[nplatform+1], nodeType = "Platform")
      currNode = progress
    }
     
  }
  
  SetGraphStyle(treeRoot, rankdir = "LR")
  SetEdgeStyle(treeRoot, fontname = 'helvetica', label = GetEdgeLabel)
  SetNodeStyle(treeRoot, style="filled", shape = GetNodeShape, fontcolor = "black", fillcolor = GetFillColour, fontname = 'helvetica')
  plot(treeRoot)
}

cohortSummary <- function(MultiPobject, sequence, platform_input, variables_of_interest){
  nperm = which(MultiPobject$TSER.summary$Platform.Sequence == sequence)[1]
  strata_perm = MultiPobject$Results.stratification[[nperm]]
  
  platform_outcome = strata_perm %>% filter(platform == platform_input) %>% select(SampleID, strata)
  clin_filtered = colData(MultiPobject$MultiAssayExperiment) %>% 
    data.frame() %>% 
    mutate(SampleID = rownames(.)) %>%
    inner_join(platform_outcome, by = "SampleID") %>%
    select(-SampleID) %>%
    select(all_of(variables_of_interest), strata)
  
  table1(~ . | strata, data=clin_filtered, rowlabelhead = platform_input)
}

featureImportance <- function(MultiPobject, feat, n_feat){
  
  selected_features  <- MultiPobject$generatedModelsFeatures[[feat]]
  n_feat = min(n_feat, nrow(selected_features))
  selected_features$Feature <- factor(selected_features$Feature, levels = selected_features$Feature)
  
  p <- selected_features %>%
    head(n_feat) %>%
    ggplot(aes(x=Feature, y=Freq)) + 
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title= paste("Top selected features:",feat), y ="Frequency", x = "Feature")
  
  return(p)
}

calculate_test_stat<- function(model){
  n1 <- model$est[[1]]$n
  xbar1 <- model$est[[1]]$xbar
  n2 <- model$est[[2]]$n
  xbar2 <- model$est[[2]]$xbar
  varp <- model$var_pool
  test_stat <- (xbar1 - xbar2)/sqrt(varp)
  return(test_stat)
}

DLDAloadings <- function(MultiPobject, platform, n_features = 20){
  n_features <- min(n_features, length(MultiPobject$generatedModels[[platform]]))
  model_list <- MultiPobject$generatedModels[[platform]]
  
  means <- sapply(model_list, calculate_test_stat) %>% apply(1, mean)
  sds <- sapply(model_list, calculate_test_stat) %>% apply(1, sd)
  features <- names(means)
  pvals <- pnorm(abs(means), lower.tail = FALSE)
  
  df_test_stats <- data.frame(features = features, means = means, sds = sds, pvals = pvals)
  df_test_stats <-  mutate(df_test_stats, signif_lvl = case_when(
    pvals < 0.001 ~ "***",
    pvals < 0.01 ~ "**",
    pvals < 0.05 ~ "*",
    pvals < 0.1 ~ "B0",
    TRUE ~ ""
  )) %>% 
    arrange(-abs(means)) %>% 
    head(n_features) %>% 
    arrange(-means) %>%
    mutate(features = factor(features, levels = features))
  
  p <- df_test_stats %>%
    ggplot() +
    geom_point(aes(x=means, y=features), position=position_dodge(width=0.5)) +
    geom_vline(aes(xintercept = 0), linetype="dotted") +
    geom_errorbarh(aes(y=features, xmin=means - 2*sds, xmax=means + 2*sds, 
                       height = 0.25), 
                   position=position_dodge(width=0.5),
                   alpha=0.75) +
    geom_text(aes(y=features, x=Inf, label=signif_lvl), 
              position=position_dodge(width=0.5),
              hjust=0, vjust=0.75,
              size=6, show.legend = FALSE) +
    theme_bw() +
    theme(strip.background = element_rect(colour = NA, fill="grey90"),
          legend.position = "bottom",
          plot.margin = unit(c(1,2,1,1), "lines")
    ) +
    coord_cartesian(clip = 'off') +
    labs(title=paste0(platform, ": Test statistics"),
         col="Adjusted for",
         y="",
         x="Test stats")
  p
  
}
