library(e1071)
library(randomForest)
library(sparsediscrim)
library(hardhat)
library(caret)
library(magrittr)
library(rms)
library(Hmisc)
library(skimr)

perm <- function(v){
  n <- length(v)
  if (n == 1) v
  else {
    X <- NULL
    for(i in 1:n) X <- rbind(X, cbind(v[i], perm(v[-i])))
    X
  }
}

balanced_accuracy <- function(fit, test){
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
    TPR <- TP/(TP+FN)
    TNR <- TN/(TN+FP)
    BA <- (TPR + TNR)/2
    return(BA)
  }
}

generateModels <- function(MultiAssayExperiment, 
                           classes,
                           platformList, 
                           crossValParams, modellingParams,
                           characteristics, verbose, seed){
    n = nrow(colData(MultiAssayExperiment))
    n_platforms = length(platformList)

    resultList = list()
    modelList = list()
    featureList = list()
    
    for(i in 1:n_platforms){
      platform = platformList[i]

      if(platform == "clinical"){
        original_selectParams <- modellingParams@selectParams
        modellingParams@selectParams <- NULL
      }
      
      use <- data.frame(assay = platform, feature = "all")
      set.seed(seed)
      DMresults <- suppressWarnings(runTests(MultiAssayExperiment, 
                                             useFeatures = use, 
                                             outcome = classes,
                                             crossValParams = crossValParams, 
                                             modellingParams = modellingParams,
                                             characteristics = characteristics, 
                                             verbose = verbose))
      
      sampleIDs = rownames(colData(MultiAssayExperiment))
      pred_table <- table(DMresults@predictions$sample, DMresults@predictions$class)
      
      class1 <- colnames(pred_table)[1]
      class2 <- colnames(pred_table)[2]
      class1_preds = pred_table[,1]
      class1_preds = class1_preds[sampleIDs]
      class2_preds = pred_table[,2]
      class2_preds = class2_preds[sampleIDs]
      
      repeats = crossValParams@permutations
      
      preds <- data.frame(sampleID = sampleIDs, class1_pred = class1_preds, class2_pred = class2_preds)
      preds <- preds %>% 
        mutate(predict_class = ifelse(class1_pred > class2_pred, class1, class2),
               true_class = as.character(colData(MultiAssayExperiment)[,classes])) %>% 
        mutate(sser = ifelse(class1 == true_class, class2_pred/repeats, class1_pred/repeats))
      
      if(platform == "clinical"){
        selected_features <- data.frame(Feature = DMresults@originalFeatures$feature, Freq = 1)
        modellingParams@selectParams <- original_selectParams
      } else{
        selected_features <- sapply(DMresults@chosenFeatures, function(x) x %>% data.frame %>% pull("feature")) %>% 
          unlist %>% 
          table %>% 
          data.frame %>% 
          arrange(-Freq) 
        
        names(selected_features)[1] <- "Feature"
        
        selected_features <- selected_features %>%
          mutate(Feature = as.character(Feature),
                 Freq = Freq/length(DMresults@chosenFeatures))
      }
      
      resultList[[platform]] = preds 
      modelList[[platform]] = DMresults@models
      featureList[[platform]] = selected_features
    }
    return(list(resultList = resultList, modelList = modelList, featureList = featureList))
}

modelsPredict = function(MultiAssayExperiment,
                         generatedModelsResults,
                         platform,
                         classes, 
                         minPlatformSize,
                         verbose){

    # Check CV has enough samples
    if(length(unlist(colnames(MultiAssayExperiment[1]))) <= minPlatformSize){
      if(verbose){
        print(length(unlist(colnames(MultiAssayExperiment[1]))))
        print("Error: Low row count, unable to use CV")
      }
      errorTable =  tibble(SampleID = character(),
                           sser = double(),
                           predict_class = character(),
                           true_class = character())
    }else{
      keep_id <- rownames(colData(MultiAssayExperiment))
      preds <- generatedModelsResults[[platform]]
      preds_filter <- preds[keep_id,]

      errorTable =  tibble(SampleID = preds_filter$sampleID,
                           sser = preds_filter$sser,
                           predict_class = preds_filter$predict_class,
                           true_class = preds_filter$true_class) 
    }
  
  return(list(removed = NULL, SSER = errorTable))
}



MultiP.algm <- function(MultiAssayExperiment,
                      platformList, fixedPlatforms = 0,
                      confidenceCutoff,platformUnitCosts = c(100, 500, 1000), platformUnitCosts2 = NA, 
                      properties = c("TSA.retained", "Cost.Total"), properties_direction = c("high", "low"), weights = c(0.5, 0.5),
                      performanceType = "Sample Error",
                      classes = NULL, crossValParams = NULL, modellingParams = NULL, 
                      characteristics = NULL, minPlatformSize = 10, 
                      seed=1, verbose=F, predict = FALSE, resultList = NULL){
  
  nplatforms = length(platformList)

  if(fixedPlatforms == 0){
    myperms = perm(1:nplatforms)}
  else{  # implicitly assumes fixedPlatforms < nplatforms
    perm_layers = perm((fixedPlatforms + 1):nplatforms)
    fixed_layers = matrix(1:fixedPlatforms, ncol = fixedPlatforms, nrow = nrow(perm_layers), byrow = TRUE)
    myperms = cbind(fixed_layers, perm_layers)
  }
  
  nperms = nrow(myperms)
  MTlist = list()  # to store results of each perm

  ## geneate models
  if(!predict){
    generatedModels <- generateModels(MultiAssayExperiment, 
                                      classes,
                                      platformList, 
                                      crossValParams, modellingParams,
                                      characteristics, verbose, seed)
  } else{
    generatedModels <- list(resultList = resultList, modelList = NA, featureList = NA)
  }
  
  #for each permutation (row)
  for(nperm in 1:nperms){
    MTunits = list()  #store results for each platform in this perm
    id.retained = list()  #create list for # retained for this permutation
    platformseq = NULL  #create sequence for current platform
    
    for(nplatform in 1:nplatforms){  #for each platform of each permutation
      platformseq = c(platformseq, platformList[[myperms[nperm, nplatform]]])  #set platformseq and platform as per myperms
    }
    
    platformseq = paste0(platformseq, collapse="-")  #create platform-sequence string
    
    for(nplatform in 1:nplatforms){  #for each platform of each perm
      #Final platform check: this will later enforce that all samples are allocated to a leaf node regardless of error
      finalPlatform = FALSE
      if(nplatform == nplatforms || (nrow(colData(MultiAssayExperiment)) - length(id.retained)) < minPlatformSize)
      {
        finalPlatform = TRUE
      }
      
      z = myperms[nperm, nplatform]  #z = which layer of the MAE is being used
      platform = platformList[[z]]
      confidenceCutoff_block = confidenceCutoff
      
      if(finalPlatform){
        confidenceCutoff_block = 1
      }
      
      #get the column index of the class from the MAE
      classIndex = grep(classes,colnames(colData(MultiAssayExperiment)))
      
      #Remove retained
      MAEwithoutRetained = MultiAssayExperiment
      retainList = unlist(id.retained)
      retainListLogical = !(rownames(colData(MAEwithoutRetained)) %in% retainList)
      
      if(length(id.retained)>0){
          MAEwithoutRetained = MAEwithoutRetained[,retainListLogical, ]
      }
      
      #Create platformstring e.g. (1) Clin-Histo-Nano: <Clin>
      print(paste0("(", nperm, ") ", platformseq, ": <", platform, ">"))

      #pass model and data to MTblock which returns retained status of current sequence layer
      MTunits[[nplatform]] = MTblock(data=MAEwithoutRetained, #id.retained=id.retained,
                                 generatedModelsResults = generatedModels$resultList,
                                 confidenceCutoff=confidenceCutoff_block, platform=platform, plotlabel=platform,
                                 seed=seed, verbose, 
                                 classes=classes, crossValParams=crossValParams, modellingParams=modellingParams, 
                                 characteristics=characteristics, performanceType=performanceType, finalPlatform=finalPlatform, 
                                 classIndex, minPlatformSize = minPlatformSize)

      #retained is now current + new retained
      id.retained = c(id.retained, MTunits[[nplatform]]$id$id.retained)
      
      retained = MTunits[[nplatform]]$id$id.retained %>% unique() %>% length()
      toprogress = MTunits[[nplatform]]$id$id.toprogress %>% unique() %>% length()
      notprocessed = MTunits[[nplatform]]$id$id.notprocessed %>% unique() %>% length()
      total = retained+toprogress+notprocessed
      processed=retained+toprogress
      
      if(verbose){
        #not enough samples - terminate early
        print(paste0("Size Check"))
        print(paste0(total - retained))
        print(paste0("    Total = ", total,  ""))
        print(paste0("    Processed =", processed, 
                     " (", retained, " retained, ",
                     toprogress, " to progress to next platform)"))
        print(paste0("    Not processed = ", notprocessed)) 
      }
    }
    
    MTlist[[nperm]] = MTunits
  }
  
  summary_results = MultiP.cost.summary(MTlist, myperms, platformList, confidenceCutoff, platformUnitCosts, platformUnitCosts2, properties, properties_direction, weights)
  
  return(MultiPobject = list(MTlist=MTlist, 
                           myperms=myperms, 
                           TSER.summary = summary_results$TSER.summary,
                           Results.TSA.retained.overall = summary_results$Results.TSA.retained.overall,
                           Results.stratification = summary_results$Results.stratification,
                           generatedModelsResults = generatedModels$resultList,
                           generatedModels = generatedModels$modelList,
                           generatedModelsFeatures = generatedModels$featureList,
                           MultiAssayExperiment = MultiAssayExperiment,
                           properties = properties, 
                           properties_direction = properties_direction, 
                           weights = weights))
}


MultiP.summary <- function(MTlist, myperms, platformList, confidenceCutoff){
  
  nplatforms=dim(myperms)[2]
  nperms=dim(myperms)[1]

  TSER.overall <- list()
  strat.overall <- list()
  id <- list()
  platformsequence <- vector()
  balancedAcc <- vector()
  TSERoverall.retained <- vector()
  TSERoverall.notretained <- vector()
  N.retained <- vector()
  N.notretained <- vector()
  
  for(nperm in 1:nperms){
    id.retained = NULL
    id.all = NULL
    id.notretained = NULL
    platformlabel = NULL
    TSER.retained = NULL
    TSER.notretained = NULL
    empty.platform = FALSE
    
    for(nplatform in 1:nplatforms){
      id.retained = c(id.retained,
                      MTlist[[nperm]][[nplatform]]$id$id.retained)
      id.retained = id.retained %>% unique()

      id.all = c(id.all,
                 MTlist[[nperm]][[nplatform]]$id$id.retained,
                 MTlist[[nperm]][[nplatform]]$id$id.toprogress,
                 MTlist[[nperm]][[nplatform]]$id$id.notprocessed)
      id.all = id.all %>% unique()
      
      id.notretained = setdiff(id.all, id.retained)
      
      if(!empty.platform){
        platformlabel = c(platformlabel, platformList[[myperms[nperm, nplatform]]])
      }
      
      if(length(id.notretained) == 0){
        empty.platform = TRUE
      }

      TSER.retained = dplyr::bind_rows(
        TSER.retained,
        MTlist[[nperm]][[nplatform]]$SSER$SSER %>% 
          dplyr::filter(SampleID %in% id.retained)
      )
      TSER.notretained = dplyr::bind_rows(
        TSER.notretained,
        MTlist[[nperm]][[nplatform]]$SSER$SSER  %>%
          dplyr::filter(SampleID %in% id.notretained)
      )
    }

    platformlabel = paste(platformlabel, collapse = '-')
    
    id[[nperm]] = list(id.retained=id.retained, 
                       id.notretained=id.notretained)
    platformsequence[[nperm]] = platformlabel

    TSER.overall[[nperm]] = dplyr::bind_rows(TSER.retained %>% 
                                              dplyr::summarise(tser=(1-mean(sser, na.rm=T))) %>% 
                                              dplyr::mutate(strata=factor("Retained", levels=c("Retained", "Not retained"))),
                                            TSER.notretained %>% 
                                              dplyr::summarise(tser=(1-mean(sser, na.rm=T))) %>% 
                                              dplyr::mutate(strata=factor("Not retained", levels=c("Retained", "Not retained"))))
    
    n.retained=0
    lvls = NULL
    stra = NULL
    
    for(nplatform in 1:nplatforms){
      n.from=n.retained + 1
      n.end = n.retained + nrow(MTlist[[nperm]][[nplatform]]$Stratification)
      
      if(nrow(MTlist[[nperm]][[nplatform]]$Stratification) == 0){
        if(verbose){print("Empty, skipping platform")}
        next
      }
      
      stra.MTunit <- MTlist[[nperm]][[nplatform]]$Stratification %>% 
        dplyr::arrange(strata, sser) %>% 
        dplyr::mutate(platform=platformList[[myperms[nperm, nplatform]]]) %>% 
        tibble::tibble(., ID = n.from:n.end)
      
      n.retained = n.retained +
        MTlist[[nperm]][[nplatform]]$Stratification %>% 
        dplyr::filter(strata=="Retained") %>% nrow()
      
      stra = dplyr::bind_rows(stra,  stra.MTunit)
      
      lvls = c(lvls,
               paste0("Retained (", platformList[[myperms[nperm, nplatform]]], ")"), 
               paste0("To progress (", platformList[[myperms[nperm, nplatform]]], ")"),
               paste0("Not processed (", platformList[[myperms[nperm, nplatform]]], ")"))
    }
    
    stra = stra %>% 
      dplyr::mutate(
        platformstrata=paste0(strata, " (", platform, ")"),
        platformstrata=factor(platformstrata, levels=lvls)) %>% 
      dplyr::arrange(strata, ID, platformstrata, sser)

    strat.overall[[nperm]] = stra
    
    TSERoverall.retained[[nperm]] = TSER.overall[[nperm]] %>% 
      dplyr::filter(strata=="Retained") %>% 
      dplyr::group_by(strata) %>% 
      dplyr::summarise(mean=mean(tser, na.rm=T)) %$% mean
    
    TSERoverall.notretained[[nperm]] = TSER.overall[[nperm]] %>% 
      dplyr::filter(strata=="Not retained") %>% 
      dplyr::group_by(strata) %>% 
      dplyr::summarise(mean=mean(tser, na.rm=T)) %$% mean
    N.retained[[nperm]] = id[[nperm]]$id.retained %>% length()
    N.notretained[[nperm]] = id[[nperm]]$id.notretained %>% length()
    
    stra_final = stra %>% dplyr::filter(strata=="Retained")
    balancedAcc[[nperm]] = balanced_accuracy(stra_final$predict_class, stra_final$true_class)
  }
  
  TSER.summary <- tibble::tibble(Platform.Sequence = platformsequence,
                                 balancedAcc = balancedAcc,
                                 TSA.retained = TSERoverall.retained,
                                 N.retained = N.retained,
                                 TSA.notretained = TSERoverall.notretained,
                                 N.notretained = N.notretained,
                                 Threshold = confidenceCutoff)
  
  return(list(TSER.summary = TSER.summary,
              IDs = id,
              Results.TSER.overall = TSER.overall,
              Results.stratification = strat.overall))
}


MultiP.cost.summary <- function(MTlist, myperms, platformList, confidenceCutoff, platformUnitCosts, platformUnitCosts2, properties, properties_direction, weights){

  if(length(platformUnitCosts2) == 1){
    platformUnitCosts2 = rep(NA, length(platformList))
  }
  
  names(platformUnitCosts) = platformList
  names(platformUnitCosts2) = platformList
  
  MultiPsummary = MultiP.summary(MTlist, myperms, platformList, confidenceCutoff)
  strat.overall = MultiPsummary$Results.stratification

  nplatforms = length(platformList)
  nperms = nrow(MultiPsummary$TSER.summary)
  Costs = NULL
  
  if(length(platformUnitCosts)!=nplatforms) return(print("length of platformUnitCost not equal to number of platforms"))
  
  Cost.byPlatform = vector()
  Cost.Total = vector()
  Cost.Total_str = vector()
  Cost.Total2 = vector()
  
  for(nperm in 1:nperms){
    platformseq = MultiPsummary$TSER.summary$Platform.Sequence[nperm] %>% str_split("-") %>% as_vector()

    costs = strat.overall[[nperm]] %>% 
      dplyr::filter(!strata=="Not processed") %>% 
      dplyr::count(platform) %>% 
      dplyr::mutate(cost=n*platformUnitCosts[platform],
                    cost2=n*platformUnitCosts2[platform],
                    platform=factor(platform, levels=platformseq)) %>% 
      dplyr::arrange(platform)
    
    Cost.byPlatform = c(Cost.byPlatform, str_c(paste0("$", costs$cost), collapse="-"))
    Cost.Total = c(Cost.Total, sum(costs$cost, na.rm=T))
    Cost.Total2 = c(Cost.Total2, sum(costs$cost2, na.rm=T))
    Cost.Total_str = c(Cost.Total_str, str_c(paste0("$", sum(costs$cost, na.rm=T))))
  }

  TSER.summary = MultiPsummary$TSER.summary %>% 
    dplyr::mutate(N.Total = N.retained + N.notretained,
                  Prop.notretained = N.notretained/N.Total) %>% 
    dplyr::select(Platform.Sequence, balancedAcc,
                  TSA.retained, TSA.notretained,
                  N.Total, N.retained, N.notretained, 
                  Prop.notretained, Threshold) %>% 
    dplyr::mutate(Cost.byPlatform = Cost.byPlatform,
                  Cost.Total = Cost.Total,
                  Cost.Total_str = Cost.Total_str,
                  Cost.Total2 = Cost.Total2)
  
  return(list(TSER.summary = TSER.summary,
              IDs = MultiPsummary$IDs,
              Results.TSA.retained.overall = MultiPsummary$Results.TSER.overall,
              Results.stratification = MultiPsummary$Results.stratification,
              properties = properties, 
              properties_direction = properties_direction, 
              weights = weights))
}

MTblock <- function(data=NULL, 
                    generatedModelsResults = NULL,
                    confidenceCutoff=0.2, platform = NULL, plotlabel = "",
                    classes = NULL, crossValParams = NULL, modellingParams = NULL, characteristics = NULL,
                    performanceType = "Sample Error", seed=1, verbose=F, finalPlatform,  classIndex, minPlatformSize = 10){
  
  SSER <- modelsPredict(data,
                        generatedModelsResults,
                        platform,
                        classes, 
                        minPlatformSize,
                        verbose)
  
  #platform specific error rate
  Stratification = tser(SSER, sser, confidenceCutoff, finalPlatform, minPlatformSize)

  id.retained = Stratification %>% 
    dplyr::filter(strata=="Retained") %>% 
    dplyr::select(SampleID) %>% 
    as_vector()
  
  id.toprogress = Stratification %>% 
    dplyr::filter(strata=="To progress") %>% 
    dplyr::select(SampleID) %>% 
    as_vector()
  
  id.notprocessed = Stratification %>% 
    dplyr::filter(strata=="Not processed") %>% 
    as_vector()

  names(id.retained) <- NULL
  names(id.toprogress) <- NULL
  names(id.notprocessed) <- NULL
  
  #Update to user defined threshold
  if(length(id.toprogress) < minPlatformSize){
    print("Not enough samples to progress, classifying at current layer!")
    id.retained = c(id.retained,id.toprogress)
    id.toprogress = NULL
  }
  
  if(finalPlatform == TRUE){
    id.retained = c(id.retained,id.toprogress)
    id.toprogress = NULL
  }
  
  return(list(SSER=SSER,
              Stratification=Stratification,
              id=list(id.retained=id.retained,
                      id.toprogress=id.toprogress,
                      id.notprocessed=id.notprocessed))
  )
}

tser <- function(SSER=NULL, col_name=sser, mycutoff=0.5, finalPlatform, minPlatformSize=10){
  Removed = as.data.frame(SSER$removed)
  SSER=SSER$SSER
  
  #If this platform was processed
  if(nrow(SSER)>0){
    if(!finalPlatform){
      id.retained = SSER %>% 
        dplyr::filter({{col_name}}<=mycutoff | {{col_name}}>= 1-mycutoff) %>% 
        dplyr::select(SampleID) %>%
        as_vector() %>% unique()      
    } else{
      id.retained = SSER %>% 
        dplyr::select(SampleID) %>%
        as_vector() %>% unique()
    }
    
    id.toprogress = setdiff(SSER %>% 
                              dplyr::select(SampleID) %>% 
                              as_vector(), 
                            id.retained)
    
    if(length(id.toprogress) < minPlatformSize){
      id.retained = c(id.retained,id.toprogress)
      id.toprogress = NULL
    }
    
    n = nrow(SSER)
    n_retained = length(id.retained)
    n_progress = n-n_retained
    
    strata = dplyr::bind_rows(
      SSER %>% 
        dplyr::filter(SampleID %in% id.retained) %>% 
        dplyr::select(SampleID, sser, predict_class, true_class) %>% 
        dplyr::mutate(strata = "Retained"),
      SSER %>% 
        dplyr::filter(SampleID %in% id.toprogress) %>% 
        dplyr::select(SampleID, sser, predict_class, true_class) %>% 
        dplyr::mutate(strata = "To progress"),
      Removed %>%
        dplyr::mutate(strata = "Not processed")
    ) %>%
      dplyr::mutate(strata=factor(strata, levels=c("Retained", "To progress", "Not processed")))

    #Else empty results object
  }else{
    strata = tibble(SampleID=character(),sser=double(), predict_class = factor(), true_class = factor(), strata=factor(),Removed=character())
  }
  
  return(strata)
}
