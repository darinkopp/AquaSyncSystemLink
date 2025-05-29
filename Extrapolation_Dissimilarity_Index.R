# Checking for extrapolation

library(readxl)
library(parallel)
library(Rfast)

# Data files downloaded from https://seafile.rlp.net/d/789f89a8d7a7402eb73a/
setwd("Extrapolation")

PredImp <- read_xlsx("predictor_importance.xlsx")
TrainData <- readRDS("actual_sites.rds")
NewData <- readRDS("theoretical_sites.rds")

# Remove NA - predictions were not made at these locations
TrainData <- TrainData[complete.cases(TrainData),]

#NewData <- NewData[complete.cases(NewData),]

#insecticides

PredImpFilter <- PredImp[PredImp$class == "insecticides" & PredImp$type == "conc",]

# predictors that were the same for training data and new data 
# in importance are not considered.
# PredImpFilter$pred[!PredImpFilter$pred%in%colnames(TrainData)]

preds <- intersect(PredImpFilter$pred, colnames(TrainData))
PredImpFilter <- PredImpFilter[PredImpFilter$pred %in% preds,]
PredImpFilter <- setNames(PredImpFilter$V1, PredImpFilter$pred)
Importance <- PredImpFilter

# function to scale and weight predictor variables
Scale_Weight <- function(TrainData, Importance){
  
  #scale parameters predictor variables
  ScaleParms <- apply(TrainData[, names(Importance)], 2, 
                      function(x) c(center = mean(x), 
                                    scale = sd(x)))
  #scale predictor variables
  TrainDataScaled <- sapply(
    colnames(ScaleParms), function(x) 
      scale(TrainData[,x],
            center = ScaleParms["center",x],
            scale = ScaleParms["scale",x]))
  
  # apply non-standardized imporance as weight
  TrainDataScaledWeighted <- t(apply(
    TrainDataScaled[,names(Importance)], 1,
    function(x) x*Importance))
  
  # calculate mean distance among all training data in 
  # predictor space
  # Rfast::Dist(df, method = "euclidean")
  meanDist <- mean(dist(TrainDataScaledWeighted))
  
  return(list(ScaleParms = ScaleParms, 
              TrainDataScaledWeighted = TrainDataScaledWeighted, 
              MeanDist = meanDist))
}

# function to calculate dissimilarity index
Test4Extrapolation <- function(TrainDataWeightedList, NewDataPt, Importance){
  
  # a reminder that the scale weight function should be run first
  if(!is.list(TrainDataWeightedList)){
    stop("TrainDataWeighted needs to be names list")
  }
  
  ScaleParms <- TrainDataWeightedList[["ScaleParms"]]
  TrainDataScaledWeighted <- TrainDataWeightedList[["TrainDataScaledWeighted"]]
  MeanDist <- TrainDataWeightedList[["MeanDist"]]
    
  #apply scale parameters to new data
  NewDataScaled <- sapply(
    1:ncol(ScaleParms),
    function(x) 
      scale(NewDataPt[x],
            center = ScaleParms["center", x],
            scale = ScaleParms["scale", x]))
  
  # weight new data as importance 
  NewDataScaledWeighted <- NewDataScaled*Importance
  
  # calculate distance between new data point and all training data
  d2 <- sapply(1:nrow(TrainDataScaledWeighted), function(n){
    #n<-1
    df <- rbind(NewDataScaledWeighted, TrainDataScaledWeighted[n,])
    #d <- vegan::vegdist(df, method = "euclidean")
    #d <- stats::dist(df)
    # Rfast runs is much faster than other options
    d <- max(Rfast::Dist(df, method = "euclidean"))
    return(d)
  })
  
  # calculate the dissimilarity index 
  DI <- data.frame(MinDist=min(d2), DI = min(d2)/MeanDist)
  
  return(DI)
}

# scale and weight training data
TrainList <- Scale_Weight(TrainData, Importance)

# calculate dissimilarity index for new data
# estimated time for 8million records ~4days
cl <- makeCluster(detectCores())
clusterExport(cl, c("TrainList", "Importance", "Test4Extrapolation"))
t1 <- Sys.time()
Results <- parApply(cl, NewData[, names(Importance)], 1, function(x){
  Test4Extrapolation(TrainDataWeightedList = TrainList, 
                     NewDataPt = x, 
                     Importance = Importance)
})
t2 <- Sys.time()
stopCluster(cl)
t2-t1

#saveRDS(Results,"Results_1mil.rds")
#saveRDS(Results,"Results_7mil.rds")

saveRDS(Results, "Results_insecticides_Conc.rds")

DIndex <- do.call(rbind, Results)
# rows are organized the same as input New Data File
# PredImp[PredImp$class == "insecticides" & PredImp$type == "conc",]
#write.csv(DIndex, "DissimilarityIndex_insecticides_conc.csv")

# rows are organized the same as input New Data File
# PredImp[PredImp$class == "herbicides" & PredImp$type == "conc",]
#write.csv(DIndex, "DissimilarityIndex_herbicides_conc.csv")

# rows are organized the same as input New Data File
# PredImp[PredImp$class == "metals" & PredImp$type == "conc",]
#write.csv(DIndex, "DissimilarityIndex_metals_conc.csv")

#file complete 3/7/2025
#PredImp[PredImp$class == "fungicides" & PredImp$type == "conc",]
#write.csv(DIndex, "DissimilarityIndex_fungicides_conc.csv")

rm(a)
rm(b)
rm(DIndex)
gc()

