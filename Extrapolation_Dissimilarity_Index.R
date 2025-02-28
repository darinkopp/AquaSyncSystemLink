# Checking for extrapolation

setwd("Extrapolation")

library(readxl)
library(parallel)
library(Rfast)

PredImp <- read_xlsx("predictor_importance.xlsx")
TrainData <- readRDS("actual_sites.rds")
NewData <- readRDS("theoretical_sites.rds")

# Remove NA - predictions were not made at these locations
TrainData <- TrainData[complete.cases(TrainData),]
NewData <- NewData[complete.cases(NewData),]

PredImpFilter <- PredImp[PredImp$class == "fungicides" & PredImp$type == "conc",]

# predictors that were the same for training data and new data 
# in importance are not considerd.
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
Results <- parApply(cl, NewData[1:1000000, names(Importance)], 1, function(x){
  Test4Extrapolation(TrainDataWeightedList = TrainList, 
                     NewDataPt = x, 
                     Importance = Importance)
})
t2 <- Sys.time()
stopCluster(cl)
t2-t1

