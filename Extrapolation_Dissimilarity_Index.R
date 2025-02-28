# Checking for extrapolation

setwd("Raw_Data/Extrapolation")

library(readxl)
library(parallel)
library(Rfast)

PredImp <- read_xlsx("predictor_importance.xlsx")
TrainData <- readRDS("actual_sites.rds")
NewData <- readRDS("theoretical_sites.rds")

# need to ask why NAs
TrainData <- TrainData[complete.cases(TrainData),]
NewData <- NewData[complete.cases(NewData),]

PredImpFilter <- PredImp[PredImp$class == "fungicides" & PredImp$type == "conc",]

# predictors in importance are not included in the training data?
# PredImpFilter$pred[!PredImpFilter$pred%in%colnames(TrainData)]
preds <- intersect(PredImpFilter$pred, colnames(TrainData))
PredImpFilter <- PredImpFilter[PredImpFilter$pred %in% preds,]
PredImpFilter <- setNames(PredImpFilter$V1,PredImpFilter$pred)
Importance <- PredImpFilter

# function to scale and weight predictor variables
Scale_Weight <- function(TrainData, Importance){
  
  #scale predictor variables
  ScaleParms <- apply(TrainData[, names(Importance)], 2, 
                      function(x) c(center = mean(x), 
                                    scale = sd(x)))
  
  TrainDataScaled <- sapply(
    colnames(ScaleParms), function(x) 
      scale(TrainData[,x],
            center = ScaleParms["center",x],
            scale = ScaleParms["scale",x]))
  
  TrainDataScaledWeighted <- t(apply(
    TrainDataScaled[,names(Importance)], 1,
    function(x) x*Importance))
  
  meanDist <- mean(vegan::vegdist(TrainDataScaledWeighted, method = "euclidean"))
  
  return(list(ScaleParms = ScaleParms, 
              TrainDataScaledWeighted = TrainDataScaledWeighted, 
              MeanDist = meanDist))
}

# function to calculate dissimilarity index
Test4Extrapolation <- function(TrainDataWeightedList, NewDataPt, Importance){
  
  if(!is.list(TrainDataWeightedList)){
    stop("TrainDataWeighted needs to be names list")
  }
  
  ScaleParms <- TrainDataWeightedList[["ScaleParms"]]
  TrainDataScaledWeighted <- TrainDataWeightedList[["TrainDataScaledWeighted"]]
  MeanDist <- TrainDataWeightedList[["MeanDist"]]
    
  NewDataScaled <- sapply(
    #colnames(ScaleParms),
    1:ncol(ScaleParms),
    function(x) 
      scale(NewDataPt[x],
            center = ScaleParms["center", x],
            scale = ScaleParms["scale", x]))
  
  
  NewDataScaledWeighted <- NewDataScaled*Importance
  t1 <- Sys.time()
  d2 <- sapply(1:nrow(TrainDataScaledWeighted), function(n){
    #n<-1
    df <- rbind(NewDataScaledWeighted, TrainDataScaledWeighted[n,])
    #d <- vegan::vegdist(df, method = "euclidean")
    #d <- stats::dist(df)
    d <- max(Rfast::Dist(df, method = "euclidean"))
    return(d)
  })
  t2 <- Sys.time()
  t2-t1
  
  DI <- min(d2)/MeanDist
  return(min(DI))
}

# scale and weight training data
TrainList <- Scale_Weight(TrainData, Importance)


cl <- makeCluster(detectCores())
clusterExport(cl, c("TrainList", "Importance", "Test4Extrapolation"))
t1 <- Sys.time()
Results <- parApply(cl, NewData[1:1000, names(Importance)], 1, function(x){
  Test4Extrapolation(TrainDataWeightedList = TrainList, 
                     NewDataPt = x, 
                     Importance = Importance)
})
t2 <- Sys.time()
stopCluster(cl)
t2-t1


#######
Test4Extrapolation(TrainDataWeightedList = TrainList, 
                   NewDataPt = x[1, names(Importance)], 
                   Importance = Importance)
install.packages("fastDist")
nrow(TrainList[["TrainDataScaledWeighted"]][,-1])

sapply(1:12, function(x){
  Test4Extrapolation(TrainDataWeightedList = TrainList, 
                     NewDataPt = NewData[x,names(Importance)], 
                     Importance = Importance)
  
})

apply(NewData[1:13, names(Importance)],1, function(x){
  Test4Extrapolation(TrainDataWeightedList = TrainList, 
                     NewDataPt = x, 
                     Importance = Importance)
})

NewDataPt<-NewData[1, names(Importance)]
for (i in 1:13){
  di <- Test4Extrapolation(TrainDataWeightedList = TrainList, 
                           NewDataPt = NewData[i,], 
                           Importance)  
}


#scale predictor variables
########################
ScaleParms <- apply(TrainData[, -1], 2, 
                    function(x) c(center = mean(x), scale = sd(x)))

TrainDataScaled <- sapply(
  colnames(ScaleParms), function(x) 
    scale(TrainData[,x],
          center = ScaleParms["center",x],
          scale = ScaleParms["scale",x]))

NewDataScaled <- sapply(
  colnames(ScaleParms), function(x) 
    scale(NewData[,x],
          center = ScaleParms["center",x],
          scale = ScaleParms["scale",x]))

PredImpFilter <- PredImp[PredImp$class=="fungicides",]
#predictors in importance are not included in the training data?
# PredImpFilter$pred[!PredImpFilter$pred%in%colnames(TrainDataScaled)]
preds <- intersect(PredImpFilter$pred,colnames(TrainDataScaled))
PredImpFilter <- PredImpFilter[PredImpFilter$pred%in%preds,]
PredImpFilter <- setNames(PredImpFilter$V1,PredImpFilter$pred)

TrainDataScaledWeighted <- t(apply(
  TrainDataScaled[,names(PredImpFilter)], 1,
  function(x) x*PredImpFilter))

NewDataScaledWeighted <- t(apply(
  NewDataScaled[,names(PredImpFilter)], 1,
  function(x) x*PredImpFilter))


for (i in 1:10){
  NewDataPt <- NewDataScaledWeighted[i,]
  d2 <- sapply(1:nrow(TrainDataScaledWeighted), function(n){
    df <- rbind(NewDataPt,
                TrainDataScaledWeighted[n,])
    d <- vegan::vegdist(df, method = "euclidean")
    return(d)
  })
  
  print(min(d2))
}
###########################################################
