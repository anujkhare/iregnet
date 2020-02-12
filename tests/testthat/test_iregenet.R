context("iregnet prediction")
library(iregnet)


testFold.dir.vec <- Sys.glob(file.path(
  "data", "*", "cv", "*", "testFolds",
  "*"))

n.pred.vec <- sapply(testFold.dir.vec, function(testFold.dir){
  path <- file.path(dirname(testFold.dir), "models", "*", "predictions.csv")
  length(Sys.glob(path))
})

testFold.dir <- testFold.dir.vec[[1]]
data.m.list <- list()

## format data

FormatData <- function( testFold.dir ){
  library(data.table)
  test.fold <- as.integer(basename(testFold.dir))
  cv.type.dir <- dirname(dirname(testFold.dir))
  data.dir <- dirname(dirname(cv.type.dir))
  folds.csv <- file.path(cv.type.dir, "folds.csv")
  folds.dt <- fread(folds.csv)
  data.list <- list()
  for(data.type in c("inputs", "outputs")){
    csv.xz <- file.path(data.dir, paste0(data.type, ".csv.xz"))
    dt <- fread(cmd=paste("xzcat", csv.xz))
    stopifnot(nrow(dt) == nrow(folds.dt))
    m <- as.matrix(dt[, -1, with=FALSE])
    rownames(m) <- dt$sequenceID
    data.list[[data.type]] <- m
  }
  rep.val.vec <- c(
    log.log.bases="log2.n",
    n.loglog="log2.n",
    "diff abs.identity.quantile.50%"="log.hall",
    log.sd="log.hall")
  for(old.name in names(rep.val.vec)){
    new.name <- rep.val.vec[[old.name]]
    colnames(data.list$inputs)[colnames(data.list$inputs) == old.name] <- new.name
  }
  keep.inputs <- apply(is.finite(data.list$inputs), 2, all)
  data.list$inputs <- data.list$inputs[, keep.inputs, drop=FALSE]
  id.list <- list(
    train=folds.dt[fold != test.fold, sequenceID],
    test=folds.dt[fold == test.fold, sequenceID])
  set.list <- list()
  for(set.name in names(id.list)){
    set.id.vec <- id.list[[set.name]]
    set.list[[set.name]] <- lapply(data.list, function(m){
      m[set.id.vec,]
    })
  }
  
  X.train <- matrix(set.list$train$inputs, nrow(set.list$train$inputs), ncol(set.list$train$inputs))
  Y.train <-  matrix( set.list$train$outputs , nrow(set.list$train$outputs) , ncol(set.list$train$outputs))
  
  data.vec <- list( X = X.train ,y = Y.train)
  return(data.vec)
}

## test all dataset folds for error and place as test case

## testFold.dir <- testFold.dir.vec[[1]]
data.m.list <- FormatData(testFold.dir.vec[[1]])
X.train <- data.m.list$X
Y.train <- data.m.list$y


set.seed(1)
### FAILED
cv.iregnet(X.train , Y.train , family= "gaussian" , scale_init= NA , estimate_scale = TRUE)
cv.iregnet(X.train , Y.train , family= "gaussian" , scale_init= 1 , estimate_scale = FALSE)
cv.iregnet(X.train , Y.train , family= "logistic" , scale_init= NA , estimate_scale = TRUE)
cv.iregnet(X.train , Y.train , family= "logistic" , scale_init= 1 , estimate_scale = FALSE)

### PASSED
#none


## testFold.dir <- testFold.dir.vec[[2]]
data.m.list <- FormatData(testFold.dir.vec[[2]])
X.train <- data.m.list$X
Y.train <- data.m.list$y

set.seed(1)
### FAILED
cv.iregnet(X.train , Y.train , family= "gaussian" , scale_init= NA , estimate_scale = TRUE)
cv.iregnet(X.train , Y.train , family= "gaussian" , scale_init= 1 , estimate_scale = FALSE)
cv.iregnet(X.train , Y.train , family= "logistic" , scale_init= NA , estimate_scale = TRUE)
cv.iregnet(X.train , Y.train , family= "logistic" , scale_init= 1 , estimate_scale = FALSE)

### PASSED


## testFold.dir <- testFold.dir.vec[[3]]
data.m.list <- FormatData(testFold.dir.vec[[3]])
X.train <- data.m.list$X
Y.train <- data.m.list$y

set.seed(1)
### FAILED
cv.iregnet(X.train , Y.train , family= "gaussian" , scale_init= NA , estimate_scale = TRUE)
cv.iregnet(X.train , Y.train , family= "gaussian" , scale_init= 1 , estimate_scale = FALSE)
cv.iregnet(X.train , Y.train , family= "logistic" , scale_init= NA , estimate_scale = TRUE)
cv.iregnet(X.train , Y.train , family= "logistic" , scale_init= 1 , estimate_scale = FALSE)

### PASSED



