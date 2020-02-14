context("iregnet prediction")
library(iregnet)

load("data/neuroblastoma.Rdata")
### load("data/neuroblastoma.Rdata")
samp.data <- cbind(neuroblastoma.data$inputs , neuroblastoma.data$outputs)
samp.data <- matrix( samp.data , nrow(samp.data) , ncol(samp.data))

## randomize data
samp.data <- samp.data[ sample(nrow(samp.data)),]

#put into folds (5)
fold.size <- nrow(samp.data) / 5
samp.fold.1 <- samp.data[c(1: 93),]
samp.fold.2 <- samp.data[c(94 : 186),]
samp.fold.3 <- samp.data[c(187 : 279),]
samp.fold.4 <- samp.data[c(280 : 372),]
samp.fold.5 <- samp.data[c(373: 465),]

### input
samp.fold.1.X <- samp.fold.1[, c(1 : 29)] 
samp.fold.2.X <- samp.fold.2[, c(1 : 29)] 
samp.fold.3.X <- samp.fold.3[, c(1 : 29)] 
samp.fold.4.X <- samp.fold.4[, c(1 : 29)] 
samp.fold.5.X <- samp.fold.5[, c(1 : 29)] 

### output
samp.fold.1.y <- samp.fold.1[, c(30 : 31)] 
samp.fold.2.y <- samp.fold.2[, c(30 : 31)] 
samp.fold.3.y <- samp.fold.3[, c(30 : 31)] 
samp.fold.4.y <- samp.fold.4[, c(30 : 31)] 
samp.fold.5.y <- samp.fold.5[, c(30 : 31)]

## test all dataset folds for error and place as test case

 # test fold = 5
X.input <- rbind(samp.fold.1.X , samp.fold.2.X)
X.input <- rbind( X.input, samp.fold.3.X)
X.input <- rbind( X.input, samp.fold.4.X)


Y <- rbind(samp.fold.1.y , samp.fold.2.y)
Y <- rbind(Y , samp.fold.3.y)
Y <- rbind(Y , samp.fold.4.y)

### Failed
iregnet( X.input , Y , family = "gaussian" )



# test fold = 4
X.input <- rbind(samp.fold.1.X , samp.fold.2.X)
X.input <- rbind( X.input, samp.fold.3.X)
X.input <- rbind( X.input, samp.fold.5.X)


Y <- rbind(samp.fold.1.y , samp.fold.2.y)
Y <- rbind(Y , samp.fold.3.y)
Y <- rbind(Y , samp.fold.5.y)

### Failed
iregnet( X.input , Y , family = "gaussian" )



# test fold = 3
X.input <- rbind(samp.fold.1.X , samp.fold.2.X)
X.input <- rbind( X.input, samp.fold.4.X)
X.input <- rbind( X.input, samp.fold.5.X)


Y <- rbind(samp.fold.1.y , samp.fold.2.y)
Y <- rbind(Y , samp.fold.4.y)
Y <- rbind(Y , samp.fold.5.y)


### Failed
iregnet( X.input , Y , family = "gaussian" )



# test fold = 2
X.input <- rbind(samp.fold.1.X , samp.fold.3.X)
X.input <- rbind( X.input, samp.fold.4.X)
X.input <- rbind( X.input, samp.fold.5.X)


Y <- rbind(samp.fold.1.y , samp.fold.3.y)
Y <- rbind(Y , samp.fold.4.y)
Y <- rbind(Y , samp.fold.5.y)

### Failed
iregnet( X.input , Y , family = "gaussian" )



# test fold = 1
X.input <- rbind(samp.fold.2.X , samp.fold.3.X)
X.input <- rbind( X.input, samp.fold.4.X)
X.input <- rbind( X.input, samp.fold.5.X)


Y <- rbind(samp.fold.1.y , samp.fold.3.y)
Y <- rbind(Y , samp.fold.4.y)
Y <- rbind(Y , samp.fold.5.y)

### Failed
iregnet( X.input , Y , family = "gaussian" )




###################

set.seed(1)
### FAILED
iregnet(X.train , Y.train , family= "gaussian" , scale_init= NA , estimate_scale = TRUE)
iregnet(X.train , Y.train , family= "gaussian" , scale_init= 1 , estimate_scale = FALSE)
iregnet(X.train , Y.train , family= "logistic" , scale_init= NA , estimate_scale = TRUE)
iregnet(X.train , Y.train , family= "logistic" , scale_init= 1 , estimate_scale = FALSE)

### split data into train and valid and test on different folds


