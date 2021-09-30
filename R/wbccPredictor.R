##' Train the cell-type composition predictor using partial least squares regression on mehtylation values
##'
##' @title Train a predictor for cell-type composition
##' @param data Matrix with DNA methylation beta-values
##' @param covar Matrix of covariates to correct for in the model
##' @param cellPer Matrix of cell percentages on which the predictor will be trained
##' @param model Formula (default: cellPer ~ covar + data)
##' @param ncomp Number of PLS components to train (default: 50)
##' @param keep.model Logical specifying whether to return the model (default: FALSE)
##' @param ... Additional parameters for plsr see ?plsr
##' @return Prediction model PLSR object
##' @author mvaniterson
##' @export
##' @importFrom pls plsr

train_wbcc <- function(data, covar, cellPer, model=formula(cellPer ~ covar + data), 
                       ncomp = 50, keep.model = FALSE, ...){
  
  if(is.null(data) | is.null(covar) | is.null(cellPer))
    stop("Data, covariates, and cell percentages must be provided.")
  if(ncol(data) != nrow(covar) | nrow(covar) != nrow(cellPer))
    stop("Data column number must match covariate and cell percentages row number.")
  if(!all.equal(colnames(data), rownames(covar)) | !all.equal(rownames(covar), rownames(cellPer)))
    stop("Data column names must match covariate and cell percentages row names.")
  
  matchID <- match(rownames(covar), colnames(data))
  if(any(is.na(matchID)))
    stop("Data column names must match covariate row names.")
  covar <- covar[matchID, ]
  matchID <- match(rownames(covar), colnames(data))
  if(any(is.na(matchID)))
    stop("Data column names must match cell percentage row names.")
  cellPer <- cellPer[matchID, ]
  
  if(any(is.na(sum(data))) | any(is.na(sum(covar))) | any(is.na(sum(cellPer))))
    stop("Missing values are not allowed when training the predictor.")
  
  ## Model training using PLSR
  predictor <- plsr(model, ncomp=ncomp, data=list(cellPer = cellPer, 
                                                  covar = covar, 
                                                  data = t(data)), ...)
  
  ## Remove model if not being kept
  if(!keep.model) {
    predictor$model <- NULL
  }
  invisible(predictor)
}

##' Predict cell percentages based on a matrix or mvr-based model
##'
##' @title Predict cell-type composition from methylation values
##' @param pred A matrix or mvr predictor, trained using train_wbcc
##' @param data Matrix of DNA methylation beta values 
##' @param covar Matrix of covariates to correct for (must match those in predictor)
##' @param transformation Transformation to apply to predicted values (default: no transformation)
##' @param ncomp Optimal number of components 
##' @param impute Whether to impute missing values (default: TRUE)
##' @param ... Additional parameters for plsr see ?plsr
##' @return Predicted cell percentages
##' @author mvaniterson, lsinke
##' @export
##' @import pls
##' @importFrom stats coef median formula

predict_wbcc <- function(pred, data, covar, transformation=function(x) x, ncomp=NULL, impute=TRUE, ...) {
  
  if(is.null(data) | is.null(covar))               # Check if data and covariates provided
    stop("Both data and covariates must be provided.")
  if(ncol(data) != nrow(covar))                    # Check dimensions
    stop("Data column number must match covariate row number.")
  if(!isTRUE(all.equal(colnames(data), rownames(covar))))  # Check same names
    stop("Data column names must match covariate row names.")
  
  matchID <- match(rownames(covar), colnames(data))
  if(any(is.na(matchID)))
    stop("Data column names must match covariate row names.")
  covar <- covar[matchID, ]
  
  if(class(pred) == "mvr")
    names <- dimnames(coef(pred))[[1]]
  else if(class(pred) == "matrix")
    names <- rownames(pred)
  else
    stop(paste("This function is not designed for a", class(pred), "class of predictor."))
  
  # covaNames <- gsub("covar", "", grep("covar", names, value=TRUE))
  # dataNames <- gsub("data", "", grep("data", names, value=TRUE))
  # 
  # matchID <- match(covaNames, colnames(covar))
  # if(any(is.na(matchID)))
  #   stop("Covariates in the same do not match those in the predictor.")
  # covar <- covar[ , matchID]
  
  if(any(is.na(covar)) & !impute) {
    stop("Missing values are not allowed in the covariates if imputation is not specified.")
  } 
  else if(any(is.na(covar)) & impute) {
    print(paste("There are", sum(is.na(covar)), "NA's in the covariate matrix.",
                "These will be median imputed."))
    covar <- apply(covar, 2, function(x) {
      x[is.na(x)] = median(x, na.rm=TRUE)
      x})
  }
  
  # matchID <- match(dataNames, rownames(data))
  # if(any(is.na(matchID)))
  #   warning("Row names of the sample do not match those of the predictor.")
  # data <- data[matchID, ]
  if(any(is.na(data)) & !impute) {
    stop("Missing values are not allowed in the data if imputation is not specified")
  } 
  else if(any(is.na(data)) & impute) {
    print(paste("There are", sum(is.na(data)), "NA's in the data matrix.",
                "These will be median imputed."))
    nas <- apply(data, 1, function(x) any(is.na(x)))
    data[nas,] <- apply(data[nas, ], 1, function(x) median(x, na.rm=TRUE)) 
    data[is.na(data)] <- median(data, na.rm=TRUE) 
  }
  
  # Prediction
  if(class(pred) == "mvr") {
    predicted <- pls:::predict.mvr(pred, newdata = list(covar = covar, data = t(data)), ncomp=ncomp, ...)
    predicted <- predicted[ , , 1]
  }
  else if(class(pred) == "matrix") {
    predicted <- cbind(1, covar, t(data)) %*% pred
  }
  
  invisible(transformation(predicted))
}

##' Plots to validate predictor a predictor trained using the train_wbcc function
##'
##' @title Validation plots
##' @param measured Measured cell percentages
##' @param predicted Cell percentages predicted by the model
##' @param ... Additional parameters for plsr see ?plsr
##' @return Correlation plots for measured and predicted cell percentages
##' @author ljsinke
##' @export
##' @importFrom ggplot2 ggplot aes geom_point geom_smooth facet_wrap
##' @importFrom reshape2 melt
##' @importFrom stats cor

plot_wbcc <- function(measured, predicted, ...)
{ 
  corrMat <- round(cor(measured,predicted),4)
  
  for (k in 1:ncol(predicted)) {
    colnames(predicted)[k] <- paste(colnames(predicted)[k], " (correlation: ", corrMat[k,k], ")", sep="")
  }
  
  predicted <- melt(predicted)
  measured <- melt(measured)
  
  ggFrame <- data.frame(predicted = predicted[, 3], measured = measured[, 3], type = predicted[, 2])
  
  ggplot(data = ggFrame, mapping = aes(x = measured, y = predicted)) + 
    geom_point(shape=1) + 
    geom_smooth(method='lm',formula=y~x, se=FALSE, color="#D67D87", size=0.5) + 
    facet_wrap(~type, scales="free")
}


