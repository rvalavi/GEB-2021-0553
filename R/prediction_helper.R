#' Orthogonal quadratic polynomials for glmnet
#'
#' A function to creat quadratic terms for glmnet functions i.e. lasso and ridge regression.
#' The output is an object of make_quadratic that can be used to predict on rasters and data.frames
#' for creating the quadratic terms.
#'
#' @param df a data.frame, typically the training data.
#' @param cols the name or index of the columns to be transformed. If NULL, all the columns will be transformed.
#' The factor columns won't be transfromed.
#'
#' @author Roozbeh Valavi
#'
#' @return an object of make_quadratic that can be used to predict on rasters and data.frames
#' @export
#'
#' @examples
make_quadratic <- function(df, cols = NULL, verbose = TRUE){
  if(is.null(cols)){
    cols <- colnames(df)
  }
  if(is.numeric(cols)){
    cols <- colnames(df)[cols]
  }
  # remove the factors
  if(any(sapply(df[,cols, drop = FALSE], is.factor))){
    if(verbose){
      message("The factor columns were removed form cols: ", cols[which(sapply(df[,cols, drop = FALSE], is.factor))])
    }
    cols <- cols[-which(sapply(df[,cols, drop = FALSE], is.factor))]
  }
  if(!all(is.element(cols, colnames(df)))){
    stop("The cols should be the same as the column names.")
  }
  xbar <- apply(df[,cols, drop = FALSE], 2, mean)
  x1 <- data.frame(mapply(`-`, df[,cols, drop = FALSE], xbar, SIMPLIFY = FALSE))
  alpha <- colSums(x1 ^ 3) / colSums(x1 ^ 2)
  # specify the output class
  finalList <- list(names = cols, xbars = xbar, alphas = alpha)
  class(finalList) <- c("make_quadratic")
  return(finalList)
}


#' @export
#' @method predict make_quadratic
predict.make_quadratic <- function(object, newdata, cols_from_obj = TRUE, ...){
  if(!methods::is(object, "make_quadratic"))
    stop("object should be a make_quadratic object.")
  # if(!all(object$names %in% names(newdata)))
  #   stop("The newdata does not have the same names as the object.")
  # ncl <- object$names
  # if(!all(names(newdata) %in% object$names))
  #   stop("The newdata does not have the same names as the object.")
  if(cols_from_obj){
    ncl <- object$names
  } else{
    ncl <- names(newdata)
  }
  if(methods::is(newdata, "Raster")){
    for(i in ncl){
      x1 <- newdata[[i]] - object$xbars[i]
      x2 <- (x1 ^ 2) - (object$alphas[i] * x1)
      if(raster::nlayers(newdata) > 1){
        newdata <- newdata[[-which(names(newdata) == i)]]
        newdata <- raster::stack(newdata, x1)
      } else{
        newdata <- x1
      }
      names(newdata)[raster::nlayers(newdata)] <- paste0(i, "_1")
      newdata <- raster::stack(newdata, x2)
      names(newdata)[raster::nlayers(newdata)] <- paste0(i, "_2")
    }
  } else if(methods::is(newdata, "data.frame")){
    for(i in ncl){
      x1 <- newdata[,i] - object$xbars[i]
      x2 <- x1 ^ 2 - object$alphas[i] * x1
      newdata[,ncol(newdata) + 1] <- x1
      names(newdata)[ncol(newdata)] <- paste0(i, "_1")
      newdata[,ncol(newdata) + 1] <- x2
      names(newdata)[ncol(newdata)] <- paste0(i, "_2")
      newdata <- newdata[,-which(names(newdata) == i)]
    }
  } else stop("newdata should be a raster or a data.frame.")
  return(newdata)
}

# extract fold ids for spatial tuning
foldsID <- function(fld){
  folds <- list()
  for(i in 1:5){
    all_data <- unlist(map(fld[-i], pluck, 1))
    fold <- fld[[i]][[1]]
    folds[[i]] <- unique(as.character(all_data[which(! all_data %in% fold)]))
  }
  return(folds)
}

# function for simultaneous tuning maxent regularisation multiplier and features
maxent_param <- function(data, kf = 5, folds = NULL, filepath){
  require(dismo)
  require(caret)
  require(precrec)
  if(is.null(folds)){
    folds <- caret::createFolds(y = as.factor(data$occ), k = kf)
  }
  ms <- c(0.5, 1, 2, 3, 4)
  grid <- expand.grid(regmult = paste0("betamultiplier=", ms),
                      splitrule = list(c("noautofeature", "nothreshold"), # LQHP
                                       c("noautofeature", "nothreshold", "noproduct"), # LQH
                                       c("noautofeature", "nothreshold", "nohinge", "noproduct"), # LQ
                                       c("noautofeature", "nothreshold", "nolinear", "noquadratic", "noproduct"), # H
                                       c("noautofeature", "nothreshold", "noquadratic", "nohinge", "noproduct")), # L
                      stringsAsFactors = FALSE)
  AUCs <- c()
  for(n in seq_along(grid[,1])){
    full_pred <- data.frame()
    for(k in seq_len(length(folds))){
      trainSet <- unlist(folds[-k])
      testSet <- unlist(folds[k])
      if(inherits(try(
        maxmod <- dismo::maxent(x = data[trainSet, 2:ncol(data)],
                             p = data[trainSet, ]$occ,
                             removeDuplicates = FALSE,
                             path = filepath,
                             args = as.character(unlist(grid[n, ]))
        )
      ), "try-error")){
        next
      }
      modpred <- predict(maxmod, data[testSet, 2:ncol(data)], args = "outputformat=cloglog")
      pred_df <- data.frame(score = modpred, label = data[testSet, ]$occ)
      full_pred <- rbind(full_pred, pred_df)
    }
    AUCs[n] <- precrec::auc(precrec::evalmod(scores = full_pred$score, 
                                             labels = full_pred$label))[1,4]
  }
  best_param <- as.character(unlist(grid[which.max(AUCs), ]))
  return(best_param)
}

# for easy saving the parameters
param_to_txt <- function(x){
  x <- x[-(1:3)]
  if(length(x) == 0){
    features <- "LQHP"
  } else if(all(c("noquadratic", "nohinge", "noproduct") %in% x)){
    features <- "L"
  } else if(all(c("nolinear", "noquadratic", "noproduct") %in% x)){
    features <- "H"
  } else if(all(c("nohinge", "noproduct") %in% x)){
    features <- "LQ"
  } else if("noproduct" == x){
    features <- "LQH"
  }
  return(features)
}
