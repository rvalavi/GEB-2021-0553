# Author: Roozbeh Valavi
# contact: valavi.r@gmail.com
# Date : August 2020
# Version 0.1
# Licence GPL v3
# Spatial blocking for train-test datasets
spatial_block <- function(train_data, # training data - sf object
                          test_data, # testing data - sf object
                          response = "occ", # the response column
                          blocks, # spatial blocks - sf object
                          k = 5L, # number of folds
                          selection = "random", # random allocation
                          iteration = 100L, # find evenly distributed folds
                          seed = 32, # random seed for reproducibility
                          progress = FALSE, # show progress bar
                          verbose = TRUE){ # print information
  
  if(!is.element(selection, c("systematic", "random"))){
    stop("The selection argument must be 'random', 'systematic'.")
  }
  # always set numLimit to 0 for searching random folds
  numLimit <- 0
  if(methods::is(blocks, "SpatialPolygons")){
    blocks <- sf::st_as_sf(blocks)
  } else if(!methods::is(blocks, "sf")){
    stop("blocks, should be a spatial or sf object")
  }
  ## check if response is a col in train_data
  if(!is.null(response)){
    if(response %in% colnames(train_data) == FALSE){
      warning("There is no match between the column names in 'train_data' and 'response' argument (response variable).\n")
      response <- NULL
    }
  }
  # subset the blocks to train and test points
  subBlocks <- blocks[rbind(train_data, test_data), ]
  iteration <- as.integer(iteration)
  if(iteration < 1 || !is.numeric(iteration)){
    iteration <- 1L
    message("The interation has been set to 1! \n", "Iteration must be a positive integer value.\n")
  }
  if(progress==TRUE && numLimit == 0){
    pb <- progress::progress_bar$new(format = " Progress [:bar] :percent in :elapsed",
                                     total=iteration, clear=FALSE, width=75) # add progress bar
  }
  ## do the intersection once and outside of the loop  
  subBlocksDF <- as.data.frame(sf::st_intersects(sf::st_geometry(train_data), sf::st_geometry(subBlocks)))
  names(subBlocksDF) <- c("records", "blocks")
  subBlocksDF2 <- as.data.frame(sf::st_intersects(sf::st_geometry(test_data), sf::st_geometry(subBlocks)))
  names(subBlocksDF2) <- c("records", "blocks")
  # randomly remove the repeated records occurred on the edges of blocks
  if(nrow(subBlocksDF) > nrow(train_data)){
    if(!is.null(seed)){
      set.seed(seed)
    }
    subBlocksDF <- subBlocksDF[sample(nrow(subBlocksDF)), ]
    subBlocksDF <- subBlocksDF[!duplicated(subBlocksDF$records), ]
  } else if(nrow(subBlocksDF) < nrow(train_data) || anyNA(subBlocksDF)){
    nonoverlap <- nrow(train_data) - nrow(subBlocksDF)
    warning("At least ", nonoverlap, " of the points are not on the train spatial blocks")
  }
  if(nrow(subBlocksDF2) > nrow(test_data)){
    if(!is.null(seed)){
      set.seed(seed)
    }
    subBlocksDF2 <- subBlocksDF2[sample(nrow(subBlocksDF2)), ]
    subBlocksDF2 <- subBlocksDF2[!duplicated(subBlocksDF2$records), ]
  } else if(nrow(subBlocksDF2) < nrow(test_data) || anyNA(subBlocksDF2)){
    nonoverlap <- nrow(test_data) - nrow(subBlocksDF2)
    warning("At least ", nonoverlap, " of the points are not on the test spatial blocks")
  }
  nrowBlocks <- nrow(subBlocks)
  maxNumRecord <- 0
  maxSD <- Inf
  if(!is.null(seed)){
    set.seed(seed)
  }
  for(i in seq_len(iteration)){
    if(k > nrowBlocks){
      stop("'k' is bigger than the number of spatial blocks\n",
           "The number of spatial blocks is: ", nrowBlocks)
    } else if(k < 2){
      stop("'k' must be 2 or higher")
    }
    if(selection=='systematic'){
      foldDF <- data.frame(blocks = seq_len(nrowBlocks), folds = systematicNum(subBlocks, k))
      subBlocksDF <- merge(x = subBlocksDF, y = foldDF, by = "blocks", all.x = TRUE)
      subBlocksDF2 <- merge(x = subBlocksDF2, y = foldDF, by = "blocks", all.x = TRUE)
    } else if(selection=='random'){
      subBlocksDF <- subBlocksDF[, c("records", "blocks")] # to avoid repetition in iterations
      subBlocksDF2 <- subBlocksDF2[, c("records", "blocks")] # to avoid repetition in iterations
      foldDF <- data.frame(blocks = seq_len(nrowBlocks), folds = 0)
      # create random folds
      num <- floor(nrowBlocks / k)
      foldDF$folds[seq_len(num * k)] <- sample(rep(seq_len(k), num), num * k)
      if(nrowBlocks %% k != 0){
        rest <- nrowBlocks %% k
        unfold <- which(foldDF$folds==0)
        foldDF$folds[unfold] <- sample(seq_len(k), rest, replace = FALSE)
      }
      # merge the folds to blocks
      subBlocksDF <- merge(x = subBlocksDF, y = foldDF, by = "blocks", all.x = TRUE)
      subBlocksDF2 <- merge(x = subBlocksDF2, y = foldDF, by = "blocks", all.x = TRUE)
    }
    # create records table
    if(is.null(response)){
      trainTestTable <- data.frame(train=rep(0, k), test=0)
    } else{
      cl <- sort(unique(train_data[, response, drop = TRUE]))
      clen <- length(cl)
      trainTestTable <- as.data.frame(matrix(0, nrow = k, ncol = clen * 2))
      names(trainTestTable) <- c(paste("train", cl, sep = "_"), paste("test", cl, sep = "_"))
    }
    # count the number of points in each fold
    foldList <- list()
    for(p in seq_len(k)){
      trainSet <- subBlocksDF$records[which(subBlocksDF$folds != p)]
      testSet <- subBlocksDF2$records[which(subBlocksDF2$folds == p)]
      foldList[[p]] <- assign(paste0("fold", p), list(trainSet, testSet))
      if(is.null(response)){
        trainTestTable$train[p] <- length(trainSet)
        trainTestTable$test[p] <- length(testSet)
      } else{
        countrain <- table(train_data[trainSet, response, drop = TRUE])
        countest <- table(test_data[testSet, response, drop = TRUE])
        trainTestTable[p, which(cl %in% names(countrain))] <- countrain
        trainTestTable[p, clen + which(cl %in% names(countest))] <- countest
      }
    }
    if(selection == "random"){
      if(is.numeric(numLimit) && numLimit > 0){
        if(any(trainTestTable < numLimit)==FALSE){ # exit the loop if meet the limit number
          break
        }
      } else if(numLimit == 0){ # find the highest minimum number in the table and store relevant objects
        if(min(trainTestTable) >= maxNumRecord && stats::sd(unlist(trainTestTable)) < maxSD){
          trainTestTable2 <- trainTestTable
          maxNumRecord <- min(trainTestTable2)
          maxSD <- stats::sd(unlist(trainTestTable))
          subBlocksDFx <- subBlocksDF
          foldList2 <- foldList
          iter <- i
        }
      } else stop("numLimit argument should be a numeric value equal or hagher than 0 or be NULL")
      if(progress == TRUE && numLimit == 0){
        pb$tick() # update progress bar
      }
    } else{
      break
    }
  }
  if(numLimit == 0 && selection == "random"){ # return the best blocks, table etc.
    # subBlocksDF <- subBlocksDFx
    trainTestTable <- trainTestTable2
    foldList <- foldList2
    if(verbose) cat(paste0("The best folds was in iteration ", iter, ":\n"))
  }
  if(verbose) print(trainTestTable)
  if(any(trainTestTable <= numLimit)){
    zerofolds <- which(apply(trainTestTable, 1, function(x) any(x == numLimit)))
    if(length(zerofolds) > 1){
      warning("The folds ", paste(zerofolds, collapse = ", "), " have class(es) with ", numLimit, " (or less) records")
    } else{
      warning("The fold ", zerofolds, " has class(es) with ", numLimit, " (or less) records")
    }
  }

  theList <- list(folds = foldList,
                  k = k,
                  response = response,
                  records = trainTestTable)
  class(theList) <- c("SpatialBlock")
  return(theList)
}

