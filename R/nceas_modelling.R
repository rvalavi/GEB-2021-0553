suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(precrec))
suppressPackageStartupMessages(library(disdat))
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(prg))

suppressPackageStartupMessages(library(dismo))
suppressPackageStartupMessages(library(gbm))
suppressPackageStartupMessages(library(ranger))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(e1071))
suppressPackageStartupMessages(library(plotmo))
suppressPackageStartupMessages(library(earth))
suppressPackageStartupMessages(library(caret))

# list of covariates to use in modelling
myvars <- list(
  AWT = c("bc04",  "bc05",  "bc06",  "bc12",  "bc15",  "slope", "topo", "tri"),
  CAN = c("alt", "asp2", "ontprec", "ontslp", "onttemp", "ontveg", "watdist"), 
  NSW = c("cti", "disturb", "mi", "rainann", "raindq", "rugged", "soildepth", "soilfert", "solrad", "tempann", "topo"), # "vegsys"), 
  NZ = c("age", "deficit", "hillshade", "mas", "mat", "r2pet", "slope", "sseas", "toxicats", "tseas", "vpd"), 
  SA = c("sabio12", "sabio15", "sabio17", "sabio18", "sabio2", "sabio4", "sabio5", "sabio6"), 
  SWI = c("bcc", "calc", "ccc", "ddeg", "nutri", "pday", "precyy", "sfroyy", "slope", "sradyy", "swb", "topo")
)

# read required functions and variables
source("R/spatial_blocking.R")
source("R/nceas_model_vars.R")
source("R/prediction_helper.R")
# provide names for regions to be modelled - here we model all 6:
regions <- c("AWT", "CAN", "NSW", "NZ", "SA", "SWI")
#regions <- c("SWI")
# specify names of all categorical variables across all regions:
categoricalvars <- c("ontveg", "toxicats", "age", "calc") # "vegsys"
# background directory
bgdir <- "data/tgb_samples"


# output directory
outdir <- "output/nceas_model_output"

dir.create(file.path("output"))
dir.create(file.path(outdir))


# save MaxEnt parameters
tune_df <- data.frame(speice = NA, fold = NA, model = NA, reg = NA, args = NA)

n <- 0
z <- 0
for(r in regions){
  # reading presence-only and background species data for this region, one file per region:
  presences <- disPo(r)
  # extract names for all species
  species <- unique(presences$spid)
  
  print(paste("*****", r, "has", length(species), "species *****"))
  
  # now for each species, prepare data for modeling and testing 
  # with randomForests, model, predict to the testing file with 
  # environmental data. 
  for(s in species){
    # subset presence records of species for this species
    sp_presence <- presences[presences$spid == s, ]
    # find the group of the species
    grp <- sp_presence$group[1]
    # load the target-group-background data
    background <- read.csv(file.path(bgdir, r, grp, "tgbsample.csv"))
    
    # add background data
    pr_bg <- rbind(sp_presence, background)
    
    # find the testing and evaluation files – for some regions this means identifying the taxonomic group
    if (r %in% c("AWT", "NSW")) {
      testing <- disEnv(r, group = grp) # env value for prediction
      evals <- disPa(r, group = grp) # pa records for evaluation
    } else {
      testing <- disEnv(r)
      evals <- disPa(r)
    }
    
    # convert categorical vars to factor in both training and testing data. 
    # We use the package forcats to ensure that the levels of the factor in the testing data match 
    # those in the training data, regardless of whether all levels are present in the testing data. 
    for(i in 1:ncol(pr_bg)){
      if(colnames(pr_bg)[i] %in% categoricalvars){
        fac_col <- colnames(pr_bg)[i]
        pr_bg[ ,fac_col] <- as.factor(pr_bg[ ,fac_col])
        testing[ ,fac_col] <- as.factor(testing[ ,fac_col])
        all_levals <- unique(c(levels(pr_bg[ ,fac_col]), levels(testing[ ,fac_col])))
        pr_bg[ ,fac_col] <- forcats::fct_expand(pr_bg[, fac_col], all_levals)
        testing[ ,fac_col] <- forcats::fct_expand(testing[, fac_col], all_levals)
        testing[ ,fac_col] <- forcats::fct_relevel(testing[,fac_col], levels(pr_bg[,fac_col]))
      }
    }
    
    # rest rownames for spatial tuning
    rownames(pr_bg) <- NULL
    rownames(testing) <- NULL
    rownames(evals) <- NULL
    
    # extract the relevant columns for modelling
    training <- pr_bg[, c("occ", myvars[[r]])]
    testing <- testing[, myvars[[r]]]
    
    ## generating spatial and random cross-validation
    train_sf <- sf::st_as_sf(pr_bg[, c("siteid", "x", "y", "occ")], coords = 2:3, crs = disCRS(r))
    test_sf <- evals[, c("siteid", "x", "y", s)] %>% 
      setNames(c("siteid", "x", "y", "occ")) %>% 
      sf::st_as_sf(coords = 2:3, crs = disCRS(r))
    
    # read the block spatial layer
    # myblock <- st_read(sprintf("data/blocks/%s.gpkg", r), quiet = TRUE)
    myblock <- st_read(sprintf("data/blocks/%s.gpkg", r), quiet = TRUE)
    # generate spatial blocks
    spblock <- spatial_block(train_data = train_sf, 
                             test_data = test_sf, 
                             response = "occ", 
                             blocks = myblock, 
                             k = 5, 
                             seed = 32,
                             iteration = 100)
    
    # skip species with zero or low records
    if(min(spblock$records) < 1) next
    if(any(spblock$records$train_1 < 5)) next
    if(sum(spblock$records$test_1) < 10) next
    
    # subset folds
    spfolds <- spblock$folds
    
    # create random folds
    prTrain <- which(training$occ == 1)
    bgTrain <- which(training$occ == 0)
    prTest <- which(evals[,s] == 1)
    abTest <- which(evals[,s] == 0)
    
    myseed <- 33
    spname <- s

    for(f in seq_len(5)){
      # take the sample for random-cv equal to spatial-cv
      # train and test sets for spatial cv
      training_scv <- training[spfolds[[f]][[1]], ]
      testing_scv <- testing[spfolds[[f]][[2]], ]
      evals_scv <- evals[spfolds[[f]][[2]], ]
      # train and test sets for random cv
      tr_pr <- sample(prTrain, spblock$records$train_1[f], replace = FALSE) # pr in fold
      tr_bg <- sample(bgTrain, spblock$records$train_0[f], replace = FALSE) # bg in fold
      training_rcv <- training[c(tr_pr, tr_bg), ]
      ts_pr <- sample(prTest, spblock$records$test_1[f], replace = FALSE) # pr in fold
      ts_ab <- sample(abTest, spblock$records$test_0[f], replace = FALSE) # ab in fold
      testing_rcv <- testing[c(ts_pr, ts_ab), ]
      evals_rcv <- evals[c(ts_pr, ts_ab), ]
      # update test pr and ab and remove the previous folds
      prTest <- prTest[!prTest %in% ts_pr]
      abTest <- abTest[!abTest %in% ts_ab]
      
      # start modelling!
      ##*******************************************##
      ##* regression methods
      prNum_scv <- as.numeric(table(training_scv$occ)["1"])
      bgNum_scv <- as.numeric(table(training_scv$occ)["0"])
      
      prNum_rcv <- as.numeric(table(training_rcv$occ)["1"])
      bgNum_rcv <- as.numeric(table(training_rcv$occ)["0"])
      
      wt_scv <- ifelse(training_scv$occ == 1, 1, prNum_scv/bgNum_scv) # down-weighting
      wt_rcv <- ifelse(training_rcv$occ == 1, 1, prNum_rcv/bgNum_rcv) # down-weighting
      #*******************************************##
      ##* BRT
      brt_scv <- NULL
      b <- 0
      set.seed(myseed)
      while(is.null(brt_scv)){
        b <- b + 1
        if(b < 11){
          ntrees <- 50
          lrate <- 0.001
        } else if(b < 21){
          lrate <- 0.0001
        } else if(b < 31){
          ntrees <- 25
          lrate <- 0.0001
        } else{
          break
        }
        ptm <- proc.time()
        brt_scv <-  dismo::gbm.step(data = training_scv,
                                    gbm.x = 2:ncol(training_scv), # column indices for covariates
                                    gbm.y = 1, # column index for response
                                    family = "bernoulli",
                                    tree.complexity = ifelse(sum(training_scv$occ) < 50, 1, 5),
                                    learning.rate = lrate,
                                    bag.fraction = 0.75,
                                    max.trees = 10000,
                                    n.trees = ntrees,
                                    n.folds = 5, # 5-fold cross-validation
                                    site.weights = wt_scv,
                                    silent = TRUE)
        t <- proc.time() - ptm
      }
      if(is.null(brt_scv)) next
      pred_brt_scv <- predict(brt_scv, testing_scv, 
                              n.trees = brt_scv$gbm.call$best.trees, 
                              type = "response")
      
      df_scv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_scv,
                           model = "BRT",
                           score = pred_brt_scv,
                           label = evals_scv[, s],
                           cv = "spatial",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      
      ##--------------------------------##
      brt_rcv <- NULL
      b <- 0
      set.seed(myseed)
      while(is.null(brt_rcv)){
        b <- b + 1
        if(b < 11){
          ntrees <- 50
          lrate <- 0.001
        } else if(b < 21){
          lrate <- 0.0001
        } else if(b < 31){
          ntrees <- 25
          lrate <- 0.0001
        } else{
          break
        }
        ptm <- proc.time()
        brt_rcv <- dismo::gbm.step(data = training_rcv,
                                   gbm.x = 2:ncol(training_rcv), # column indices for covariates
                                   gbm.y = 1, # column index for response
                                   family = "bernoulli",
                                   tree.complexity = ifelse(sum(training_rcv$occ) < 50, 1, 5),
                                   learning.rate = lrate,
                                   bag.fraction = 0.75,
                                   max.trees = 10000,
                                   n.trees = ntrees,
                                   n.folds = 5, # 5-fold cross-validation
                                   site.weights = wt_rcv,
                                   silent = TRUE)
        t <- proc.time() - ptm
      }
      if(is.null(brt_rcv)) next
      pred_brt_rcv <- predict(brt_rcv, testing_rcv, 
                              n.trees = brt_rcv$gbm.call$best.trees,
                              type = "response")
      
      df_rcv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_rcv,
                           model = "BRT",
                           score = pred_brt_rcv,
                           label = evals_rcv[, s],
                           cv = "random",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      
      final_df <- rbind(df_scv, df_rcv)
      
      write.csv(final_df, 
                sprintf("%s/%s_%s_%s_fold%s.csv", outdir, r, spname, "BRT", f), 
                row.names = FALSE)
      
      rm(final_df, df_scv, df_rcv)
      
      cat("*")
      #*******************************************##
      ##* Lasso with CV
      # generating the quadratic terms for all continuous variables
      # function to creat quadratic terms for lasso and ridge
      quad_obj <- make_quadratic(training_scv, cols = 2:ncol(training_scv), verbose = FALSE)
      # now we can predict this quadratic object on the training and testing data
      # this make two columns for each covariates used in the transformation
      training_quad <- predict.make_quadratic(quad_obj, newdata = training_scv)
      testing_quad <- predict.make_quadratic(quad_obj, newdata = testing_scv)
      # convert the data.frames to sparse matrices
      # select all quadratic (and non-quadratic) columns, except the y (occ)
      new_vars <- names(training_quad)[names(training_quad) != "occ"]
      training_sparse <- sparse.model.matrix(~. -1, training_quad[, new_vars])
      testing_sparse <- sparse.model.matrix( ~. -1, testing_quad[, new_vars])
      
      set.seed(myseed)
      ptm <- proc.time()
      # fitting glmnet with cv
      lasso_scv <- glmnet::cv.glmnet(x = training_sparse,
                                     y = training_quad$occ,
                                     family = "binomial",
                                     alpha = 1, # fitting lasso
                                     standardize = TRUE,
                                     weights = wt_scv,
                                     nfolds = 10)
      t <- proc.time() - ptm
      pred_lasso_scv <- predict(lasso_scv, testing_sparse, s = "lambda.min", type = "response")
      
      df_scv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_scv,
                           model = "GLM-lasso",
                           score = as.numeric(pred_lasso_scv),
                           label = evals_scv[, s],
                           cv = "spatial",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      rm(testing_sparse)
      ##--------------------------------##
      # generating the quadratic terms for all continuous variables
      # function to creat quadratic terms for lasso and ridge
      quad_obj <- make_quadratic(training_rcv, cols = 2:ncol(training_rcv), verbose = FALSE)
      # now we can predict this quadratic object on the training and testing data
      # this make two columns for each covariates used in the transformation
      training_quad <- predict.make_quadratic(quad_obj, newdata = training_rcv)
      testing_quad <- predict.make_quadratic(quad_obj, newdata = testing_rcv)
      # convert the data.frames to sparse matrices
      # select all quadratic (and non-quadratic) columns, except the y (occ)
      new_vars <- names(training_quad)[names(training_quad) != "occ"]
      training_sparse <- sparse.model.matrix(~. -1, training_quad[, new_vars])
      testing_sparse <- sparse.model.matrix( ~. -1, testing_quad[, new_vars])
      
      set.seed(myseed)
      ptm <- proc.time()
      # fitting glmnet with cv
      lasso_rcv <- glmnet::cv.glmnet(x = training_sparse,
                                     y = training_quad$occ,
                                     family = "binomial",
                                     alpha = 1, # fitting lasso
                                     standardize = TRUE,
                                     weights = wt_rcv,
                                     nfolds = 10)
      
      t <- proc.time() - ptm
      pred_lasso_rcv <- predict(lasso_rcv, testing_sparse, s = "lambda.min", type = "response")
      
      df_rcv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_rcv,
                           model = "GLM-lasso",
                           score = as.numeric(pred_lasso_rcv),
                           label = evals_rcv[, s],
                           cv = "random",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      
      final_df <- rbind(df_scv, df_rcv)
      
      write.csv(final_df, 
                sprintf("%s/%s_%s_%s_fold%s.csv", outdir, r, spname, "Lasso", f), 
                row.names = FALSE)
      
      rm(final_df, df_scv, df_rcv)
      cat("*")
      #*******************************************##
      ##* MaxEnt
      occurrences <- training_scv$occ # presnece (1s) and background (0s) points
      covariates <- training_scv[, 2:ncol(training_scv)] # predictor covariates
      set.seed(myseed)
      ptm <- proc.time()
      mxnt_scv <- dismo::maxent(x = covariates,
                                p = occurrences,
                                removeDuplicates = FALSE,
                                path = "output/maxent_files2", # path to save maxent files
                                args = c("nothreshold"))
      t <- proc.time() - ptm
      pred_mxnt_scv <- predict(mxnt_scv, testing_scv, args = "outputformat=cloglog")
      
      df_scv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_scv,
                           model = "MaxEnt",
                           score = pred_mxnt_scv,
                           label = evals_scv[, s],
                           cv = "spatial",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      ##--------------------------------##
      occurrences <- training_rcv$occ # presnece (1s) and background (0s) points
      covariates <- training_rcv[, 2:ncol(training_rcv)] # predictor covariates
      set.seed(myseed)
      ptm <- proc.time()
      mxnt_rcv <- dismo::maxent(x = covariates,
                                p = occurrences,
                                removeDuplicates = FALSE,
                                path = "output/maxent_files2", # path to save maxent files
                                args = c("nothreshold"))
      t <- proc.time() - ptm
      pred_mxnt_rcv <- predict(mxnt_rcv, testing_rcv, args = "outputformat=cloglog")
      
      df_rcv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_rcv,
                           model = "MaxEnt",
                           score = pred_mxnt_rcv,
                           label = evals_rcv[, s],
                           cv = "random",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      
      final_df <- rbind(df_scv, df_rcv)
      
      write.csv(final_df, 
                sprintf("%s/%s_%s_%s_fold%s.csv", outdir, r, spname, "MaxEnt", f), 
                row.names = FALSE)
      
      rm(final_df, df_scv, df_rcv)
      cat("*")
      #*******************************************##
      ##* MaxEnt no-clamping
      occurrences <- training_scv$occ # presnece (1s) and background (0s) points
      covariates <- training_scv[, 2:ncol(training_scv)] # predictor covariates
      set.seed(myseed)
      ptm <- proc.time()
      mxnt_nocmp_scv <- dismo::maxent(x = covariates,
                                      p = occurrences,
                                      removeDuplicates = FALSE,
                                      path = "output/maxent_files2", # path to save maxent files
                                      args = c("nothreshold"))
      t <- proc.time() - ptm
      pred_mxnt_nocmp_scv <- predict(mxnt_nocmp_scv, testing_scv, args = "doclamp=false")
      
      df_scv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_scv,
                           model = "MaxEnt-noclamp",
                           score = pred_mxnt_nocmp_scv,
                           label = evals_scv[, s],
                           cv = "spatial",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      ##--------------------------------##
      occurrences <- training_rcv$occ # presnece (1s) and background (0s) points
      covariates <- training_rcv[, 2:ncol(training_rcv)] # predictor covariates
      set.seed(myseed)
      ptm <- proc.time()
      mxnt_nocmp_rcv <- dismo::maxent(x = covariates,
                                      p = occurrences,
                                      removeDuplicates = FALSE,
                                      path = "output/maxent_files2", # path to save maxent files
                                      args = c("nothreshold"))
      t <- proc.time() - ptm
      pred_mxnt_nocmp_rcv <- predict(mxnt_nocmp_rcv, testing_rcv, args = "doclamp=false")
      
      df_rcv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_rcv,
                           model = "MaxEnt-noclamp",
                           score = pred_mxnt_nocmp_rcv,
                           label = evals_rcv[, s],
                           cv = "random",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      
      final_df <- rbind(df_scv, df_rcv)
      
      write.csv(final_df, 
                sprintf("%s/%s_%s_%s_fold%s.csv", outdir, r, spname, "MaxEnt_noclamp", f), 
                row.names = FALSE)
      
      rm(final_df, df_scv, df_rcv, pred_mxnt_nocmp_rcv, pred_mxnt_nocmp_scv)
      cat("*")
      #*******************************************##
      ##* MaxEnt-tuned
      occurrences <- training_scv$occ # presnece (1s) and background (0s) points
      covariates <- training_scv[, 2:ncol(training_scv)] # predictor covariates
      
      ############### random tuned
      set.seed(myseed)
      ptm <- proc.time()
      
      param_optim <- NULL
      param_optim <- maxent_param(data = training_scv, 
                                  kf = 4,
                                  filepath = "output/maxent_files2")
      
      # save the tuned parameters
      tune_df[z <- z + 1, ] <- c(spname, f, "MaxEnt-tuned", param_optim[1], param_to_txt(param_optim))
      
      mxtund_scv <- dismo::maxent(x = covariates,
                                  p = occurrences,
                                  removeDuplicates = FALSE,
                                  path = "output/maxent_files2", # path to save maxent files
                                  args = param_optim)
      t <- proc.time() - ptm
      pred_mxtund_scv <- predict(mxtund_scv, testing_scv, args = "outputformat=cloglog")
      
      df_scv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_scv,
                           model = "MaxEnt-tuned",
                           score = pred_mxtund_scv,
                           label = evals_scv[, s],
                           cv = "spatial",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      
      ################ spatially tuned
      # extract folds
      myspfolds <- foldsID(fld = spfolds)
      myspfolds <- myspfolds[-f]
      
      set.seed(myseed)
      ptm <- proc.time()
      
      param_optim <- NULL
      param_optim <- maxent_param(data = training_scv, 
                                  folds = myspfolds,
                                  filepath = "output/maxent_files2")
      
      # save the tuned parameters
      tune_df[z <- z + 1, ] <- c(spname, f, "MaxEnt-spatial-tuned", param_optim[1], param_to_txt(param_optim))
      
      mxtund_scv <- NULL
      mxtund_scv <- dismo::maxent(x = covariates,
                                  p = occurrences,
                                  removeDuplicates = FALSE,
                                  path = "output/maxent_files2", # path to save maxent files
                                  args = param_optim)
      t <- proc.time() - ptm
      pred_mxtund_scv <- predict(mxtund_scv, testing_scv, args = "outputformat=cloglog")
      
      df_scv2 <- data.frame(species = spname,
                            region = r,
                            occ = prNum_scv,
                            model = "MaxEnt-spatial-tuned",
                            score = pred_mxtund_scv,
                            label = evals_scv[, s],
                            cv = "spatial",
                            fold = paste0("fols", f),
                            time = t[3], 
                            stringsAsFactors = FALSE)
      
      ##--------------------------------##
      occurrences <- training_rcv$occ # presenece (1s) and background (0s) points
      covariates <- training_rcv[, 2:ncol(training_rcv)] # predictor covariates
      set.seed(myseed)
      ptm <- proc.time()
      
      param_optim <- NULL
      param_optim <- maxent_param(data = training_rcv, 
                                  kf = 4, 
                                  filepath = "output/maxent_files2")

      # save the tuned parameters
      tune_df[z <- z + 1, ] <- c(spname, f, "MaxEnt-random-tuned", param_optim[1], param_to_txt(param_optim))
      
      mxtund_rcv <- dismo::maxent(x = covariates,
                                  p = occurrences,
                                  removeDuplicates = FALSE,
                                  path = "output/maxent_files2", # path to save maxent files
                                  args = param_optim)
      t <- proc.time() - ptm
      pred_mxtund_rcv <- predict(mxtund_rcv, testing_rcv, args = "outputformat=cloglog")
      
      df_rcv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_rcv,
                           model = "MaxEnt-tuned",
                           score = pred_mxtund_rcv,
                           label = evals_rcv[, s],
                           cv = "random",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      
      final_df <- rbind(df_scv, df_scv2, df_rcv)
      
      write.csv(final_df, 
                sprintf("%s/%s_%s_%s_fold%s.csv", outdir, r, spname, "MaxEnt_tuned", f), 
                row.names = FALSE)
      
      rm(final_df, df_scv, df_scv2, df_rcv, pred_mxtund_rcv, pred_mxtund_scv)
      cat("*")
      #*******************************************##
      ##* MaxEnt Simple
      occurrences <- training_scv$occ # presnece (1s) and background (0s) points
      covariates <- training_scv[, 2:ncol(training_scv)] # predictor covariates
      set.seed(myseed)
      ptm <- proc.time()
      mxt_simp_scv <- dismo::maxent(x = covariates,
                                    p = occurrences,
                                    removeDuplicates = FALSE,
                                    path = "output/maxent_files2", # path to save maxent files
                                    args = c("noautofeature", "nothreshold", 
                                             "nohinge", "noproduct", 
                                             "betamultiplier=1"))
      t <- proc.time() - ptm
      pred_mxt_simp_scv <- predict(mxt_simp_scv, testing_scv, args = "outputformat=cloglog")
      
      df_scv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_scv,
                           model = "MaxEnt-LQ",
                           score = pred_mxt_simp_scv,
                           label = evals_scv[, s],
                           cv = "spatial",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      ##--------------------------------##
      occurrences <- training_rcv$occ # presnece (1s) and background (0s) points
      covariates <- training_rcv[, 2:ncol(training_rcv)] # predictor covariates
      set.seed(myseed)
      ptm <- proc.time()
      mxt_simp_rcv <- dismo::maxent(x = covariates,
                                    p = occurrences,
                                    removeDuplicates = FALSE,
                                    path = "output/maxent_files2", # path to save maxent files
                                    args = c("noautofeature", "nothreshold", 
                                             "nohinge", "noproduct", 
                                             "betamultiplier=1"))
      t <- proc.time() - ptm
      pred_mxt_simp_rcv <- predict(mxt_simp_rcv, testing_rcv, args = "outputformat=cloglog")
      
      df_rcv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_rcv,
                           model = "MaxEnt-LQ",
                           score = pred_mxt_simp_rcv,
                           label = evals_rcv[, s],
                           cv = "random",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      
      final_df <- rbind(df_scv, df_rcv)
      
      write.csv(final_df, 
                sprintf("%s/%s_%s_%s_fold%s.csv", outdir, r, spname, "MaxEnt_LQ", f), 
                row.names = FALSE)
      
      rm(final_df, df_scv, df_rcv, pred_mxt_simp_rcv, pred_mxt_simp_scv)
      cat("*")
      #*******************************************##
      ##*******************************************##
      ##* Classification methods
      training_scv$occ <- as.factor(training_scv$occ)
      training_rcv$occ <- as.factor(training_rcv$occ)
      
      cwt_scv <- 1 / prop.table(table(training_scv$occ)) # inverse weighting
      prNum_scv <- as.numeric(table(training_scv$occ)["1"])
      samsize_scv <- c("0" = prNum_scv, "1" = prNum_scv)
      cwt_rcv <- 1 / prop.table(table(training_rcv$occ)) # inverse weighting
      prNum_rcv <- as.numeric(table(training_rcv$occ)["1"])
      samsize_rcv <- c("0" = prNum_rcv, "1" = prNum_rcv)
      ##*******************************************##
      ##* RF down-sampled
      set.seed(myseed)
      ptm <- proc.time()
      rf_dws_scv <- randomForest::randomForest(occ ~ .,
                                               data = training_scv,
                                               ntree = 1000,
                                               sampsize = samsize_scv,
                                               replace = TRUE)
      t <- proc.time() - ptm
      pred_rf_scv <- predict(rf_dws_scv, testing_scv, type = "prob")[,"1"]
      
      df_scv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_scv,
                           model = "RF down-sample",
                           score = pred_rf_scv,
                           label = evals_scv[, s],
                           cv = "spatial",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      ##--------------------------------##
      set.seed(myseed)
      ptm <- proc.time()
      rf_dws_rcv <- randomForest::randomForest(occ ~ .,
                                               data = training_rcv,
                                               ntree = 1000,
                                               sampsize = samsize_rcv,
                                               replace = TRUE)
      t <- proc.time() - ptm
      pred_rf_rcv <- predict(rf_dws_rcv, testing_rcv, type = "prob")[,"1"]
      
      df_rcv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_rcv,
                           model = "RF down-sample",
                           score = pred_rf_rcv,
                           label = evals_rcv[, s],
                           cv = "random",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      
      final_df <- rbind(df_scv, df_rcv)
      
      write.csv(final_df, 
                sprintf("%s/%s_%s_%s_fold%s.csv", outdir, r, spname, "RF", f), 
                row.names = FALSE)
      
      rm(final_df, df_scv, df_rcv)
      cat("*")
      ##*******************************************##
      ##* ranger with shallow trees
      set.seed(myseed)
      ptm <- proc.time()
      rf_ranger_scv <- ranger::ranger(formula = occ ~ .,
                                      data = training_scv, 
                                      num.trees = 2000,
                                      probability = TRUE,
                                      splitrule = "hellinger", 
                                      max.depth = 2,
                                      num.threads = 6,
                                      replace = TRUE)
      t <- proc.time() - ptm
      pred_ranger_scv <- predict(rf_ranger_scv, testing_scv, type = "response")$predictions[,"1"]
      
      df_scv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_scv,
                           model = "RF-shallow",
                           score = pred_ranger_scv,
                           label = evals_scv[, s],
                           cv = "spatial",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      ##--------------------------------##
      set.seed(myseed)
      ptm <- proc.time()
      rf_ranger_rcv <- ranger::ranger(formula = occ ~ .,
                                      data = training_rcv, 
                                      num.trees = 2000,
                                      probability = TRUE,
                                      splitrule = "hellinger", 
                                      max.depth = 2,
                                      num.threads = 6,
                                      replace = TRUE)
      t <- proc.time() - ptm
      pred_ranger_rcv <- predict(rf_ranger_rcv, testing_rcv, type = "response")$predictions[,"1"]
      
      df_rcv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_rcv,
                           model = "RF-shallow",
                           score = pred_ranger_rcv,
                           label = evals_rcv[, s],
                           cv = "random",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      
      final_df <- rbind(df_scv, df_rcv)
      
      write.csv(final_df, 
                sprintf("%s/%s_%s_%s_fold%s.csv", outdir, r, spname, "Ranger", f), 
                row.names = FALSE)
      
      rm(final_df, df_scv, df_rcv)
      cat("*")
      #*******************************************##
      ##* SVM
      set.seed(myseed)
      ptm <- proc.time()
      svm_scv <-  e1071::svm(formula = occ ~ .,
                             data = training_scv,
                             kernel = "radial",
                             scale = TRUE,
                             class.weights = cwt_scv,
                             probability = TRUE)
      t <- proc.time() - ptm
      pred_svm0 <- predict(svm_scv, testing_scv, probability = TRUE)
      pred_svm_scv <- attr(pred_svm0, "probabilities")[, "1"]
      
      df_scv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_scv,
                           model = "SVM",
                           score = pred_svm_scv,
                           label = evals_scv[, s],
                           cv = "spatial",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      ##--------------------------------##
      set.seed(myseed)
      ptm <- proc.time()
      svm_rcv <- e1071::svm(formula = occ ~ .,
                            data = training_rcv,
                            kernel = "radial",
                            scale = TRUE,
                            class.weights = cwt_rcv,
                            probability = TRUE)
      t <- proc.time() - ptm
      pred_svm0 <- predict(svm_rcv, testing_rcv, probability = TRUE)
      pred_svm_rcv <- attr(pred_svm0, "probabilities")[, "1"]
      
      df_rcv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_rcv,
                           model = "SVM",
                           score = pred_svm_rcv,
                           label = evals_rcv[, s],
                           cv = "random",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      
      final_df <- rbind(df_scv, df_rcv)
      
      write.csv(final_df, 
                sprintf("%s/%s_%s_%s_fold%s.csv", outdir, r, spname, "SVM", f), 
                row.names = FALSE)
      
      rm(final_df, df_scv, df_rcv, pred_svm_rcv, pred_svm_scv)
      cat("*")
      #*******************************************##
      ##* MARS un-weighted
      # change the response to factor variables with non-numeric levels
      levels(training_scv$occ) <- c("C0", "C1")
      levels(training_rcv$occ) <- c("C0", "C1")
      mytuneGrid <- expand.grid(nprune = 2:20,
                                degree = 1) # no interaction
      mytrControl <- trainControl(method = "cv",
                                  number = 5, # 5-fold cross-validation
                                  classProbs = TRUE,
                                  summaryFunction = twoClassSummary,
                                  allowParallel = TRUE)
      
      set.seed(myseed)
      ptm <- proc.time()
      cluster <- makeCluster(6, type = "FORK") # you can use all cores of your machine instead e.g. 8
      registerDoParallel(cluster)
      mars_scv <- train(form = occ ~ .,
                        data = training_scv,
                        method = "earth",
                        metric = "ROC",
                        trControl = mytrControl,
                        tuneGrid = mytuneGrid,
                        thresh = 0.00001)
      stopCluster(cluster)
      registerDoSEQ()
      t <- proc.time() - ptm
      pred_mars_scv <- predict(mars_scv, testing_scv, type = "prob")[,"C1"]
      
      df_scv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_scv,
                           model = "MARS",
                           score = pred_mars_scv,
                           label = evals_scv[, s],
                           cv = "spatial",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      ##--------------------------------##
      set.seed(myseed)
      ptm <- proc.time()
      cluster <- makeCluster(6, type = "FORK") # you can use all cores of your machine instead e.g. 8
      registerDoParallel(cluster)
      mars_rcv <- train(form = occ ~ .,
                        data = training_rcv,
                        method = "earth",
                        metric = "ROC",
                        trControl = mytrControl,
                        tuneGrid = mytuneGrid,
                        thresh = 0.00001)
      stopCluster(cluster)
      registerDoSEQ()
      t <- proc.time() - ptm
      pred_mars_rcv <- predict(mars_rcv, testing_rcv, type = "prob")[,"C1"]
      
      df_rcv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_rcv,
                           model = "MARS",
                           score = pred_mars_rcv,
                           label = evals_rcv[, s],
                           cv = "random",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      
      final_df <- rbind(df_scv, df_rcv)
      
      write.csv(final_df, 
                sprintf("%s/%s_%s_%s_fold%s.csv", outdir, r, spname, "MARS", f), 
                row.names = FALSE)
      
      rm(final_df, df_scv, df_rcv, pred_mars_rcv, pred_mars_scv)
      cat("*")
      #*******************************************##
      
      # normalising covariates
      for(v in myvars[[r]]){
        if(v %in% categoricalvars == FALSE){
          # spatial cv covariates
          meaanv <- mean(training_scv[,v])
          sdv <- sd(training_scv[,v])
          training_scv[,v] <- (training_scv[,v] - meaanv) / sdv
          testing_scv[,v] <- (testing_scv[,v] - meaanv) / sdv
          # random cv covariates
          meaanv2 <- mean(training_rcv[,v])
          sdv2 <- sd(training_rcv[,v])
          training_rcv[,v] <- (training_rcv[,v] - meaanv2) / sdv2
          testing_rcv[,v] <- (testing_rcv[,v] - meaanv2) / sdv2
        }
      }
      
      #*******************************************##
      ##* GAM
      set.seed(myseed)
      ptm <- proc.time()
      gam_scv <-  mgcv::gam(formula = as.formula(myform[[r]]), 
                            data = training_scv,
                            family = binomial(link = "logit"),
                            weights = wt_scv,
                            method = "REML")
      t <- proc.time() - ptm
      pred_gam_scv <- predict(gam_scv, testing_scv, type = "response")
      
      df_scv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_scv,
                           model = "GAM",
                           score = pred_gam_scv,
                           label = evals_scv[, s],
                           cv = "spatial",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      
      ##--------------------------------##
      set.seed(myseed)
      ptm <- proc.time()
      gam_rcv <- mgcv::gam(formula = as.formula(myform[[r]]), 
                           data = training_rcv,
                           family = binomial(link = "logit"),
                           weights = wt_rcv,
                           method = "REML")
      t <- proc.time() - ptm
      pred_gam_rcv <- predict(gam_rcv, testing_rcv, type = "response")
      
      df_rcv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_rcv,
                           model = "GAM",
                           score = pred_gam_rcv,
                           label = evals_rcv[, s],
                           cv = "random",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      
      final_df <- rbind(df_scv, df_rcv)
      
      write.csv(final_df, 
                sprintf("%s/%s_%s_%s_fold%s.csv", outdir, r, spname, "GAM", f), 
                row.names = FALSE)
      
      rm(final_df, df_scv, df_rcv)
      cat("*")
      #*******************************************##
      ##* GLM
      set.seed(myseed)
      ptm <- proc.time()
      # forward-backward model selection
      lm1 <- glm(occ ~., 
                 data = training_scv, 
                 weights = wt_scv, 
                 family = binomial(link = "logit"))
      lm_scv <- gam::step.Gam(lm1, 
                              scope = myscope[[r]], 
                              direction = "both", 
                              data = training_scv,
                              trace = FALSE)
      t <- proc.time() - ptm
      pred_glm_w <- predict(lm_scv, testing_scv, type = "response")
      
      df_scv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_scv,
                           model = "GLM-step",
                           score = pred_glm_w,
                           label = evals_scv[, s],
                           cv = "spatial",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      rm(pred_glm_w)
      ##--------------------------------##
      set.seed(myseed)
      ptm <- proc.time()
      # forward-backward model selection
      lm1 <- glm(occ ~., 
                 data = training_rcv, 
                 weights = wt_rcv, 
                 family = binomial(link = "logit"))
      lm_rcv <- gam::step.Gam(lm1, 
                              scope = myscope[[r]], 
                              direction = "both", 
                              data = training_rcv,
                              trace = FALSE)
      t <- proc.time() - ptm
      pred_glm_w <- predict(lm_rcv, testing_rcv, type = "response")
      
      df_rcv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_rcv,
                           model = "GLM-step",
                           score = pred_glm_w,
                           label = evals_rcv[, s],
                           cv = "random",
                           fold = paste0("fols", f),
                           time = t[3], 
                           stringsAsFactors = FALSE)
      
      final_df <- rbind(df_scv, df_rcv)
      
      write.csv(final_df, 
                sprintf("%s/%s_%s_%s_fold%s.csv", outdir, r, spname, "GLM", f), 
                row.names = FALSE)
      
      rm(final_df, df_scv, df_rcv, pred_glm_w)
      cat("*")
      ##*******************************************##
      ##* Ensemble
      pred_ens_scv <- data.frame(lasso = scales::rescale(pred_lasso_scv, to = c(0,1)),
                                 gam = scales::rescale(pred_gam_scv, to = c(0,1)),
                                 maxent = scales::rescale(pred_mxnt_scv, to = c(0,1)),
                                 brt = scales::rescale(pred_brt_scv, to = c(0,1)),
                                 rf = scales::rescale(pred_rf_scv, to = c(0,1))) %>% 
        rowMeans()
      
      df_scv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_scv,
                           model = "Ensemble",
                           score = pred_ens_scv,
                           label = evals_scv[, s],
                           cv = "spatial",
                           fold = paste0("fols", f),
                           time = NA, 
                           stringsAsFactors = FALSE)
      ##--------------------------------##
      pred_ens_rcv <- data.frame(lasso = scales::rescale(pred_lasso_rcv, to = c(0,1)),
                                 gam = scales::rescale(pred_gam_rcv, to = c(0,1)),
                                 maxent = scales::rescale(pred_mxnt_rcv, to = c(0,1)),
                                 brt = scales::rescale(pred_brt_rcv, to = c(0,1)),
                                 rf = scales::rescale(pred_rf_rcv, to = c(0,1))) %>% 
        rowMeans()
      
      df_rcv <- data.frame(species = spname,
                           region = r,
                           occ = prNum_rcv,
                           model = "Ensemble",
                           score = pred_ens_rcv,
                           label = evals_rcv[, s],
                           cv = "random",
                           fold = paste0("fols", f),
                           time = NA, 
                           stringsAsFactors = FALSE)
      
      final_df <- rbind(df_scv, df_rcv)
      
      write.csv(final_df, 
                sprintf("%s/%s_%s_%s_fold%s.csv", outdir, r, spname, "Ensemble", f), 
                row.names = FALSE)
      
      rm(final_df, df_scv, df_rcv)
      cat("*\n")
      #*******************************************##
    }
    n <- n + 1
    cat("Species", n, "is done!\n")
  }
  write.csv(tune_df, "output/maxexnt_param_total.csv", row.names = FALSE)
}

sessionInfo()
