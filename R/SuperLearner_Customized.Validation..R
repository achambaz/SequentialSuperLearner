SuperLearner_Customized.Validation. <- function (Y, X, newX = NULL, family = stat::gaussian(),
                                                 SL.library, 
                                                 train.valid.Rows, ## only new argument
                                                 method = "method.NNLS", id = NULL, verbose = FALSE, 
                                                 control = list(), cvControl = list(), obsWeights = NULL, 
                                                 env = parent.frame()) {
  time_start = proc.time()
  if (is.character(method)) {
    if (exists(method, mode = "list")) {
      method <- get(method, mode = "list")
    }
    else if (exists(method, mode = "function")) {
      method <- get(method, mode = "function")()
    }
  }
  else if (is.function(method)) {
    method <- method()
  }
  if (!is.list(method)) {
    stop("method is not in the appropriate format. Check out help('method.template')")
  }
  if (!is.null(method$require)) {
    sapply(method$require, function(x) require(force(x), 
                                               character.only = TRUE))
  }
  control <- do.call(SuperLearner::SuperLearner.control, control)
  cvControl <- do.call(SuperLearner::SuperLearner.CV.control, cvControl)
  library <- SuperLearner:::.createLibrary(SL.library)
  SuperLearner:::.check.SL.library(library = c(unique(library$library$predAlgorithm), 
                                               library$screenAlgorithm))
  call <- match.call(expand.dots = TRUE)
  if (!inherits(X, "data.frame")) 
    message("X is not a data frame. Check the algorithms in SL.library to make sure they are compatible with non data.frame inputs")
  varNames <- colnames(X)
  N <- dim(X)[1L]
  NZ <- sum(
    sapply(train.valid.Rows, function(ll){length(ll[[2]])})
  )
  p <- dim(X)[2L]
  k <- nrow(library$library)
  kScreen <- length(library$screenAlgorithm)
  Z <- matrix(NA, NZ, k)
  libraryNames <- paste(library$library$predAlgorithm, library$screenAlgorithm[library$library$rowScreen], 
                        sep = "_")
  if (p < 2 & !identical(library$screenAlgorithm, "All")) {
    warning("Screening algorithms specified in combination with single-column X.")
  }
  fitLibEnv <- new.env()
  assign("fitLibrary", vector("list", length = k), 
         envir = fitLibEnv)
  assign("libraryNames", libraryNames, envir = fitLibEnv)
  evalq(names(fitLibrary) <- libraryNames, envir = fitLibEnv)
  errorsInCVLibrary <- rep(0, k)
  errorsInLibrary <- rep(0, k)
  if (is.null(newX)) {
    newX <- X
  }
  if (!identical(colnames(X), colnames(newX))) {
    stop("The variable names and order in newX must be identical to the variable names and order in X")
  }
  if (sum(is.na(X)) > 0 | sum(is.na(newX)) > 0 | sum(is.na(Y)) > 
      0) {
    stop("missing data is currently not supported. Check Y, X, and newX for missing values")
  }
  if (!is.numeric(Y)) {
    stop("the outcome Y must be a numeric vector")
  }
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (family$family != "binomial" & isTRUE("cvAUC" %in% 
                                           method$require)) {
    stop("'method.AUC' is designed for the 'binomial' family only")
  }
  # validRows <- CVFolds(N = N, id = id, Y = Y, cvControl = cvControl)
  validRows <- sapply(train.valid.Rows, function(ii) ii[[2]])
  if (is.null(id)) {
    id <- seq(N)
  }
  if (!identical(length(id), N)) {
    stop("id vector must have the same dimension as Y")
  }
  if (is.null(obsWeights)) {
    obsWeights <- rep(1, N)
  }
  if (!identical(length(obsWeights), N)) {
    stop("obsWeights vector must have the same dimension as Y")
  }
  if (length(train.valid.Rows) != cvControl$V) {
    stop("train.valid.Rows must be of length V")
  }
  .crossValFUN <- function(train.valid, Y, dataX, id, obsWeights, 
                           library, kScreen, k, p, libraryNames, saveCVFitLibrary) {
    tempLearn <- dataX[train.valid[[1]], , drop = FALSE]
    tempOutcome <- Y[train.valid[[1]]]
    tempValid <- dataX[train.valid[[2]], , drop = FALSE]
    tempWhichScreen <- matrix(NA, nrow = kScreen, ncol = p)
    tempId <- id[train.valid[[1]]]
    tempObsWeights <- obsWeights[train.valid[[1]]]
    for (s in seq(kScreen)) {
      screen_fn = get(library$screenAlgorithm[s], envir = env)
      testScreen <- try(do.call(screen_fn, list(Y = tempOutcome, 
                                                X = tempLearn, family = family, id = tempId, 
                                                obsWeights = tempObsWeights)))
      if (inherits(testScreen, "try-error")) {
        warning(paste("replacing failed screening algorithm,", 
                      library$screenAlgorithm[s], ", with All()", 
                      "\n "))
        tempWhichScreen[s, ] <- TRUE
      }
      else {
        tempWhichScreen[s, ] <- testScreen
      }
      if (verbose) {
        message(paste("Number of covariates in ", 
                      library$screenAlgorithm[s], " is: ", 
                      sum(tempWhichScreen[s, ]), sep = ""))
      }
    }
    out <- matrix(NA, nrow = nrow(tempValid), ncol = k)
    if (saveCVFitLibrary) {
      model_out <- vector(mode = "list", length = k)
    }
    else {
      model_out <- NULL
    }
    for (s in seq(k)) {
      pred_fn = get(library$library$predAlgorithm[s], envir = env)
      testAlg <- try(do.call(pred_fn, list(Y = tempOutcome, 
                                           X = subset(tempLearn, select = tempWhichScreen[library$library$rowScreen[s], 
                                           ], drop = FALSE), newX = subset(tempValid, 
                                                                           select = tempWhichScreen[library$library$rowScreen[s], 
                                                                           ], drop = FALSE), family = family, id = tempId, 
                                           obsWeights = tempObsWeights)))
      if (inherits(testAlg, "try-error")) {
        warning(paste("Error in algorithm", library$library$predAlgorithm[s], 
                      "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n"))
      }
      else {
        out[, s] <- testAlg$pred
        if (saveCVFitLibrary) {
          model_out[[s]] <- testAlg$fit
        }
      }
      if (verbose) 
        message(paste("CV", libraryNames[s]))
    }
    if (saveCVFitLibrary) {
      names(model_out) <- libraryNames
    }
    invisible(list(out = out, model_out = model_out))
  }

  time_train_start = proc.time()
  crossValFUN_out <- lapply(train.valid.Rows, FUN = .crossValFUN, 
                            Y = Y, dataX = X, id = id, obsWeights = obsWeights, library = library, 
                            kScreen = kScreen, k = k, p = p, libraryNames = libraryNames, 
                            saveCVFitLibrary = control$saveCVFitLibrary)
  # Z[unlist(validRows, use.names = FALSE), ] <- do.call("rbind", 
  Z <- do.call("rbind", lapply(crossValFUN_out, "[[", "out"))
  
  if (control$saveCVFitLibrary) {
    cvFitLibrary <- lapply(crossValFUN_out, "[[", "model_out")
  }
  else {
    cvFitLibrary <- NULL
  }
  errorsInCVLibrary <- apply(Z, 2, function(x) anyNA(x))
  if (sum(errorsInCVLibrary) > 0) {
    Z[, as.logical(errorsInCVLibrary)] <- 0
  }
  if (all(Z == 0)) {
    stop("All algorithms dropped from library")
  }
  getCoef <- method$computeCoef(Z = Z, Y = Y[unlist(validRows, use.names = FALSE)],
                                libraryNames = libraryNames,
                                obsWeights = obsWeights[unlist(validRows, use.names = FALSE)],
                                control = control, verbose = verbose,
                                errorsInLibrary = errorsInCVLibrary)
  coef <- getCoef$coef
  names(coef) <- libraryNames
  time_train = proc.time() - time_train_start
  if (!("optimizer" %in% names(getCoef))) {
    getCoef["optimizer"] <- NA
  }
  m <- dim(newX)[1L]
  predY <- matrix(NA, nrow = m, ncol = k)
  .screenFun <- function(fun, list) {
    screen_fn = get(fun, envir = env)
    testScreen <- try(do.call(screen_fn, list))
    if (inherits(testScreen, "try-error")) {
      warning(paste("replacing failed screening algorithm,", 
                    fun, ", with All() in full data", "\n "))
      out <- rep(TRUE, ncol(list$X))
    }
    else {
      out <- testScreen
    }
    return(out)
  }
  time_predict_start = proc.time()
  whichScreen <- sapply(library$screenAlgorithm, FUN = .screenFun, 
                        list = list(Y = Y, X = X, family = family, id = id, obsWeights = obsWeights), 
                        simplify = FALSE)
  whichScreen <- do.call(rbind, whichScreen)
  .predFun <- function(index, lib, Y, dataX, newX, whichScreen, 
                       family, id, obsWeights, verbose, control, libraryNames) {
    pred_fn = get(lib$predAlgorithm[index], envir = env)
    testAlg <- try(do.call(pred_fn, list(Y = Y, X = subset(dataX, 
                                                           select = whichScreen[lib$rowScreen[index], ], drop = FALSE), 
                                         newX = subset(newX, select = whichScreen[lib$rowScreen[index], 
                                         ], drop = FALSE), family = family, id = id, obsWeights = obsWeights)))
    if (inherits(testAlg, "try-error")) {
      warning(paste("Error in algorithm", lib$predAlgorithm[index], 
                    " on full data", "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n"))
      out <- rep.int(NA, times = nrow(newX))
    }
    else {
      out <- testAlg$pred
      if (control$saveFitLibrary) {
        eval(bquote(fitLibrary[[.(index)]] <- .(testAlg$fit)), 
             envir = fitLibEnv)
      }
    }
    if (verbose) {
      message(paste("full", libraryNames[index]))
    }
    invisible(out)
  }
  predY <- do.call("cbind", lapply(seq(k), FUN = .predFun, 
                                   lib = library$library, Y = Y, dataX = X, newX = newX, 
                                   whichScreen = whichScreen, family = family, id = id, 
                                   obsWeights = obsWeights, verbose = verbose, control = control, 
                                   libraryNames = libraryNames))
  errorsInLibrary <- apply(predY, 2, function(algorithm) anyNA(algorithm))
  if (sum(errorsInLibrary) > 0) {
    if (sum(coef[as.logical(errorsInLibrary)]) > 0) {
      warning(paste0("Re-running estimation of coefficients removing failed algorithm(s)\n", 
                     "Original coefficients are: \n", paste(coef, 
                                                            collapse = ", "), "\n"))
      Z[, as.logical(errorsInLibrary)] <- 0
      if (all(Z == 0)) {
        stop("All algorithms dropped from library")
      }
      getCoef <- method$computeCoef(Z = Z, Y = Y, libraryNames = libraryNames, 
                                    obsWeights = obsWeights, control = control, verbose = verbose, 
                                    errorsInLibrary = errorsInLibrary)
      coef <- getCoef$coef
      names(coef) <- libraryNames
    }
    else {
      warning("Coefficients already 0 for all failed algorithm(s)")
    }
  }
  getPred <- method$computePred(predY = predY, coef = coef, 
                                control = control)
  time_predict = proc.time() - time_predict_start
  colnames(predY) <- libraryNames
  if (sum(errorsInCVLibrary) > 0) {
    getCoef$cvRisk[as.logical(errorsInCVLibrary)] <- NA
  }
  time_end = proc.time()
  times = list(everything = time_end - time_start, train = time_train, 
               predict = time_predict)
  out <- list(call = call, libraryNames = libraryNames, SL.library = library, 
              SL.predict = getPred, coef = coef, library.predict = predY, 
              newX = newX,
              Z = Z, cvRisk = getCoef$cvRisk, family = family, fitLibrary = get("fitLibrary", 
                                                                                envir = fitLibEnv), cvFitLibrary = cvFitLibrary, 
              varNames = varNames, validRows = validRows, method = method, 
              whichScreen = whichScreen, control = control, cvControl = cvControl, 
              errorsInCVLibrary = errorsInCVLibrary, errorsInLibrary = errorsInLibrary, 
              metaOptimizer = getCoef$optimizer, env = env, times = times)
  class(out) <- c("SuperLearner")
  return(out)
}

environment(SuperLearner_Customized.Validation.) <- asNamespace("SuperLearner")
