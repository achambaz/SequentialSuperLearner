#' SequentialSuperLearner: Implements the overarching sequential super learning algorithm
#'
#' The overarching sequential super learning algorithm is a variant of the super 
#' learning algorithm. Designed to learn from times series, it sequentially identifies 
#' the best algorithm in a library, or the best combination of algorithms in the library,
#' where the said library consists of several super learners.
#' 
#' 
#' @docType package
#' @name SequentialSuperLearner
NULL
#> NULL


#' Trains the overarching Super Learner
#'
#' This man  page relies heavily  on the  man page of  the \code{SuperLearner}
#' function.
#' 
#' 
#'
#' The \code{overarching_SuperLearner} function takes  a pair \code{(X,Y)} and trains
#' the  overarching  Super  Learner.   The   weights  for  all  algorithms  in the
#' \code{SL.library} entry of the \code{meta_learning} argument are estimated,  along 
#' with the  fits  of all  base and meta algorithms.
#' The \code{overarching_SuperLearner} function can also return  the predicted values based
#' on a validation set.
#' 
#' The prescreen  algorithms.  These  algorithms first  rank the  variables in
#' \code{X}  based   on  either  a   univariate  regression  p-value   of  the
#' \code{randomForest}  variable importance.   A  subset of  the variables  in
#' \code{X} is selected  based on a pre-defined cut-off.  With  this subset of
#' the X variables, the algorithms in \code{SL.library} are then fit.
#' 
#' The SuperLearner package contains a few prediction and screening algorithm
#' wrappers. The full list of wrappers can be viewed with
#' \code{listWrappers()}. The design of the SuperLearner package is such that
#' the user can easily add their own wrappers. A website is also maintained with
#' additional examples of wrapper functions at
#' \url{https://github.com/ecpolley/SuperLearnerExtra}.
#'
#' 
#' 
#' @param Y The outcome in the training data set. Must be a numeric vector.
#'
#' @param X The predictor variables in the training data set, usually a
#' data.frame.
#'
#' @param newX The predictor variables in the validation data set. The
#' structure should match X. If missing, uses X for newX.
#'
#' @param  family  Currently  allows \code{gaussian}  or  \code{binomial}  to
#'   describe  the  error  distribution.  Link function  information  will  be
#'   ignored and should be contained in the method argument below.
#'
#' @param base_learning A list  with three entries: \itemize{\item SL.library:
#'   same  as the  'SL.library' argument  of the  original \code{SuperLearner}
#'   function, i.e., either  a character vector of prediction  algorithms or a
#'   list containing character vectors. See  details below for examples on the
#'   structure.  A list of functions  included in the SuperLearner package can
#'   be  found  with  \code{listWrappers()}.    By  default,  'SL.library'  is
#'   \code{c("SL.mean",  "SL.glm")}. \item  train.valid.Rows: a  list with  as
#'   many  entries as  folds,  each  entry a  list  itself  whose first  entry
#'   specifies which  data can  be used  for training  and whose  second entry
#'   specifies which data can be used for validation. \item cvControl: same as
#'   the 'cvControl'  argument of  the original  \code{SuperLearner} function,
#'   i.e.,  a list  of  parameters to  control  the cross-validation  process.
#'   Parameters  include   \code{V},  \code{stratifyCV},   \code{shuffle}  and
#'   \code{validRows}.     See    \code{\link{SuperLearner.CV.control}}    for
#'   details.} The argument 'meta_learning' has the exact same structure.
#'
#' @param meta_learning A list  with three entries: \itemize{\item SL.library:
#'   same  as the  'SL.library' argument  of the  original \code{SuperLearner}
#'   function, i.e., either  a character vector of prediction  algorithms or a
#'   list containing character vectors. See  details below for examples on the
#'   structure.  A list of functions  included in the SuperLearner package can
#'   be  found  with  \code{listWrappers()}.    By  default,  'SL.library'  is
#'   \code{c("SL.mean",  "SL.glm")}. \item  train.valid.Rows: a  list with  as
#'   many  entries as  folds,  each  entry a  list  itself  whose first  entry
#'   specifies which  data can  be used  for training  and whose  second entry
#'   specifies which data can be used for validation. \item cvControl: same as
#'   the 'cvControl'  argument of  the original  \code{SuperLearner} function,
#'   i.e.,  a list  of  parameters to  control  the cross-validation  process.
#'   Parameters  include   \code{V},  \code{stratifyCV},   \code{shuffle}  and
#'   \code{validRows}.     See    \code{\link{SuperLearner.CV.control}}    for
#'   details.} The argument 'base_learning' has the exact same structure.
#'
#' @param  overarching   A  list  with  one   single  entry:  \itemize{\item
#'   train.valid.Rows: a list with as many entries as folds, each entry a list
#'   itself whose  first entry specifies which  data can be used  for training
#'   and whose second entry specifies which data can be used for validation.}
#' 
#' @param method A list (or a function to create a list) containing details on
#'   estimating the  coefficients for  the \strong{overarching  super learner}
#'   and   the   overarching-meta-algorithm    to   combine   the   compteting
#'   meta-algorithms in the library.  See \code{?method.template} for details.
#'   Currently, the built  in options are either  "method.NNLS" (the default),
#'   "method.NNLS2",             "method.NNloglik",            "method.CC_LS",
#'   "method.CC_nloglik",                                                   or
#'   "method.AUC".  NNLS and NNLS2 are non-negative least squares based on the
#'   Lawson-Hanson  algorithm and  the  dual method  of  Goldfarb and  Idnani,
#'   respectively.  NNLS  and NNLS2 will  work for both gaussian  and binomial
#'   outcomes.  NNloglik  is a  non-negative binomial  likelihood maximization
#'   using  the  BFGS  quasi-Newton   optimization  method.  NN*  methods  are
#'   normalized  so weights  sum to  one.   CC_LS uses  Goldfarb and  Idnani's
#'   quadratic programming algorithm to  calculate the best convex combination
#'   of weights to minimize the squared error loss.  CC_nloglik calculates the
#'   convex combination  of weights  that minimize  the negative  binomial log
#'   likelihood  on   the  logistic  scale  using   the  sequential  quadratic
#'   programming algorithm.  AUC,  which only works for  binary outcomes, uses
#'   the  Nelder-Mead method  via the  optim  function to  minimize rank  loss
#'   (equivalent to maximizing AUC).
#'
#' @param   id   Optional   cluster  identification   variable.    For   the
#'   cross-validation  splits,  \code{id}  forces  observations  in  the  same
#'   cluster to  be in the same  validation fold.  \code{id} is  passed to the
#'   prediction and screening  algorithms in SL.library, but be  sure to check
#'   the individual wrappers as many of them ignore the information.
#'
#' @param verbose logical; TRUE for printing progress during the computation
#'   (helpful for debugging).
#'
#' @param  control A  list of  parameters to  control the  estimation process.
#'   Parameters  include   \code{saveFitLibrary}  and   \code{trimLogit}.  See
#'   \code{\link{SuperLearner.control}} for details.
#'
#' @param obsWeights Optional observation  weights variable. As with \code{id}
#'   above,  \code{obsWeights}  is  passed  to the  prediction  and  screening
#'   algorithms, but many  of the built in wrappers ignore  (or can't use) the
#'   information. If you are using  observation weights, make sure the library
#'   you specify uses the information. 
#' 
#' @param env Environment containing the learner functions. Defaults to the
#' calling environment.
#' 
#' @return  The function  returns  a list  with five  entries:\itemize{\item
#'   base_learners: a complete summary of the successive trainings of the base
#'   learners, of which the structure is the  same as that of an output of the
#'   \code{SuperLearner} function. \item meta_learners:  a complete summary of
#'   the successive trainings of the  meta-learners, of which the structure is
#'   the same as that of an  output of the \code{SuperLearner} function. \item
#'   coef_overarching:  the   weights  sequence   of  the   overarching  Super
#'   Learner. \item  predictions_overarching_training: the  so-called Z-matrix
#'   for the  overarching Super Learner.   \item predictions_overarching_newX:
#'   the overarching Super Learner's predictions for the 'newX' data.}
#'
#' @author  Geoffrey   Ecoto  \email{geoffrey.ecoto@@gmail.com}  and  Antoine
#'   Chambaz \email{antoine.chambaz@@u-paris.fr}
#'
#' @references van  der Laan, M. J.,  Polley, E. C. and Hubbard,  A. E. (2008)
#'   Super Learner,  \emph{Statistical Applications of Genetics  and Molecular
#'   Biology}, \bold{6}, article 25.
#'
#'   Ecoto, G.,  Bibaut, A. and  Chambaz, A. (2021) One-step  ahead sequential
#'   Super Learning from  short times series of many  slightly dependent data,
#'   and     anticipating      the     cost     of      natural     disasters,
#'   \url{https://arxiv.org/abs/2107.13291}.
#' 
#' @examples
#' 
#' X1 <- rnorm(1500, 0, 1)
#' X2 <- rexp(1500, 0.8)
#' epoch <- c(rep(1990, 200),
#'            rep(1991, 400),
#'            rep(1992, 200),
#'            rep(1993, 200),
#'            rep(1994, 300),
#'            rep(1995, 200))
#' Y <- 2*X1 + X2 + rnorm(2*X1 + X2 , 2)
#'
#' dat <- data.frame(X1, X2, epoch, Y)
#'
#' train.valid.Rows_base <- list(
#'   list(which(epoch <= 1990), which(epoch == 1991)),
#'   list(which(epoch <= 1991), which(epoch == 1992)),
#'   list(which(epoch <= 1992), which(epoch == 1993)),
#'   list(which(epoch <= 1993), which(epoch == 1994))
#' )
#'
#' epoch_meta <- epoch[epoch >= 1991] 
#' train.valid.Rows_meta <- list(
#'   list(which(epoch_meta <= 1991),
#'        which(epoch_meta == 1992)),
#'   list(which(epoch_meta <= 1992),
#'        which(epoch_meta == 1993)),
#'   list(which(epoch_meta <= 1993),
#'        which(epoch_meta == 1994))
#' )
#'
#' epoch_overarching <- epoch[epoch >= 1992]
#'
#' train.valid.Rows_overarching <- list(
#'    list(which(epoch_overarching <= 1992),
#'         which(epoch_overarching == 1993)),
#'    list(which(epoch_overarching <= 1993),
#'         which(epoch_overarching == 1994))
#' )
#'
#' base_learning <- list(SL.library = c("SL.mean", "SL.rpart", "SL.lm"), 
#'                      train.valid.Rows = train.valid.Rows_base,
#'                      cvControl = SuperLearner::SuperLearner.CV.control(V = 4))
#'
#' meta_learning <- list(SL.library = list("SL.mean", "SL.rpart", "SL.lm"), 
#'                       train.valid.Rows = train.valid.Rows_meta,
#'                       cvControl = SuperLearner::SuperLearner.CV.control(V = 3))
#'
#' overarching <- list(train.valid.Rows = train.valid.Rows_overarching)
#' 
#' SL <- overarching_SuperLearner(Y = dat$Y[epoch <= 1994],
#'                                X = dat[epoch <= 1994, -4],
#'                                obsWeights = rep(1, nrow(dat[epoch <= 1994, ])),
#'                                base_learning = base_learning,
#'                                meta_learning = meta_learning,
#'                                overarching = overarching,
#'                                family = "gaussian")
#'
#' preds <- predict(SL, dat[epoch == 1995, -4]) 
#' 
#' 
#' @export
overarching_SuperLearner <- function(Y, X, newX = NULL, family = stat::gaussian(),
                                     base_learning = list(SL.library = c("SL.mean", "SL.glm"), 
                                                          train.valid.Rows = list(), 
                                                          cvControl = list()),
                                     meta_learning = list(SL.library = c("SL.mean", "SL.glm"), 
                                                          train.valid.Rows = list(),
                                                          cvControl = list()), 
                                     overarching = list(train.valid.Rows = list()), 
                                     method = "method.NNLS", id = NULL, verbose = FALSE, 
                                     control = list(), obsWeights, # obsWeights = rep(1, length(Y)) 
                                     env = parent.frame()) {  
  ## base_learners
  
  base_learners_args <- as.list(sys.call())[-1]
  base_learners_args$base_learning <- base_learners_args$meta_learning <- NULL
  base_learners_args$overarching <- NULL
  base_learners_args <- c(base_learners_args, base_learning)
  
  base_learners <- do.call(SequentialSuperLearner:::SuperLearner_Customized.Validation., base_learners_args)
  
  ## meta_learners
  
  meta_learners_args <- as.list(sys.call())[-1]
  
  meta_learners_args$base_learning <- meta_learners_args$meta_learning <- NULL
  meta_learners_args$overarching <- NULL
  
  meta_learners_args$newX <- as.data.frame(base_learners$library.predict)
  names(meta_learners_args$newX) <- paste0(base_learners$libraryNames, "_BLpreds")
  meta_learners_args$newX <- cbind(meta_learners_args$newX, base_learners$newX)

  meta_learners_args$obsWeights <- obsWeights[unlist(base_learners$validRows, use.names = FALSE)]
  
  meta_learners_args <- c(meta_learners_args, meta_learning)
  meta_learners_args$Y <- Y[unlist(base_learners$validRows, use.names = FALSE)]
  meta_learners_args$X <- as.data.frame(base_learners$Z)
  names(meta_learners_args$X) <- paste0(base_learners$libraryNames, "_BLpreds")
  meta_learners_args$X <- cbind(meta_learners_args$X,
                                X[unlist(base_learners$validRows, use.names = FALSE), ])
  meta_learners <- do.call(SequentialSuperLearner:::SuperLearner_Customized.Validation., meta_learners_args)
  
  ## overarching super learner
  
  N_epoch <- length(overarching$train.valid.Rows)
  coef_overarching <- matrix(NA, N_epoch, 1 + dim(meta_learners$Z)[2])
  predictions_overarching <- c()
  
  for(epoch in 1:N_epoch){
    
    Y_meta <- Y[unlist(base_learners$validRows, use.names = FALSE)]
    Y_meta <- Y_meta[unlist(meta_learners$validRows, use.names = FALSE)]

    obsWeights_meta <- obsWeights[unlist(base_learners$validRows, use.names = FALSE)]
    obsWeights_meta <- obsWeights[unlist(meta_learners$validRows, use.names = FALSE)]
        
    idx_epoch_train <- overarching$train.valid.Rows[[epoch]][[1]]
    idx_epoch_validation <- overarching$train.valid.Rows[[epoch]][[2]]
    
    ## Compute weights for each algorithm in library.
    getCoef <- meta_learners$method$computeCoef(Z = meta_learners$Z[idx_epoch_train, ],
                                  Y = Y_meta[idx_epoch_train],
                                  libraryNames = meta_learners$libraryNames,
                                  obsWeights = obsWeights_meta[idx_epoch_train],
                                  control = control, verbose = verbose,
                                  errorsInLibrary = meta_learners$errorsInCVLibrary)
    coef <- getCoef$coef
    names(coef) <- meta_learners$libraryNames
    
    ## Compute super learner predictions at each epoch.
    getPred <- meta_learners$method$computePred(predY = meta_learners$Z[idx_epoch_validation, ],
                                                coef = coef, control = control)
    
    coef_overarching[epoch, ] <- c(epoch, coef)
    predictions_overarching <- c(predictions_overarching, getPred)
    
  }
  
  colnames(coef_overarching) <- c("Epoch", names(coef))

  ## Compute super learner predictions on newX.
  getPred_NewX <- meta_learners$method$computePred(predY = meta_learners$library.predict,
                                                   coef = coef, control = control)
  out <- list(base_learners = base_learners,
              meta_learners = meta_learners,
              coef_overarching = coef_overarching,
              predictions_overarching_training = predictions_overarching,
              predictions_overarching_newX = getPred_NewX)
  class(out) <- "overarching_SuperLearner"
  return(out)
  
}

environment(overarching_SuperLearner) <- asNamespace("SuperLearner")

#' @method print overarching_SuperLearner
#' @export
print.overarching_SuperLearner <- function(x, ...) {
  cat("A fit object from the 'overarching_SuperLearner' function.\n\n")
  cat("A list five entries (see the function's manual):\n\n")
  cat("* base_learners\n\n")
  cat("\n\n* meta_learners\n\n")
  cat("\n\n* coef_overarching\n\n")
  utils::str(x$coef_overarching)
  cat("\n\n* predictions_overarching_training\n\n")
  utils::str(x$predictions_overarching_training)
  cat("\n\n* predictions_overarching_newX\n\n")
  utils::str(x$predictions_overarching_newX)
  
  return(invisible())  
}
