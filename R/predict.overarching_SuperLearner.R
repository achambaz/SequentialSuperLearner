#' Predict method for SequentialSuperLearner object
#'
#' Obtains predictions on a new data set from a SequentialSuperLearner fit.
#'
#' Each algorithm in the Super Learner library needs to have a
#' corresponding prediction function with the ``predict.'' prefixed onto the
#' algorithm name (e.g. \code{predict.SL.glm} for \code{SL.glm}).
#'
#' @param object Fitted object from \code{overarching_SuperLearner}
#' @param newdata New X values for prediction
#'
#' @return
#' \item{base_learners_predictions}{Predicted values for each base algorithm in library}
#' \item{meta_learners_predictions}{Predicted values for each meta-algorithm in library}
#' \item{overarching_predictions}{ Predicted values from the overarching super learner fit}
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
#' @seealso \code{\link{SuperLearner}}
#'
#' @keywords models
#' @method predict overarching_SuperLearner
#' @export
predict.overarching_SuperLearner <- function(object, newdata){
  
  base_learners_predictions <- SuperLearner::predict.SuperLearner(object = object$base_learners,
                                                                  newdata = newdata)
  
  base_learners_predictions <- as.data.frame(base_learners_predictions$library.predict)
  names(base_learners_predictions) <- paste0(names(base_learners_predictions),
                                             "_BLpreds")
  base_learners_predictions <- cbind(base_learners_predictions, newdata)
  
  meta_learners_predictions <- SuperLearner::predict.SuperLearner(object = object$meta_learners,
                                                                  newdata = base_learners_predictions)
  
  overarching_predictions <- object$meta_learners$method$computePred(predY = meta_learners_predictions$library.predict,
                                                              coef = as.numeric(object$coef_overarching[nrow(object$coef_overarching), -1]),
                                                              control = object$meta_learners$control)
  out <-     list(base_learners_predictions = base_learners_predictions,
                  meta_learners_predictions = meta_learners_predictions,
                  overarching_predictions = overarching_predictions
  )
  class(out) <- "predict.overarching_SuperLearner"
  return(out)
}

#' @method print predict.overarching_SuperLearner
#' @export
print.predict.overarching_SuperLearner <- function(object, n = 7) {
  cat("A prediction object based on an 'overarching_SuperLearner' fit.\n\n")
  cat("A list consisting of three objects:\n\n")
  cat("* a data.frame that notably contains the predictions for each base algorithm in the library\n\n")
  str(object$base_learners_predictions)
  cat("\n\n* a list with two entries that notably contains the predictions for each meta algorithm in the library\n\n")
  str(object$meta_learners_predictions)
  cat("\n\n* a vector containing the predictions made by the overarching super learner\n\n")
  str(object$overarching_predictions)
  return(invisible())  
}
