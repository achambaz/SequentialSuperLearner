% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overarching_SuperLearner.R
\name{overarching_SuperLearner}
\alias{overarching_SuperLearner}
\title{Trains the overarching Super Learner}
\usage{
overarching_SuperLearner(
  Y,
  X,
  newX = NULL,
  family = stat::gaussian(),
  base_learning = list(SL.library = c("SL.mean", "SL.glm"), train.valid.Rows = list(),
    cvControl = list()),
  meta_learning = list(SL.library = c("SL.mean", "SL.glm"), train.valid.Rows = list(),
    cvControl = list()),
  overarching = list(train.valid.Rows = list()),
  method = "method.NNLS",
  id = NULL,
  verbose = FALSE,
  control = list(),
  obsWeights,
  env = parent.frame()
)
}
\arguments{
\item{Y}{The outcome in the training data set. Must be a numeric vector.}

\item{X}{The predictor variables in the training data set, usually a
data.frame.}

\item{newX}{The predictor variables in the validation data set. The
structure should match X. If missing, uses X for newX.}

\item{family}{Currently  allows \code{gaussian}  or  \code{binomial}  to
describe  the  error  distribution.  Link function  information  will  be
ignored and should be contained in the method argument below.}

\item{base_learning}{A list  with three entries: \itemize{\item SL.library:
same  as the  'SL.library' argument  of the  original \code{SuperLearner}
function, i.e., either  a character vector of prediction  algorithms or a
list containing character vectors. See  details below for examples on the
structure.  A list of functions  included in the SuperLearner package can
be  found  with  \code{listWrappers()}.    By  default,  'SL.library'  is
\code{c("SL.mean",  "SL.glm")}. \item  train.valid.Rows: a  list with  as
many  entries as  folds,  each  entry a  list  itself  whose first  entry
specifies which  data can  be used  for training  and whose  second entry
specifies which data can be used for validation. \item cvControl: same as
the 'cvControl'  argument of  the original  \code{SuperLearner} function,
i.e.,  a list  of  parameters to  control  the cross-validation  process.
Parameters  include   \code{V},  \code{stratifyCV},   \code{shuffle}  and
\code{validRows}.     See    \code{\link{SuperLearner.CV.control}}    for
details.} The argument 'meta_learning' has the exact same structure.}

\item{meta_learning}{A list  with three entries: \itemize{\item SL.library:
same  as the  'SL.library' argument  of the  original \code{SuperLearner}
function, i.e., either  a character vector of prediction  algorithms or a
list containing character vectors. See  details below for examples on the
structure.  A list of functions  included in the SuperLearner package can
be  found  with  \code{listWrappers()}.    By  default,  'SL.library'  is
\code{c("SL.mean",  "SL.glm")}. \item  train.valid.Rows: a  list with  as
many  entries as  folds,  each  entry a  list  itself  whose first  entry
specifies which  data can  be used  for training  and whose  second entry
specifies which data can be used for validation. \item cvControl: same as
the 'cvControl'  argument of  the original  \code{SuperLearner} function,
i.e.,  a list  of  parameters to  control  the cross-validation  process.
Parameters  include   \code{V},  \code{stratifyCV},   \code{shuffle}  and
\code{validRows}.     See    \code{\link{SuperLearner.CV.control}}    for
details.} The argument 'base_learning' has the exact same structure.}

\item{overarching}{A  list  with  one   single  entry:  \itemize{\item
train.valid.Rows: a list with as many entries as folds, each entry a list
itself whose  first entry specifies which  data can be used  for training
and whose second entry specifies which data can be used for validation.}}

\item{method}{A list (or a function to create a list) containing details on
estimating the  coefficients for  the \strong{overarching  super learner}
and   the   overarching-meta-algorithm    to   combine   the   compteting
meta-algorithms in the library.  See \code{?method.template} for details.
Currently, the built  in options are either  "method.NNLS" (the default),
"method.NNLS2",             "method.NNloglik",            "method.CC_LS",
"method.CC_nloglik",                                                   or
"method.AUC".  NNLS and NNLS2 are non-negative least squares based on the
Lawson-Hanson  algorithm and  the  dual method  of  Goldfarb and  Idnani,
respectively.  NNLS  and NNLS2 will  work for both gaussian  and binomial
outcomes.  NNloglik  is a  non-negative binomial  likelihood maximization
using  the  BFGS  quasi-Newton   optimization  method.  NN*  methods  are
normalized  so weights  sum to  one.   CC_LS uses  Goldfarb and  Idnani's
quadratic programming algorithm to  calculate the best convex combination
of weights to minimize the squared error loss.  CC_nloglik calculates the
convex combination  of weights  that minimize  the negative  binomial log
likelihood  on   the  logistic  scale  using   the  sequential  quadratic
programming algorithm.  AUC,  which only works for  binary outcomes, uses
the  Nelder-Mead method  via the  optim  function to  minimize rank  loss
(equivalent to maximizing AUC).}

\item{id}{Optional   cluster  identification   variable.    For   the
cross-validation  splits,  \code{id}  forces  observations  in  the  same
cluster to  be in the same  validation fold.  \code{id} is  passed to the
prediction and screening  algorithms in SL.library, but be  sure to check
the individual wrappers as many of them ignore the information.}

\item{verbose}{logical; TRUE for printing progress during the computation
(helpful for debugging).}

\item{control}{A  list of  parameters to  control the  estimation process.
Parameters  include   \code{saveFitLibrary}  and   \code{trimLogit}.  See
\code{\link{SuperLearner.control}} for details.}

\item{obsWeights}{Optional observation  weights variable. As with \code{id}
above,  \code{obsWeights}  is  passed  to the  prediction  and  screening
algorithms, but many  of the built in wrappers ignore  (or can't use) the
information. If you are using  observation weights, make sure the library
you specify uses the information.}

\item{env}{Environment containing the learner functions. Defaults to the
calling environment.}
}
\value{
The function  returns  a list  with five  entries:\itemize{\item
base_learners: a complete summary of the successive trainings of the base
learners, of which the structure is the  same as that of an output of the
\code{SuperLearner} function. \item meta_learners:  a complete summary of
the successive trainings of the  meta-learners, of which the structure is
the same as that of an  output of the \code{SuperLearner} function. \item
coef_overarching:  the   weights  sequence   of  the   overarching  Super
Learner. \item  predictions_overarching_training: the  so-called Z-matrix
for the  overarching Super Learner.   \item predictions_overarching_newX:
the overarching Super Learner's predictions for the 'newX' data.}
}
\description{
This man  page relies heavily  on the  man page of  the \code{SuperLearner}
function.
}
\details{
The \code{overarching_SuperLearner} function takes  a pair \code{(X,Y)} and trains
the  overarching  Super  Learner.   The   weights  for  all  algorithms  in the
\code{SL.library} entry of the \code{meta_learning} argument are estimated,  along
with the  fits  of all  base and meta algorithms.
The \code{overarching_SuperLearner} function can also return  the predicted values based
on a validation set.

The prescreen  algorithms.  These  algorithms first  rank the  variables in
\code{X}  based   on  either  a   univariate  regression  p-value   of  the
\code{randomForest}  variable importance.   A  subset of  the variables  in
\code{X} is selected  based on a pre-defined cut-off.  With  this subset of
the X variables, the algorithms in \code{SL.library} are then fit.

The SuperLearner package contains a few prediction and screening algorithm
wrappers. The full list of wrappers can be viewed with
\code{listWrappers()}. The design of the SuperLearner package is such that
the user can easily add their own wrappers. A website is also maintained with
additional examples of wrapper functions at
\url{https://github.com/ecpolley/SuperLearnerExtra}.
}
\examples{

X1 <- rnorm(1500, 0, 1)
X2 <- rexp(1500, 0.8)
epoch <- c(rep(1990, 200),
           rep(1991, 400),
           rep(1992, 200),
           rep(1993, 200),
           rep(1994, 300),
           rep(1995, 200))
Y <- 2*X1 + X2 + rnorm(2*X1 + X2 , 2)

dat <- data.frame(X1, X2, epoch, Y)

train.valid.Rows_base <- list(
  list(which(epoch <= 1990), which(epoch == 1991)),
  list(which(epoch <= 1991), which(epoch == 1992)),
  list(which(epoch <= 1992), which(epoch == 1993)),
  list(which(epoch <= 1993), which(epoch == 1994))
)

epoch_meta <- epoch[epoch >= 1991] 
train.valid.Rows_meta <- list(
  list(which(epoch_meta <= 1991),
       which(epoch_meta == 1992)),
  list(which(epoch_meta <= 1992),
       which(epoch_meta == 1993)),
  list(which(epoch_meta <= 1993),
       which(epoch_meta == 1994))
)

epoch_overarching <- epoch[epoch >= 1992]

train.valid.Rows_overarching <- list(
   list(which(epoch_overarching <= 1992),
        which(epoch_overarching == 1993)),
   list(which(epoch_overarching <= 1993),
        which(epoch_overarching == 1994))
)

base_learning <- list(SL.library = c("SL.mean", "SL.rpart", "SL.lm"), 
                     train.valid.Rows = train.valid.Rows_base,
                     cvControl = SuperLearner::SuperLearner.CV.control(V = 4))

meta_learning <- list(SL.library = list("SL.mean", "SL.rpart", "SL.lm"), 
                      train.valid.Rows = train.valid.Rows_meta,
                      cvControl = SuperLearner::SuperLearner.CV.control(V = 3))

overarching <- list(train.valid.Rows = train.valid.Rows_overarching)

SL <- overarching_SuperLearner(Y = dat$Y[epoch <= 1994],
                               X = dat[epoch <= 1994, -4],
                               obsWeights = rep(1, nrow(dat[epoch <= 1994, ])),
                               base_learning = base_learning,
                               meta_learning = meta_learning,
                               overarching = overarching,
                               family = "gaussian")

preds <- predict(SL, dat[epoch == 1995, -4]) 


}
\references{
van  der Laan, M. J.,  Polley, E. C. and Hubbard,  A. E. (2008)
Super Learner,  \emph{Statistical Applications of Genetics  and Molecular
Biology}, \bold{6}, article 25.

Ecoto, G.,  Bibaut, A. and  Chambaz, A. (2021) One-step  ahead sequential
Super Learning from  short times series of many  slightly dependent data,
and     anticipating      the     cost     of      natural     disasters,
\url{https://arxiv.org/abs/2107.13291}.
}
\author{
Geoffrey   Ecoto  \email{geoffrey.ecoto@gmail.com}  and  Antoine
Chambaz \email{antoine.chambaz@u-paris.fr}
}
