# SequentialSuperLearner


The `SequentialSuperLearner` `R` package implements the so-called *overarching
sequential  super  learning  algorithm*,  a  variant  of  the  super  learning
algorithm.  Designed  to learn from  times series, it  sequentially identifies
the best algorithm in a library, or  the best combination of algorithms in the
library, where the said library consists of several super learners. 

Please refer  to the preprint  [One-step ahead sequential Super  Learning from
short times series of many slightly  dependent data, and anticipating the cost
of  natural  disasters](https://arxiv.org/abs/2107.13291)  for  a  theoretical
analysis and an application to the prediction of the cost of natural disasters
in France. The abstract of the preprint reads:
> Suppose that we observe a short time series where each time-t-specific data-structure consists of many slightly dependent data indexed by a and that we want to estimate a feature of the law of the experiment that depends neither on t nor on a. We develop and study an algorithm to learn sequentially which base algorithm in a user-supplied collection best carries out the estimation task in terms of excess risk and oracular inequalities. The analysis, which uses dependency graph to model the amount of conditional independence within each t-specific data-structure and a concentration inequality by Janson [2004], leverages a large ratio of the number of distinct a's to the degree of the dependency graph in the face of a small number of t-specific data-structures. The so-called one-step ahead Super Learner is applied to the motivating example where the challenge is to anticipate the cost of natural disasters in France. 

The package is developed by G. Ecoto (CCR and MAP5, UMR CNRS 8145, Université
de Paris) and A.  Chambaz (MAP5, UMR CNRS 8145, Université de Paris). 


## Introduction

The main function of the package is `overarching_SuperLearner()`.

```r
> library("SequentialSuperLearner")
> example(overarching_SuperLearner)
```

## Citation

To cite the package, see 

```r
> citation("SequentialSuperLearner")
> toBibtex(citation("SequentialSuperLearner"))
```

## Installation 

Both  a   stable  version  and   a  development  version  are   available  via
[GitHub](https://github.com/achambaz/SequentialSuperLearner)   and    can   be
installed in R as:

```r 
devtools::install_github("achambaz/SequentialSuperLearner", ref = "main")
```

or 

```r 
devtools::install_github("achambaz/SequentialSuperLearner", ref = "develop")
```

