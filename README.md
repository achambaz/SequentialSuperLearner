# SequentialSuperLearner


The `SequentialSuperLearner` `R` package implements the so-called *overarching
sequential  super  learning  algorithm*,  a  variant  of  the  super  learning
algorithm.  Designed  to learn from  times series, it  sequentially identifies
the best algorithm in a library, or  the best combination of algorithms in the
library, where the said library consists of several super learners. 

The package is developed by G. Ecoto  (CCR and MAP5, UMR CNRS 8145, Université
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

