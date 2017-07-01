# BayesNetBP

This package implements belief propagation methods in Bayesian Networks to propagate evidence through the network. This package supports reasoning or inference in discrete, continuous and hybrid networks under the framework of Conditional Gaussian Bayesian networks.

This package is also available on [CRAN](https://cran.r-project.org/package=BayesNetBP), but will be more frequently updated on GitHub. To install the package from GitHub, please use

```{r, eval=FALSE}
library("devtools")
install_github("hyu-ub/BayesNetBP")
```

The [vignette](https://github.com/hyu-ub/BayesNetBP/blob/master/inst/doc/BayesNetBP_intro.pdf) has instructions on how to use this package to perform probabilistic reasoning and inferece in Bayesian networks. To view the [vignette](https://github.com/hyu-ub/BayesNetBP/blob/master/inst/doc/BayesNetBP_intro.pdf), please run

```{r, eval=FALSE}
library("BayesNetBP")
browseVignettes("BayesNetBP")
```
