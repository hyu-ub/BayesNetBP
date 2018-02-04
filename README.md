# BayesNetBP

This package implements belief propagation methods in Bayesian Networks to propagate evidence through the network. This package supports reasoning or inference in discrete, continuous and hybrid networks under the framework of Conditional Gaussian Bayesian networks. To cite this package, please use

>Han Yu, Moharil Janhavi, Rachael Hageman Blair. "BayesNetBP: An R package for probabilistic reasoning in Bayesian Networks". Submitted.

The belief propagation methods is implemented through interfacing the work by [Cowell, 2005](http://www.jmlr.org/papers/volume6/cowell05a/cowell05a.pdf) and the sum-product algorithm as described by Daphne Koller and Nir Friedman. Probabilistic graphical models: principles and techniques. MIT press, 2009.

This package is also available on [CRAN](https://cran.r-project.org/package=BayesNetBP), but will be more frequently updated on GitHub. To install the package from GitHub, please use

```{r, eval=FALSE}
library("devtools")
install_github("hyu-ub/BayesNetBP/BayesNetBP")
```
