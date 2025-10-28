
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MRGNgeneral

MRGNgeneral builds on and extends the capability of MRGN (for a trio of
a genotype node and two phenotype nodes) and infers a causal network
from genomic data. The network includes genotype nodes (as instrumental
variables for Mendelian randomization), phenotype nodes, and different
types of confounding variables. Directed edges are likely to have a
causal interpretation under the assumption of Mendelian randomization.
There may be undirected edges when the data are not sufficiently
informative.

This package depends on the MRGN package:
<https://github.com/audreyfulab/MRGN>, as well as a few other packages.
You need to install MRGN before installing MRGNgeneral.

## Installation

You can install the development version of MRGNgeneral from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("audreyfulab/MRGNgeneral")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MRGNgeneral)
?MRGN
# generate some data
seed = 123
simdata <- sample.graph.data(
  n_t = 100, n_v.t = 1, family.n_v = NULL, conf.num.vec = c(W = 50, Z = 50, U = 200, K = 0, I = 10), graph_type = "scale-free", 
  degree = 3, theta = 0.4, b0 = 0, b.snp = c(-0.5, 0.5), b.med = c(-0.8, 0.8), sigma = 0.1, neg.freq = 0.5, 
  conf.coef.ranges = list(W = c(0.4, 0.5), Z = c(1, 1.5), U = c(0.4, 0.5), K = c(0.15, 0.5)), 
  scale = FALSE, sample.size = 500, seed = seed
)
# model with MRGNgeneral
#install.packages("bnlearn", repos = "https://cran.r-project.org") if you set bn.methods, otherwise run.methods() only run MRGNgeneral method
library(bnlearn)
result <- run.methods(
  simdata = simdata,    
  bn.methods = c("tabu", "hc"), #bn.methods = c("none",  "tabu", "hc", "pc.stable", "mmhc")
  selection_FDRcontrol = "qvalue", selection_fdr = 0.05, nb.cl = 8,
  verbose = 1
)
result_8$fits$MRGN$adjacency
```

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
