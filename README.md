
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
devtools::install_github("audreyfulab/MRGN", dependencies = TRUE)
devtools::install_github("audreyfulab/MRGNgeneral")
```

## Example

``` r
library(MRGNgeneral)
## Simulate some data with 20 phenotypes, 10 intermediate and 10 common children confounders, 40 confounders
set.seed(1123)
net20data <- sample.graph.data (n_t = 20,
                                n_v.t = 1,
                                family.n_v = NULL,
                                conf.num.vec = c(W = 10, Z = 10,
                                                 U = 40, K = 0, I = 20),
                                graph_type = "scale-free",
                                degree = 3,
                                theta = .4,
                                b0 = 0,
                                b.snp = c(-0.5, 0.5),
                                b.med = c(-0.8, 0.8),
                                sigma = 0.1,
                                neg.freq = 0.5,
                                conf.coef.ranges = list(W = c(0.4, 0.5),
                                                        Z = c(1, 1.5),
                                                        U = c(0.4, 0.5),
                                                        K = c(0.01, 0.1)),
                                scale = FALSE,
                                sample.size = 100)

## Performing confounding variable selection
confnet20 <- get.conf.sets(data = net20data$data,
                           n_v = net20data$dims$n_v,
                           n_t = net20data$dims$n_t,
                           n_c = NCOL(net20data$data) - net20data$dims$n_v - net20data$dims$n_t,
                           blocksize = 10,
                           T.measure = 'partial',
                           C.measure = 'partial',
                           FDRcontrol = 'qvalue',
                           adjust_by = 'individual',
                           alpha = 0.01,
                           fdr = 0.05,
                           lambda = 0.05,
                           pi0.method = 'smoother')

## Recall and precision of the selection procedure
Perf <- assess.conf.selection (confnet20,
                               adjacency = net20data$adjacency,
                               n_v = net20data$dims$n_v,
                               n_t = net20data$dims$n_t,
                               n_w = net20data$dims$n_w,
                               n_z = net20data$dims$n_z,
                               n_u = net20data$dims$n_u)
Perf$recall
Perf$precision

## MRGN example
# a small network
# Build a graph skeleton
Adjacency0 <- get.initial.skeleton (data = net20data$data,
                                    n_v = net20data$dims$n_v,
                                    n_t = net20data$dims$n_t,
                                    threshold_v = 0.2,
                                    threshold_m = 0.05,
                                    conf.sets = confnet20)
# run MRGN:
MRGNfit <- MRGN(data = net20data$data,
                n_v = net20data$dims$n_v,
                n_t = net20data$dims$n_t,
                Qlabels = confnet20$WZindices,
                n_q = length(confnet20$WZindices),
                n_u = net20data$dims$n_w + net20data$dims$n_z +
                  net20data$dims$n_u + net20data$dims$n_k +
                  net20data$dims$n_i - length(confnet20$WZindices),
                adjacency = Adjacency0,
                confounders = confnet20$confounders,
                alpha = 0.01,
                FDRcontrol = 'bonferroni',
                fdr = 0.05,
                verbose = TRUE)


#a large network
data ('networkA11')  # .rda structure with simulated data
data ('confsetsA11') #  obtained from get.conf.sets() on 'networkA11'

# Build a graph skeleton
Adjacency0 <- get.initial.skeleton (data = networkA11$data,
                                    n_v = networkA11$dims$n_v,
                                    n_t = networkA11$dims$n_t,
                                    threshold_v = 0.2,
                                    threshold_m = 0.05,
                                    conf.sets = confsetsA11)


# Run MRGN: 
MRGNfit <- MRGN(data = networkA11$data,
                n_v = networkA11$dims$n_v,
                n_t = networkA11$dims$n_t,
                Qlabels = confsetsA11$WZindices,
                n_q = length(confsetsA11$WZindices),
                n_u = networkA11$dims$n_w + networkA11$dims$n_z +
                  networkA11$dims$n_u + networkA11$dims$n_k +
                  networkA11$dims$n_i - length(confsetsA11$WZindices),
                adjacency = Adjacency0,
                confounders = confsetsA11$confounders,
                alpha = 0.01,
                FDRcontrol = 'bonferroni',
                fdr = 0.05,
                verbose = TRUE)

# run.methods() is a wrap-up function that can run MRGN and other bn.learn methods
# generate some data
seed = 123
simdata <- sample.graph.data(
  n_t = 100, n_v.t = 1, family.n_v = NULL, conf.num.vec = c(W = 50, Z = 50, U = 200, K = 0, I = 10), graph_type = "scale-free", 
  degree = 3, theta = 0.4, b0 = 0, b.snp = c(-0.5, 0.5), b.med = c(-0.8, 0.8), sigma = 0.1, neg.freq = 0.5, 
  conf.coef.ranges = list(W = c(0.4, 0.5), Z = c(1, 1.5), U = c(0.4, 0.5), K = c(0.15, 0.5)), 
  scale = FALSE, sample.size = 500, seed = seed
)
# model with MRGNgeneral
#install.packages("bnlearn") if you set bn.methods, otherwise run.methods() only run MRGN method
library(bnlearn)
result <- run.methods(
  simdata = simdata,    
  bn.methods = c("tabu", "hc"), #bn.methods = c("none",  "tabu", "hc", "pc.stable", "mmhc")
  selection_FDRcontrol = "qvalue", selection_fdr = 0.05, nb.cl = 8,
  verbose = 1
)
result$fits$MRGN$adjacency

```

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
