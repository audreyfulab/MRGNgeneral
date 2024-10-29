
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MRGNgeneral

<img src="C:/Users/Bruin/Desktop/Research Assistantship/MRGN_logo/MRGN-logos.png" width="100%" />
<!-- badges: start --> <!-- badges: end -->

MRGNgeneral builds on and extends the capability of MRGN (for a trio of a genotype node and two phenotype nodes) and infers a causal network from genomic data.  The network includes genotype nodes (as instrumental variables for Mendelian randomization), phenotype nodes, and different types of confounding variables.  Directed edges are likely to have a causal interpretation under the assumption of Mendelian randomization.  There may be undirected edges when the data are not sufficiently informative.

This package depends on the MRGN package: https://github.com/Jarred6068/MRGN, as well as a few other packages.  You need to install MRGN before installing MRGNgeneral.

## Installation

You can install the development version of MRGNgeneral from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("audreyqyfu/MRGNgeneral")
```

## Example

This is an example in the documentation for the main function MRGN():

``` r
library(MRGNgeneral)
?MRGN
```

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
