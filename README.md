# unbiased enigm RNA-binding protein prediction

This package provides functions for identifying RNA-binding proteins (RBPs) from their amino acid.
In contrast to common Hidden-Markov-Model-based procedures RBDpredict does not use prior knowledge such as known RNA-binding domains.
RBDpredict uses supervised learning to identify sequence structures which are abundant in RBPs but not in non-RBPs.
Thus, it is capable of identifying RBPs which bind to RNA due to unstructured amino acid sequences -- enigm RBPs.
These enigm RBPs have been found to map to disordered regions and are thought to make one half of all RBPs 
[Castello et al](http://www.cell.com/molecular-cell/abstract/S1097-2765(16)30287-8).

## Getting started

The package can be installed from `R` with `devtools::install_github("mRcSchwering/RBDpredict", build_vignettes = TRUE)`.
This will also create vignettes which explain the package functions `browseVignettes("RBDpredict")`.
Building the vignette can take a while.
If no vignette is needed use `devtools::install_github("mRcSchwering/RBDpredict")` to install the package without vignettes.

