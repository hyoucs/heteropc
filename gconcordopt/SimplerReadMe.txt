##########
Edited version that takes mask matrix as one of the input
Thus has different penalties at each entry
##########

To install:
[in terminal][**outside the gconcordopt folder**]
R CMD build gconcordopt
R CMD INSTALL gconcordopt_1.0.tar.gz

May need:
[in R shell]
install.packages('Rcpp')
install.packages('RcppEigen')

To import:
[in .r script]
library(gconcordopt)
library(MASS)

Enjoy & check bugs!
--Sikun