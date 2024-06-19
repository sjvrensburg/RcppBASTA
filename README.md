# RcppBASTA

`Rcpp` implementation of [Fryzlewicz and Subba Rao's (2014)](https://doi.org/10.1111/rssb.12054) binary segmentation for transformed auto-regressive conditional heteroscedasticity (BASTA) algorithm to detect multiple-change-point in auto-regressive conditional heteroscedastic processes.

Many thanks to Prof. Piotr Fryzlewicz who made their R implementation of BASTA freely available. You may view their original code at (http://stats.lse.ac.uk/fryzlewicz/basta/basta.html)[http://stats.lse.ac.uk/fryzlewicz/basta/basta.html].

## Roadmap

At this point, I have only implemented `stat.resid`, `bin.segm` and `inner.prod.iter` in C++. 