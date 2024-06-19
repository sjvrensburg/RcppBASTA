# RcppBASTA

`Rcpp` implementation of [Fryzlewicz and Subba Rao's (2014)](https://doi.org/10.1111/rssb.12054) binary segmentation for transformed auto-regressive conditional heteroscedasticity (BASTA) algorithm to detect multiple-change-point in auto-regressive conditional heteroscedastic processes.

Many thanks to [Prof. Piotr Fryzlewicz](https://stats.lse.ac.uk/fryzlewicz/) who made their R implementation of BASTA freely available. You may view the original code at (http://stats.lse.ac.uk/fryzlewicz/basta/basta.html)[http://stats.lse.ac.uk/fryzlewicz/basta/basta.html].

## Installation

If you have the `remotes` package installed, then you can use the following to install `RcppBASTA`:

```r
# Uncomment the next line if you don't have the remotes package installed.
# install.packages("remotes")
remotes::install_github("sjvrensburg/RcppBASTA")
```

## Roadmap

At this point, I only implemented `stat.resid`, `bin.segm` and `inner.prod.iter` in C++. The rest of the code is mostly unchanged from that of Prof. Piotr Fryzlewicz. I intend to implement more of the code in C++ at some stage.
