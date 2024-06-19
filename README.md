# RcppBASTA

`Rcpp` implementation of [Fryzlewicz and Subba Rao's (2014)](https://doi.org/10.1111/rssb.12054) binary segmentation for transformed auto-regressive conditional heteroscedasticity (BASTA) algorithm to detect multiple-change-point in auto-regressive conditional heteroscedastic processes.

## Installation

If you have the `remotes` package installed, then you can use the following to install `RcppBASTA`:

```r
# Uncomment the next line if you don't have the remotes package installed.
# install.packages("remotes")
remotes::install_github("sjvrensburg/RcppBASTA")
```

## Ackowledgements

Many thanks to [Prof. Piotr Fryzlewicz](https://stats.lse.ac.uk/fryzlewicz/) who made their R implementation of BASTA freely available. You may view the original code at (http://stats.lse.ac.uk/fryzlewicz/basta/basta.html)[http://stats.lse.ac.uk/fryzlewicz/basta/basta.html].

**PLEASE CITE** the paper [_"Multiple-change-point detection for auto-regressive conditional heteroscedastic processes"_](https://doi.org/10.1111/rssb.12054) by Fryzlewicz and Subba Rao (2014) if you use this package.

```bibtex
@article{10.1111/rssb.12054,
    author = {Fryzlewicz, P. and Subba Rao, S.},
    title = "{Multiple-Change-Point Detection for Auto-Regressive Conditional Heteroscedastic Processes}",
    journal = {Journal of the Royal Statistical Society Series B: Statistical Methodology},
    volume = {76},
    number = {5},
    pages = {903-924},
    year = {2013},
    month = {11},
    issn = {1369-7412},
    doi = {10.1111/rssb.12054},
    url = {https://doi.org/10.1111/rssb.12054},
    eprint = {https://academic.oup.com/jrsssb/article-pdf/76/5/903/49506908/jrsssb\_76\_5\_903.pdf},
}
```

## Roadmap

At this point, I only implemented `stat.resid`, `bin.segm` and `inner.prod.iter` in C++. The rest of the code is mostly unchanged from that of Prof. Piotr Fryzlewicz. I intend to implement more of the code in C++ at some stage.
