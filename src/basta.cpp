#define ARMA_DONT_USE_OPENMP
#include "RcppArmadillo.h"

#include <algorithm>
#include <vector>
#include <boost/math/special_functions/sign.hpp>

using Rcpp::_;
using Rcpp::as;
using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::LogicalVector;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::Range;
using Rcpp::Rcout;
using Rcpp::seq;

using arma::colvec;
using arma::uword;

// [[Rcpp::export(name = "stat.resid")]]
arma::colvec stat_resid(const arma::colvec& x,
    Rcpp::Nullable<NumericVector> a__ = R_NilValue,
    unsigned int order = 1,
    const double factor = 8,
    const double epsilon = 0.0)
{
    NumericVector a_;
    auto const n = x.n_elem;

    // Determine order if a_ is not NULL
    if (a__.isNotNull())
    {
        a_ = a__;
        order = a_.length() - 1;
    }

    if (a__.isNull())
    {
        Rcpp::Function f("stat.est");
        a_ = f(x, _["order"] = order);

        for (size_t i = 1; i < a_.length(); ++i)
        {
            a_[i] = a_[i] / factor;
        }
    }

    colvec a = as<colvec>(a_);
    colvec y(n-order);
    const uword n_y = y.n_elem;
    double sum_a, numerator, denominator;

    for (uword i = 0; i < n_y; ++i)
    {
        numerator = (x(n - i - 1)) * (x(n - i - 1));
        colvec tmp = arma::join_vert(
            arma::ones<colvec>(1), arma::square(
                arma::reverse(x.subvec(n - 1 - i - order, n - 2 - i))));
        denominator = arma::sum(a % tmp);
        denominator += epsilon * numerator;
        y(n_y - 1 - i) = numerator / denominator;
    }

    return y;
}

// [[Rcpp::export(name = "bin.segm")]]
Rcpp::List bin_segm(Rcpp::List buh, double th)
{
    using boost::math::sign;
    using std::abs;
    using std::ceil;

    List trees = buh["tree"];
    size_t J = trees.size();
    NumericMatrix tree_j;
    NumericMatrix tree_j_min_1;
    size_t K;

    // First loop: modify based on threshold
    for (size_t j = 0; j < J; j++)
    {
        tree_j = as<NumericMatrix>(trees[j]);
        K = tree_j.ncol();

        for (size_t k = 0; k < K; k++)
        {
            tree_j(1, k) *= abs(tree_j(1, k)) > th ? 1 : 0;
        }

        trees[j] = tree_j;
    }

    // Second loop: adjust values based on parent nodes
    for (size_t j = 1; j < J; j++)
    {
        tree_j = as<NumericMatrix>(trees[j]);
        tree_j_min_1 = as<NumericMatrix>(trees[j - 1]);
        K = tree_j.ncol();
        int par_ind = 0;

        for (size_t k = 0; k < K; k++)
        {
            while (ceil(0.5 * tree_j(0, k)) > tree_j_min_1(0, par_ind))
            {
                par_ind++;
            }

            tree_j(1, k) = tree_j(1, k) * abs(sign(tree_j_min_1(1, par_ind)));
        }

        trees[j] = tree_j;
    }

    buh["tree"] = trees;

    return buh;
}

// [[Rcpp::export(name = "inner.prod.iter")]]
arma::colvec inner_prod_iter(const arma::colvec &x)
{
    using std::sqrt;

    uword n = x.n_elem;
    colvec I_plus = arma::zeros(n - 1);
    colvec I_minus = arma::zeros(n - 1);

    I_plus(0) = sqrt(1.0 - 1.0 / n) * x(0);
    I_minus(0) = arma::sum(x.tail(n - 1)) / sqrt(std::pow(n, 2.0) - n);

    double factor;

    if (n != 2)
    {
        for (uword m = 0; m < n - 2; ++m)
        {
            factor = sqrt((n - m - 2.0) * (m + 1.0) / (m + 2.0)/(n- m - 1.0));
            I_plus(m + 1) = I_plus(m) * factor + x(m + 1) * sqrt( 1.0 / (m + 2.0) - 1.0 / n);
            I_minus(m + 1) = I_minus(m) / factor - x(m + 1) / sqrt(n * n / (m + 2.0) - n);
        }        
    }

    return I_plus - I_minus;
}
