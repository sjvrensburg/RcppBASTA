#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <boost/math/special_functions/sign.hpp>

using Rcpp::_;
using Rcpp::as;
using Rcpp::IntegerVector;
using Rcpp::wrap;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

// [[Rcpp::export(name = "stat.resid")]]
Rcpp::NumericVector stat_resid(const Rcpp::NumericVector& x,
    Rcpp::Nullable<NumericVector> a_ = R_NilValue,
    unsigned int order = 1,
    const double factor = 8,
    const double epsilon = 0.0)
{
    NumericVector a;
    auto const n = x.length();

    // Determine order if a_ is not NULL
    if (a_.isNotNull())
    {
        a = a_;
        order = a.length() - 1;
    }

    if (a_.isNull())
    {
        Rcpp::Function f("stat.est");
        a = f(x, _["order"] = order);

        for (R_xlen_t i = 1; i < a.length(); ++i)
        {
            a(i) = a(i) / factor;
        }
    }

    NumericVector y(n-order);
    R_xlen_t n_y = y.length();

    double numerator, denominator;

    for (R_xlen_t i = 0; i < n_y; ++i)
    {
        numerator = (x(n - i - 1)) * (x(n - i - 1));
        
        IntegerVector r = Rcpp::rev(Rcpp::seq(n - 1 - i - order, n - 1 - i));
        NumericVector segment(order);
        for (R_xlen_t j = 0; j < order; ++j)
        {
            segment(order - 1 - j) = x(j) * x(j);
        }
        
        NumericVector tmp(order + 1);
        tmp(0) = 1;
        for (R_xlen_t j = 1; j < order + 1; ++j)
        {
            tmp(j) = segment(j-1);
        }

        denominator = Rcpp::sum(a * tmp);
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
            tree_j(1, k) *= abs(tree_j(1, k)) > th ? 1.0 : 0.0;
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
Rcpp::NumericVector inner_prod_iter(Rcpp::NumericVector x)
{
    using std::sqrt;
    using Rcpp::tail;
    using Rcpp::sum;

    R_xlen_t n = x.length();
    NumericVector I_plus(n-1); I_plus.fill(0.0);
    NumericVector I_minus(n-1); I_minus.fill(0.0);

    I_plus(0) = sqrt(1.0 - 1.0 / n) * x(0);
    I_minus(0) = sum(tail(x, n - 1)) / sqrt(std::pow(n, 2.0) - n);

    double factor;

    if (n != 2)
    {
        for (R_xlen_t m = 0; m < n - 2; ++m)
        {
            factor = sqrt((n - m - 2.0) * (m + 1.0) / (m + 2.0)/(n- m - 1.0));
            I_plus(m + 1) = I_plus(m) * factor + x(m + 1) * sqrt( 1.0 / (m + 2.0) - 1.0 / n);
            I_minus(m + 1) = I_minus(m) / factor - x(m + 1) / sqrt(n * n / (m + 2.0) - n);
        }        
    }

    return I_plus - I_minus;
}

// [[Rcpp::export]]
double med(Rcpp::NumericVector x)
{
    if (x.length() == 1) return x(0);

    NumericVector y = Rcpp::clone(x);
    std::sort(y.begin(), y.end());

    R_xlen_t j = std::floor(0.5 * (x.length() - 1));
    const double g = 0.5 * (x.length() - 1.0) - j;

    if ((std::fabs(g) <= 1E-8) and (j % 2 == 0)) return y(j + 1);

    return y(j);
}