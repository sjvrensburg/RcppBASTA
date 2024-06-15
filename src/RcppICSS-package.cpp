#include <Rcpp.h>

using Rcpp::_;
using Rcpp::as;
using Rcpp::List;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;

//' Centered Cusum of Values
//'
//' @param y Vector of centered, cumulative sum of squares
//' @returns Vector of centered, cumulative sum of squares
// [[Rcpp::export]]
Rcpp::NumericVector CenteredCusumValues(const Rcpp::NumericVector& y)
{
    const int T = y.length();
    NumericVector Ck = Rcpp::cumsum(y * y);
    double CT = Ck(T-1);
    IntegerVector ks = Rcpp::seq_len(T);
    NumericVector t1 = Ck / CT;
    NumericVector t2 = as<NumericVector>(ks) / T;

    return t1 - t2;
}

//' Check Critical Value
//'
//' @param Dk Vector of centered, cumulative sum of squares
//' @returns List with position and an indicator of significance
// [[Rcpp::export]]
Rcpp::List check_critical_value(const Rcpp::NumericVector& Dk)
{
    using std::fabs;
    using std::sqrt;
    using Rcpp::abs;
    using Rcpp::max;
    using Rcpp::which_max;

    const NumericVector Dk_abs = abs(Dk);
    const double value = max(Dk_abs);
    const unsigned int position = which_max(Dk_abs) + 1;

    Rcpp::LogicalVector comparison = Dk_abs >= value;
    if (Rcpp::sum(comparison) > 1)
    {
        return List::create(
            _["position"] = NA_INTEGER,
            _["exceeds"] = NA_LOGICAL);
    }
    
    const double M = sqrt(0.5 * Dk.length()) * fabs(value);
    const bool exceeds = M > 1.358;

    return List::create(
        _["position"] = position,
        _["exceeds"] = exceeds);
}

//' Check Convergence
//'
//' @param oldVec Integer vector of old positions.
//' @param newVec Integer vector of new positions.
//' @returns Boolean that indicates if ICSS algorithm converged.
// [[Rcpp::export]]
bool is_converged(const Rcpp::IntegerVector& oldVec, const Rcpp::IntegerVector& newVec)
{
    int low, high;

    if (oldVec.length() == newVec.length())
    {
        for (size_t i = 0; i < newVec.length(); i++)
        {
            low = oldVec(i) < newVec(i) ? oldVec(i) : newVec(i);
            high = oldVec(i) > newVec(i) ? oldVec(i) : newVec(i);
            if ((high - low) > 2) return false;       
        }
    } else {
        return false;
    }
    
    return true;
}

// [[Rcpp::export]]
Rcpp::IntegerVector ICSS_step_1_and_2(const Rcpp::NumericVector &x)
{
    IntegerVector change_points(2);

    NumericVector Dk = CenteredCusumValues(x);
    Rcpp::List tmp = check_critical_value(Dk);
    bool exceeds = tmp["exceeds"];
    int position_step1 = tmp["position"];

    int t1, t2, k_first, k_last, position_step2;
    IntegerVector r;
    NumericVector Dk_step2a;
    NumericVector Dk_step2b;
    
    int position = position_step1;

    k_first = 1;
    k_last = position;

    Rcpp::Rcout << "exceeds = (" << exceeds << ")\tposition = (" << position << ")\n";
    
    if (exceeds)
    {
        while (exceeds)
        {
            Rcpp::checkUserInterrupt();
            t2 = position;
            r = Rcpp::Range(0, t2-1);
            Dk_step2a = CenteredCusumValues(x[r]);
            tmp = check_critical_value(Dk_step2a);

            exceeds = tmp["exceeds"];
            position = tmp["position"];

            Rcpp::Rcout << "r = (" << r(0) << ", ..., " << r(r.length() - 1) << ")\n";
            Rcpp::Rcout << "exceeds = (" << exceeds << ")\tposition = (" << position << ")\n";

            if (R_IsNA(position)) Rcpp::stop("check_critical_value returned NA");
            if (position <= 0) Rcpp::stop("check_critical_value returned non-positive position");
        }

        k_first = t2 + 1;
        position = position_step1 + 1;
        exceeds = true;

        while (exceeds)
        {
            Rcpp::checkUserInterrupt();
            t1 = position;
            r = Rcpp::Range(t1-1, x.length()-1);

            Dk_step2b = CenteredCusumValues(x[r]);
            tmp = check_critical_value(Dk_step2b);

            exceeds = tmp["exceeds"];
            position_step2 = tmp["position"];

            Rcpp::Rcout << "r = (" << r(0) << ", ..., " << r(r.length() - 1) << ")\n";
            Rcpp::Rcout << "exceeds = (" << exceeds << ")\tposition_step2 = (" << position_step2 << ")\n";

            if (R_IsNA(position_step2)) Rcpp::stop("check_critical_value returned NA");
            if (position_step2 <= 0) Rcpp::stop("check_critical_value returned non-positive position");

            position = position_step2 + position;
            Rcpp::Rcout << "position = " << position << "\n";
        }

        k_last = t1;
    }

    change_points[0] = k_first;
    change_points[1] = k_last;

    return change_points;
}
