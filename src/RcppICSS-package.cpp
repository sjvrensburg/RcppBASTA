#include "RcppArmadillo.h"

using arma::colvec;
using arma::uvec;
using Rcpp::List;

//' Centered Cusum of Values
//'
//' @param y Vector of centered, cumulative sum of squares
//' @returns Vector of centered, cumulative sum of squares
// [[Rcpp::export]]
arma::colvec CenteredCusumValues(const arma::colvec& y)
{
    auto const T = y.n_elem;
    colvec Ck = arma::cumsum(arma::square(y));
    double CT = Ck(T-1);
    colvec ks = arma::linspace(1, T, T);

    return Ck / CT - ks / T;
}

//' Check Critical Value
//'
//' @param Dk Vector of centered, cumulative sum of squares
//' @returns List with position and an indicator of significance
// [[Rcpp::export]]
Rcpp::List check_critical_value(const arma::colvec& Dk)
{
    using std::fabs;
    using std::sqrt;
    using Rcpp::_;
    const uvec positions = arma::sort_index(arma::abs(Dk), "descend");
    if (fabs(Dk(positions(0))) == fabs(Dk(positions(1))))
    {
        return List::create(
            _["position"] = NA_INTEGER,
            _["exceeds"] = NA_LOGICAL);
    }
    
    const double M = sqrt(0.5 * Dk.n_elem) * fabs(Dk(positions(0)));
    const bool exceeds = M > 1.358;

    return List::create(
        _["position"] = positions(0),
        _["exceeds"] = exceeds);
}

//' Check Convergence
//'
//' @param oldVec Integer vector of old positions.
//' @param newVec Integer vector of new positions.
//' @returns Boolean that indicates if ICSS algorithm converged.
// [[Rcpp::export]]
bool is_converged(const arma::uvec& oldVec, const arma::uvec& newVec)
{
    int low, high;

    if (oldVec.n_elem == newVec.n_elem)
    {
        for (size_t i = 0; i < newVec.n_elem; i++)
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

arma::colvec ICSS_step_1_and_2(const arma::colvec& x)
{
    colvec change_points();

    auto const Dk = CenteredCusumValues(x);
    auto tmp = check_critical_value(Dk);
    int position = tmp[0];
    
    if (Rcpp::traits::is_nan<INTSXP>(position))
    {
        /* code */
    }
    
}

/*
ICSS_step_1_and_2 <- function(x){

  change_points <- vector(mode = "numeric")

  if (is.null(x)) return(NA)

  Dk <- CenteredCusumValues(x);
  tmp <- check_critical_value(Dk);
  if(length(tmp$position)==0) return(NA)
  exceeds <- tmp$exceeds
  position_step1 <- tmp$position

  if(exceeds){
    position <- position_step1

    while(exceeds){
      t2 <- position
      Dk_step2a = CenteredCusumValues(x[1:t2])
      tmp <- check_critical_value(Dk_step2a)
      if(length(tmp$position)==0) return(NA)
      exceeds <- ifelse(is.na(tmp$exceeds), FALSE, tmp$exceeds)
      position <- tmp$position

    }

    k_first <- t2

    position <- position_step1 + 1
    exceeds <- TRUE

    while(exceeds){
      t1 <- position
      Dk_step2b <- CenteredCusumValues(x[t1:length(x)])
      tmp <- check_critical_value(Dk_step2b)
      if(length(tmp$position)==0) return(NA)
      exceeds <- ifelse(is.na(tmp$exceeds), FALSE, tmp$exceeds)
      position2 <- tmp$position
      position <- position2 + position
    }

    k_last <- t1 - 1

    if (k_first == k_last){
      change_points <- k_first
    }else{
      change_points = c(k_first, k_last)
    }
  }

  return(change_points)
}
*/