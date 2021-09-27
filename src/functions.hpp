#include <armadillo>

arma::mat create_tridiagonal(int n, const vec &a, const vec &d, const vec &e);

arma::mat create_tridiagonal(int n, double a, double d, double e);

arma::double max_offdiag_symmetric(mat A, int &k, int &l);
