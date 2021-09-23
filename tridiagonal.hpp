#include <armadillo>

using namespace arma;

mat create_tridiagonal(int n, const vec &a, const vec &d, const vec &e);

mat create_tridiagonal(int n, double a, double d, double e);