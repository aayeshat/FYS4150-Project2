#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>

using namespace arma;
using namespace std;

arma::mat create_tridiagonal(int n, const arma::vec &a, const arma::vec &d, const arma::vec &e)
{
  // Start from identity matrix
  arma::mat A = arma::mat(n, n, fill::eye);

  for (int i = 0; i < n - 1; i++)
  {
    A(i, i) = d(i);
    A(i, i + 1) = e(i);
    A(i + 1, i) = a(i);
  }
  A(n - 1, n - 1) = d(n - 1);
  return A;
}

arma::mat create_tridiagonal(int n, double a, double d, double e)
{
  arma::vec a_vec = arma::vec(n - 1, arma::fill::ones) * a;
  arma::vec d_vec = arma::vec(n, arma::fill::ones) * d;
  arma::vec e_vec = arma::vec(n - 1, arma::fill::ones) * e;
  return create_tridiagonal(n, a_vec, d_vec, e_vec);
}

int main()
{

  int n = 6;

  // - all n-1 elements on the subdiagonal have value a
  // - all n elements on the diagonal have value d
  // - all n-1 elements on the superdiagonal have value e
  arma::mat A = create_tridiagonal(n, -1, 2., -1);
  A.print("A = ");

  mat B =
  vec eigval;
  vec eigvec;
  eig_sym(eig_val, eigvec, B)

  return 0;
}

//g++ -c tridiag.cpp -std=c++11 && g++ -o tridiag.out tridiag.o -larmadillo && ./tridiag.out
