#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>

using namespace arma;
using namespace std;

// mat A(6,6).fill.(0.));
// mat B = A.t()*A; //generate symmetrix matrix

int n = 6;

// Create tridiagonal matrix from vectors.
// - lower diagonal: vector a, lenght n-1
// - main diagonal:  vector d, lenght n
// - upper diagonal: vector e, lenght n-1

vec a = vec(n - 1).fill(-1.);
vec d = vec(n).fill(2.);
vec e = vec(n - 1).fill(-1.);
arma::mat create_tridiagonal(const arma::vec& a, const arma::vec& d, const arma::vec& e)
{
  // Start from identity matrix
  arma::mat A = arma::mat(n, n, fill::eye);

  // Fill first row (row index 0)

  // Loop that fills rows 2 to n-1 (row indices 1 to n-2)

  // Fill last row (row index n-1)

  return A;
}

// Create a tridiagonal matrix tridiag(a,d,e) of size n*n,
// from scalar input a, d, and e. That is, create a matrix where
// - all n-1 elements on the subdiagonal have value a
// - all n elements on the diagonal have value d
// - all n-1 elements on the superdiagonal have value e
arma::mat create_tridiagonal(int n, double a, double d, double e)
{
  // Start from identity matrix
  arma::mat A = arma::mat(n, n, arma::fill::eye);

  // Fill the first row (row index 0), e.g.
  A(0,0) = d;
  A(0,1) = e;
  A(1,0) = a;

  // Loop that fills rows 2 to n-1 (row indices 1 to n-2)
  for (int i = 1; i < n - 1; i++){
    double d = A(i, i);
  }

  // Fill last row (row index n-1)
  for (int i = 1; i < n - 2; i++){
    double e = A(i, i + 1);
    double a = (A(i + 1, i));
  }

  return A;
  cout << A << endl;
}

// // Create a symmetric tridiagonal matrix tridiag(a,d,a) of size n*n
// // from scalar input a and d.
// arma::mat create_symmetric_tridiagonal(int n, double a, double d)
// {
//   // Call create_tridiagonal and return the result
//   return create_tridiagonal(n, a, d, a);
// }


// g++ -c tridiag.cpp -std=c++11 && g++ -o tridiag.out tridiag.o -larmadillo && ./tridiag.out
