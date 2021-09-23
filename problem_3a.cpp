#include "tridiagonal.hpp"
#include <iostream>
#include <string>
#include <iomanip>
#include <armadillo>


using namespace arma;
using namespace std;


int main()
{





  int n = 6;
double h = 1/double(n);
double h_2 = h*h;


  vec a = vec(n - 1).fill(-1.)*(1/h_2);
  vec d = vec(n).fill(2.)*(1/h_2);
  vec e = vec(n - 1).fill(-1.)*(1/h_2);

  // - all n-1 elements on the subdiagonal have value a
  // - all n elements on the diagonal have value d
  // - all n-1 elements on the superdiagonal have value e

  mat A = create_tridiagonal(n, a, d, e);
  A.print("A = ");

  mat I = A.t() * A; // generate a symmetric matrix

  //I.print("I = ");

  vec eigval;
  mat eigvec;

  eig_sym(eigval, eigvec, A);

  eigval.print("eigval = ");
  eigvec.print("eigvec = ");

  return 0;
}
