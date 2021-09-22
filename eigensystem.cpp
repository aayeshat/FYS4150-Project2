#include "tridiagonal.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cmath>

using namespace arma;
using namespace std;


int main()
{

  int n = 6;

  vec a = vec(n - 1).fill(-1.);
  vec d = vec(n).fill(2.);
  vec e = vec(n - 1).fill(-1.);

  // - all n-1 elements on the subdiagonal have value a
  // - all n elements on the diagonal have value d
  // - all n-1 elements on the superdiagonal have value e

  mat A = create_tridiagonal(n, a, d, e);
  A.print("A = ");

  vec eigval;
  mat eigvec;

  eig_sym(eigval, eigvec, A);

  eigval.print("eigenvalues = ");
  eigvec.print("eigvectors = ");
  //norm_eigenvec.print("normalised eigenvectors =");

  //Analytical solutions to eigenvalues
  double pi = 3.14159265358979323846;

  //Eigenvectors
  for (int i = 0; i < n ; ++i){
    vec eigvals_analytical = vec(n);
    eigvals_analytical(i) = d(i) + 2*a(i)*cos((i)*pi/(n + 1));
    i += 1;
    cout << eigvals_analytical << endl;
  }

  return 0;
}
