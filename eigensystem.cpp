#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cmath>

using namespace arma;
using namespace std;

mat create_tridiagonal(int n, const vec &a, const vec &d, const vec &e)
{
  // Start from identity matrix
  mat A = mat(n, n, fill::eye);

  for (int i = 0; i < n - 1; i++)
  {
    A(i, i) = d(i);
    A(i, i + 1) = e(i);
    A(i + 1, i) = a(i);
  }
  A(n - 1, n - 1) = d(n - 1);
  return A;
}

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
  double cos_arg = pi/(n + 1.0);

  //Eigenvectors
  for (int i = 0; i < n ; ++i){
    vec eigvals_analytical = vec(n);
    eigvals_analytical(i) = d(i) + 2*a(i)*cos((cos_arg)*(i + 1));
  cout << "analytical =" << endl << eigvals_analytical << endl;
  //eigvals_analytical.print("Analytical eigenvalues = ");
  }

  return 0;
}

//g++ -c eigensystem.cpp -std=c++11 && g++ -o eigensystem.out eigensystem.o -larmadillo && ./eigensystem.out
