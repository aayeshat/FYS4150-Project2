#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>

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

  mat I = A.t() * A; // generate a symmetric matrix

  I.print("I = ");

  vec eigval;
  mat eigvec;

  eig_sym(eigval, eigvec, A);

  eigval.print("eigval = ");
  eigvec.print("eigvec = ");

  return 0;
}
