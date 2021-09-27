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
  double h = 1/double(n);
  double h_2 = h*h;
  double a = -1./h_2;
  double d = 2./h_2;
  double e = -1./h_2;

  // - all n-1 elements on the subdiagonal have value a
  // - all n elements on the diagonal have value d
  // - all n-1 elements on the superdiagonal have value e

  mat A = create_tridiagonal(n, a, d, e);
  A.print("A = ");

  vec eigval;
  mat eigvec;

  eig_sym(eigval, eigvec, A);

  eigval.print("Armadillo eigenvalues = ");
  eigvec.print("Armadillo eigenvectors = ");
  //norm_eigenvec.print("normalised eigenvectors =");

  //Analytical solutions to eigenvalues
  double pi = 3.14159265358979323846;
  double arg = pi/(n + 1.0);

  int width = 12;
  int prec = 4;


  vec eigvals_analytical = vec(n);

  //Eigenvalues
  for (int i = 0; i < n ; ++i){
    eigvals_analytical(i) = d + 2*a*cos((arg)*(i + 1));
    //cout << setw(width) << setprecision(prec) << std::scientific << eigvals_analytical(i) << endl;
  }
  eigvals_analytical.print("Analytical eigenvalues = ");

  //Eigenvectors
  mat eigvec_analytical = mat(n,n);
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      eigvec_analytical(i,j) = sin(arg*(i+1)*(j+1));
    }
  }
  normalise(eigvec_analytical).print("Analytical eigenvectors = ");
//
//
return 0;
}

//g++ -c eigensystem_cana.cpp -std=c++11 && g++ -o eigensystem_cana.out eigensystem_cana.o -larmadillo && ./eigensystem_cana.out
