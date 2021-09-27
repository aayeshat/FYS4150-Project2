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

  // - all n-1 elements on the subdiagonal have value a
  // - all n elements on the diagonal have value d
  // - all n-1 elements on the superdiagonal have value e
  mat A = create_tridiagonal(n, -1, 2., -1);
  A.print("A = ");

  return 0;
}
