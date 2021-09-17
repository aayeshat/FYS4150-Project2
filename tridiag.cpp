#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>

using namespace arma;
using namespace std;

arma::mat create_tridiagonal(int n, const arma::vec &a, const arma::vec &d, const arma::vec &e){
    // Start from identity matrix
    arma::mat A = arma::mat(n, n, fill::eye);

    for (int i = 0; i < n-1; i++){
      A(i, i) = d(i);
      A(i, i+1) = e(i);
      A(i+1, i) = a(i);
    }
    A(n-1, n-1) = d(n-1);
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
  arma::mat A = create_tridiagonal(n, a, d, e);
  A.print("A = ");
}
