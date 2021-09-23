#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <assert.h>

using namespace arma;
using namespace std;


//create tridiagonal matrix using create_tridiagonal(n, a, d, e)

//while a_kl > epsilon:
//find max off-diagonal element using max_offdiag_symmetric(arma::mat &A, int &k, int &l)
//compute a rotation using jacobi_rotate()
//input matrix A and R
//A is our tridiagonal, R = I identity matrix
//assumes symmetric matrix, so we only loop through upper tridiagonal
//modifies input matrixes A and R


void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l)
{
  assert(A.is_square());
  int n = sqrt(A.size());
  //compute tau, t, c, s
  //chose solution of tau that gives smallest t
  double s, c;
  if (A(k, l) != 0.0)
  {
      double t, tau;
      tau = (A(l, l) - A(k, k)) / (2 * A(k, l));
      if (tau >= 0)
      {
          t = 1.0 / (tau + sqrt(1.0 + tau * tau));
      }
      else
      {
          t = -1.0 / (-tau + sqrt(1.0 + tau * tau));
      }
      c = 1 / sqrt(1 + t * t);
      s = c * t;
  }
  else
  {
      c = 1.0;
      s = 0.0;
  }
  double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
  a_kk = A(k,k)
  a_ll = A(l,l)
  //Change matrix elements with indices kk and ll
  A(k, k) = c * c * a_kk - 2.0 * c * s * A(k, l) + s * s * a_ll;
  A(l, l) = s * s * a_kk + 2.0 * c * s * A(k, l) + c * c * a_ll;
  A(k, l) = 0.0; // hard-coding non-diagonal elements by hand
  A(l, k) = 0.0; // I have written not needed in my lecture notes

  //For all i != k, l
  for (int i = 0; i < n; i++)
  {
      if (i != k && i != l)
      {
          a_ik = A(i, k); //Temporary variable to extract old number
          a_il = A(i, l);
          A(i, k) = c * a_ik - s * a_il;
          A(k, i) = A(i, k);
          A(i, l) = c * a_il + s * a_ik;
          A(l, i) = A(i, l);
      }
      // Updating overall rotation matrix R
      r_ik = R(i, k);
      r_il = R(i, l);
      R(i, k) = c * r_ik - s * r_il;
      R(i, l) = c * r_il + s * r_ik;
  }
  return;
}

// Jacobi method eigensolver:
// - Runs jacobo_rotate until max off-diagonal element < eps
// - Writes the eigenvalues as entries in the vector "eigenvalues"
// - Writes the eigenvectors as columns in the matrix "eigenvectors"
//   (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
// - Stops if it the number of iterations reaches "maxiter"
// - Writes the number of iterations to the integer "iterations"
// - Sets the bool reference "converged" to true if convergence was reached before hitting
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors,
                        const int maxiter, int& iterations, bool& converged)
{



}



//compute tau
