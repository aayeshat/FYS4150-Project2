#include "tridiagonal.hpp"
#include "max_offdiag.hpp"
#include "jacobi_rotate.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <assert.h>

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


mat create_tridiagonal(int n, double a, double d, double e)
{
  vec a_vec = vec(n - 1, fill::ones) * a;
  vec d_vec = vec(n, fill::ones) * d;
  vec e_vec = vec(n - 1, fill::ones) * e;
  return create_tridiagonal(n, a_vec, d_vec, e_vec);
}

double max_offdiag_symmetric(mat A, int &k, int &l)
{

    assert(A.is_square());

    int n = sqrt(A.size());

    //assigning minimum value to variable
    double maxval = -1;
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            double aij = fabs(A(i, j));
            if (aij > maxval)
            {
                maxval = aij;
                k = i;
                l = j;
            }
        }
    }

    return maxval;
}

void jacobi_rotate(arma::mat &A, arma::mat &R, int k, int l)
{
    assert(A.is_square());
    int n = sqrt(A.size());

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
    a_kk = A(k, k);
    a_ll = A(l, l);
    A(k, k) = c * c * a_kk - 2.0 * c * s * A(k, l) + s * s * a_ll;
    A(l, l) = s * s * a_kk + 2.0 * c * s * A(k, l) + c * c * a_ll;
    A(k, l) = 0.0; // hard-coding non-diagonal elements by hand
    A(l, k) = 0.0; // same here
    for (int i = 0; i < n; i++)
    {
        if (i != k && i != l)
        {
            a_ik = A(i, k);
            a_il = A(i, l);
            A(i, k) = c * a_ik - s * a_il;
            A(k, i) = A(i, k);
            A(i, l) = c * a_il + s * a_ik;
            A(l, i) = A(i, l);
        }
        //  eigenvectors
        r_ik = R(i, k);
        r_il = R(i, l);
        R(i, k) = c * r_ik - s * r_il;
        R(i, l) = c * r_il + s * r_ik;
    }
    return;
} // end of function jacobi_rotate

int main()
{

    int n = 6;
    double h = 1 / double(n);
    double h_2 = h * h;

    double a = -1./h_2;
    double d = 2./h_2;
    double e = -1./h_2;
    mat A = create_tridiagonal(n, a, d, e);
    vec eigval;
    mat R;
    eig_sym(eigval, R, A);

    A.print("A = ");
    normalise(R).print("R_eigsym = ");

    int iteration = 0;
    int max_iterations = 1000;
    int i, j;
    double max_offdiag = max_offdiag_symmetric(A, i, j);

    double tolerance = 1E-30;
    //cout << "tolerance = " << iteration << endl;
    R = mat(n,n,fill::eye);
    for (int iteration = 0; iteration < max_iterations; iteration++)
    {

        jacobi_rotate(A, R, i, j);
        max_offdiag = max_offdiag_symmetric(A, i, j);
        // cout << "max_offdiag ("
        //      << "i" << i << "j" << j << ") = " << max_offdiag << endl;

        // jacobi_rotate(A, R, i, j);
        //cout << "jacobi_rotate iteration = " << iteration << endl;

        if (max_offdiag < tolerance)
        {
            break;
        }
    }

    A.print("A = ");
    normalise(R,2,1).print("R = ");
    vec eigenvals;
    eigenvals = diagvec(A);
    eigenvals.print("Eigenvectors = ");
    return 0;
}
