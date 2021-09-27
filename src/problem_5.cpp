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
    eigval.print("Eigenvalues_eigsym = ");

    int iteration = 0;
    int max_iterations = 1000;
    int i, j;
    double max_offdiag = max_offdiag_symmetric(A, i, j);

    double tolerance = 1E-30;
    R = mat(n,n,fill::eye);
    for (int iteration = 0; iteration < max_iterations; iteration++)
    {

        jacobi_rotate(A, R, i, j);
        max_offdiag = max_offdiag_symmetric(A, i, j);
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
