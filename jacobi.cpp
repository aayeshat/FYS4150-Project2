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

    vec a = vec(n - 1).fill(-1.);
    vec d = vec(n).fill(2.);
    vec e = vec(n - 1).fill(-1.);
    mat A = create_tridiagonal(n, a, d, e);
    mat I = A.t() * A;
    vec eigval;
    mat R;
    eig_sym(eigval, R, I);

    A.print("A = ");
    R.print("R = ");

    int iteration = 0;
    int max_iterations = 5;

    double max_offdiag = -1;

    double tolerance = 1;
    cout << "tolerance = " << iteration << endl;

    for (int iteration = 0; iteration < max_iterations; iteration++)
    {

        cout << "jacobi_rotate iteration = " << iteration << endl;

        int i, j;

        max_offdiag = max_offdiag_symmetric(A, i, j);
        cout << "max_offdiag ("
             << "i" << i << "j" << j << ") = " << max_offdiag << endl;

        jacobi_rotate(A, R, i, j);

        if (max_offdiag > tolerance)
        {
            break;
        }
    }

    A.print("A = ");
    R.print("R = ");
    return 0;
}