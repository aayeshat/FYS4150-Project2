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
    double h;
    double h_2;

    double a;
    double d;
    double e;
    mat A;
    vec eigval;
    mat R;
    //eig_sym(eigval, R, A);

    // A.print("A = ");
    // normalise(R).print("R_eigsym = ");

    int iteration = 0;
    int max_iterations = 10000;
    int i, j;
    double max_offdiag = max_offdiag_symmetric(A, i, j);

    double tolerance = 1E-30;
    cout << "tolerance = " << iteration << endl;
    //R = mat(n,n,fill::eye);

    for (int n = 1; n < 100; ++n)
    {
        double h = 1 / double(n);
        double h_2 = h * h;
        double a = -1. / h_2;
        double d = 2. / h_2;
        double e = -1. / h_2;
        mat A = create_tridiagonal(n, a, d, e);
        R = mat(n, n, fill::eye);
        for (int iteration = 0; iteration < max_iterations; iteration++)
        {
            jacobi_rotate(A, R, i, j);
            max_offdiag = max_offdiag_symmetric(A, i, j);
            // cout << "max_offdiag ("
            //      << "i" << i << "j" << j << ") = " << max_offdiag << endl;

            // jacobi_rotate(A, R, i, j);
            if (max_offdiag < tolerance)
            {
                break;
            }
            //cout << "jacobi_rotate iteration = " << iteration << endl;
        }
        n += 10 - 1;
        cout << "n = " << n << endl;
    }

    // A.print("A = ");
    //normalise(R,2,1).print("R = ");
    // vec eigenvals;
    // eigenvals = diagvec(A);
    // eigenvals.print("Eigenvectors = ");
    return 0;
}

//g++ -c problem_6.cpp -std=c++11 && g++ -o problem_6.out problem_6.o -larmadillo && ./problem_6.out
//g++ -c problem5_cana.cpp -std=c++11 && g++ -o problem5_cana.out problem5_cana.o -larmadillo && ./problem5_cana.out
