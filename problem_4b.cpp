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

arma::mat create_tridiagonal(int n, double a, double d, double e)
{
    arma::vec a_vec = arma::vec(n - 1, arma::fill::ones) * a;
    arma::vec d_vec = arma::vec(n, arma::fill::ones) * d;
    arma::vec e_vec = arma::vec(n - 1, arma::fill::ones) * e;
    return create_tridiagonal(n, a_vec, d_vec, e_vec);
}

double max_offdiag_symmetric(arma::mat A, int &k, int &l)
{

    assert(A.is_square());

    int n = sqrt(A.size());
    cout << "n = " << n << endl;

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

int main()
{

    int n = 4;

    arma::mat A = arma::mat(n, n, fill::eye);
    A(0, 0) = 1;
    A(0, 1) = 0;
    A(0, 3) = 0.5;

    A(1, 2) = -0.7;

    A(2, 1) = -0.7;
    A(2, 3) = 0;

    A(3, 0) = 0.5;

    A.print("A = ");

    int i = 0;
    int j = 0;
    double maxval = max_offdiag_symmetric(A, i, j);

    cout << "maxval ("
         << "i" << i << "j" << j << ") = " << maxval << endl;

    return 0;
}