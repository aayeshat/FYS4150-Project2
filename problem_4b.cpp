#include "tridiagonal.hpp"
#include "max_offdiag.hpp"
#include <iostream>
#include <string>
#include <iomanip>
#include <armadillo>

using namespace arma;
using namespace std;

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
         << "i_" << i << "j_" << j << ") = " << maxval << endl;

    return 0;
}