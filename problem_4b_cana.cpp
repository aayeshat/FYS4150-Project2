#include "tridiagonal.hpp"
#include "max_offdiag.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>

using namespace arma;
using namespace std;

int main()
{
    //int n = 4;
    mat A = mat(4, 4, fill::eye);
    A(0, 3) = 0.5;
    A(3, 0) = 0.5;
    A(2, 1) = -0.7;
    A(1, 2) = -0.7;

    cout << A << endl;
    int k, l;

    double result = max_offdiag_symmetric(A, k, l);
    cout << "maxwal = " << result << endl;

    return 0;
}
