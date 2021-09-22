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

    int n = 6;

    vec a = vec(n - 1).fill(-1.);
    vec d = vec(n).fill(2.);
    vec e = vec(n - 1).fill(-1.);
    mat A = create_tridiagonal(n, a, d, e);

   

    int i = 0;
    int j = 0;
    double maxval = max_offdiag_symmetric(A, i, j);

    A.print("A = ");

    cout << "maxval ("
         << "i" << i << "j" << j << ") = " << maxval << endl;

    return 0;
}