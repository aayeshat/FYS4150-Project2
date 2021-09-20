#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>

using namespace arma;
using namespace std;

mat create_tridiagonal(int n)
{
    // Start from identity matrix
    mat A = mat(4, 4, fill::eye);

    A(0, n) = 0.5;
    A(n, 0) = 0.5;
    A(2, 1) = -0.7;
    A(1, 2) = -0.7;

    return A;
}

double max_offdiag_symmetric(arma::mat A, int &k, int &l)

{
    int n = A.size();
    int k, l;
    double maxval = 0.0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            double aij = fabs(A(i, j));
            if (aij > maxval)
            {
                maxval = aij;
                l = i;
                k = j;
            }
        }
    }

    return maxval;
}

int main()
{
    //int n = 4;
    mat A = mat(4, 4, fill::eye);
    A(0, 3) = 0.5;
    A(3, 0) = 0.5;
    A(2, 1) = -0.7;
    A(1, 2) = -0.7;
    //vec a = vec(n - 1)*a;
    //vec d = vec(n) * d;
    //vec e = vec(n - 1)*e;
    // mat A = create_tridiagonal(n, a, d, e);
    //mat I = A.t() * A;


    double max_offdiag = 0;
    cout << A << endl;

    // A.print("A = ");
  
    double max_offdiag = max_offdiag_symmetric(A, int &k, int &l);
    //cout << "n " << n << endl;

    cout << "maxval = " << max_offdiag_symmetric(A, int &k, int &l) << endl;

    return 0;
}