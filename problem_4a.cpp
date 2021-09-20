#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>

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

double max_offdiag_symmetric(arma::mat A, int &k, int &l)

{
    int n = A.size();
    double max;
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            double aij = fabs(A(i, j));
            if (aij > max)
            {
                max = aij;
                k = i;
                l = j;
            }
        }
    }

    return max;
}

int main()
{

    int n = 6;

    vec a = vec(n - 1).fill(-1.);
    vec d = vec(n).fill(2.);
    vec e = vec(n - 1).fill(-1.);
    mat A = create_tridiagonal(n, a, d, e);
    mat I = A.t() * A;

    A.print("A = ");

    cout<<"max = "<<max_offdiag_symmetric<<endl;

    return 0;
}