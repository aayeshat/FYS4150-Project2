#include "max_offdiag.hpp"
#include <armadillo>
#include <assert.h>

using namespace arma;

double max_offdiag_symmetric(mat A, int &k, int &l)
{

    assert(A.is_square());

    int n = sqrt(A.size());
    //cout << "n = " << n << endl;

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