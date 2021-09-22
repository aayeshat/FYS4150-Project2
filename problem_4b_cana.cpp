#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>

using namespace arma;
using namespace std;

double max_offdiag_symmetric(const arma::mat& A, int& k, int& l)
{
    int n = A.size();
    //int k, l;
    double maxval;
    for (int i = 1; i < n - 1 ; ++i)
    {
        for (int j = i + 1; j < n - 1; ++j)
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

    cout << A << endl;
    int k, l;

    double result = max_offdiag_symmetric(A, k, l);
    cout << "maxwal = " << result << endl;

    return 0;
}

// g++ -c problem_4b.cpp -std=c++11 && g++ -o problem_4b.out problem_4b.o -larmadillo && ./problem_4b.out
