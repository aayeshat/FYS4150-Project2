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
<<<<<<< HEAD

// g++ -c problem_4b.cpp -std=c++11 && g++ -o problem_4b.out problem_4b.o -larmadillo && ./problem_4b.out
//g++ -c problem_4a.cpp -std=c++11 && g++ -o problem_4a.out problem_4a.o -larmadillo && ./problem_4a.out
=======
>>>>>>> e61816766f769bc7228fc3a02462c8dae527d5d5
