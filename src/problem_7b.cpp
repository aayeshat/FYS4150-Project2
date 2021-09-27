#include "problem_7.hpp"
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
    int n = 100;
    int max_iterations = 2500;
    problem_7(n, max_iterations);
    return 0;
}
