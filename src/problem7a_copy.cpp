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
    int n = 10;

    double h = 1 / double(n);
    double h_2 = h * h;

    vec x = linspace(0, 1, n + 2);
    // vec x = arma::vec(n);
    // x(0) = 0;
    // for (int i = 0; i < n; i++)
    // {
    //     x(i) = x(0) + i * h;
    // }
    x.print("X =");

    double a = -1./h_2;
    double d = 2./h_2;
    double e = -1./h_2;
    mat A = create_tridiagonal(n, a, d, e);
    vec eigval;
    mat R;
    eig_sym(eigval, R, A);

    A.print("A = ");
    normalise(R).print("R_eigsym = ");
    eigval.print("Eigenvalues_eigsym = ");

    int iteration = 0;
    int max_iterations = 25000;
    int i, j;
    double max_offdiag = max_offdiag_symmetric(A, i, j);

    double tolerance = 1E-30;
    cout << "tolerance = " << iteration << endl;
    R = mat(n, n, fill::eye);
    for (int iteration = 0; iteration < max_iterations; iteration++)
    {

        jacobi_rotate(A, R, i, j);
        max_offdiag = max_offdiag_symmetric(A, i, j);

        if (max_offdiag < tolerance)
        {
            break;
        }
    }

    A.print("A = ");
    vec eigenvals = diagvec(A);
    eigenvals.print("Eigenvalues = ");

    //normalise(R, 2, 1).print("Eigen vector = ");
    R.print("Eigen vector = ");

    cout << endl
         << "Minimum eigen values" << endl
         << endl;
    uvec sort_indexes = sort_index(eigenvals);

    vec v[3];
    for (int s = n - 1; s >= n - 3; s--)
    {

        int index = sort_indexes(s);
        double eigenval = eigenvals(index);

        cout << "index " << index << " --- eigen value " << scientific << eigenval << endl;
        vec indexV = matric_to_vector(n, R, index);
        v[n - s - 1] = indexV;
        //v.print("V =");

        cout << endl;
    }

    int width = 18;
    int prec = 4;
    string filename = "../out/problem7_n" + to_string(n) + ".txt";
    ofstream ofile;
    ofile.open(filename);

    double x_n = 1.0;
    double v_n = 0.0;
    ofile
        << std::setw(width) << std::setprecision(prec) << std::scientific << v_n
        << std::setw(width) << std::setprecision(prec) << std::scientific << v_n
        << std::setw(width) << std::setprecision(prec) << std::scientific << v_n
        << std::endl;
    for (int i = 0; i < n; i++)
    {
      cout
            << std::setw(width) << std::setprecision(prec) << std::scientific << v[0](i)
             << std::setw(width) << std::setprecision(prec) << std::scientific << v[1](i)
             << std::setw(width) << std::setprecision(prec) << std::scientific << v[2](i)
             << std::endl;

        ofile
              << std::setw(width) << std::setprecision(prec) << std::scientific << v[0](i) << 0
              << std::setw(width) << std::setprecision(prec) << std::scientific << v[1](i) << 0
              << std::setw(width) << std::setprecision(prec) << std::scientific << v[2](i) << 0
              << std::endl;
    }
    ofile
        << std::setw(width) << std::setprecision(prec) << std::scientific << v_n
        << std::setw(width) << std::setprecision(prec) << std::scientific << v_n
        << std::setw(width) << std::setprecision(prec) << std::scientific << v_n
        << std::endl;

    ofile.close();


    string filename2 = "../out/problem7_x_n" + to_string(n) + ".txt";
    ofile.open(filename2);
    ofile
          << std::setw(width) << std::setprecision(prec) << std::scientific << x << endl;
    ofile.close();

    return 0;
}
