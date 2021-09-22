#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>

using namespace arma;
using namespace std;

mat create_jacobional(int n, const vec &a, const vec &d, const vec &e)
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

double max_offdiag_symmetric(arma::mat &A, int &k, int &l)
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

void jacobi_rotate(arma::mat &A, arma::mat &R, int k, int l)
{
    int n = A.size();

    double s, c;
    if (A(k, l) != 0.0)
    {
        double t, tau;
        tau = (A(l, l) - A(k, k)) / (2 * A(k, l));
        if (tau >= 0)
        {
            t = 1.0 / (tau + sqrt(1.0 + tau * tau));
        }
        else
        {
            t = -1.0 / (-tau + sqrt(1.0 + tau * tau));
        }
        c = 1 / sqrt(1 + t * t);
        s = c * t;
    }
    else
    {
        c = 1.0;
        s = 0.0;
    }
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k, k);
    a_ll = A(l, l);
    A(k, k) = c * c * a_kk - 2.0 * c * s * A(k, l) + s * s * a_ll;
    A(l, l) = s * s * a_kk + 2.0 * c * s * A(k, l) + c * c * a_ll;
    A(k, l) = 0.0; // hard-coding non-diagonal elements by hand
    A(l, k) = 0.0; // same here
    for (int i = 0; i < n; i++)
    {
        if (i != k && i != l)
        {
            a_ik = A(i, k);
            a_il = A(i, l);
            A(i, k) = c * a_ik - s * a_il;
            A(k, i) = A(i, k);
            A(i, l) = c * a_il + s * a_ik;
            A(l, i) = A(i, l);
        }
        // And finally the new eigenvectors
        r_ik = R(i, k);
        r_il = R(i, l);
        R(i, k) = c * r_ik - s * r_il;
        R(i, l) = c * r_il + s * r_ik;
    }
    return;
} // end of function jacobi_rotate

int main()
{

    int n = 6;

    vec a = vec(n - 1).fill(-1.);
    vec d = vec(n).fill(2.);
    vec e = vec(n - 1).fill(-1.);
    mat A = create_jacobional(n, a, d, e);
    mat I = A.t() * A;
    vec eigval;
    mat R;
    eig_sym(eigval, R, I);

   // A.print("A = ");
   // R.print("R = ");
    double tolerance = 1.0E-8;

    int max_iterations = 5;
    double max_offdiag = 0;
    int iterations = 0;
    while (max_offdiag > tolerance && iterations <= max_iterations)
    {
        int k, l;

        max_offdiag = max_offdiag_symmetric(A, k, l);

        jacobi_rotate(A, R, k, l);

        iterations++;
        cout<<"jacobi_rotate"<<jacobi_rotate<<endl;
    }

    A.print("A = ");
    R.print("R = ");
    return 0;
}


//g++ -c jacobi.cpp -std=c++11 && g++ -o jacobi.out jacobi.o -larmadillo && ./jacobi.out
