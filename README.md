# FYS4150 Project2

### problem tridiag

`g++ tridiag.cpp tridiagonal.cpp -larmadillo && ./a.out`

### problem eigensystem

`g++ eigensystem.cpp tridiagonal.cpp -larmadillo && ./a.out`

### problem 3

`g++ problem_3a.cpp tridiagonal.cpp -larmadillo && ./a.out`

### problem 4a

`g++ problem_4a.cpp tridiagonal.cpp max_offdiag.cpp -larmadillo && ./a.out`
### problem 4b

`g++ problem_4b.cpp tridiagonal.cpp max_offdiag.cpp -larmadillo && ./a.out`

### problem 5

`g++ problem_5.cpp tridiagonal.cpp max_offdiag.cpp jacobi_rotate.cpp -larmadillo && ./a.out`

### problem 6

`g++ problem_6.cpp tridiagonal.cpp max_offdiag.cpp jacobi_rotate.cpp -larmadillo && ./a.out`

`python3 problem6_plot.py`

### problem 7

`g++ problem_7a.cpp tridiagonal.cpp max_offdiag.cpp jacobi_rotate.cpp problem_7.cpp -larmadillo && ./a.out`

`g++ problem_7b.cpp tridiagonal.cpp max_offdiag.cpp jacobi_rotate.cpp problem_7.cpp -larmadillo && ./a.out`

`python3 problem6_plot.py`

#### mac compiler 

`g++ -c tridiag.cpp tridiagonal.cpp -std=c++11 && g++ -o tridiag.out tridiag.o -larmadillo && ./tridiag.out`
