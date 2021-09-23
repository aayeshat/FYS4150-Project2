# project2

### problem tridiag

`g++ tridiag.cpp tridiagonal.cpp -larmadillo`

### problem eigensystem

`g++ eigensystem.cpp tridiagonal.cpp -larmadillo`

### problem 3

`g++ problem_3a.cpp tridiagonal.cpp -larmadillo`

### problem 4a

`g++ problem_4a.cpp tridiagonal.cpp max_offdiag.cpp -larmadillo`

### problem 4b

`g++ problem_4b.cpp tridiagonal.cpp max_offdiag.cpp -larmadillo`

### problem jacobi_rotate

`g++ jacobi.cpp tridiagonal.cpp max_offdiag.cpp jacobi_rotate.cpp -larmadillo`

#### mac compiler 

`g++ -c tridiag.cpp tridiagonal.cpp -std=c++11 && g++ -o tridiag.out tridiag.o -larmadillo && ./tridiag.out`