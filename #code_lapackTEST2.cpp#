/* Description
   Testing function dgesv for solving Ax=B
   Matrix 9x9
 */

// ==>   g++ code_lapackTEST2.cpp -o code_lapackTEST2.exe -L/usr/lib/x86_64-linux-gnu/ -llapacke

#include <iostream>
#include <lapacke.h>
using namespace std;
int main() {
  // k =1. ; dx = 1.; dt = 0.01;
  // Ax = b; A is a tridiagonal matrix: -a b -c
  // a, c - k*dt/(dx*dx) = 1*0.01/1 = 0.01
  // b = 1 + 2*k*dt/(dx*dx) = 1.02

  const int N = 9, NRHS = 1, LDA = N, LDB = N;
  int ipiv[N], info;
  double A[LDA*N] = {
		     1.02, -0.01, 0.,  0.,  0.,  0.,  0.,  0.,  0., 
		     -0.01, 1.02, -0.01,  0.,  0.,  0.,  0.,  0.,  0.,
		     0., -0.01, 1.02, -0.01,  0.,  0.,  0.,  0.,  0.,
		     0.,  0., -0.01, 1.02, -0.01,  0.,  0.,  0.,  0.,
		     0.,  0., 0.,  -0.01, 1.02, -0.01,  0.,  0.,  0.,
		     0.,  0., 0.,  0.,  -0.01, 1.02, -0.01,  0.,  0.,
		     0.,  0., 0.,  0.,  0.,  -0.01, 1.02, -0.01,  0.,
		     0.,  0., 0.,  0.,  0.,  0., -0.01, 1.02, -0.01,
		     0.,  0., 0.,  0.,  0.,  0.,  0.,  -0.01, 1.02
  };
  
  double b[LDB*NRHS] = {
			0.,0.,0.,0.,1.,0.,0.,0.,0.
  };
  info = LAPACKE_dgesv( LAPACK_COL_MAJOR, N, NRHS, A, LDA, ipiv, b, LDB);
  for(int i=0; i<= N-1; i++) {
    cout << b[i] << "\t";
  }
  cout << "\n";
  
}
