/* Description
   Testing function dgesv for solving Ax=B
   Matrix 3x3
 */

#include <iostream>
#include <lapacke.h>
using namespace std;
int main() {
  const int N = 3, NRHS = 1, LDA = N, LDB = N;
  int ipiv[N], info;
  // Test 1: Identity matrix A
  // double A[LDA*N] = {
  // 		     1, 0, 0,
  // 		     0, 1, 0,
  // 		     0, 0, 1
  // };

  // Test 2: arbitrary matrix B
  double A[LDA*N] = {
		     -2, 4, 0,
		     0, 1, 0,
		     0, 0, 1
  };
  
  double b[LDB*NRHS] = {
			1.7, 0.5, -3.2
  };
  info = LAPACKE_dgesv( LAPACK_COL_MAJOR, N, NRHS, A, LDA, ipiv, b, LDB);
  for(int i=0; i<= N-1; i++) {
    cout << b[i] << "\t";
  }
  cout << "\n";
  
}
