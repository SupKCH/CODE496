
// Practice creating variable (phi) array

// g++ practice_init.cpp -o practice_init.exe -L/usr/lib/x86_64-linux-gpu/ -llapack

#include <iostream>
#include <lapacke.h>
#include <fstream>
using namespace std;



/* =================== Subroutine ===================== */
void node_initialize(double* var, int nx_func) {

  // initialize a 1D heat steel rod
  
  for (int i = 0; i <= nx_func-1; i++) {
    var[i] = 0.;
  }
}

void node_1DHeatCondition(double* var) {

  // Let's set temperature at one end of the rod to be constant at 100degC
  var[0] = 100.;

  // No fixed condition at the other end, just 0 as the other nodes
  
}

void visualize(double* var, int nx_func) {
  
  for (int i = 0; i <= nx_func-1; i++) {
    cout << var[i] << ' ';
  }
  cout << "\n";

}

void explicit_sim(double* var, double* var_new, int nx_func, double coeff, double dt, double dx) {
  // Central finite difference except two ends of the rod
  for (int i = 1; i <= nx_func - 2; i++) {
    var_new[i] = var[i] + coeff*dt*(var[i+1] - 2*var[i] + var[i-1])/(dx*dx);
  }

  var_new[0] = var[0];
  var_new[nx_func-1] = var[nx_func-1] + coeff*dt*(var[nx_func-3] - 2*var[nx_func-2] + var[nx_func-1])/(dx*dx);

  for (int i = 0; i <= nx_func - 1; i++) {
    var[i] = var_new[i];
  }
}

void implicit_sim(double* var, double* var_new, int nx_func, double coeff, double dt, double dx, double* b_func) {

  // We won't solve on boundary condition node
  const int N = nx_func-1, NRHS = 1, LDA = N, LDB = N;
  int ipiv[N], info;

  double a1 = -coeff*dt/(dx*dx);
  double a2 = 1 + 2*(-a1);

  double A[LDA*N];
  for (int i = 0; i < N*N; i++) {
    A[i] = 0.0;
  }
  
  for (int i = 0; i <= N*N-1; i++) {
    if (i == 0) {
      A[0] = a2;
      A[1] = a1;
    }
    else if (i % N == 0 && i != N*(N-1)) {
      A[i + int(i/N)-1] = a1;
      A[i+1 + int(i/N)-1] = a2;
      A[i+2 + int(i/N)-1] = a1;
    }
    else if (i == N*(N-1)) {
      A[N*N -3] = a1;
      A[N*N -2] = 2*(-a1);
      A[N*N -1] = (1 + a1);
    }
  }

  /*
  double b[LDB*NRHS] = {
			-a1*100,
    			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			0.0
  };
  */

  /*
  for (int i = 0; i < N*N; i++) {
    if (i % N == 0) {
      cout << "\n";
    }
    cout << A[i] << " ";
  }
  cout << "\n";
  */

  b_func[0] = -a1*100;
  info = LAPACKE_dgesv( LAPACK_COL_MAJOR, N, NRHS, A, LDA, ipiv, b_func, LDB);
  /*
  for(int i=0; i<= N-1; i++) {
    cout << b_func[i] << "\t";
  }
  cout << "\n";
  */

}

/*
void to_paraview() {
  
  ofstream myfile;
  myfile.open("vtk_files/practice/phi_practice.vtk");

  // Paraview header
  myfile << "# vtk DataFile Version 2.0\n";
  myfile << "FlowField\n";
  myfile << "ASCII\n";

  // Grid
  myfile << "DATASET STRUCTURED_GRID\n";
  myfile << "DIMENSIONS " << nx << " " << 1 << " " << 1 << "\n";
  myfile << "POINTS " << nx*1*1-1 << " float\n";
  for (int i = 0; i <= nx-2; i++) { 
    myfile << i << " " << "0" << " 0\n";
  }
  
  // Data
  myfile << "\n";
  myfile << "POINT_DATA";
  myfile << " " << nx*1-1 << "\n";

  myfile << "\n";
  myfile << "SCALARS PHI float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for (int i = 0; i <= nx-2; i++) { 
    myfile << b[i] << "\n";
  }
  
  myfile.close();
}
*/  


/* =================== MAIN ======================= */
int main() {
  
  int nx = 21;
  double k = 1.;
  double dt = 0.05;
  double dx = 1.;
  // if dx = 0.1 and dt = 0.01 --> simulation is blown up!! --> increase dx, decrease dt will help
  
  double temperature[nx];
  double temperature_new[nx];

  
  node_initialize(temperature, nx);
  node_1DHeatCondition(temperature);
  visualize(temperature, nx);
  
  /*
    // (1) Explicit scheme
  for (int i = 0; i <= 10; i++) { 
    explicit_sim(temperature, temperature_new, nx, k, dt, dx);
    visualize(temperature_new, nx);
  }
  */


  // (2) Implicit scheme
  /*
  double b[nx-1] = {
		    -(-k*dt/(dx*dx))*100,
		    0.0,
		    0.0,
		    0.0,
		    0.0,
		    0.0,
		    0.0,
		    0.0,
		    0.0,
		    0.0,
		    0.0,
		    0.0,
		    0.0,
		    0.0,
		    0.0,
		    0.0,
		    0.0,
		    0.0,
		    0.0
		    }; */
  /*
  for (int iteration = 0; iteration <= 500; iteration++) {
    implicit_sim(temperature, temperature_new, nx, k, dt, dx, b);
    visualize(b, nx-1);
  }
  */
  
  
}
