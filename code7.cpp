// Solving 1D Heat Equation

// d(phi)/dt = k(d^2 phi)/dx^2

// dt = time step size
// dx = space between 2 nodes
// -------------------------------------

// Central finite different
// d^2 phi / dx^2 ==> approx. (phi[i+1] - 2*phi[i] + phi[i-1]) / dx^2 ===> RHS

// Explicit Euler
// d(phi)/dt = RHS ==> { phi^(n+1) - phi^(n) } / dt = RHS^2

// bc[0], bc[10] = 0
// phi[0], phi[nx-1] = 0

//            - - - - - - - - - -

#include <iostream> // cout
#include <cmath> // math lib
#include <fstream> // file I/O
#include <stdlib.h> // malloc
#include <unistd.h> // sleep
using namespace std;

void initialize(double *var, double *var_new, int nx) {
  // initialize var
  for (int i = 0; i <= nx-1; i++) {
    var[i] = 0.;
    var_new[i] = 0.;
  }
  var[(nx-1)/2] = 1.;
}


void visualize(double *var, int nx) {
  //visualization
  for (int i=0; i <= nx-1; i++) {
    cout << var[i] << " ";
  }
  cout << "\n";
}


void simulation(double *var, double *var_new, int n, double c, double d, double dt) {
  double RHS;
  // simulation
  for (int i = 1; i <= n-2; i++) {
    RHS = c*((var[i+1] - 2.*var[i] + var[i-1])/(d*d));
    var_new[i] = var[i] + dt*RHS;
  }
}

void update(double *var, double *var_new, int n) {
  for (int i = 0; i <= n-1; i++) {
    var[i] = var_new[i];
  }
}

int main() {
  int nx = 11;
  double k = 1.;
  double dx = 1.;
  double dt = 0.01;

  // ---------------------------------------------
  double *phi;
  phi = (double *) malloc (nx * sizeof(double));
  double *phi_new;
  phi_new = (double *) malloc (nx * sizeof (double));
  // --------------------------------------------

  initialize(phi, phi_new, nx);
  visualize(phi, nx);

  // begin time loop
  for (int n = 1; n <= 1000; n++) {
    simulation(phi, phi_new, nx, k, dx, dt);
    update(phi, phi_new, nx);
    visualize(phi, nx);
    //sleep(1);
  }
  
  
}

