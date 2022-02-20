// Solving 2D Heat Equation

// d(phi)/dt = k(d^2 phi)/dx^2

// dt = time step size
// dx = space between 2 nodes
// dy = space between 2 nodes
// -------------------------------------

// Central finite different
// d^2 phi / dx^2 ==> approx. (phi[i+1][j] - 2*phi[i][j] + phi[i-1][j]) / dx^2 ===> RHS
// d^2 phi / dy^2 ==> approx. (phi[i][j+1] - 2*phi[i][j] + phi[i][j-1]) / dy^2 ===> RHS

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
#include <string> // string
using namespace std;

void initialize(double **var, double **var_new, int nx, int ny) {
  // initialize var
  for (int i = 0; i <= nx-1; i++) {
    for (int j = 0; j <= ny-1; j++) {
    var[i][j] = 0.;
    var_new[i][j] = 0.;
    }
  }
  var[(nx-1)/2][(ny-1)/2] = 1.;
  var[(nx-1)/2-1][(ny-1)/2] = 1.;
  var[(nx-1)/2+1][(ny-1)/2] = 1.;
  var[(nx-1)/2][(ny-1)/2-1] = 1.;
  var[(nx-1)/2][(ny-1)/2+1] = 1.;
  var[(nx-1)/2-1][(ny-1)/2-1] = 1.;
  var[(nx-1)/2-1][(ny-1)/2+1] = 1.;
  var[(nx-1)/2+1][(ny-1)/2-1] = 1.;
  var[(nx-1)/2+1][(ny-1)/2+1] = 1.;
}


void visualize(double **var, int nx, int ny) {
  //visualization
  for (int i=0; i <= nx-1; i++) {
    for (int j=0; j <= ny-1; j++) {   
      cout << var[i][j] << "\t";
    }
    cout << "\n";
  }
  cout << "\n";
}


void simulation(double **var, double **var_new, int nx, int ny, double c, double dx, double dy, double dt) {
  double RHS;
  // simulation
  for (int i = 1; i <= nx-2; i++) {
    for (int j = 1;j <= ny-2; j++) {
      RHS = c*((var[i+1][j] - 2.*var[i][j] + var[i-1][j])/(dx*dx)) +
            c*((var[i][j+1] - 2.*var[i][j] + var[i][j-1])/(dy*dy));
      var_new[i][j] = var[i][j] + dt*RHS;
    }
  }
}

void update(double **var, double **var_new, int nx, int ny) {
  for (int i = 0; i <= nx-1; i++) {
    for (int j = 0; j <= ny-1; j++) {
      var[i][j] = var_new[i][j];
    }
  }
}

void paraview(string fileName, double **var, int nx, int ny, double dx, double dy) {
  ofstream myfile;
  myfile.open(fileName);
  // -----------------------------------------------------------------//
  // Paraview header
  myfile << " vtk DataFile Version 2.0\n";
  myfile << "FlowField\n";
  myfile << "ASCII\n";

  // Grid
  myfile << "DATASET STRUCTURED_GRID\n";
  myfile << "DIMENSIONS " << nx << " " << 1 << " " << ny << "\n";
  myfile << "POINTS" << nx*1*ny << " float\n";
  for (int j = 0; j <= ny-1; j++) {
    for (int i = 0; i <= nx-1; i++) {
      myfile << dx*i << " " << dy*j << "0\n";
    }
  }
  
  // Data
  myfile << "\n";
  myfile << "POINT_DATA";
  myfile << " " << nx*ny << "\n";

  myfile << "\n";
  myfile << "SCALARS PHI float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for (int j = 0; j <= ny-1; j++) {
    for (int i = 0; i <= nx-1; i++) {
      myfile << var[i][j] << "\n";
    }
  }
  myfile.close();
}

int main() {
  int nx = 11;
  int ny = 11;
  double k = 1.;
  double dx = 1.;
  double dy = 1.;
  double dt = 0.001;
  string fileName;

  // ---------------------------------------------
  double **phi;
  phi = (double **) malloc (nx * sizeof(double));
  for (int i = 0; i < nx; i++) {
    phi[i] = (double *) malloc (ny * sizeof(double));
  }

  double **phi_new;
  phi_new = (double **) malloc (nx * sizeof(double));
  for (int i = 0; i < nx; i++) {
    phi_new[i] = (double *) malloc (ny * sizeof(double));
  }
  // --------------------------------------------

  initialize(phi, phi_new, nx, ny);
  //visualize(phi, nx, ny);
  //paraview(phi, nx, ny, dx,  dy);

  // begin time loop
  
  for (int n = 1; n <= 1000; n++) {
    simulation(phi, phi_new, nx, ny, k, dx, dy, dt);
    update(phi, phi_new, nx, ny);
    //visualize(phi, nx, ny);
    //sleep(0.25);

    if (n%100 == 0) {
      fileName = "var_" + to_string(n) + ".vtk";
      paraview(fileName, phi, nx, ny, dx,  dy);
    }
  }
      
  //fileName = "var.vtk";
  //paraview(fileName, phi, nx, ny, dx,  dy);
  
  
}

