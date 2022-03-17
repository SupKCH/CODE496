#include <iostream>
#include <cmath>
#include <fstream> // save & load restart files
#include <string>
#include <iomanip> // std::setprecision()
using namespace std;
// Validation_6.vtk (recheck) gives results identical to Validation_4.vtk (old data)
//    ^ x_max = 32                                          ^ x_max = 32
// note: _5.vtk ---> x_max = 64 (too much), presents different convergence (that is usual for  sure!)

ofstream myfileO;  // output file stream
ifstream myfileI;  // input file stream

double du2_dx(double** u, int i, int j, double dx) {
  return (pow((u[i][j] + u[i+1][j])/2.0, 2) - pow((u[i-1][j] + u[i][j])/2.0, 2))/dx;
}

double duv_dy(double** u, double** v, int i, int j, double dy) {
  return ((v[i][j] + v[i+1][j])*(u[i][j] + u[i][j+1])  -  (v[i][j-1] + v[i+1][j-1])*(u[i][j-1] + u[i][j]))/(4.0*dy);
}

double d2u_dx2(double** u, int i, int j, double dx) {
  return (u[i+1][j] - 2*u[i][j] + u[i-1][j])/pow(dx, 2);
}

double d2u_dy2(double** u, int i, int j, double dy) {
  return (u[i][j+1] - 2*u[i][j] + u[i][j-1])/pow(dy, 2);
}

double duv_dx(double** u, double** v, int i, int j, double dx) {
  return ((u[i][j] + u[i][j+1])*(v[i][j] + v[i+1][j])  -  (u[i-1][j] + u[i-1][j+1])*(v[i-1][j] + v[i][j]))/(4.0*dx);
}

double dv2_dy(double** v, int i, int j, double dy) {
  return (pow( v[i][j] + v[i][j+1] ,2)  -  pow( v[i][j-1] + v[i][j] ,2))/(4*dy);
}

double d2v_dx2(double** v, int i, int j, double dx) {
  return (v[i+1][j] - 2*v[i][j] + v[i-1][j])/pow(dx,2);
}

double d2v_dy2(double** v, int i, int j, double dy) {
  return (v[i][j+1] - 2*v[i][j] + v[i][j-1])/pow(dy,2);
}

void initialize_u_3D(double*** var, int nx, int ny, int nz) {
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      for (int k=0; k < nz; k++) {
	var[i][j][k] = 0.0;
      }
    }
  }

  for (int i=1; i <= nx-2; i++) {
    for (int j=1; j <= ny-2; j++) {
      var[i][j][0] = 1.0;
    }
  }
  /* // No-slip condition
  for (int i=0; i <= nx-1; i++) {
    var[i][0] = -1.0*var[i][1];
    var[i][ny-1] = -1.0*var[i][ny-2];
  }
  */
}

void initialize_zero(double** var, int nx, int ny) {
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      var[i][j] = 0.0;
    }
  }
}

void initialize_phi(double** var, int nx, int ny) {
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      var[i][j] = 0.0;
    }
  }

  for (int j=0; j < int(ny/2); j++) {
    var[0][j] = 1.0;
  }
}

void initialize_pressure(double** var, int nx, int ny) {
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      if ((int)(i+j) % 2 == 0) {	
	var[i][j] = 1.0;
      }
      else var[i][j] = 0.5;
    }
  }
}

void quick_visualize_3D(double*** var, int nx, int ny, int nz) {
  for (int k=0; k <= 3; k++) {
    for (int j=ny-1; j >= 0; j-- ) {
      for (int i=0; i <= nx-1; i++) {
	cout << var[i][j][k] << " ";
      }
      cout << "\n";
    }
    cout << "\n";
  }
  cout << "\n";
}

/*
double g_value(int VALUE) {
  if (VALUE == 1) {return 0.0;}
  else if (VALUE == 2) {return 0.0;}
  else if (VALUE == 3) {return 1.0;}
  else return 0.0;
  }*/

void F_calculation(double** F_n, double** u, double** v, int nx, int ny, double dt, double Re, double dx, double dy) {
  for (int i=1; i <= nx-3; i++) {
    for (int j=1; j <= ny-2; j++) {
      F_n[i][j] = u[i][j] + dt*( (d2u_dx2(u,i,j,dx) + d2u_dy2(u,i,j,dy))/Re - du2_dx(u,i,j,dx) - duv_dy(u,v,i,j,dy) + 0.0 );
    }
  }
  for (int j=1; j <= ny-2; j++) {
    F_n[0][j] = u[0][j];
    F_n[nx-2][j] = u[nx-2][j];
  }
}

void G_calculation(double** G_n, double** u, double** v, int nx, int ny, double dt, double Re, double dx, double dy) {
  for (int i=1; i <= nx-2; i++) {
    for (int j=1; j <= ny-3; j++) {
      G_n[i][j] = v[i][j] + dt*( (d2v_dx2(v,i,j,dx) + d2v_dy2(v,i,j,dy))/Re - duv_dx(u,v,i,j,dx) - dv2_dy(v,i,j,dy) + 0.0 );
    }
  }
  for (int i=1; i <= nx-2; i++) {
    G_n[i][0] = v[i][0];
    G_n[i][ny-2] = v[i][ny-2];
  }
}

double eW(int i) {
  if (i == 1) return 0.0;
  else if (i > 1) return 1.0;
  return 1.0;
}

double eE(int i, int nx) {
  if (i < nx-2) return 1.0;
  else if (i == nx-2) return 0.0;
  return 0.0;
}

double eS(int j) {
  if (j == 1) return 0.0;
  else if (j > 1) return 1.0;
  return 1.0;
}

double eN(int j, int ny) {
  if (j < ny-2) return 1.0;
  else if (j == ny-2) return 0.0;
  return 0.0;
}

void pressure_solver(double** p, double SOR,  double** F_n, double** G_n, int nx, int ny, double dt, double dx, double dy) {
  for (int iteration=1; iteration <= 30; iteration++) {
    for (int i=1; i <= nx-2; i++) {
      for (int j=1; j <= ny-2; j++) {
	double rhs_ij = ((F_n[i][j] - F_n[i-1][j])/dx + (G_n[i][j] - G_n[i][j-1])/dy)/dt;
	p[i][j] = (1-SOR) * p[i][j] + SOR/((eE(i, nx) + eW(i))/(dx*dx) + (eN(j, ny) + eS(j))/(dy*dy)) * ((eE(i, nx)*p[i+1][j] + eW(i)*p[i-1][j])/(dx*dx) + (eN(j, ny)*p[i][j+1] + eS(j)*p[i][j-1])/(dy*dy) - rhs_ij);
      }
    }
    
    for (int j=1; j <= ny-2; j++) {
      p[0][j] = p[1][j];
      p[nx-1][j] = p[nx-2][j];
    }
    for (int i=1; i <= nx-2; i++) {
      p[i][0] = p[i][1];
      p[i][ny-1] = p[i][ny-2];
    }
  }
}

void u_calculation(double** u_new, double** u, double** F_n, double** p, double dt, double dx, int nx, int ny) {
  for (int i=1; i <= nx-3; i++) {
    for (int j=1; j <= ny-2; j++) {
      u_new[i][j] = F_n[i][j] - (dt/dx)*(p[i+1][j] - p[i][j]);
    }
  }

  // confirm inflow
  for (int j=1; j <= ny-2; j++) {
    u_new[0][j] = 1.0;
  }

  // outflow
  for(int j=1; j <= ny-2; j++) {
    u_new[nx-2][j] = u_new[nx-3][j];
    u_new[nx-1][j] = u_new[nx-2][j];
  }

  // No-slip -- update all u-nodes at outlet (nx-1) still works like (nx-2)
  for (int i=1; i <= nx-1; i++) {
    // TOP
    u_new[i][ny-1] = -1.0*u_new[i][ny-2];
    // BOTTOM
    u_new[i][0] = -1.0*u_new[i][1];
  }

  // update
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      u[i][j] = u_new[i][j];
    }
  }
}

void v_calculation(double** v_new, double** v, double** G_n, double** p, double dt, double dy, int nx, int ny) {
  for (int i=1; i <= nx-2; i++) {
    for (int j=1; j <= ny-3; j++) {
      v_new[i][j] = G_n[i][j] - (dt/dy)*(p[i][j+1] - p[i][j]);
    }
  }

  // inflow
  for (int j=1; j <= ny-3; j++) {
    v_new[0][j] = -1.0*v_new[1][j];
  }
  
  // outflow
  for(int j=1; j <= ny-3; j++) {
    v_new[nx-1][j] = v_new[nx-2][j];
  }

  // No-Slip
  for (int i=1; i <= nx-1; i++) {
    // TOP
    v_new[i][ny-2] = 0.0;
    // BOTTOM
    v_new[i][0] = 0.0;
  }
  
  // update
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      v[i][j] = v_new[i][j];
    }
  }
}


void paraview(string fileName, double **var, int nx, int ny, double dx, double dy, string parameter) {
  ofstream myfile;
  myfile.open(fileName);
  // -----------------------------------------------------------------//
  // Paraview header
  myfile << "# vtk DataFile Version 2.0\n";
  myfile << "FlowField\n";
  myfile << "ASCII\n";

  // Grid
  myfile << "DATASET STRUCTURED_GRID\n";
  myfile << "DIMENSIONS " << nx << " " << 1 << " " << ny << "\n";
  myfile << "POINTS " << nx*1*ny << " float\n";
  for (int j = 0; j <= ny-1; j++) {
    for (int i = 0; i <= nx-1; i++) {
      myfile << dx*i << " " << dy*j << " 0\n";
    }
  }
  
  // Data
  myfile << "\n";
  myfile << "POINT_DATA";
  myfile << " " << nx*ny << "\n";

  myfile << "\n";
  myfile << "SCALARS " << parameter << " float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for (int j = 0; j <= ny-1; j++) {
    for (int i = 0; i <= nx-1; i++) {
      myfile << setprecision(15) << var[i][j] << "\n";
    }
  }
  myfile.close();
}

void paraview2D(string fileName, double** u_new, double** v_new, double** p, double** phi_new, int nx, int ny, double dx, double dy, int precision) {
  ofstream myfile;
  myfile.open(fileName);

  myfile << "# vtk DataFile Version 2.0\n";
  myfile << "FlowField\n";
  myfile << "ASCII\n";

  myfile << "DATASET STRUCTURED_GRID\n";
  myfile << "DIMENSIONS " << nx << " " << 1 << " " << ny << "\n";
  myfile << "POINTS " << nx*1*ny << " float\n";
  for (int j = 0; j <= ny-1; j++) {
    for (int i = 0; i <= nx-1; i++) {
      myfile << dx*i << " " << dy*j << " 0\n";
    }
  }
  
  // Data: (u, v), p, phi (4 variables)
  myfile << "\n";
  myfile << "POINT_DATA";
  myfile << " " << nx*ny << "\n";

  // Pressure (p) [SCALAR]
  myfile << "\n";
  myfile << "SCALARS " << "p" << " float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for (int j = 0; j <= ny-1; j++) {
    for (int i = 0; i <= nx-1; i++) {
      myfile << setprecision(precision) << p[i][j] << "\n";
    }
  }

  // Pollution (phi) [SCALAR]
  myfile << "\n";
  myfile << "SCALARS " << "phi" << " float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for (int j = 0; j <= ny-1; j++) {
    for (int i = 0; i <= nx-1; i++) {
      myfile << setprecision(precision) << phi_new[i][j] << "\n";
    }
  }
  
  // Velocity vector (u, v) [VECTOR]
  myfile << "\n";
  myfile << "VECTORS " << "Velocity" << " float\n";
  for (int j = 0; j <= ny-1; j++) {
    for (int i = 0; i <= nx-1; i++) {
      myfile << setprecision(precision) << u_new[i][j] << " " << v_new[i][j] << " 0\n";
    }
  }
  myfile.close();
}

void paraview3D(string fileName, double** u_new, double** v_new, double** w_new, double** p, double** phi_new, int nx, int ny, int nz, double dx, double dy, double dz) {
  
  
  
}

void explicit_passiveScalar(double** phi, double** phi_new, double** u_new, double** v_new, int nx, int ny, double dx, double dy, double dt, double Re) {
  for (int i=1; i <= nx-2; i++) {
    for (int j=1; j <= ny-2; j++) {
      phi_new[i][j] = phi[i][j] + dt*(((phi[i+1][j] - 2*phi[i][j] + phi[i-1][j])/(dx*dx) + (phi[i][j+1] - 2*phi[i][j] + phi[i][j-1])/(dy*dy))/Re - u_new[i][j]*(phi[i+1][j] - phi[i-1][j])/(2*dx) - v_new[i][j]*(phi[i][j+1] - phi[i][j-1])/(2*dy));
    }
  }
  for (int j=0; j <= ny-1; j++) {
    phi_new[nx-1][j] = phi_new[nx-2][j];
  }
  for (int i=0; i <= nx-1; i++) {
    phi_new[i][0] = phi_new[i][1];
    phi_new[i][ny-1] = phi_new[i][ny-2];
  }

  // Update "phi" with "phi_new"
  for (int i=0; i <= nx-1; i++) {
    for (int j=0; j <= ny-1; j++) {
      phi[i][j] = phi_new[i][j];
    }
  }
}

void save_restartfile(double **var, string name_prefix, string variable_name, int iteration, int nx, int ny, int precision){
  ofstream myfileO;  // output file stream
  myfileO.open(name_prefix + variable_name + "_" + to_string(iteration) + ".dat");
  for (int i = 0; i <= nx-1; i++) {
     for (int j = 0; j <= ny-1; j++) {
       myfileO << setprecision(precision) << var[i][j] << " ";
    }
     myfileO << "\n";
  }
  myfileO.close();
}

void read_restartfile(double **var, string name_prefix, string variable_name, int iteration, int nx, int ny){
  ifstream myfileI;  // input file stream
  myfileI.open(name_prefix + to_string(nx) + "x" + to_string(ny) + "_" + variable_name + "_" + to_string(iteration) + ".dat");
  for (int i = 0; i <= nx-1; i++) {
     for (int j = 0; j <= ny-1; j++) {
       myfileI >> var[i][j];
    }
  }
  myfileI.close();
}

int main() {
  const int nx = 10;
  const int ny = 20;
  const int nz = 15;
  const double dx = 20.0/(double)ny;
  const double dy = 10.0/(double)ny;
  const double dz = 15.0/(double)nz;
  const double dt = 0.001; // Previously 0.005
  const double diameter = 1.0;
  const double Re = 1000.0;
  const double SOR = 1.7;

  const double centerX = 5.0;
  const double centerY = 5.0;
  
  // ---------------------------------------------
  double ***u = (double ***) malloc (nx * sizeof(double**));
  double ***v = (double ***) malloc (nx * sizeof(double**));
  double ***w = (double ***) malloc (nx * sizeof(double**));
  double ***phi = (double ***) malloc (nx * sizeof(double**));
  double ***p = (double ***) malloc (nx * sizeof(double**));
  double ***F_n = (double ***) malloc (nx * sizeof(double**));
  double ***G_n = (double ***) malloc (nx * sizeof(double**));
  double ***H_n = (double ***) malloc (nx * sizeof(double**));
  double ***u_new = (double ***) malloc (nx * sizeof(double**));
  double ***v_new = (double ***) malloc (nx * sizeof(double**));
  double ***w_new = (double ***) malloc (nx * sizeof(double**));
  double ***phi_new = (double ***) malloc (nx * sizeof(double**));
  double ***phi_half = (double ***) malloc (nx * sizeof(double**));
  
  for (int i = 0; i < nx; i++) {
    u[i] = (double **) malloc (ny * sizeof(double*));
    v[i] = (double **) malloc (ny * sizeof(double*));
    w[i] = (double **) malloc (ny * sizeof(double*));
    phi[i] = (double **) malloc (ny * sizeof(double*));
    p[i] = (double **) malloc (ny * sizeof(double*));
    F_n[i] = (double **) malloc (ny * sizeof(double*));
    G_n[i] = (double **) malloc (ny * sizeof(double*));
    H_n[i] = (double **) malloc (ny * sizeof(double*));
    u_new[i] = (double **) malloc (ny * sizeof(double*));
    v_new[i] = (double **) malloc (ny * sizeof(double*));
    w_new[i] = (double **) malloc (ny * sizeof(double*));
    phi_new[i] = (double **) malloc (ny * sizeof(double*));
    phi_half[i] = (double **) malloc (ny * sizeof(double*));
    
    for (int ii = 0; ii < ny; ii++) {
      u[i][ii] = (double *) malloc (nz * sizeof(double));
      v[i][ii] = (double *) malloc (nz * sizeof(double));
      w[i][ii] = (double *) malloc (nz * sizeof(double));
      phi[i][ii] = (double *) malloc (nz * sizeof(double));
      p[i][ii] = (double *) malloc (nz * sizeof(double));
      F_n[i][ii] = (double *) malloc (nz * sizeof(double));
      G_n[i][ii] = (double *) malloc (nz * sizeof(double));
      H_n[i][ii] = (double *) malloc (nz * sizeof(double));
      u_new[i][ii] = (double *) malloc (nz * sizeof(double));
      v_new[i][ii] = (double *) malloc (nz * sizeof(double));
      w_new[i][ii] = (double *) malloc (nz * sizeof(double));
      phi_new[i][ii] = (double *) malloc (nz * sizeof(double));
      phi_half[i][ii] = (double *) malloc (nz * sizeof(double));
    }
  }
  
  // --------------------------------------------
  initialize_u_3D(u, nx, ny, nz);
  quick_visualize_3D(u, nx, ny, nz);
  
}
