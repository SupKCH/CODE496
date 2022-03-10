#include <iostream>
#include <cmath>
#include <fstream> // save & load restart files
#include <string>
using namespace std;

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

/*
double dvar_2Bydspace_1(double** var, double dspace, int i, int j) {
  double cal_value = (pow((var[i][j] + var[i+1][j])/2.0, 2) - pow((var[i-1][j] + var[i][j])/2.0, 2))/dspace;
  return cal_value;
}

double dvar1var2Bydspace_1(double** var1, double** var2, double dspace, int i, int j) {
  double cal_value = ((var2[i][j] + var2[i+1][j])*(var1[i][j] + var1[i][j+1])/4.0 + (var2[i][j-1] + var2[i+1][j-1])*(var1[i][j-1] + var1[i][j])/4.0)/dspace;
  return cal_value;
}

double dvar_2Bydspace_2(double** var, double dspace, int i, int j) {
  double cal_value = (var[i+1][j] - 2*var[i][j] + var[i-1][j])/(dspace*dspace);
  return cal_value;
}
*/

void initialize_u(double** var, int nx, int ny) {
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      var[i][j] = 0.0;
    }
  }

  for (int j=1; j < ny-1; j++) {
    var[0][j] = 1.0;
  }
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

void quick_visualize(double** var, int nx, int ny) {
  for (int j=ny-1; j >= 0; j-- ) {
    for (int i=0; i < nx; i++) {
      cout << var[i][j] << " ";
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
	//cout << rhs_ij << "\n";
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
      myfile << var[i][j] << "\n";
    }
  }
  myfile.close();
}


int main() {
  const int nx = 400;
  const int ny = 200; // increases by 4 times, previously =50
  const double dy = (1.0/(double)ny); // becomes smaller by 4 times
  const double dx = 16.0*dy; // increases by 4 times, previously =4.0*dy
  const double dt = 0.001; // Previously 0.005
  const double Re = 300.0;
  const double SOR = 1.7;
  
  // ---------------------------------------------
  double **u;
  double **v;
  double **phi;
  double **p;
  double **F_n;
  double **G_n;

  double **u_new;
  double **v_new;
  double **phi_new;
  double **p_new;
  
  u = (double **) malloc (nx * sizeof(double));
  for (int i = 0; i < nx; i++) {
    u[i] = (double *) malloc (ny * sizeof(double));
  }
  
  v = (double **) malloc (nx * sizeof(double));
  for (int i = 0; i < nx; i++) {
    v[i] = (double *) malloc (ny * sizeof(double));
  }
  
  phi = (double **) malloc (nx * sizeof(double));
  for (int i = 0; i < nx; i++) {
    phi[i] = (double *) malloc (ny * sizeof(double));
  }
  
  p = (double **) malloc (nx * sizeof(double));
  for (int i = 0; i < nx; i++) {
    p[i] = (double *) malloc (ny * sizeof(double));
  }


  u_new = (double **) malloc (nx * sizeof(double));
  for (int i = 0; i < nx; i++) {
    u_new[i] = (double *) malloc (ny * sizeof(double));
  }
  
  v_new = (double **) malloc (nx * sizeof(double));
  for (int i = 0; i < nx; i++) {
    v_new[i] = (double *) malloc (ny * sizeof(double));
  }
  
  phi_new = (double **) malloc (nx * sizeof(double));
  for (int i = 0; i < nx; i++) {
    phi_new[i] = (double *) malloc (ny * sizeof(double));
  }
  
  p_new = (double **) malloc (nx * sizeof(double));
  for (int i = 0; i < nx; i++) {
    p_new[i] = (double *) malloc (ny * sizeof(double));
  }
  
  F_n = (double **) malloc (nx * sizeof(double));
  for (int i = 0; i < nx; i++) {
    F_n[i] = (double *) malloc (ny * sizeof(double));
  }
  
  G_n = (double **) malloc (nx * sizeof(double));
  for (int i = 0; i < nx; i++) {
    G_n[i] = (double *) malloc (ny * sizeof(double));
  }

  // --------------------------------------------

  
  initialize_u(u, nx, ny);
  initialize_u(u_new, nx, ny);
  initialize_zero(v, nx, ny);
  initialize_zero(v_new, nx, ny);
  //initialize_phi(phi, nx, ny);
  initialize_zero(F_n, nx, ny);
  initialize_zero(G_n, nx, ny);
  initialize_pressure(p, nx, ny);

  string filename_prefix = "vtk_files/valid2/valid2_";
  paraview(filename_prefix + "u_" + to_string(1) + ".vtk", u_new, nx, ny, dx, dy, "u");
  paraview(filename_prefix + "v_" + to_string(1) + ".vtk", v_new, nx, ny, dx, dy, "v");
  paraview(filename_prefix + "p_" + to_string(1) + ".vtk", p, nx, ny, dx, dy, "p");
   

  for (int i=2; i<=25000; i++) {
    // ############ STEP 1: Calculate F_n & G_n from u_n, v_n ############
    F_calculation(F_n, u, v, nx, ny, dt, Re, dx, dy);
    G_calculation(G_n, u, v, nx, ny, dt, Re, dx, dy);
  
    // ############ STEP 2: Solve Poisson Eqn. for pressure @ n+1, p_n+1  ############ 
    pressure_solver(p, SOR, F_n, G_n, nx, ny, dt, dx, dy);
  

    // ############ STEP 3: compute new velocity field (u_n+1, v_n+1) with p_n+1  ############
    u_calculation(u_new, u, F_n, p, dt, dx, nx, ny);
    v_calculation(v_new, v, G_n, p, dt, dy, nx, ny);
    if (i % 500 == 0) {
      paraview(filename_prefix + "u_" + to_string(i) + ".vtk", u_new, nx, ny, dx, dy, "u");
      paraview(filename_prefix + "v_" + to_string(i) + ".vtk", v_new, nx, ny, dx, dy, "v");
      paraview(filename_prefix + "p_" + to_string(i) + ".vtk", p, nx, ny, dx, dy, "p");
      cout << i << "\n";
    }
    else if (i >= 24910 && i <= 24990) {
      paraview(filename_prefix + "u_" + to_string(i) + ".vtk", u_new, nx, ny, dx, dy, "u");
      paraview(filename_prefix + "v_" + to_string(i) + ".vtk", v_new, nx, ny, dx, dy, "v");
      paraview(filename_prefix + "p_" + to_string(i) + ".vtk", p, nx, ny, dx, dy, "p");
      cout << i << "\n";
    }
  }
}
