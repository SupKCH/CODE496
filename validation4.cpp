#include <iostream>
#include <cmath>
#include <fstream> // save & load restart files
#include <string>
#include <iomanip> // std::setprecision()
#include <cstring>
#include <algorithm> // std::max
#include <lapacke.h>
// ==>   g++ code_lapackTEST2.cpp -o code_lapackTEST2.exe -L/usr/lib/x86_64-linux-gnu/ -llapacke

using namespace std;

ofstream myfileO;  // output file stream
ifstream myfileI;  // input file stream



double gamma(double** u, double** v, int i, int j, double dx, double dy, double dt) {
  return max(abs(u[i][j]*dt/dx),abs(v[i][j]*dt/dy));
}

double du2_dx(double** u, double** v, int i, int j, double dx, double dy, double dt) {
  return (pow((u[i][j] + u[i+1][j])/2.0, 2) - pow((u[i-1][j] + u[i][j])/2.0, 2))/dx + gamma(u,v,i,j,dx,dy,dt)*(abs(u[i][j] + u[i+1][j])*(u[i][j] - u[i+1][j]) - abs(u[i-1][j] + u[i][j])*(u[i-1][j] - u[i][j]))/(4*dx);
}

double duv_dy(double** u, double** v, int i, int j, double dx, double dy, double dt) {
  return ((v[i][j] + v[i+1][j])*(u[i][j] + u[i][j+1])  -  (v[i][j-1] + v[i+1][j-1])*(u[i][j-1] + u[i][j]))/(4.0*dy) + gamma(u,v,i,j,dx,dy,dt)*(abs(v[i][j] + v[i+1][j])*(u[i][j] - u[i][j+1]) - abs(v[i][j-1] + v[i+1][j-1])*(u[i][j-1] - u[i][j]))/(4.0*dy);
}

/*
double duv_dy(double** u, double** v, int i, int j, double dx, double dy, double dt) {
  return v[i][j]*((u[i][j] + u[i][j+1])  -  (u[i][j-1] + u[i][j]))/(2.0*dy);
}
*/

double d2u_dx2(double** u, int i, int j, double dx) {
  return (u[i+1][j] - 2*u[i][j] + u[i-1][j])/pow(dx, 2);
}

double d2u_dy2(double** u, int i, int j, double dy) {
  return (u[i][j+1] - 2*u[i][j] + u[i][j-1])/pow(dy, 2);
}

double duv_dx(double** u, double** v, int i, int j, double dx, double dy, double dt) {
  return ((u[i][j] + u[i][j+1])*(v[i][j] + v[i+1][j])  -  (u[i-1][j] + u[i-1][j+1])*(v[i-1][j] + v[i][j]))/(4.0*dx) + gamma(u,v,i,j,dx,dy,dt)*(abs(u[i][j] + u[i][j+1])*(v[i][j] - v[i+1][j]) - abs(u[i-1][j] + u[i-1][j+1])*(v[i-1][j] - v[i][j]))/(4.0*dx);
}

/*
double duv_dx(double** u, double** v, int i, int j, double dx, double dy, double dt) {
  return u[i][j]*((v[i][j] + v[i+1][j])  -  (v[i-1][j] + v[i][j]))/(2.0*dx);
}
*/

double dv2_dy(double** u, double** v, int i, int j, double dx, double dy, double dt) {
  return (pow( v[i][j] + v[i][j+1] ,2)  -  pow( v[i][j-1] + v[i][j] ,2))/(4*dy) + gamma(u,v,i,j,dx,dy,dt)*(abs(v[i][j] + v[i][j+1])*(v[i][j] - v[i][j+1]) - abs(v[i][j-1] + v[i][j])*(v[i][j-1] - v[i][j]))/(4.0*dy);
}

double d2v_dx2(double** v, int i, int j, double dx) {
  return (v[i+1][j] - 2*v[i][j] + v[i-1][j])/pow(dx,2);
}

double d2v_dy2(double** v, int i, int j, double dy) {
  return (v[i][j+1] - 2*v[i][j] + v[i][j-1])/pow(dy,2);
}

void initialize_u(double** var, int nx, int ny) {
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      var[i][j] = 1.0;
    }
  }

  for (int j=1; j < ny-1; j++) {
    var[0][j] = 1.0;
  }
  for (int i=0; i <= nx-1; i++) {
    var[i][0] = -1.0*var[i][1];
    var[i][ny-1] = -1.0*var[i][ny-2];
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
  for (int i=0; i <= nx-1; i++) {
    for (int j=0; j <= ny-1; j++) {
      var[i][j] = 1.0;
    }
  }
  for (int j=0; j <= ny-1; j++) {
    var[0][j] = 1.0;
  }
}

void quick_visualize(double** var, int nx, int ny) {
  for (int j=ny-1; j >= 0; j-- ) {
    for (int i=0; i <= nx-1; i++) {
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
      F_n[i][j] = u[i][j] + dt*( (d2u_dx2(u,i,j,dx) + d2u_dy2(u,i,j,dy))/Re - du2_dx(u,v,i,j,dx,dy,dt) - duv_dy(u,v,i,j,dx,dy,dt) + 0.0 );
    }
  }
  for (int j=1; j <= ny-2; j++) {
    //F_n[0][j] = u[0][j];
    F_n[0][j] = 1.0;
    //F_n[nx-2][j] = u[nx-2][j];
    F_n[nx-2][j] = F_n[nx-3][j];
    F_n[nx-1][j] = F_n[nx-2][j];
  }
}

void G_calculation(double** G_n, double** u, double** v, int nx, int ny, double dt, double Re, double dx, double dy) {
  for (int i=1; i <= nx-2; i++) {
    for (int j=1; j <= ny-3; j++) {
      G_n[i][j] = v[i][j] + dt*( (d2v_dx2(v,i,j,dx) + d2v_dy2(v,i,j,dy))/Re - duv_dx(u,v,i,j,dx,dy,dt) - dv2_dy(u,v,i,j,dx,dy,dt) + 0.0 );
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
  for (int iteration=1; iteration <= 5; iteration++) {
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
  //memcpy(u, u_new, sizeof(u));
  
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
    //v_new[0][j] = -1.0*v_new[1][j];
    v_new[0][j] = 0.0;
    v_new[1][j] = 0.0;
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
  //memcpy(v, v_new, sizeof(v));
  
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
      myfile << setprecision(6) << var[i][j] << "\n";
    }
  }
  myfile.close();
}

void implicit_passiveScalar(double** phi, double** phi_half, double** phi_new, double** u_new, double** v_new, int nx, int ny, double dx, double dy, double dt, double Re, bool verbose) {
  int N = nx-2, NRHS = 1, LDA = N, LDB = N;
  int ipiv[N], info;
  double A[LDA*N];

  double a = 1/(Re*(dx*dx));
  double b = 1/(Re*(dy*dy));
  double c = 2/dt;

  for (int i = 0; i <= N*N-1; i++) {
    A[i] = 0.0;
  }

  // ----------------- Start half-step X cycles (fixed Y)  -----------------
  for (int j = 1; j <= ny-2; j++) {
    for (int i = 0; i <= N*N-1; i++) {
      if (i == 0) {
	A[0] = c - 2*a;
	A[1] = u_new[i+1][j] + a;
      }
      else if (i % N == 0 && i != N*(N-1)) {
	int real_i = int(i/N)-1;
	A[i + real_i] = -1.0*u_new[real_i+1][j] + a;
	A[i+1 + real_i] = c - 2*a;
	A[i+2 + real_i] = u_new[real_i+1][j] + a;
      }
      else if (i == N*(N-1)) {
	int real_i = int(i/N)-1;
	A[N*N -2] = -1.0*u_new[real_i+1][j] + a;
	A[N*N -1] = c - 2*a + (u_new[real_i+1][j] + a);
      }
    }

    if (verbose) {
      for (int i = 0; i < N*N; i++) {
	if (i % N == 0) {
	  cout << "\n";
	}
	cout << A[i] << " ";
      }
      cout << "\n";
    }
  
    double bx[LDB*NRHS];
    // column matrix which "phi" is at single "j"
    for (int i = 0; i <= LDB*NRHS -1; i++) {
      bx[i] = c*phi[i+1][j] - v_new[i+1][j]*(phi[i+1][j+1] - phi[i+1][j-1])/(2*dy) + (phi[i+1][j+1] - 2*phi[i+1][j] + phi[i+1][j-1])/(Re*(dy*dy));
      if (i == 0) {
	bx[i] -= (-1.0*u_new[i+1][j] + a)*phi[i-1+1][j];
      }
      else if (i == LDB*NRHS-1) {
	bx[i] -= 0.0;
      }
    }
    if (verbose) {
      for(int i=0; i<= LDB*NRHS-1; i++) {
	cout << bx[i] << "\t";
      }
      cout << "\n";
    }

    info = LAPACKE_dgesv(LAPACK_COL_MAJOR, N, NRHS, A, LDA, ipiv, bx, LDB);
    // solution comes out as "bx"
    if (verbose) {
      for(int i=0; i<= LDB*NRHS-1; i++) {
	cout << bx[i] << "\t";
      }
      cout << "\n";
    }
    
    // Update to "phi_half"
    for (int i = 0; i <= LDB*NRHS-1; i++) {
      phi_half[i+1][j] = bx[i];
    }
   }//----------------- End of half-step X cycle -----------------

  // outflow (half-step)
  for (int j=1; j <= ny-2; j++) {
    phi_half[nx-1][j] = phi_half[nx-2][j];
  }
  // Top-Bottom no flux flow (half-step)
  for (int i = 0; i <= nx-1; i++) {
    phi_half[i][0] = phi_half[i][1];
    phi_half[i][ny-1] = phi_half[i][ny-2];
  }

  // checking X-half step
  //quick_visualize(phi_half, nx, ny);


  N = ny-2;
  NRHS = 1;
  LDA = N;
  LDB = N;
  int ipiv2[N], info2;
  double A2[LDA*N];

  for (int i = 0; i <= N*N-1; i++) {
    A2[i] = 0.0;
  }

  //----------------- Start half-step Y cycles -----------------
  for (int i = 1; i <= nx-2; i++) {
    for (int j = 0; j <= N*N-1; j++) {
      if (j == 0) {
	int real_j = int(j/N)-1;
	A2[0] = c - 2*b + (v_new[i][real_j+1] + b);
	A2[1] = v_new[i][j+1] + b;
      }
      else if (j % N == 0 && j != N*(N-1)) {
	int real_j = int(j/N)-1;
	A2[j + real_j] = -1.0*v_new[i][real_j+1] + b;
	A2[j+1 + real_j] = c - 2*b;
	A2[j+2 + real_j] = v_new[i][real_j+1] + b;
      }
      else if (j == N*(N-1)) {
	int real_j = int(j/N)-1;
	A2[N*N -2] = -1.0*v_new[i][real_j+1] + b;
	A2[N*N -1] = c - 2*b + (v_new[i][real_j+1] + b);
      }
    }

    if (verbose) {
      for (int j = 0; j < N*N; j++) {
	if (j % N == 0) {
	  cout << "\n";
	}
	cout << A2[j] << " ";
      }
      cout << "\n";
    }
  
    double by[LDB*NRHS];
    // column matrix which "phi" is at single "i"
    for (int j = 0; j <= LDB*NRHS -1; j++) {
      by[j] = c*phi_half[i][j+1] - u_new[i][j+1]*(phi_half[i+1][j+1] - phi_half[i-1][j+1])/(2*dx) + (phi_half[i+1][j+1] - 2*phi_half[i][j+1] + phi_half[i-1][j+1])/(Re*(dx*dx));
      if (j == 0) {
	by[j] -= 0.0;
      }
      else if (j == LDB*NRHS-1) {
	by[j] -= 0.0;
      }
    }

    if (verbose) {
      for(int j=0; j<= LDB*NRHS-1; j++) {
	cout << by[j] << "\t";
      }
      cout << "\n";
    }

    info2 = LAPACKE_dgesv(LAPACK_COL_MAJOR, N, NRHS, A2, LDA, ipiv2, by, LDB);
    // solution comes out as "bx"
    if (verbose) {
      for(int j=0; j<= LDB*NRHS-1; j++) {
	cout << by[j] << "\t";
      }
      cout << "\n";
    }
    
    // Update to "phi_new"
    for (int j = 0; j <= LDB*NRHS-1; j++) {
      phi_new[i][j+1] = by[j];
    }
  }//----------------- end of half-step Y cycles -----------------

  // outflow (full-step)
  for (int j=1; j <= ny-2; j++) {
    phi_new[nx-1][j] = phi_new[nx-2][j];
  }
  
  // Top-Bottom no flux flow (full-step)
  for (int i = 0; i <= nx-1; i++) {
    phi_new[i][0] = phi_new[i][1];
    phi_new[i][ny-1] = phi_new[i][ny-2];
  }
  
  // checking Y-half step
  //quick_visualize(phi_new, nx, ny);
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

int main() {
  const int nx = 800; //400// increases from 400 to 800, to reach 63.84 in X-dimensionless distance, instead of 31.92
                      // changed back to 400 (reduced load)
  const int ny = 200; //800// increases by 4 times, previously =50
  const double dy = (1.0/(double)ny); // becomes smaller by 4 times
  const double dx = 16.0*dy; // increases by 4 times, previously =4.0*dy // changed back to 4.0 from 16.0
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
  double **phi_half;
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

  phi_half = (double **) malloc (nx * sizeof(double));
  for (int i = 0; i < nx; i++) {
    phi_half[i] = (double *) malloc (ny * sizeof(double));
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
  
  initialize_phi(phi, nx, ny);
  initialize_phi(phi_half, nx, ny);
  initialize_phi(phi_new, nx, ny);

  //quick_visualize(phi, nx, ny);
  
  initialize_zero(F_n, nx, ny);
  initialize_zero(G_n, nx, ny);
  initialize_pressure(p, nx, ny);
  
  for (int it=1; it <= 100; it++) {
    F_calculation(F_n, u, v, nx, ny, dt, Re, dx, dy);
    G_calculation(G_n, u, v, nx, ny, dt, Re, dx, dy);
    pressure_solver(p, SOR, F_n, G_n, nx, ny, dt, dx, dy);  
    u_calculation(u_new, u, F_n, p, dt, dx, nx, ny);
    v_calculation(v_new, v, G_n, p, dt, dy, nx, ny);
    if (it % 5 == 0) {
      cout << it << "\n";
    }
  }
  /*
  for (int it2=1; it2 <= 1; it2++) {
    //implicit_passiveScalar(phi, phi_half, phi_new, u_new, v_new, nx, ny, dx,  dy, dt, Re, false);
    if (it2 % 1 == 0) {
      cout << it2 << ".\n";
    }
  }
  */
  
  for (int it = 1; it <= 10000; it++) {
    explicit_passiveScalar(phi, phi_new, u_new, v_new, nx, ny, dx, dy, dt, Re);
    if (it % 100 == 0) {
      cout << it << ".\n";
    }
  }
  //quick_visualize(phi_new, nx, ny);
  //quick_visualize(u_new, nx, ny);
  //quick_visualize(v_new, nx, ny);

  string filename_prefix = "vtk_files/valid_2_1/valid_2_1_";
  paraview(filename_prefix + "phi_" + to_string(1) + ".vtk", phi_new, nx, ny, dx, dy, "phi");
  paraview(filename_prefix + "u_" + to_string(1) + ".vtk", u_new, nx, ny, dx, dy, "u");
  paraview(filename_prefix + "v_" + to_string(1) + ".vtk", v_new, nx, ny, dx, dy, "v");
  paraview(filename_prefix + "p_" + to_string(1) + ".vtk", p, nx, ny, dx, dy, "p");
  
  /*
  string filename_prefix = "vtk_files/valid_2_1/valid_2_1_";
  paraview(filename_prefix + "u_" + to_string(1) + ".vtk", u_new, nx, ny, dx, dy, "u");
  paraview(filename_prefix + "v_" + to_string(1) + ".vtk", v_new, nx, ny, dx, dy, "v");
  paraview(filename_prefix + "p_" + to_string(1) + ".vtk", p, nx, ny, dx, dy, "p");
  paraview(filename_prefix + "Fn_" + to_string(1) + ".vtk", F_n, nx, ny, dx, dy, "F_n");
  paraview(filename_prefix + "Gn_" + to_string(1) + ".vtk", G_n, nx, ny, dx, dy, "G_n");
 

  for (int i=2; i<=5000; i++) {
    // ############ STEP 1: Calculate F_n & G_n from u_n, v_n ############
    F_calculation(F_n, u, v, nx, ny, dt, Re, dx, dy);
    G_calculation(G_n, u, v, nx, ny, dt, Re, dx, dy);
  
    // ############ STEP 2: Solve Poisson Eqn. for pressure @ n+1, p_n+1  ############ 
    pressure_solver(p, SOR, F_n, G_n, nx, ny, dt, dx, dy);
  

    // ############ STEP 3: compute new velocity field (u_n+1, v_n+1) with p_n+1  ############
    u_calculation(u_new, u, F_n, p, dt, dx, nx, ny);
    v_calculation(v_new, v, G_n, p, dt, dy, nx, ny);
    if (i % 1 == 0) {
      paraview(filename_prefix + "u_" + to_string(i) + ".vtk", u_new, nx, ny, dx, dy, "u");
      paraview(filename_prefix + "v_" + to_string(i) + ".vtk", v_new, nx, ny, dx, dy, "v");
      paraview(filename_prefix + "p_" + to_string(i) + ".vtk", p, nx, ny, dx, dy, "p");
      paraview(filename_prefix + "Fn_" + to_string(i) + ".vtk", F_n, nx, ny, dx, dy, "u");
      paraview(filename_prefix + "Gn_" + to_string(i) + ".vtk", G_n, nx, ny, dx, dy, "v");
      cout << i << "\n";
    }
    else if (i >= 4910 && i <= 4990 && i % 10 == 0) {
      paraview(filename_prefix + "u_" + to_string(i) + ".vtk", u_new, nx, ny, dx, dy, "u");
      paraview(filename_prefix + "v_" + to_string(i) + ".vtk", v_new, nx, ny, dx, dy, "v");
      paraview(filename_prefix + "p_" + to_string(i) + ".vtk", p, nx, ny, dx, dy, "p");
      cout << i << "\n";
    }
  }


  */
  // RESTART FILE
  /*
  string restart_filename = "res_files/valid_2_1.dat";
  save_restartfile(phi, restart_filename, nx, ny);
  read_restartfile(phi, restart_filename, nx, ny);
  */
}
