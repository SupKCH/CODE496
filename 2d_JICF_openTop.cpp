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

void initialize_v(double** var, int nx, int ny, double dx) {
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      var[i][j] = 0.0;
    }
  }

  for (int i=int((5-1)/dx); i <= int((5+1)/dx) ; i++) {
    var[i][0] = 2.0;
  } 
}

void initialize_zero(double** var, int nx, int ny) {
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      var[i][j] = 0.0;
    }
  }
}

void initialize_phi(double** var, int nx, int ny, double dx) {
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      var[i][j] = 0.0;
    }
  }

  for (int i=int((5-1)/dx); i <= int((5+1)/dx) ; i++) {
    var[i][0] = 1.0;
  } 
}

void initialize_pressure(double** var, int nx, int ny) {
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      if ((int)(i+j) % 2 == 0) {	
	var[i][j] = 1.0;
      }
    }
  }
}
// remove 0.5 of pressure section

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

void pressure_cal(double** p, double SOR,  double** F_n, double** G_n, int nx, int ny, double dt, double dx, double dy, int i, int j) {
  double rhs_ij = ((F_n[i][j] - F_n[i-1][j])/dx + (G_n[i][j] - G_n[i][j-1])/dy)/dt;
  p[i][j] = \
    (1-SOR) * p[i][j] +							\
    SOR/((eE(i, nx) + eW(i))/(dx*dx) + \
	 (eN(j, ny) + eS(j))/(dy*dy)) *					\
    ((eE(i, nx)*p[i+1][j] + eW(i)*p[i-1][j])/(dx*dx) + \
     (eN(j, ny)*p[i][j+1] + eS(j)*p[i][j-1])/(dy*dy) - rhs_ij);
}

void pressure_solver(double** p, double SOR,  double** F_n, double** G_n, int nx, int ny, double dt, double dx, double dy, int ext_it) {
  for (int iteration=1; iteration <= 30; iteration++) {
    if (ext_it % 2 == 0) {
      for (int i=1; i <= nx-2; i++) {
	for (int j=1; j <= ny-2; j++) {
	  pressure_cal(p, SOR, F_n, G_n, nx, ny, dt, dx, dy, i, j);
	}
      }
    }
    else {
      for (int i=nx-2; i >= 1; i--) {
	for (int j=1; j <= ny-2; j++) {
	  pressure_cal(p, SOR, F_n, G_n, nx, ny, dt, dx, dy, i, j);
	}
      }
    }
    
    for (int j=1; j <= ny-2; j++) {
      p[0][j] = p[1][j];
      p[nx-1][j] = p[nx-2][j];
    }
    for (int i=0; i <= nx-1; i++) {
      p[i][0] = p[i][1];
      p[i][ny-1] = p[i][ny-2];
    }
    // L2-Norm + Max-Norm
    if (iteration % 5 == 0) {
      double tmp = 0.0;
      double max_norm = 0.0;
      for (int i=1; i <= nx-2; i++) {
	for (int j=1; j <= ny-2; j++) {
	  double rhs_ij = ((F_n[i][j] - F_n[i-1][j])/dx + \
			   (G_n[i][j] - G_n[i][j-1])/dy)/dt;
	  double residual_ij = \
	    (eE(i,nx)*(p[i+1][j] - p[i][j]) - eW(i)*(p[i][j] - p[i-1][j]))/(dx*dx) + \
	    (eN(j,ny)*(p[i][j+1] - p[i][j]) - eS(j)*(p[i][j] - p[i][j-1]))/(dy*dy) - rhs_ij;
	  tmp += pow(residual_ij, 2);
	  if (residual_ij > max_norm) {
	    max_norm = residual_ij;
	  }	  
	}
      }
      cout << iteration << "\t" << sqrt(tmp/double((nx-2)*(ny-2))) << "\t" << max_norm << "\n";
    }
  }
  cout << "------------------\n";
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

void v_calculation(double** v_new, double** v, double** G_n, double** p, double dt, double dy, int nx, int ny, double dx) {
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

    if (i >= int((5-1)/dx) && i <= int((5+1)/dx)) {
      v_new[i][0] = 2.0;
    }
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

void explicit_passiveScalar(double** phi, double** phi_new, double** u_new, double** v_new, int nx, int ny, double dx, double dy, double dt, double Re) {
  for (int i=1; i <= nx-2; i++) {
    for (int j=1; j <= ny-2; j++) {
      phi_new[i][j] = phi[i][j] + dt*(\
				      ((phi[i+1][j] - 2*phi[i][j] + phi[i-1][j])/(dx*dx) + \
				       (phi[i][j+1] - 2*phi[i][j] + phi[i][j-1])/(dy*dy))/Re - \
				      (u_new[i][j] + u_new[i-1][j])*(phi[i+1][j] - phi[i-1][j])/(4*dx) - \
				      (v_new[i][j] + v_new[i][j-1])*(phi[i][j+1] - phi[i][j-1])/(4*dy)\
				      );
    }
  }
  for (int j=0; j <= ny-1; j++) {
    phi_new[nx-1][j] = phi_new[nx-2][j];
  }
  for (int i=0; i <= nx-1; i++) {
    phi_new[i][0] = phi_new[i][1];
    phi_new[i][ny-1] = phi_new[i][ny-2];
  }

  for (int i=int((5-1)/dx); i <= int((5+1)/dx) ; i++) {
    phi_new[i][0] = 1.0;
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
  const int nx = 400;//20*10*1; // increases from 400 to 800, to reach 63.84 in X-dimensionless distance, instead of 31.92
  const int ny = 200;//15*10*2; // increases by 4 times, previously =50
  //const double dy = 15.0/(double)ny; 
  const double dy = (15.0/(double)ny); // becomes smaller by 4 times
  //const double dx = 20.0/(double)nx; 
  const double dx = 1.0*dy; // increases by 4 times, previously =4.0*dy
  const double dt = 0.001; // Previously 0.005
  const double Re = 300.0;
  const double SOR = 1.7;
  
  // ---------------------------------------------
  double **u = (double **) malloc (nx * sizeof(double*));
  double **v = (double **) malloc (nx * sizeof(double*));
  double **phi = (double **) malloc (nx * sizeof(double*));
  double **p = (double **) malloc (nx * sizeof(double*));
  double **F_n = (double **) malloc (nx * sizeof(double*));
  double **G_n = (double **) malloc (nx * sizeof(double*));

  double **u_new = (double **) malloc (nx * sizeof(double*));
  double **v_new = (double **) malloc (nx * sizeof(double*));
  double **phi_new = (double **) malloc (nx * sizeof(double*));
  double **phi_half = (double **) malloc (nx * sizeof(double*));

  for (int i = 0; i < nx; i++) {
    u[i] = (double *) malloc (ny * sizeof(double));
    v[i] = (double *) malloc (ny * sizeof(double));
    phi[i] = (double *) malloc (ny * sizeof(double));
    p[i] = (double *) malloc (ny * sizeof(double));
    F_n[i] = (double *) malloc (ny * sizeof(double));
    G_n[i] = (double *) malloc (ny * sizeof(double));
    u_new[i] = (double *) malloc (ny * sizeof(double));
    v_new[i] = (double *) malloc (ny * sizeof(double));
    phi_new[i] = (double *) malloc (ny * sizeof(double));
    phi_half[i] = (double *) malloc (ny * sizeof(double));
  }

  // --------------------------------------------
  initialize_u(u, nx, ny);
  initialize_u(u_new, nx, ny);
  // initialize_zero(v, nx, ny);
  // initialize_zero(v_new, nx, ny);
  
  // initialize_phi(phi, nx, ny);
  // initialize_phi(phi_half, nx, ny);
  // initialize_phi(phi_new, nx, ny);
 
  initialize_v(v, nx, ny, dx);
  initialize_v(v_new, nx, ny, dx);
  
  initialize_phi(phi, nx, ny, dx);
  initialize_phi(phi_half, nx, ny, dx);
  initialize_phi(phi_new, nx, ny, dx);
  //quick_visualize(phi, nx, ny);
  
  initialize_zero(F_n, nx, ny);
  initialize_zero(G_n, nx, ny);
  initialize_pressure(p, nx, ny);
  
  string f_name = "vtk_files/2dJICF_openTop/2dJICF_openTop_";
  //string name_prefix = "checkpoints/Validation_2DChannelFlow/";
  const int precision = 10;
  const int save_precision = 6;
  paraview2D(f_name + to_string(1) + ".vtk", u_new, v_new, p, phi_new, nx, ny, dx, dy, precision);
  
  for (int it=1; it <= 100000; it++) { //100000
    F_calculation(F_n, u, v, nx, ny, dt, Re, dx, dy);
    G_calculation(G_n, u, v, nx, ny, dt, Re, dx, dy);
    pressure_solver(p, SOR, F_n, G_n, nx, ny, dt, dx, dy, it);  
    u_calculation(u_new, u, F_n, p, dt, dx, nx, ny);
    v_calculation(v_new, v, G_n, p, dt, dy, nx, ny, dx);
    explicit_passiveScalar(phi, phi_new, u_new, v_new, nx, ny, dx, dy, dt, Re);
    if (it % 5 == 0) { //20
      cout << it << "\n";
      if (it % 100 == 0) { //1000
	paraview2D(f_name + to_string(it) + ".vtk", u_new, v_new, p, phi_new, nx, ny, dx, dy, precision);
      }
      /*
      if (it % 200 == 0) {
	save_restartfile(u, name_prefix, "u", it, nx, ny, save_precision);
	save_restartfile(v, name_prefix, "v", it, nx, ny, save_precision);
	save_restartfile(F_n, name_prefix, "F_n", it, nx, ny, save_precision);
	save_restartfile(G_n, name_prefix, "G_n", it, nx, ny, save_precision);
	save_restartfile(u_new, name_prefix, "u_new", it, nx, ny, save_precision);
	save_restartfile(v_new, name_prefix, "v_new", it, nx, ny, save_precision);
	save_restartfile(p, name_prefix, "p", it, nx, ny, save_precision);
	save_restartfile(phi, name_prefix, "phi", it, nx, ny, save_precision);
	save_restartfile(phi_half, name_prefix, "phi_half", it, nx, ny, save_precision);
	save_restartfile(phi_new, name_prefix, "phi_new", it, nx, ny, save_precision);
      }
      */
    }
  }
  // RESTART FILE  
  // read_restartfile(u, name_prefix, "u", it, nx, ny);  
}
