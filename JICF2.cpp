#include <iostream>
#include <cmath>
#include <fstream> // save & load restart files
#include <string>
#include <iomanip> // std::setprecision()
#include <algorithm> // std:min
using namespace std;

ofstream myfileO;  // output file stream
ifstream myfileI;  // input file stream

double du2_dx(double*** u, int i, int j, int k, double dx) {
  return (pow((u[i][j][k] + u[i+1][j][k])/2.0, 2) - pow((u[i-1][j][k] + u[i][j][k])/2.0, 2))/dx;
}

double duv_dy(double*** u, double*** v, int i, int j, int k, double dy) {
  return ((v[i][j][k] + v[i+1][j][k])*(u[i][j][k] + u[i][j+1][k])  -  (v[i][j-1][k] + v[i+1][j-1][k])*(u[i][j-1][k] + u[i][j][k]))/(4.0*dy);
}

double d2u_dx2(double*** u, int i, int j, int k, double dx) {
  return (u[i+1][j][k] - 2*u[i][j][k] + u[i-1][j][k])/pow(dx, 2);
}

double d2u_dy2(double*** u, int i, int j, int k, double dy) {
  return (u[i][j+1][k] - 2*u[i][j][k] + u[i][j-1][k])/pow(dy, 2);
}

double duv_dx(double*** u, double*** v, int i, int j, int k, double dx) {
  return ((u[i][j][k] + u[i][j+1][k])*(v[i][j][k] + v[i+1][j][k])  -  (u[i-1][j][k] + u[i-1][j+1][k])*(v[i-1][j][k] + v[i][j][k]))/(4.0*dx);
}

double dv2_dy(double*** v, int i, int j, int k, double dy) {
  return (pow( v[i][j][k] + v[i][j+1][k] ,2)  -  pow( v[i][j-1][k] + v[i][j][k] ,2))/(4*dy);
}

double d2v_dx2(double*** v, int i, int j, int k, double dx) {
  return (v[i+1][j][k] - 2*v[i][j][k] + v[i-1][j][k])/pow(dx,2);
}

double d2v_dy2(double*** v, int i, int j, int k, double dy) {
  return (v[i][j+1][k] - 2*v[i][j][k] + v[i][j-1][k])/pow(dy,2);
}


double d2u_dz2(double*** u, int i, int j, int k, double dz) {
  return (u[i][j][k+1] - 2*u[i][j][k] + u[i][j][k-1])/pow(dz,2);
}

double duw_dz(double*** u, double*** w, int i, int j, int k, double dz) {
  return ((w[i][j][k] + w[i+1][j][k])*(u[i][j][k] + u[i][j][k+1]) - (w[i][j][k-1] + w[i+1][j][k-1])*(u[i][j][k-1] + u[i][j][k]))/(4*dz);
}

double d2v_dz2(double*** v, int i, int j, int k, double dz) {
  return (v[i][j][k+1] - 2*v[i][j][k] + v[i][j][k-1])/pow(dz,2);
}

double dvw_dz(double*** v, double*** w, int i, int j, int k, double dz) {
  return ((w[i][j][k] + w[i][j+1][k])*(v[i][j][k] + v[i][j][k+1]) - (w[i][j][k-1] + w[i][j+1][k-1])*(v[i][j][k-1] + v[i][j][k]))/(4*dz);
}

double d2w_dx2(double*** w, int i, int j, int k, double dx) {
  return (w[i+1][j][k] - 2*w[i][j][k] + w[i-1][j][k])/pow(dx,2);
}

double d2w_dy2(double*** w, int i, int j, int k, double dy) {
  return (w[i][j+1][k] - 2*w[i][j][k] + w[i][j-1][k])/pow(dy,2);
}

double d2w_dz2(double*** w, int i, int j, int k, double dz) {
  return (w[i][j][k+1] - 2*w[i][j][k] + w[i][j][k-1])/pow(dz,2);
}

double duw_dx(double*** u, double*** w, int i, int j, int k, double dx) {
  return ((u[i][j][k] + u[i][j][k+1])*(w[i][j][k] + w[i+1][j][k]) - (u[i-1][j][k] + u[i-1][j][k+1])*(w[i-1][j][k] + w[i][j][k]))/(4*dx);
}

double dvw_dy(double*** v, double*** w, int i, int j, int k, double dy) {
  return ((v[i][j][k] + v[i][j][k+1])*(w[i][j][k] + w[i][j+1][k]) - (v[i][j-1][k] + v[i][j-1][k+1])*(w[i][j-1][k] + w[i][j][k]))/(4*dy);
}

double dw2_dz(double*** w, int i, int j, int k, double dz) {
  return (pow(w[i][j][k] + w[i][j][k+1],2) - pow(w[i][j][k-1] + w[i][j][k],2))/(4*dz);
}

void initialize_u_3D(double*** var, int nx, int ny, int nz, double dx, double dy, double radius, double centerX, double centerY) {
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      for (int k=0; k < nz; k++) {
	var[i][j][k] = 0.0;
      }
    }
  }
  
  for (int k=1; k <= nz-1; k++) {
    for (int j=0; j <= ny-1; j++) {
      var[0][j][k] = 0.5;
    }
  }

  for (int j=0; j <= ny-1; j++) {
    var[0][j][0] = -1.0*var[0][j][1];
  }
}

void initialize_v_3D(double*** var, int nx, int ny, int nz, double dx, double dy, double radius, double centerX, double centerY) {
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      for (int k=0; k < nz; k++) {
	var[i][j][k] = 0.0;
      }
    }
  }
}

void initialize_w_3D(double*** var, int nx, int ny, int nz, double dx, double dy, double radius, double centerX, double centerY) {
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      for (int k=0; k < nz; k++) {
	var[i][j][k] = 0.0;
      }
    }
  }
  
  for (int i=1; i <= nx-2; i++) {
    for (int j=1; j <= ny-2; j++) {
      if (sqrt(pow(i*dx - centerX ,2) + pow(j*dy - centerY ,2)) < radius + sqrt(pow(dx,2) + pow(dy,2))*0.001) {
	var[i][j][0] = 1.0; // 2.0
      }
    }
  }
}


void initialize_zero_3D(double*** var, int nx, int ny, int nz) {
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      for (int k=0; k < nz; k++) {
	var[i][j][k] = 0.0;
      }
    }
  }
}

void initialize_phi_3D(double*** var, int nx, int ny, int nz, double dx, double dy, double radius, double centerX, double centerY) {
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      for (int k=0; k < nz; k++) {
	var[i][j][k] = 0.0;
      }
    }
  }
  for (int i=1; i <= nx-2; i++) {
    for (int j=1; j <= ny-2; j++) {
      if (sqrt(pow(i*dx - centerX ,2) + pow(j*dy - centerY ,2)) < radius + sqrt(pow(dx,2) + pow(dy,2))*0.001) {
	var[i][j][0] = 1.0;
      }
    }
  }
}

void initialize_pressure_3D(double*** var, int nx, int ny, int nz) {
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      for (int k=0; k < nz; k++) {
	var[i][j][k] = 0.0;
      }
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

void F_calculation(double*** F_n, double*** u, double*** v, double*** w, int nx, int ny, int nz, double dt, double Re, double dx, double dy, double dz) {
  for (int i=1; i <= nx-3; i++) {
    for (int j=1; j <= ny-2; j++) {
      for (int k=1; k <= nz-2; k++) {
	F_n[i][j][k] = u[i][j][k] + dt*( \
					(d2u_dx2(u,i,j,k,dx) + d2u_dy2(u,i,j,k,dy) + d2u_dz2(u,i,j,k,dz))/Re - \
					du2_dx(u,i,j,k,dx) - duv_dy(u,v,i,j,k,dy) - duw_dz(u,w,i,j,k,dz) + 0.0 \
					 );
      }
    }
  }
  for (int k=1; k <= nz-2; k++) {
    for (int j=1; j <= ny-2; j++) {
      F_n[0][j][k] = u[0][j][k];
      F_n[nx-2][j][k] = u[nx-2][j][k];
    }
  }
}

void G_calculation(double*** G_n, double*** u, double*** v, double*** w, int nx, int ny, int nz, double dt, double Re, double dx, double dy, double dz) {
  for (int i=1; i <= nx-2; i++) {
    for (int j=1; j <= ny-3; j++) {
      for (int k=1; k <= nz-2; k++) {
	G_n[i][j][k] = v[i][j][k] + dt*( \
					(d2v_dx2(v,i,j,k,dx) + d2v_dy2(v,i,j,k,dy) + d2v_dz2(v,i,j,k,dz))/Re - \
					duv_dx(u,v,i,j,k,dx) - dv2_dy(v,i,j,k,dy) - dvw_dz(v,w,i,j,k,dz) + 0.0 \
					 );
      }
    }
  }
  for (int k=1; k <= nz-2; k++) {
    for (int i=1; i <= nx-2; i++) {
      G_n[i][0][k] = v[i][0][k];
      G_n[i][ny-2][k] = v[i][ny-2][k];
    }
  }
}

void H_calculation(double*** H_n, double*** u, double*** v, double*** w, int nx, int ny, int nz, double dt, double Re, double dx, double dy, double dz) {
  for (int i=1; i <= nx-2; i++) {
    for (int j=1; j <= ny-2; j++) {
      for (int k=1; k <= nz-3; k++) {
	H_n[i][j][k] = w[i][j][k] + dt*( \
					(d2w_dx2(w,i,j,k,dx) + d2w_dy2(w,i,j,k,dy) + d2w_dz2(w,i,j,k,dz))/Re - \
					duw_dx(u,w,i,j,k,dx) - dvw_dy(v,w,i,j,k,dy) - dw2_dz(w,i,j,k,dz) + 0.0 \
					 );
      }
    }
  }
  for (int j=1; j <= ny-2; j++) {
    for (int i=1; i <= nx-2; i++) {
      H_n[i][j][0] = w[i][j][0];
      H_n[i][j][nz-2] = w[i][j][nz-2];
    }
  }
}

double eW(int i) {
  if (i == 1) return 1.0;
  else if (i > 1) return 1.0;
  return 1.0;
}
// 0 1 1
double eE(int i, int nx) {
  if (i < nx-2) return 1.0;
  else if (i == nx-2) return 1.0;
  return 1.0;
}
// 1 0 0
double eS(int j) {
  if (j == 1) return 1.0;
  else if (j > 1) return 1.0;
  return 1.0;
}
// 0 1 1
double eN(int j, int ny) {
  if (j < ny-2) return 1.0;
  else if (j == ny-2) return 1.0;
  return 1.0;
}
// 1 0 0
double eB(int k) {
  if (k == 1) return 1.0;
  else if (k > 1) return 1.0;
  return 1.0;
}
// 0 1 1
double eT(int k, int nz) {
  if (k < nz-2) return 1.0;
  else if (k == nz-2) return 1.0;
  return 1.0;
}
// 1 0 0
// return all to 1.0; since we replace pressure boundary

void pressure_cal(double*** p, double SOR,  double*** F_n, double*** G_n, double*** H_n, int nx, int ny, int nz, double dt, double dx, double dy, double dz, int i, int j, int k) {
  double rhs_ijk = ((F_n[i][j][k] - F_n[i-1][j][k])/dx + \
		    (G_n[i][j][k] - G_n[i][j-1][k])/dy + \
		    (H_n[i][j][k] - H_n[i][j][k-1])/dz)/dt;
	  
  p[i][j][k] = (1-SOR) * p[i][j][k] + \
    (SOR/((eE(i, nx) + eW(i))/(dx*dx) + (eN(j, ny) + eS(j))/(dy*dy) + (eT(k, nz) + eB(k))/(dz*dz))) * \
    ((eE(i, nx)*p[i+1][j][k] + eW(i)*p[i-1][j][k])/(dx*dx) + \
     (eN(j, ny)*p[i][j+1][k] + eS(j)*p[i][j-1][k])/(dy*dy) + \
     (eT(k, nz)*p[i][j][k+1] + eB(k)*p[i][j][k-1])/(dz*dz) - \
     rhs_ijk);
}

void pressure_solver(double*** p, double SOR,  double*** F_n, double*** G_n, double*** H_n, int nx, int ny, int nz, double dt, double dx, double dy, double dz, int ext_it) {  
  int iteration_limit;
  if (ext_it == 1) iteration_limit = 100;
  else iteration_limit = 30;
  
  for (int iteration=1; iteration <= iteration_limit; iteration++) {
    for (int k=0; k <= nz-1; k++) {
      for (int j=0; j <= ny-1; j++) {
	p[0][j][k] = p[1][j][k];
	p[nx-1][j][k] = p[nx-2][j][k];
      }
    }
    for (int k=0; k <= nz-1; k++) {
      for (int i=0; i <= nx-1; i++) {
	p[i][0][k] = p[i][1][k];
	p[i][ny-1][k] = p[i][ny-2][k];
      }
    }
    for (int j=0; j <= ny-1; j++) {
      for (int i=0; i <= nx-1; i++) {
	p[i][j][0] = p[i][j][1];
	p[i][j][nz-1] = p[i][j][nz-2];
      }
    }

    if (ext_it % 8 == 1) {
      for (int i=1; i <= nx-2; i++) {
	for (int j=1; j <= ny-2; j++) {
	  for (int k=1; k <= nz-2; k++) {
	    pressure_cal(p, SOR, F_n, G_n, H_n, nx, ny, nz, dt, dx, dy, dz, i, j, k);
	  }
	}
      }
    }
    else if (ext_it % 8 == 2) {
      for (int i=nx-2; i >= 1; i--) {
	for (int j=ny-2; j >= 1; j--) {
	  for (int k=nz-2; k >= 1; k--) {
	    pressure_cal(p, SOR, F_n, G_n, H_n, nx, ny, nz, dt, dx, dy, dz, i, j, k);
	  }
	}
      }
    }    
    else if (ext_it % 8 == 3) {
      for (int i=1; i <= nx-2; i++) {
	for (int j=ny-2; j >= 1; j--) {
	  for (int k=nz-2; k >= 1; k--) {
	    pressure_cal(p, SOR, F_n, G_n, H_n, nx, ny, nz, dt, dx, dy, dz, i, j, k);
	  }
	}
      }
    }
    else if (ext_it % 8 == 4) {
      for (int i=nx-2; i >= 1; i--) {
	for (int j=1; j <= ny-2; j++) {
	  for (int k=1; k <= nz-2; k++) {
	    pressure_cal(p, SOR, F_n, G_n, H_n, nx, ny, nz, dt, dx, dy, dz, i, j, k);
	  }
	}
      }
    }
    else if (ext_it % 8 == 5) {
      for (int i=1; i <= nx-2; i++) {
	for (int j=ny-2; j >= 1; j--) {
	  for (int k=1; k <= nz-2; k++) {
	    pressure_cal(p, SOR, F_n, G_n, H_n, nx, ny, nz, dt, dx, dy, dz, i, j, k);
	  }
	}
      }
    }
    else if (ext_it % 8 == 6) {
      for (int i=nx-2; i >= 1; i--) {
	for (int j=1; j <= ny-2; j++) {
	  for (int k=nz-2; k >= 1; k--) {
	    pressure_cal(p, SOR, F_n, G_n, H_n, nx, ny, nz, dt, dx, dy, dz, i, j, k);
	  }
	}
      }
    }
    else if (ext_it % 8 == 7) {
      for (int i=1; i <= nx-2; i++) {
	for (int j=1; j <= ny-2; j++) {
	  for (int k=nz-2; k >= 1; k--) {
	    pressure_cal(p, SOR, F_n, G_n, H_n, nx, ny, nz, dt, dx, dy, dz, i, j, k);
	  }
	}
      }
    }
    else if (ext_it % 8 == 0) {
      for (int i=nx-2; i >= 1; i--) {
	for (int j=ny-2; j >= 1; j--) {
	  for (int k=1; k <= nz-2; k++) {
	    pressure_cal(p, SOR, F_n, G_n, H_n, nx, ny, nz, dt, dx, dy, dz, i, j, k);
	  }
	}
      }
    }
        
    // L2-Norm + Max-Norm
    if (iteration % 5 == 0) {
      double tmp = 0.0;
      double max_norm = 0.0;
      for (int i=1; i <= nx-2; i++) {
	for (int j=1; j <= ny-2; j++) {
	  for (int k=1; k <= nz-2; k++) {
	    double rhs_ijk = ((F_n[i][j][k] - F_n[i-1][j][k])/dx + \
			      (G_n[i][j][k] - G_n[i][j-1][k])/dy + \
			      (H_n[i][j][k] - H_n[i][j][k-1])/dz)/dt;
	    double residual_ijk = \
	      (eE(i,nx)*(p[i+1][j][k] - p[i][j][k]) - eW(i)*(p[i][j][k] - p[i-1][j][k]))/(dx*dx) + \
	      (eN(j,ny)*(p[i][j+1][k] - p[i][j][k]) - eS(j)*(p[i][j][k] - p[i][j-1][k]))/(dy*dy) + \
	      (eT(k,nz)*(p[i][j][k+1] - p[i][j][k]) - eB(k)*(p[i][j][k] - p[i][j][k-1]))/(dz*dz) - rhs_ijk;
	    tmp += pow(residual_ijk, 2);
	    if (residual_ijk > max_norm) {
	      max_norm = residual_ijk;
	    }
	  }
	}
      }
      double norm = sqrt(tmp/double((nx-2)*(ny-2)*(nz-2)));
      cout << iteration << "\t" << norm << "\t" << max_norm << "\n";
      if (iteration >= iteration_limit && norm <= 1.0) break;
    }
  }
  cout << "------------------\n";
}

void u_calculation(double*** u_new, double*** u, double*** F_n, double*** p, double dt, double dx, int nx, int ny, int nz) {
  for (int i=1; i <= nx-3; i++) {
    for (int j=1; j <= ny-2; j++) {
      for (int k=1; k <= nz-2; k++) {
	u_new[i][j][k] = F_n[i][j][k] - (dt/dx)*(p[i+1][j][k] - p[i][j][k]);
      }
    }
  }

  // confirm inflow x=West (YZ)
  for (int k=1; k <= nz-2; k++) {
    for (int j=1; j <= ny-2; j++) {
      u_new[0][j][k] = 0.5;
    }
  }

  // outflow x=East (YZ)
  for (int k=1; k <= nz-2; k++) {
    for(int j=1; j <= ny-2; j++) {
      u_new[nx-2][j][k] = u_new[nx-3][j][k];
      u_new[nx-1][j][k] = u_new[nx-2][j][k];
    }
  }

  // No-slip -- update all u-nodes at outlet (nx-1) still works like (nx-2) ********* <----- previously nx-1
  for (int j=1; j <= ny-2; j++) {
    for (int i=1; i <= nx-2; i++) {
      // z=Bottom (XY)
      u_new[i][j][0] = -1.0*u_new[i][j][1];
    }
  }

  // Free-slip
  // y=South,North (XZ)
  for (int i=1; i <= nx-2; i++) {
    for (int k=1; k <= nz-2; k++) {
      u_new[i][0][k] = u_new[i][1][k];
      u_new[i][ny-1][k] = u_new[i][ny-2][k];
    }
  }
  // z=Top (XY)
  for (int j=1; j <= ny-2; j++) {
    for (int i=1; i <= nx-2; i++) {
      u_new[i][j][nz-1] = u_new[i][j][nz-2];
      //u_new[i][j][nz-1] = -1.0*u_new[i][j][nz-2]; // no-slip
    }
  }

  // update
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      for (int k=0; k < nz; k++) {
	u[i][j][k] = u_new[i][j][k];
      }
    }
  }
}

void v_calculation(double*** v_new, double*** v, double*** G_n, double*** p, double dt, double dy, int nx, int ny, int nz) {
  for (int i=1; i <= nx-2; i++) {
    for (int j=1; j <= ny-3; j++) {
      for (int k=1; k <= nz-2; k++) {
	v_new[i][j][k] = G_n[i][j][k] - (dt/dy)*(p[i][j+1][k] - p[i][j][k]);
      }
    }
  }

  // inflow x=West/Free-slip x=East (YZ)
  for (int k=1; k <= nz-2; k++) {
    for (int j=1; j <= ny-2; j++) {
      v_new[0][j][k] = -1.0*v_new[1][j][k]; //  v_avg = 0.0 at inlet
      v_new[nx-1][j][k] = v_new[nx-2][j][k]; // free-slip
    }
  }
  
  // outflow y=South,North (XZ)
  for (int i=1; i <= nx-2; i++) {
    for (int k=1; k <= nz-2; k++) {
      v_new[i][0][k] = v_new[i][1][k];
      v_new[i][ny-2][k] = v_new[i][ny-3][k];
      v_new[i][ny-1][k] = v_new[i][ny-2][k];
    }
  }
  
  // Free-slip k=Top/ No-slip k=Bottom (XY)
  for (int i=1; i <= nx-1; i++) {
    for (int j=0; j <= ny-1; j++) {
      v_new[i][j][nz-1] = v_new[i][j][nz-2]; // open top
      //v_new[i][j][nz-1] = -1.0*v_new[i][j][nz-2]; // no-slip
      v_new[i][j][0] = -1.0*v_new[i][j][1]; // ground at bottom
    }
  }

  // update
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      for (int k=0; k < nz; k++) {
	v[i][j][k] = v_new[i][j][k];
      }
    }
  }
}

void w_calculation(double*** w_new, double*** w, double*** H_n, double*** p, double dt, double dx, double dy, double dz, int nx, int ny, int nz, double centerX, double centerY, double radius) {
  for (int i=1; i <= nx-2; i++) {
    for (int j=1; j <= ny-2; j++) {
      for (int k=1; k <= nz-3; k++) {
	w_new[i][j][k] = H_n[i][j][k] - (dt/dz)*(p[i][j][k+1] - p[i][j][k]);
      }
    }
  }
  
  // No-flow (except hole) k=Bottom
  for (int i=1; i <= nx-2; i++) {
    for (int j=1; j <= ny-2; j++) {
      if (sqrt(pow(i*dx - centerX ,2) + pow(j*dy - centerY ,2)) >= radius + sqrt(pow(dx,2) + pow(dy,2))*0.001) { // not hole
	w_new[i][j][0] = 0.0;
      }
      else { // confirm jet flow
	w_new[i][j][0] = 1.0; //2.0
      }
    }
  }
  
  // outflow k=Top (XY)
  for (int i=1; i <= nx-2; i++) {
    for (int j=1; j <= ny-2; j++) {
      w_new[i][j][nz-2] = w_new[i][j][nz-3];
      w_new[i][j][nz-1] = w_new[i][j][nz-2];
      //w_new[i][j][nz-2] = 0.0;
      //w_new[i][j][nz-1] = 0.0;
    }
  }
  
  // Free-slip x=West,East (YZ) 
  for (int k=1; k <= nz-1; k++) {
    for (int j=0; j <= ny-1; j++) {
      w_new[0][j][k] = w_new[1][j][k];
      w_new[nx-1][j][k] = w_new[nx-2][j][k];
    }
  }

  // Free-slip y=South,North (XZ)
  for (int i=0; i <= nx-1; i++) {
    for (int k=1; k <= nz-1; k++) {
      w_new[i][0][k] = w_new[i][1][k];
      w_new[i][ny-1][k] = w_new[i][ny-2][k];
    }
  }
  // update
  for (int i=0; i < nx; i++) {
    for (int j=0; j < ny; j++) {
      for (int k=0; k < nz; k++) {
	w[i][j][k] = w_new[i][j][k];
      }
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

void paraview3D(string fileName, double*** u_new, double*** v_new, double*** w_new, double*** p, double*** phi_new, int nx, int ny, int nz, double dx, double dy, double dz, int precision, double t) {
  ofstream myfile;
  myfile.open(fileName);

  myfile << "# vtk DataFile Version 2.0\n";
  myfile << "FlowField\n";
  myfile << "ASCII\n";

  myfile << "DATASET STRUCTURED_GRID\n";
  myfile << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
  myfile << "POINTS " << nx*ny*nz << " float\n";
  for (int k = 0; k <= nz-1; k++) {
    for (int j = 0; j <= ny-1; j++) {
      for (int i = 0; i <= nx-1; i++) {
	myfile << dx*i << " " << dy*j << " " << dz*k << "\n";
      }
    }
  }
  
  // Data: (u, v, w), p, phi (5 variables) + k
  myfile << "\n";
  myfile << "POINT_DATA";
  myfile << " " << nx*ny*nz << "\n";

  // Pressure (p) [SCALAR]
  myfile << "\n";
  myfile << "SCALARS " << "p" << " float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for (int k = 0; k <= nz-1; k++) {
    for (int j = 0; j <= ny-1; j++) {
      for (int i = 0; i <= nx-1; i++) {
	myfile << setprecision(precision) << p[i][j][k] << "\n";
      }
    }
  }

  // Pollution (phi) [SCALAR]
  myfile << "\n";
  myfile << "SCALARS " << "phi" << " float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for (int k = 0; k <= nz-1; k++) {
    for (int j = 0; j <= ny-1; j++) {
      for (int i = 0; i <= nx-1; i++) {
	myfile << setprecision(precision) << phi_new[i][j][k] << "\n";
      }
    }
  }

  // Kinetic Energy (k) [SCALAR]
  myfile << "\n";
  myfile << "SCALARS " << "k" << " float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  for (int k = 0; k <= nz-1; k++) {
    for (int j = 0; j <= ny-1; j++) {
      for (int i = 0; i <= nx-1; i++) {
	myfile << setprecision(precision) << pow(u_new[i][j][k], 2) + pow(v_new[i][j][k], 2) + pow(w_new[i][j][k], 2) << "\n";
      }
    }
  }
  
  // Velocity vector (u, v, w) [VECTOR]
  myfile << "\n";
  myfile << "VECTORS " << "Velocity" << " float\n";
  for (int k = 0; k <= nz-1; k++) {
    for (int j = 0; j <= ny-1; j++) {
      for (int i = 0; i <= nx-1; i++) {
	myfile << setprecision(precision) << u_new[i][j][k] << " " << v_new[i][j][k] << " " << w_new[i][j][k] << "\n";
      }
    }
  }
  myfile.close();
}

void explicit_passiveScalar(double*** phi, double*** phi_new, double*** u_new, double*** v_new, double*** w_new, int nx, int ny, int nz, double dx, double dy, double dz, double dt, double Re, double centerX, double centerY, double radius) {
  for (int i=1; i <= nx-2; i++) {
    for (int j=1; j <= ny-2; j++) {
      for (int k=1; k <= nz-2; k++) {
	phi_new[i][j][k] = phi[i][j][k] +\
	  dt*(((phi[i+1][j][k] - 2*phi[i][j][k] + phi[i-1][j][k])/(dx*dx) + \
	       (phi[i][j+1][k] - 2*phi[i][j][k] + phi[i][j-1][k])/(dy*dy) + \
	       (phi[i][j][k+1] - 2*phi[i][j][k] + phi[i][j][k-1])/(dz*dz))/Re - \
	      (u_new[i][j][k] + u_new[i-1][j][k])*(phi[i+1][j][k] - phi[i-1][j][k])/(4*dx) - \
	      (v_new[i][j][k] + v_new[i][j-1][k])*(phi[i][j+1][k] - phi[i][j-1][k])/(4*dy) - \
	      (w_new[i][j][k] + w_new[i][j][k-1])*(phi[i][j][k+1] - phi[i][j][k-1])/(4*dz));
      }
    }
  }
  // no gradient on all planes except the hole area
  // k=Top/Bottom (XY)
  for (int i=0; i <= nx-1; i++) {
    for (int j=0; j <= ny-1; j++) {
      phi_new[i][j][nz-1] = phi_new[i][j][nz-2];
      //phi_new[i][j][nz-1] = 0.0;
      if (sqrt(pow(i*dx - centerX ,2) + pow(j*dy - centerY ,2)) >= radius + sqrt(pow(dx,2) + pow(dy,2))*0.001) { // not hole
	phi_new[i][j][0] = phi_new[i][j][1];
	//phi_new[i][j][0] = 0.0;
      }
      else { // at hole
      	phi_new[i][j][0] = phi[i][j][0];
      }
    }
  }
  // x=West,East (YZ)
  for (int j=0; j <= ny-1; j++) {
    for (int k=0; k <= nz-1; k++) {
      phi_new[0][j][k] = phi_new[1][j][k];
      phi_new[nx-1][j][k] = phi_new[nx-2][j][k];
      //phi_new[0][j][k] = 0.0;
      //phi_new[nx-1][j][k] = 0.0;
    }
  }
  
  // y=South,North (XZ)
  for (int i=0; i <= nx-1; i++) {
    for (int k=0; k <= nz-1; k++) {
      phi_new[i][0][k] = phi_new[i][1][k];
      phi_new[i][ny-1][k] = phi_new[i][ny-2][k];
      //phi_new[i][0][k] = 0.0;
      //phi_new[i][ny-1][k] = 0.0;
    }
  }

  // Update "phi" with "phi_new"
  for (int i=0; i <= nx-1; i++) {
    for (int j=0; j <= ny-1; j++) {
      for (int k=0; k <= nz-1; k++) {
	if (phi_new[i][j][k] >= 1.0) phi_new[i][j][k] = 1.0;
	else if (phi_new[i][j][k] <= 0.0) phi_new[i][j][k] = 0.0;
	phi[i][j][k] = phi_new[i][j][k];
      }
    }
  }
}

void save_restartfile(double ***var, string name_prefix, string variable_name, int iteration, int nx, int ny, int nz, int save_precision){
  ofstream myfileO;
  myfileO.open(name_prefix + to_string(nx) + "x" + to_string(ny) + "x" + to_string(nz) + "_" + variable_name + "_" + to_string(iteration) + ".dat");
  for (int i = 0; i <= nx-1; i++) {
    for (int j = 0; j <= ny-1; j++) {
      for (int k = 0; k <= nz-1; k++) {
	myfileO << setprecision(save_precision) << var[i][j][k] << " ";
      }
    }
    myfileO << "\n";
  }
  myfileO.close();
}

void read_restartfile(double ***var, string name_prefix, string variable_name, int readAt_iteration, int nx, int ny, int nz){
  ifstream myfileI;
  myfileI.open(name_prefix + to_string(nx) + "x" + to_string(ny) + "x" + to_string(nz) + "_" + variable_name + "_" + to_string(readAt_iteration) + ".dat");
  for (int i = 0; i <= nx-1; i++) {
    for (int j = 0; j <= ny-1; j++) {
      for (int k = 0; k <= nz-1; k++) {
	myfileI >> var[i][j][k];
      }
    }
  }
  myfileI.close();
}

void save_time(double t, string name_prefix, int iteration, int nx, int ny, int nz) {
  ofstream myfileO;
  myfileO.open(name_prefix + to_string(nx) + "x" + to_string(ny) + "x" + to_string(nz) + "_t_" + to_string(iteration) + ".dat");
  myfileO << t << endl;
  myfileO.close();
}

double load_time(string name_prefix, int readAt_iteration, int nx, int ny, int nz) {
  ifstream myfileI;
  double t;
  myfileI.open(name_prefix + to_string(nx) + "x" + to_string(ny) + "x" + to_string(nz) + "_t_" + to_string(readAt_iteration) + ".dat");
  myfileI >> t;
  myfileI.close();
  return t;
}

double variable_dt(double tau, int nx, int ny, int nz, double Re, double dx, double dy, double dz, double*** u, double*** v, double*** w) {
  double u_max = 0.0;
  double v_max = 0.0;
  double w_max = 0.0;
  for (int i=0; i <= nx-1; i++) {
    for (int j=0; j <= ny-1; j++) {
      for (int k=0; k <= nz-1; k++) {
	if (u[i][j][k] > u_max) u_max = u[i][j][k];
	if (v[i][j][k] > v_max) v_max = v[i][j][k];
	if (w[i][j][k] > w_max) w_max = w[i][j][k];
      }
    }
  }
  return tau * min( (Re/2)/(1/pow(dx,2) + 1/pow(dy,2) + 1/pow(dz,2)), \
		    min(dx/u_max, \
			min(dy/v_max, dz/w_max)));
}

int main() {
  const int nx = 200;
  const int ny = 100;
  const int nz = 80;//150;
  const double dx = 20.0/(double)nx;
  const double dy = 10.0/(double)ny;
  //const double dz = 15.0/(double)nz;
  const double dz = 8.0/(double)nz;
  
  const double radius = 1.0/2;
  const double Re = 1000.0;
  const double SOR = 1.7;

  const double centerX = 5.0;
  const double centerY = 5.0;

  // **** ------------------------ ****
  const bool checkpoint = false;
  int readAt_iteration = 13000;
  // **** ------------------------ ****

  string f_name = "../vtk/JICF/JICF_";
  string name_prefix = "checkpoints/small_JICF/";
  const int precision = 4;
  const int save_precision = 8;
  
  double dt = 0.0001;
  double t = 0.0;
  double tau = 0.05;

  printf("Total nodes = %d \n", nx*ny*nz);
  cout << "Allocating varaibles ...\n";
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

  if (!checkpoint) {
    cout << "Initialize initial condition ...\n";
    // --------------------------------------------
    initialize_u_3D(u, nx, ny, nz, dx, dy, radius, centerX, centerY);
    initialize_v_3D(v, nx, ny, nz, dx, dy, radius, centerX, centerY);
    initialize_w_3D(w, nx, ny, nz, dx, dy, radius, centerX, centerY);
    initialize_u_3D(u_new, nx, ny, nz, dx, dy, radius, centerX, centerY);
    initialize_v_3D(v_new, nx, ny, nz, dx, dy, radius, centerX, centerY);
    initialize_w_3D(w_new, nx, ny, nz, dx, dy, radius, centerX, centerY);
    initialize_zero_3D(F_n, nx, ny, nz);
    initialize_zero_3D(G_n, nx, ny, nz);
    initialize_zero_3D(H_n, nx, ny, nz);
    initialize_phi_3D(phi, nx, ny, nz, dx, dy, radius, centerX, centerY);
    initialize_phi_3D(phi_new, nx, ny, nz, dx, dy, radius, centerX, centerY);
    initialize_pressure_3D(p, nx, ny, nz);

    cout << "Writing the initial .vtk file ...\n";
    paraview3D(f_name + to_string(1) + ".vtk", u_new, v_new, w_new, p, phi_new, nx, ny, nz, dx, dy, dz, precision, 0.0);
    readAt_iteration = 0;
  }
  else {
    cout << "Load checkpoint. Continue from timestep#" + to_string(readAt_iteration) << endl;
    read_restartfile(u, name_prefix, "u", readAt_iteration, nx, ny, nz);
    read_restartfile(v, name_prefix, "v", readAt_iteration, nx, ny, nz);
    read_restartfile(w, name_prefix, "w", readAt_iteration, nx, ny, nz);
    read_restartfile(p, name_prefix, "p", readAt_iteration, nx, ny, nz);
    read_restartfile(phi, name_prefix, "phi", readAt_iteration, nx, ny, nz);
    t = load_time(name_prefix, readAt_iteration, nx, ny, nz);
  }

  cout << "Start solving ...\n";
  for (int it = readAt_iteration + 1; it <= 15000; it++) {
    if (t >= 51) break;
    dt = variable_dt(tau, nx, ny, nz, Re, dx, dy, dz, u, v, w);
    cout << it << "\t" << dt << "\t" << t << "\t";
    printf("(tau=%1.2f, Re=%4.1f) \n", tau, Re);
    F_calculation(F_n, u, v, w, nx, ny, nz, dt, Re, dx, dy, dz);
    G_calculation(G_n, u, v, w, nx, ny, nz, dt, Re, dx, dy, dz);
    H_calculation(H_n, u, v, w, nx, ny, nz, dt, Re, dx, dy, dz);
    pressure_solver(p, SOR, F_n, G_n, H_n, nx, ny, nz, dt, dx, dy, dz, it);
    u_calculation(u_new, u, F_n, p, dt, dx, nx, ny, nz);
    v_calculation(v_new, v, G_n, p, dt, dy, nx, ny, nz);
    w_calculation(w_new, w, H_n, p, dt, dx, dy, dz, nx, ny, nz, centerX, centerY, radius);
    explicit_passiveScalar(phi, phi_new, u_new, v_new, w_new, nx, ny, nz, dx, dy, dz, dt, Re, centerX, centerY, radius);
    t += dt;
    if (it % 200 == 0) {
      paraview3D(f_name + to_string(it) + ".vtk", u_new, v_new, w_new, p, phi_new, nx, ny, nz, dx, dy, dz, precision, t);
    }
    if (it % 500 == 0) {
      save_restartfile(u, name_prefix, "u", it, nx, ny, nz, save_precision);
      save_restartfile(v, name_prefix, "v", it, nx, ny, nz, save_precision);
      save_restartfile(w, name_prefix, "w", it, nx, ny, nz, save_precision);
      save_restartfile(p, name_prefix, "p", it, nx, ny, nz, save_precision);
      save_restartfile(phi, name_prefix, "phi", it, nx, ny, nz, save_precision);
      save_time(t, name_prefix, it, nx, ny, nz);
    }
  }
}
