#include <iostream>
#include <cmath>
#include <fstream> // save & load restart files
#include <string>
#include <iomanip> // std::setprecision()
#include <algorithm> // std:min
using namespace std;

ofstream myfileO;  // output file stream
ifstream myfileI;  // input file stream

void save_restartfile2(double ***var, string name_prefix, string variable_name, int iteration, int nx, int ny, int nz, int save_precision){
  ofstream myfileO;
  myfileO.open(name_prefix + to_string(nx) + "x" + to_string(ny) + "x" + to_string(nz) + "_" + variable_name + "_" + to_string(iteration) + ".dat");
  for (int i = 0; i <= nx-1; i++) {
    for (int j = 0; j <= ny-1; j++) {
      for (int k = 0; k <= nz-1; k++) {
	myfileO << setprecision(save_precision) << var[i][j][k] << " ";
      }
      myfileO << "\n";
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

double load_time(string name_prefix, int readAt_iteration, int nx, int ny, int nz) {
  ifstream myfileI;
  double t;
  myfileI.open(name_prefix + to_string(nx) + "x" + to_string(ny) + "x" + to_string(nz) + "_t_" + to_string(readAt_iteration) + ".dat");
  myfileI >> t;
  myfileI.close();
  return t;
}

int main() {
  const int nx = 200;
  const int ny = 100;
  const int nz = 150;//150;
  const double dx = 20.0/(double)nx;
  const double dy = 10.0/(double)ny;
  const double dz = 15.0/(double)nz;

  //  string f_name = "/media/sf_Shared/shared_vtk/3DJICF/JICF_";
  string name_prefix = "/media/sf_Shared/shared_checkpoints/checkpoints/3DJICF/";
  string name_prefix_2 = "/media/sf_Shared/shared_checkpoints/checkpoints/3DJICF/post_processed/";
  const int save_precision = 10;
  
  double dt = 0.0001;
  double t = 0.0;
  double last_t = 0.0;
  double at = 0.0;

  const int start_steady_at = 17800; // 20000 just picked a number
  const int end_steady_at = 25800;
  const int step = 100;
  double first_t = load_time(name_prefix, start_steady_at - step, nx, ny, nz);

  cout << "Allocating varaibles ...\n";
  // ---------------------------------------------
  double ***u = (double ***) malloc (nx * sizeof(double**));
  double ***v = (double ***) malloc (nx * sizeof(double**));
  double ***w = (double ***) malloc (nx * sizeof(double**));
  double ***phi = (double ***) malloc (nx * sizeof(double**));
  double ***p = (double ***) malloc (nx * sizeof(double**));
  double ***ke = (double ***) malloc (nx * sizeof(double**));
  
  double ***u_buff = (double ***) malloc (nx * sizeof(double**));
  double ***v_buff = (double ***) malloc (nx * sizeof(double**));
  double ***w_buff = (double ***) malloc (nx * sizeof(double**));
  double ***p_buff = (double ***) malloc (nx * sizeof(double**));
  double ***phi_buff = (double ***) malloc (nx * sizeof(double**));
  double ***ke_buff = (double ***) malloc (nx * sizeof(double**));

  double ***u_tmp = (double ***) malloc (nx * sizeof(double**));
  double ***v_tmp = (double ***) malloc (nx * sizeof(double**));
  double ***w_tmp = (double ***) malloc (nx * sizeof(double**));
  double ***p_tmp = (double ***) malloc (nx * sizeof(double**));
  double ***phi_tmp = (double ***) malloc (nx * sizeof(double**));
  double ***ke_tmp = (double ***) malloc (nx * sizeof(double**));
  cout << "Finished \n";
  
  for (int i = 0; i < nx; i++) {
    u[i] = (double **) malloc (ny * sizeof(double*));
    v[i] = (double **) malloc (ny * sizeof(double*));
    w[i] = (double **) malloc (ny * sizeof(double*));
    phi[i] = (double **) malloc (ny * sizeof(double*));
    p[i] = (double **) malloc (ny * sizeof(double*));
    ke[i] = (double **) malloc (ny * sizeof(double*));

    u_buff[i] = (double **) malloc (ny * sizeof(double*));
    v_buff[i] = (double **) malloc (ny * sizeof(double*));
    w_buff[i] = (double **) malloc (ny * sizeof(double*));
    p_buff[i] = (double **) malloc (ny * sizeof(double*));
    phi_buff[i] = (double **) malloc (ny * sizeof(double*));
    ke_buff[i] = (double **) malloc (ny * sizeof(double*));

    u_tmp[i] = (double **) malloc (ny * sizeof(double*));
    v_tmp[i] = (double **) malloc (ny * sizeof(double*));
    w_tmp[i] = (double **) malloc (ny * sizeof(double*));
    p_tmp[i] = (double **) malloc (ny * sizeof(double*));
    phi_tmp[i] = (double **) malloc (ny * sizeof(double*));
    ke_tmp[i] = (double **) malloc (ny * sizeof(double*));
    
    for (int ii = 0; ii < ny; ii++) {
      u[i][ii] = (double *) malloc (nz * sizeof(double));
      v[i][ii] = (double *) malloc (nz * sizeof(double));
      w[i][ii] = (double *) malloc (nz * sizeof(double));
      phi[i][ii] = (double *) malloc (nz * sizeof(double));
      p[i][ii] = (double *) malloc (nz * sizeof(double));
      ke[i][ii] = (double *) malloc (nz * sizeof(double));

      u_buff[i][ii] = (double *) malloc (nz * sizeof(double));
      v_buff[i][ii] = (double *) malloc (nz * sizeof(double));
      w_buff[i][ii] = (double *) malloc (nz * sizeof(double));
      p_buff[i][ii] = (double *) malloc (nz * sizeof(double));
      phi_buff[i][ii] = (double *) malloc (nz * sizeof(double));
      ke_buff[i][ii] = (double *) malloc (nz * sizeof(double));

      u_tmp[i][ii] = (double *) malloc (nz * sizeof(double));
      v_tmp[i][ii] = (double *) malloc (nz * sizeof(double));
      w_tmp[i][ii] = (double *) malloc (nz * sizeof(double));
      p_tmp[i][ii] = (double *) malloc (nz * sizeof(double));
      phi_tmp[i][ii] = (double *) malloc (nz * sizeof(double));
      ke_tmp[i][ii] = (double *) malloc (nz * sizeof(double));
    }
  }

  last_t = first_t;
  for (int it = start_steady_at; it <= end_steady_at; it += step) {
    cout << "Reading timestep#" + to_string(it) << endl;
    read_restartfile(u, name_prefix, "u", it, nx, ny, nz);
    read_restartfile(v, name_prefix, "v", it, nx, ny, nz);
    read_restartfile(w, name_prefix, "w", it, nx, ny, nz);
    read_restartfile(p, name_prefix, "p", it, nx, ny, nz);
    read_restartfile(phi, name_prefix, "phi", it, nx, ny, nz);
    t = load_time(name_prefix, it, nx, ny, nz);
    dt = t - last_t;
    at += dt;
    cout << "dt = " + to_string(dt) << endl;
    cout << to_string(at) << endl;
    
    for (int i = 0; i <= nx-1; i++) {
      for (int j = 0; j <= ny-1; j++) {
	for (int k = 0; k <= nz-1; k++) {
	  u_buff[i][j][k] += u[i][j][k]*dt;
	  v_buff[i][j][k] += v[i][j][k]*dt;
	  w_buff[i][j][k] += w[i][j][k]*dt;
	  p_buff[i][j][k] += p[i][j][k]*dt;
	  phi_buff[i][j][k] += phi[i][j][k]*dt;
	  ke_buff[i][j][k] += (\
			       pow(u[i][j][k],2) + \
			       pow(v[i][j][k],2) + \
			       pow(w[i][j][k],2) )*dt;
	}
      }
    }
    last_t = t;
  }

  cout << "Averaging ...\n";
  for (int i = 0; i <= nx-1; i++) {
    for (int j = 0; j <= ny-1; j++) {
      for (int k = 0; k <= nz-1; k++) {
	u_tmp[i][j][k] = u_buff[i][j][k]/at;
	v_tmp[i][j][k] = v_buff[i][j][k]/at;
	w_tmp[i][j][k] = w_buff[i][j][k]/at;
	p_tmp[i][j][k] = p_buff[i][j][k]/at;
	phi_tmp[i][j][k] = phi_buff[i][j][k]/at;
	ke_tmp[i][j][k] = ke_buff[i][j][k]/at;
      }
    }
  }

  cout << "Saving Temporal-Means ...";
  save_restartfile2(u_tmp, name_prefix_2, "u_tmp", end_steady_at, nx, ny, nz, save_precision);
  save_restartfile2(v_tmp, name_prefix_2, "v_tmp", end_steady_at, nx, ny, nz, save_precision);
  save_restartfile2(w_tmp, name_prefix_2, "w_tmp", end_steady_at, nx, ny, nz, save_precision);
  save_restartfile2(p_tmp, name_prefix_2, "p_tmp", end_steady_at, nx, ny, nz, save_precision);
  save_restartfile2(phi_tmp, name_prefix_2, "phi_tmp", end_steady_at, nx, ny, nz, save_precision);
  save_restartfile2(ke_tmp, name_prefix_2, "ke_tmp", end_steady_at, nx, ny, nz, save_precision);
  cout << "Finished! \n";

  // Clean buffers
  for (int i = 0; i <= nx-1; i++) {
    for (int j = 0; j <= ny-1; j++) {
      for (int k = 0; k <= nz-1; k++) {
	u_buff[i][j][k] = 0.0;
	v_buff[i][j][k] = 0.0;
	w_buff[i][j][k] = 0.0;
	p_buff[i][j][k] = 0.0;
	phi_buff[i][j][k] = 0.0;
	ke_buff[i][j][k] = 0.0;
      }
    }
  }
  

    // ################  RMS-Fluctuation ################
  last_t = first_t;
  for (int it = start_steady_at; it <= end_steady_at; it += step) {
    cout << "Reading timestep#" + to_string(it) << endl;
    read_restartfile(u, name_prefix, "u", it, nx, ny, nz);
    read_restartfile(v, name_prefix, "v", it, nx, ny, nz);
    read_restartfile(w, name_prefix, "w", it, nx, ny, nz);
    t = load_time(name_prefix, it, nx, ny, nz);
    dt = t - last_t;
    cout << "dt = " + to_string(dt) << endl;

    for (int i = 0; i <= nx-1; i++) {
      for (int j = 0; j <= ny-1; j++) {
	for (int k = 0; k <= nz-1; k++) {
	  u_buff[i][j][k] += pow(u[i][j][k] - u_tmp[i][j][k], 2)*dt;
	  v_buff[i][j][k] += pow(v[i][j][k] - v_tmp[i][j][k], 2)*dt;
	  w_buff[i][j][k] += pow(w[i][j][k] - w_tmp[i][j][k], 2)*dt;
  	  // ke_buff[i][j][k] += (pow(		   \
  	  // 			   pow(u[i][j][k],2) + \
   	  // 			   pow(v[i][j][k],2) + \
   	  // 			   pow(w[i][j][k],2) - \
   	  // 			   ke_tmp[i][k][k], 2)\
   	  // 		       )*dt;
  	}
      }
    }
    last_t = t;
  }

  cout << "Performing Root-Mean-Square-Fluc  ...\n";
  for (int i = 0; i <= nx-1; i++) {
    for (int j = 0; j <= ny-1; j++) {
      for (int k = 0; k <= nz-1; k++) {
  	u_tmp[i][j][k] = sqrt(u_buff[i][j][k]/at);
  	v_tmp[i][j][k] = sqrt(v_buff[i][j][k]/at);
  	w_tmp[i][j][k] = sqrt(w_buff[i][j][k]/at);
	// ke_tmp[i][j][k] = sqrt(ke_buff[i][j][k]/at);
      }
    }
  }

  cout << "Saving RMS-Fluctuation ...";
  save_restartfile2(u_tmp, name_prefix_2, "u_rms", end_steady_at, nx, ny, nz, save_precision);
  save_restartfile2(v_tmp, name_prefix_2, "v_rms", end_steady_at, nx, ny, nz, save_precision);
  save_restartfile2(w_tmp, name_prefix_2, "w_rms", end_steady_at, nx, ny, nz, save_precision);
  // save_restartfile2(p_tmp, name_prefix_2, "p_rms", end_steady_at, nx, ny, nz, save_precision);
  // save_restartfile2(phi_tmp, name_prefix_2, "phi_rms", end_steady_at, nx, ny, nz, save_precision);
  // save_restartfile2(ke_tmp, name_prefix_2, "ke_rms", end_steady_at, nx, ny, nz, save_precision);
  cout << "Finished! \n";
  
}
