#include <iostream>
#include <cmath>
#include <fstream> // save & load restart files
#include <string>
#include <iomanip> // std::setprecision()
#include <algorithm> // std:min
using namespace std;

ofstream myfileO;  // output file stream
ifstream myfileI;  // input file stream

void save_restartfile2(double **var, string name_prefix, string variable_name, int iteration, int nx, int ny, int save_precision){
  ofstream myfileO;
  myfileO.open(name_prefix + to_string(nx) + "x" + to_string(ny) + "_" + variable_name + "_" + to_string(iteration) + ".dat");
  for (int i = 0; i <= nx-1; i++) {
    for (int j = 0; j <= ny-1; j++) {
      myfileO << setprecision(save_precision) << var[i][j] << " ";
    }
    myfileO << "\n";
  }
  myfileO.close();
}

void read_restartfile(double **var, string name_prefix, string variable_name, int readAt_iteration, int nx, int ny){
  ifstream myfileI;
  myfileI.open(name_prefix + to_string(nx) + "x" + to_string(ny) + "_" + variable_name + "_" + to_string(readAt_iteration) + ".dat");
  for (int i = 0; i <= nx-1; i++) {
    for (int j = 0; j <= ny-1; j++) {
      myfileI >> var[i][j];
    }
  }
  myfileI.close();
}

double load_time(string name_prefix, int readAt_iteration, int nx, int ny) {
  ifstream myfileI;
  double t;
  myfileI.open(name_prefix + to_string(nx) + "x" + to_string(ny) + "_t_" + to_string(readAt_iteration) + ".dat");
  myfileI >> t;
  myfileI.close();
  return t;
}

int main() {
  const int nx = 400;
  const int ny = 200;
  const double dy = (1.0/(double)ny);
  const double dx = 16.0*dy;
  
  string name_prefix = "/media/sf_Shared/shared_checkpoints/checkpoints/Validation_2DChannelFlow/";
  string name_prefix_2 = "/media/sf_Shared/shared_checkpoints/checkpoints/Validation_2DChannelFlow/post_processed/";
  const int save_precision = 10;
  
  double dt = 0.0001;
  double t = 0.0;
  double last_t = 0.0;
  double at = 0.0;

  const int start_steady_at = 24000;
  const int end_steady_at = 34000;
  const int step = 500;
  double first_t = load_time(name_prefix, start_steady_at - step, nx, ny);

  cout << "Allocating varaibles ...\n";
  // ---------------------------------------------
  double **u = (double **) malloc (nx * sizeof(double*));
  double **v = (double **) malloc (nx * sizeof(double*));
  double **phi = (double **) malloc (nx * sizeof(double*));
  double **p = (double **) malloc (nx * sizeof(double*));
 
  double **u_buff = (double **) malloc (nx * sizeof(double*));
  double **v_buff = (double **) malloc (nx * sizeof(double*));
  double **p_buff = (double **) malloc (nx * sizeof(double*));
  double **phi_buff = (double **) malloc (nx * sizeof(double*));

  double **u_tmp = (double **) malloc (nx * sizeof(double*));
  double **v_tmp = (double **) malloc (nx * sizeof(double*));
  double **p_tmp = (double **) malloc (nx * sizeof(double*));
  double **phi_tmp = (double **) malloc (nx * sizeof(double*));
  cout << "Finished \n";
  
  for (int i = 0; i < nx; i++) {
    u[i] = (double *) malloc (ny * sizeof(double));
    v[i] = (double *) malloc (ny * sizeof(double));
    phi[i] = (double *) malloc (ny * sizeof(double));
    p[i] = (double *) malloc (ny * sizeof(double));

    u_buff[i] = (double *) malloc (ny * sizeof(double));
    v_buff[i] = (double *) malloc (ny * sizeof(double));
    p_buff[i] = (double *) malloc (ny * sizeof(double));
    phi_buff[i] = (double *) malloc (ny * sizeof(double));

    u_tmp[i] = (double *) malloc (ny * sizeof(double));
    v_tmp[i] = (double *) malloc (ny * sizeof(double));
    p_tmp[i] = (double *) malloc (ny * sizeof(double));
    phi_tmp[i] = (double *) malloc (ny * sizeof(double));
  }

  last_t = first_t;
  for (int it = start_steady_at; it <= end_steady_at; it += step) {
    cout << "Reading timestep#" + to_string(it) << endl;
    read_restartfile(u, name_prefix, "u", it, nx, ny);
    read_restartfile(v, name_prefix, "v", it, nx, ny);
    read_restartfile(p, name_prefix, "p", it, nx, ny);
    read_restartfile(phi, name_prefix, "phi", it, nx, ny);
    t = load_time(name_prefix, it, nx, ny);
    dt = t - last_t;
    at += dt;
    cout << "dt = " + to_string(dt) << endl;
    cout << to_string(at) << endl;
    
    for (int i = 0; i <= nx-1; i++) {
      for (int j = 0; j <= ny-1; j++) {   
	u_buff[i][j] += u[i][j]*dt;
	v_buff[i][j] += v[i][j]*dt;
	p_buff[i][j] += p[i][j]*dt;
	phi_buff[i][j] += phi[i][j]*dt;
      }
    }
    last_t = t;
  }

  cout << "Averaging ...\n";
  for (int i = 0; i <= nx-1; i++) {
    for (int j = 0; j <= ny-1; j++) {
      u_tmp[i][j] = u_buff[i][j]/at;
      v_tmp[i][j] = v_buff[i][j]/at;
      p_tmp[i][j] = p_buff[i][j]/at;
      phi_tmp[i][j] = phi_buff[i][j]/at;
    }
  }

  cout << "Saving Temporal-Means ...";
  save_restartfile2(u_tmp, name_prefix_2, "u_tmp", end_steady_at, nx, ny, save_precision);
  save_restartfile2(v_tmp, name_prefix_2, "v_tmp", end_steady_at, nx, ny, save_precision);
  save_restartfile2(p_tmp, name_prefix_2, "p_tmp", end_steady_at, nx, ny, save_precision);
  save_restartfile2(phi_tmp, name_prefix_2, "phi_tmp", end_steady_at, nx, ny, save_precision);
  cout << "Finished! \n";

  // Clean buffers
  for (int i = 0; i <= nx-1; i++) {
    for (int j = 0; j <= ny-1; j++) {
      u_buff[i][j] = 0.0;
      v_buff[i][j] = 0.0;
      p_buff[i][j] = 0.0;
      phi_buff[i][j] = 0.0;
    }
  }
  

  // ################  RMS-Fluctuation ################
  last_t = first_t;
  for (int it = start_steady_at; it <= end_steady_at; it += step) {
    cout << "Reading timestep#" + to_string(it) << endl;
    read_restartfile(u, name_prefix, "u", it, nx, ny);
    read_restartfile(v, name_prefix, "v", it, nx, ny);
    t = load_time(name_prefix, it, nx, ny);
    dt = t - last_t;
    cout << "dt = " + to_string(dt) << endl;

    for (int i = 0; i <= nx-1; i++) {
      for (int j = 0; j <= ny-1; j++) {
	u_buff[i][j] += pow(u[i][j] - u_tmp[i][j], 2)*dt;
	v_buff[i][j] += pow(v[i][j] - v_tmp[i][j], 2)*dt;
      }
    }
    last_t = t;
  }

  cout << "Performing Root-Mean-Square-Fluc  ...\n";
  for (int i = 0; i <= nx-1; i++) {
    for (int j = 0; j <= ny-1; j++) {
      u_tmp[i][j] = sqrt(u_buff[i][j]/at);
      v_tmp[i][j] = sqrt(v_buff[i][j]/at);
    }
  }

  cout << "Saving RMS-Fluctuation ...";
  save_restartfile2(u_tmp, name_prefix_2, "u_rms", end_steady_at, nx, ny, save_precision);
  save_restartfile2(v_tmp, name_prefix_2, "v_rms", end_steady_at, nx, ny, save_precision);
  cout << "Finished! \n";
}