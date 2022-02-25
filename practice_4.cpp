
// Practice creating variable (phi) array

// g++ practice_init.cpp -o practice_init.exe -L/usr/lib/x86_64-linux-gpu/ -llapack

#include <iostream>
//#include <lapacke.h>
#include <fstream>
using namespace std;



/* =================== Subroutine ===================== */
void node_initialize(double* var, double* var_new, int nx) {
  for (int i=0; i <= nx-1; i++) {
    var[i] = 0.0;
    var_new[i] = 0.0;
  }
}

void node_1DHeatCondition(double* var, int nx) {
  var[int((nx+1)/2-1)] = 100.;
}

void visualize(double* var, int nx) {
  for (int i=0; i <= nx-1; i++) {
    cout << var[i] << ' ';
  }
  cout << "\n";
}

void rungekutta2(double* var, double* var_new, int nx, double k, double dt, double dx) {
  for (int i=2; i <= nx-3; i++) {
    double rhs_nowright = k/(dx*dx)*(var[i+2] - 2*var[i+1] + var[i]);
    //cout << rhs_nowright << "\t";
    double rhs_nowmid = k/(dx*dx)*(var[i+1] - 2*var[i] + var[i-1]);
    //cout << rhs_nowmid << "\t";
    double rhs_nowleft = k/(dx*dx)*(var[i] - 2*var[i-1] + var[i-2]);
    //    cout << rhs_nowleft << "\n";
    double rhs_new = k/(dx*dx)*( var[i+1] - 2*var[i] + var[i-1] + dt*(rhs_nowright - 2*rhs_nowmid + rhs_nowleft));
    // cout << rhs_nowmid << "\t";
    // cout << rhs_new << "\t";
    // cout << dt*(rhs_nowmid + rhs_new)/2 << "\t";
    // cout << var[i] << "\t";
    // cout << var[i] + dt*(rhs_nowmid + rhs_new)/2 << "\n";
    //var_new[i] = var[i] + dt*(1/2*rhs_nowmid + 1/2*rhs_new);
    var_new[i] = var[i] + dt*(rhs_nowmid + rhs_new)/2;
    //cout << var_new[i] << "\n";
    //cout << dt*( 1/2*rhs_nowmid + 1/2*rhs_new ) << "\n";
  }
  for (int i=2; i <= nx-3; i++) {
    var[i] = var_new[i];
  }
}

double k1_cal(double* var, int i, double coeff) {
  double k1 = coeff * (var[i+1] - 2*var[i] + var[i-1]);
  return k1;
}

double k2_cal(double* var, int i, double coeff,  double dt) {
  double k2 = coeff * (var[i+1] - 2*var[i] + var[i-1] + dt*(k1_cal(var, i+1, coeff) - 2*k1_cal(var, i, coeff) + k1_cal(var, i-1, coeff))/2);
  return k2;
}

double k3_cal(double* var, int i, double coeff,  double dt) {
  double k3 = coeff * (var[i+1] - 2*var[i] + var[i-1] + dt*(k2_cal(var, i+1, coeff, dt) - 2*k2_cal(var, i, coeff, dt) + k2_cal(var, i-1, coeff, dt))/2);
  return k3;
}

void rungekutta4(double* var, double* var_new, int nx, double k, double dt, double dx) {
  // calculate k1, k2, k3, k4, = different styles of RHS
  double coeff = k/(dx*dx);
  for (int i = 4; i <= nx-1-4; i++) { 
    double k1 = k1_cal(var, i, coeff);
    double k2 = k2_cal(var, i, coeff, dt);
    double k3 = k3_cal(var, i, coeff, dt);
    double k4 = coeff * (var[i+1] - 2*var[i] + var[i-1] + dt*(k3_cal(var, i+1, coeff, dt) - 2*k3_cal(var, i, coeff, dt) + k3_cal(var, i-1, coeff, dt)));
    var_new[i] = var[i] + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4);
  }
  for (int i=0; i <= nx-1; i++) {
    var[i] = var_new[i];
  }
}


int main() {
  int nx = 25;
  double k = 1.0;
  double dt = 0.001;
  double dx = 0.1;
  double phi[nx];
  double phi_new[nx];
  node_initialize(phi, phi_new,  nx);
  node_1DHeatCondition(phi, nx);
  visualize(phi, nx);
  /* for (int i=1; i <= 20; i++) {
     rungekutta2(phi, phi_new, nx, k, dt, dx);
     visualize(phi_new, nx);
     }*/

  for (int iteration=0; iteration <= 20; iteration++) {
    rungekutta4(phi, phi_new, nx, k, dt, dx);
    visualize(phi_new, nx);
  }
}
