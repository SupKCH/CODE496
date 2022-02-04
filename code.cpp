#include <iostream> // cout
#include <cmath>
using namespace std; // std = standard

int i,j;
const int nx = 31;
const int ny = 31;
double phi[nx][ny];

void initialize(){
  // initialize variable phi
  for (i = 0; i <= nx-1; i++) {
    for (j = 0; j <= ny-1; j++) {
      phi[i][j] = 0.0;
    }
  }
}

void visualize(){
   for (i = 0; i <= nx-1; i++) {
     for (j = 0; j <= ny-1; j++) {
      cout << phi[i][j] << " ";
    }
     cout << "\n";
  } 
}

void set_phi() {
  // Assign phi = 1 inside a circle around center of radius = 15
  // letting dx = dy = 1
  // as e.g., x_spaceing = dx*diff(index), not diff(index) itself

  int i_c = (nx-1)/2;
  int j_c = (ny-1)/2;
  double radius = 10.;

  for (i = 0; i <= nx-1; i++) {
     for (j = 0; j <= ny-1; j++) {

       if ( sqrt(pow(i-i_c,2) + pow(j-j_c,2)) < radius ) {
	 phi[i][j] = 1.;
       }
       else {phi[i][j] = 0.;}
    }   
  }
}

int main() {

  // initiallize variable phi
  initialize();
  set_phi();
  visualize();
  
} // end main
