#include <iostream> // cout
#include <cmath>  // use pow & sqrt
#include <fstream>  // file I/O
using namespace std; // std = standard

int i, j;
int i_c, j_c;
const int nx = 31;
const int ny = 31;
double phi[nx][ny];

ofstream myfileO;  // output file stream
ifstream myfileI;  // input file stream

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

  i_c = (nx-1)/2;
  j_c = (ny-1)/2;
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

// routine, function
void save_restartfile(){
  myfileO.open("file.dat");
  for (i = 0; i <= nx-1; i++) {
     for (j = 0; j <= ny-1; j++) {
      myfileO << phi[i][j] << " ";
    }
     myfileO << "\n";
  }
  myfileO.close();
}

void read_restartfile(){
  myfileI.open("file.dat");
  for (i = 0; i <= nx-1; i++) {
     for (j = 0; j <= ny-1; j++) {
       myfileI >> phi[i][j];
    }
  }
  myfileI.close();
}

int main() {

  initialize();
  set_phi();
  // visualize();
  // save_restartfile();
  //read_restartfile();
  visualize();

  ofstream myfile;
  myfile.open("phi_paraview.vtk");

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
      myfile << i << " " << j << " 0\n";
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
      myfile << phi[i][j] << "\n";
    }
  }
  
  myfile.close();

  
  
} // end main
