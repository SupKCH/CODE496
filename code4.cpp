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

void test2D(int nx, int ny, double **p) {
  for (int i = 0 ;i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      cout << p[i][j];
    }
  }
}

int main() {

  const int nx = 10;
  const int ny = 10;
  //double phi[nx][ny];

  double **phi;  // pointer to pointer
  phi = (double **) malloc (nx * sizeof(double));
  for (int i=0; i<nx; i++) {
    phi[i] = (double *) malloc (ny * sizeof(double));
  }
  test2D(nx, ny, phi);

  cout << "\n";

  // const int nx = 10;
  // // double phi[nx];

  // // Alternatively, ...
  // double *phi;
  // phi = (double *) malloc (nx * sizeof(double));

  // for (int i = 0; i < nx; i++) {
  //   phi[i] = double(i);
  // }
  // test(nx, phi);
  
  // int a = 10;
  // cout << " a = " << a << "\n";
  // cout << "&a = " << &a << "\n";   // referencing
  // cout << "\n";

  // int array[10];
  // for (int i = 0; i < 10; i++) {
  //   array[i] = double(i);
  //   cout << "array[" << i << "] = " << array[i] << ", stored at " << &array[i] << "\n";
  // }

  // cout << "\n";
  // int *p;
  // p = &a;
  // cout << "p = " << p << "\n";
  // cout << "*p = " << *p << "\n";   // de-referencing

  // cout << "&p = " << &p << "\n";
  
  
  // const int nx = 10;
  // double phi[nx];
  // for (int i = 0; i < nx; i++) {
  //   phi[i] = double(i);
  //   cout << "phi[" << i << "]" << phi[i] << "]n";
  // }
  // cout << "\n";
  // test(nx, phi);
  // cout << "phi[0] = " << phi[0] << "\n";
  // cout << "\n";
  // cout << "nx = " << nx << "\n";

  
  // initiallize variable phi
  // initialize();
  // set_phi();
  // visualize();
  // save_restartfile();
  //read_restartfile();
  //visualize();
  // summation(nx, ny);
  
} // end main
