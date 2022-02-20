#include <fstream> // file I/O
using namespace std;

int main() {
  ofstream myfile;
  myfile.open("test.vtk");

  // Paraview header
  myfile << "# vtk DataFile Version 2.0\n";
  myfile << "FlowField\n";
  myfile << "ASCII\n";

  // Grid
  myfile << "DATASET STRUCTURED_GRID\n";
  myfile << "DIMENSIONS " << 2 << " " << 1 << " " << 2 << "\n";
  myfile << "POINTS " << 2*1*2 << " float\n";
  myfile << 1.0 << " " << 1.0 << " 0\n";
  myfile << 1.0 << " " << 2.0 << " 0\n";
  myfile << 2.0 << " " << 1.0 << " 0\n";
  myfile << 2.0 << " " << 2.0 << " 0\n";

  // Data
  myfile << "\n";
  myfile << "POINT_DATA";
  myfile << " " << 2*2 << "\n";

  myfile << "\n";
  myfile << "SCALARS PHI float 1\n";
  myfile << "LOOKUP_TABLE default\n";
  myfile << 1.0 << "\n";
  myfile << 2.0 << "\n";
  myfile << 3.0 << "\n";
  myfile << 4.0 << "\n";

  myfile.close();

}
