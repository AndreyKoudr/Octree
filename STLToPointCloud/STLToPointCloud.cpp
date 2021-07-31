#include "../CloudPoints.h"

#include <iostream>

/**
  STLToPointCloud
  ===============
  Converter from STL to point cloud XYZ file with space separator in lines. 

  HOW THIS PROJECT WORKS

  - Purpose : generate a dense point cloud from a parametric geometry in STL file. Points
to be distributed unformly along triangluar surfaces at a specified mean distance. This distance is NOT 
absolute, this is a fraction of the whole model size. Input STL file may be both text and binary.

  - This project is mostly useful for generation of sample point clouds for debugging given STL files.

  - Generated XYZ file is a text file with three X,Y,Z point components in every line. They are
separated by spaces. Change line 188 in CloudPoints.h from
  
    file << v.X << " " << v.Y << " " << v.Z << "\n"; 

to

    file << v.X << "\t" << v.Y << "\t" << v.Z << "\n"; 

if you need tabs instead of spaces.
*/

int main(int argc, char* argv[])
{
  printf("\n  Converter from STL to point cloud XYZ file with space separator in lines.\n");
  printf("Call : >STLToPointCloud STL_file_name point_density\n");
  printf("e.g.\n");
  printf(">STLToPointCloud shuttle.stl 0.0005\n");
  printf("will convert shuttle.stl into shuttle.0.0005.xyz\n");
  printf("with mean distance between points 0.0005 * model_size\n\n");

  if (argc != 3) 
    return 1;

  std::string infilename = argv[1];
  double density = atof(argv[2]);
  std::string outfilename = forceExtension(infilename,"") + "." + argv[2] + ".xyz";

  printf("Output : %s\n",outfilename.c_str());

  CloudPoints<double> points;

  if (points.readSTL1(infilename,density,true))
  {
    points.writeXYZ(outfilename);
  }

  return 0;
}
