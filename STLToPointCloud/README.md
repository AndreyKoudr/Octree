  STLToPointCloud
  ===============
  Converter from STL to point cloud XYZ file with space separator in lines.<br />

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
