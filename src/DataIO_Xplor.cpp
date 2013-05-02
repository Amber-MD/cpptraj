#include "DataIO_Xplor.h"
#include "CpptrajStdio.h"
#include "DataSet_3D.h"

int DataIO_Xplor::WriteData3D(std::string const& fname, DataSet const& setIn,
                              DimArray const& Dim)
{
  if (setIn.Ndim() != 3) {
    mprinterr("Internal Error: DataSet %s in DataFile %s has %zu dimensions, expected 3.\n",
              setIn.Legend().c_str(), fname.c_str(), setIn.Ndim());
    return 1;
  }
  DataSet_3D const& set = static_cast<DataSet_3D const&>( setIn );
  // Open output file
  CpptrajFile outfile;
  if (outfile.OpenWrite( fname )) {
    mprinterr("Error: Could not open Xplor output file.\n");
    return 1;
  }
  // Title
  outfile.Printf("%s\n", title_.c_str());
  // Remarks - only 1
  outfile.Printf("%8i\n%s\n",1,remark_.c_str());
  // Header
  outfile.Printf("%8i%8i%8i",   set.NX(), -set.NX()/2 + 1, set.NX()/2 );
  outfile.Printf("%8i%8i%8i",   set.NY(), -set.NY()/2 + 1, set.NY()/2 );
  outfile.Printf("%8i%8i%8i\n", set.NZ(), -set.NZ()/2 + 1, set.NZ()/2 );
  outfile.Printf("%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n",
                 (double)set.NX() * set.DX(), 
                 (double)set.NY() * set.DY(), 
                 (double)set.NZ() * set.DZ(),
                 90.0, 90.0, 90.0);
  outfile.Printf("ZYX\n");
  // Print grid bins
  size_t NZ2 = set.NZ() / 2;
  for (size_t k = 0; k < set.NZ(); ++k) {
    outfile.Printf("%8i\n", k - NZ2 + 1);
    for (size_t j = 0; j < set.NY(); ++j) {
      int col = 1;
      for (size_t i = 0; i < set.NX(); ++i) {
        outfile.Printf("%12.5f", set.GetElement(i, j, k));
        if ( (col % 6)==0 )
          outfile.Printf("\n");
        ++col;
      }
      if ( (col-1) % 6 != 0 )
        outfile.Printf("\n");
    }
  }
  outfile.CloseFile();
  return 0;
}
