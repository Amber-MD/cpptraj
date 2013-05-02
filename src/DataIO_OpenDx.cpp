#include "DataIO_OpenDx.h"
#include "CpptrajStdio.h"
#include "DataSet_3D.h"

int DataIO_OpenDx::WriteData3D(std::string const& fname, DataSet const& setIn,
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
  if (outfile.OpenWrite(fname)) {
    mprinterr("Error: Could not open OpenDX output file.\n");
    return 1;
  }
  // Print the OpenDX header
  size_t gridsize = set.Size();
  outfile.Printf("object 1 class gridpositions counts %d %d %d\n",
                 set.NX(), set.NY(), set.NZ());
  outfile.Printf("origin %lg %lg %lg\n", set.OX(), set.OY(), set.OZ());
  outfile.Printf("delta %lg 0 0\n", set.DX());
  outfile.Printf("delta 0 %lg 0\n", set.DY());
  outfile.Printf("delta 0 0 %lg\n", set.DZ());
  outfile.Printf("object 2 class gridconnections counts %d %d %d\n",
                 set.NX(), set.NY(), set.NZ());
  outfile.Printf(
    "object 3 class array type double rank 0 items %d data follows\n",
    gridsize);

  // Now print out the data. It is already in row-major form (z-axis changes
  // fastest), so no need to do any kind of data adjustment
  for (size_t i = 0UL; i < gridsize - 2UL; i += 3UL)
    outfile.Printf("%g %g %g\n", set[i], set[i+1], set[i+2]);
  // Print out any points we may have missed
  switch (gridsize % 3) {
    case 2: outfile.Printf("%g %g\n", set[gridsize-2], set[gridsize-1]); break;
    case 1: outfile.Printf("%g\n", set[gridsize-1]); break;
  }

  // Print tail
  // TODO: Make this an option
  //if (mode_ == CENTER)
  //  outfile.Printf("\nobject \"density (%s) [A^-3]\" class field\n",
  //                 centerMask_.MaskString());
  //else
    outfile.Printf("\nobject \"density [A^-3]\" class field\n");
  outfile.CloseFile();
  return 0;
}

