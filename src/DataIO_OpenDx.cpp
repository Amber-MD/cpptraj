#include <cstdio> // sscanf
#include <cstdlib> // atof
#include "DataIO_OpenDx.h"
#include "CpptrajStdio.h"
#include "DataSet_GridFlt.h"
#include "BufferedLine.h"

// DataIO_OpenDx::ReadData()
int DataIO_OpenDx::ReadData(std::string const& fname, DataSetList& datasetlist) {
  // Add grid data set. Default to float for now.
  DataSet* ds = datasetlist.AddSet( DataSet::GRID_FLT, fname, "GRID" );
  if (LoadGrid(fname.c_str(), *ds)) {
    // Load failed. Erase grid data set.
    DataSetList::const_iterator last = datasetlist.end();
    --last;
    datasetlist.erase( last );
    return 1;
  }
  return 0;
}

// DataIO_OpenDx::LoadGrid()
int DataIO_OpenDx::LoadGrid(const char* filename, DataSet& ds)
{
  DataSet_GridFlt& grid = static_cast<DataSet_GridFlt&>( ds );
  // Open file
  BufferedLine infile;
  if (infile.OpenRead(filename)) return 1;
  // Skip comments
  std::string line = infile.GetLine();
  while (!line.empty() && line[0] == '#') line = infile.GetLine();
  if (line.empty()) {
    mprinterr("Error: Unexpected EOF in DX file %s\n", filename);
    return 1;
  }
  // object 1 class gridpositions counts nx ny nz
  int nx, ny, nz;
  if (sscanf(line.c_str(), "object 1 class gridpositions counts %d %d %d",
             &nx, &ny, &nz) != 3)
  {
    mprinterr("Error: Reading grid counts from DX file %s\n", filename);
    return 1;
  }
  // origin xmin ymin zmin 
  int oxyz[3];
  line = infile.GetLine();
  if (sscanf(line.c_str(), "origin %i %i %i", oxyz, oxyz+1, oxyz+2) != 3) {
    mprinterr("Error: Reading origin line from DX file %s\n", filename);
    return 1;
  }
  // 3x 'delta hx hy hz'
  double dxyz[3];
  double dx, dy, dz;
  dx = dy = dz = 0.0;
  for (int i = 0; i < 3; i++) {
    line = infile.GetLine();
    if (sscanf(line.c_str(), "delta %lg %lg %lg", dxyz, dxyz+1, dxyz+2) != 3) {
      mprinterr("Error: Reading delta line from DX file %s\n", filename);
      return 1;
    }
    // Check that only 1 of the 3 values is non-zero
    if (dxyz[i] != (dxyz[0] + dxyz[1] + dxyz[2])) {
      mprinterr("Error: Rotated basis in DX file %s not yet supported.\n", filename);
      return 1;
    }
    switch (i) {
      case 0: dx = dxyz[0]; break;
      case 1: dy = dxyz[1]; break;
      case 2: dz = dxyz[2]; break;
    }
  }
  // object 2 class gridconnections counts nx ny nz
  int nxyz[3];
  line = infile.GetLine();
  if (sscanf(line.c_str(), "object 2 class gridconnections counts %d %d %d",
             nxyz, nxyz+1, nxyz+2) != 3)
  {
    mprintf("Error: Reading grid connections from DX file %s\n", filename);
    return 1;
  }
  // Sanity check for conflicting grid dimensions
  if (nxyz[0] != nx || nxyz[1] != ny || nxyz[2] != nz) {
    mprinterr("Error: Conflicting grid dimensions in input DX density file %s.\n",
              filename);
    mprinterr("Error: Grid positions: %d %d %d\n", grid.NX(), grid.NY(), grid.NZ());
    mprinterr("Error: Grid connections: %d %d %d\n", nxyz[0], nxyz[1], nxyz[2]);
    return 1;
  }
  // object 3 class array type <type> rank <r> times <i>
  // This line describes whether data will be in binary or ascii format.
  line = infile.GetLine();
  if (line.compare(0, 8, "object 3") != 0) {
    mprinterr("Error: DX file %s; expected 'object 3 ...', got [%s]\n",
              filename, line.c_str());
    return 1;
  }
  if (line.find("binary") != std::string::npos) {
    mprinterr("Error: DX file %s; binary DX files not yet supported.\n", filename);
    return 1;
  }
  // Allocate Grid from dims, origin, and spacing
  if (grid.Allocate_N_O_D(nx,ny,nz, Vec3(oxyz), Vec3(dx,dy,dz))) {
    mprinterr("Error: Could not allocate grid.\n");
    return 1;
  }
  // Read in data
  size_t gridsize = grid.Size();
  size_t ndata = 0;
  infile.SetupBuffer();
  while (ndata < gridsize) {
    if (infile.Line() == 0) {
      mprinterr("Error: Unexpected EOF hit in %s\n", filename);
      infile.CloseFile();
      return 1;
    }
    int nTokens = infile.TokenizeLine(" \t");
    for (int j = 0; j < nTokens; j++) {
      if (ndata >= gridsize) {
        mprinterr("Error: Too many grid points found!\n");
        infile.CloseFile();
        return 1;
      }
      grid[ndata++] = (float)atof(infile.NextToken());
    }
  }
  // Finished
  infile.CloseFile();
  return 0;
}

// DataIO_OpenDx::WriteData3D()
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
