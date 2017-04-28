#include <cstdio> // sscanf
#include <cstdlib> // atof
#include "DataIO_OpenDx.h"
#include "CpptrajStdio.h"
#include "DataSet_GridFlt.h"
#include "BufferedLine.h"
#include "ProgressBar.h"

bool DataIO_OpenDx::ID_DataFormat( CpptrajFile& infile ) {
  bool isDX = false;
  if (!infile.OpenFile()) {
    std::string firstLine = infile.GetLine();
    if (!firstLine.empty())
      isDX = (firstLine.compare(0, 28, "object 1 class gridpositions") == 0);
    infile.CloseFile();
  }
  return isDX;
}

// DataIO_OpenDx::ReadData()
int DataIO_OpenDx::ReadData(FileName const& fname, 
                            DataSetList& datasetlist, std::string const& dsname)
{
  // Add grid data set. Default to float for now.
  DataSet* ds = datasetlist.AddSet( DataSet::GRID_FLT, dsname, "GRID" );
  if (ds==0) return 1;
  if (LoadGrid(fname.full(), *ds)) {
    // Load failed. Erase grid data set.
    datasetlist.RemoveSet( ds );
    return 1;
  }
  return 0;
}

// DataIO_OpenDx::LoadGrid()
int DataIO_OpenDx::LoadGrid(const char* filename, DataSet& ds)
{
  // TODO: This may need to be changed if new 3D types introduced.
  DataSet_GridFlt& grid = static_cast<DataSet_GridFlt&>( ds );
  // Open file
  BufferedLine infile;
  if (infile.OpenFileRead(filename)) return 1;
  // Skip comments
  std::string line = infile.GetLine();
  while (!line.empty() && line[0] == '#') {
    mprintf("\t%s", line.c_str());
    line = infile.GetLine();
  }
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
  double oxyz[3];
  line = infile.GetLine();
  if (sscanf(line.c_str(), "origin %lg %lg %lg", oxyz, oxyz+1, oxyz+2) != 3) {
    mprinterr("Error: Reading origin line from DX file %s\n", filename);
    return 1;
  }
  // 3x 'delta hx hy hz'
  double dxyz[3];
  Matrix_3x3 delta(0.0);
  bool isNonortho = false;
  int midx = 0;
  for (int i = 0; i < 3; i++, midx += 3) {
    line = infile.GetLine();
    if (sscanf(line.c_str(), "delta %lg %lg %lg", dxyz, dxyz+1, dxyz+2) != 3) {
      mprinterr("Error: Reading delta line from DX file %s\n", filename);
      return 1;
    }
    // Check that only 1 of the 3 values is non-zero. Otherwise non-ortho.
    if (dxyz[i] != (dxyz[0] + dxyz[1] + dxyz[2]))
      isNonortho = true;
    delta[midx  ] = dxyz[0];
    delta[midx+1] = dxyz[1];
    delta[midx+2] = dxyz[2];
  }
  // object 2 class gridconnections counts nx ny nz
  int nxyz[3];
  line = infile.GetLine();
  if (sscanf(line.c_str(), "object 2 class gridconnections counts %d %d %d",
             nxyz, nxyz+1, nxyz+2) != 3)
  {
    mprinterr("Error: Reading grid connections from DX file %s\n", filename);
    return 1;
  }
  // Sanity check for conflicting grid dimensions
  if (nxyz[0] != nx || nxyz[1] != ny || nxyz[2] != nz) {
    mprinterr("Error: Conflicting grid dimensions in input DX density file %s.\n",
              filename);
    mprinterr("Error: Grid positions: %d %d %d\n", nx, ny, nz);
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
  int err = 0;
  if (isNonortho) {
    // Create unit cell from delta and bins.
    delta[0] *= (double)nx; delta[1] *= (double)nx; delta[2] *= (double)nx;
    delta[3] *= (double)ny; delta[4] *= (double)ny; delta[5] *= (double)ny;
    delta[6] *= (double)nz; delta[7] *= (double)nz; delta[8] *= (double)nz;
    err = grid.Allocate_N_O_Box(nx,ny,nz, Vec3(oxyz), Box(delta));
  } else
    err = grid.Allocate_N_O_D(nx,ny,nz, Vec3(oxyz), Vec3(delta[0],delta[4],delta[8]));
  if (err != 0) { 
    mprinterr("Error: Could not allocate grid.\n");
    return 1;
  }
  grid.GridInfo();
  // Read in data
  size_t gridsize = grid.Size();
  mprintf("\tReading in %zu data elements from DX file.\n", gridsize); 
  size_t ndata = 0;
  ProgressBar progress( gridsize );
  while (ndata < gridsize) {
    if (infile.Line() == 0) {
      mprinterr("Error: Unexpected EOF hit in %s\n", filename);
      return 1;
    }
    int nTokens = infile.TokenizeLine(" \t");
    for (int j = 0; j < nTokens; j++) {
      if (ndata >= gridsize) {
        mprintf("Warning: Too many grid points found. Only reading %zu grid points.\n", gridsize);
        mprintf("Warning: Check that data region ends with a newline.\n");
        break;
      }
      grid[ndata++] = (float)atof(infile.NextToken());
    }
    progress.Update( ndata );
  }
  // Set dimensions
  // FIXME: This should be integrated with allocation
  //grid.SetDim(Dimension::X, Dimension(oxyz[0], dx, nx, "X"));
  //grid.SetDim(Dimension::Y, Dimension(oxyz[1], dy, ny, "Y"));
  //grid.SetDim(Dimension::Z, Dimension(oxyz[2], dz, nz, "Z"));
  return 0;
}

// -----------------------------------------------------------------------------
void DataIO_OpenDx::WriteHelp() {
  mprintf("\tbincenter: Center grid points on bin centers.\n"
          "\tgridwrap:  Like 'bincenter', but also wrap grid density.\n"
          "\t           Useful when grid encompasses unit cell.\n"
          "\tgridext:   Like 'bincenter', but also print extra layer of empty bins.\n");
}

int DataIO_OpenDx::processWriteArgs(ArgList& argIn) {
  if (argIn.hasKey("bincenter")) gridWriteMode_ = BIN_CENTER;
  else if (argIn.hasKey("gridwrap")) gridWriteMode_ = WRAP;
  else if (argIn.hasKey("gridext")) gridWriteMode_ = EXTENDED;
  if (gridWriteMode_ == BIN_CORNER)
    mprintf("\tOpenDx: Grid will be created using bin corners.\n");
  else if (gridWriteMode_ == BIN_CENTER)
    mprintf("\tOpenDx: Grid will be created using bin centers.\n");
  else if (gridWriteMode_ == WRAP)
    mprintf("\tOpenDx: Grid will be created using bin centers and wrapped.\n");
  else if (gridWriteMode_ == EXTENDED)
    mprintf("\tOpenDx: Grid will be created using bin centers and surrounded with empty bins.\n");
  return 0;
}

// DataIO_OpenDx::WriteData()
int DataIO_OpenDx::WriteData(FileName const& fname, DataSetList const& setList)
{
  // Open output file
  CpptrajFile outfile;
  if (outfile.OpenWrite(fname)) {
    mprinterr("Error: Could not open OpenDX output file.\n");
    return 1;
  }
  // Warn about writing multiple sets
  if (setList.size() > 1)
    mprintf("Warning: %s: Writing multiple 3D sets in OpenDX format may result in unexpected behavior\n", fname.full());
  int err = 0;
  for (DataSetList::const_iterator set = setList.begin(); set != setList.end(); ++set)
    err += WriteSet3D( *(*set), outfile );
  return err;
}

// DataIO_OpenDx::WriteSet3D()
int DataIO_OpenDx::WriteSet3D(DataSet const& setIn, CpptrajFile& outfile) const {
  if (setIn.Ndim() != 3) {
    mprinterr("Internal Error: DataSet %s in DataFile %s has %zu dimensions, expected 3.\n",
              setIn.legend(), outfile.Filename().full(), setIn.Ndim());
    return 1;
  }
  int err = 0;
  switch ( gridWriteMode_ ) {
    case BIN_CORNER:
    case BIN_CENTER: err = WriteGrid( setIn, outfile ); break;
    case WRAP:
    case EXTENDED  : err = WriteGridWrap( setIn, outfile ); break;
  }
  // Print tail
  if (err == 0) {
    // TODO: Make this an option
    //if (mode_ == CENTER)
    //  outfile.Printf("\nobject \"density (%s) [A^-3]\" class field\n",
    //                 centerMask_.MaskString());
    //else
      outfile.Printf("\nobject \"density [A^-3]\" class field\n");
  }
  return err;
}

void DataIO_OpenDx::WriteDxHeader(CpptrajFile& outfile,
                                  size_t NX, size_t NY, size_t NZ,
                                  double LX, double LY, double LZ,
                                  Matrix_3x3 const& ucell, Vec3 const& oxyz) const
{
  outfile.Printf("object 1 class gridpositions counts %zu %zu %zu\n"
                 "origin %g %g %g\ndelta %g %g %g\ndelta %g %g %g\ndelta %g %g %g\n"
                 "object 2 class gridconnections counts %zu %zu %zu\n"
                 "object 3 class array type double rank 0 items %zu data follows\n",
                 NX, NY, NZ, oxyz[0], oxyz[1], oxyz[2],
                 ucell[0]/LX, ucell[1]/LX, ucell[2]/LX,
                 ucell[3]/LY, ucell[4]/LY, ucell[5]/LY,
                 ucell[6]/LZ, ucell[7]/LZ, ucell[8]/LZ,
                 NX, NY, NZ, NX*NY*NZ);
}

int DataIO_OpenDx::WriteGridWrap(DataSet const& setIn, CpptrajFile& outfile) const {
  DataSet_3D const& set = static_cast<DataSet_3D const&>( setIn );
  // Need to construct a grid mesh around bins, with points centered on the bins.
  int mesh_x = set.NX();
  int mesh_y = set.NY();
  int mesh_z = set.NZ();
  // Origin needs to be shifted half grid spacing, i.e. it is the center of the
  // bin located at -1, -1, -1.
  Vec3 oxyz = set.BinCenter(-1, -1, -1);
  // Print the OpenDX header
  WriteDxHeader(outfile, mesh_x+2, mesh_y+2, mesh_z+2, mesh_x, mesh_y, mesh_z,
                set.Ucell(), oxyz);
  // Print out the data. Start at bin -1, end on bin N.
  int nvals = 0; // Keep track of how many values printed on current line.
  if (gridWriteMode_ == WRAP) {
    int bi, bj, bk;
    for (int ii = -1; ii <= mesh_x; ++ii) {
      if      (ii < 0      ) bi = mesh_x - 1;
      else if (ii == mesh_x) bi = 0;
      else                   bi = ii;
      for (int ij = -1; ij <= mesh_y; ++ij) {
        if      (ij < 0      ) bj = mesh_y - 1;
        else if (ij == mesh_y) bj = 0;
        else                   bj = ij;
        for (int ik = -1; ik <= mesh_z; ++ik) {
          if      (ik < 0      ) bk = mesh_z - 1;
          else if (ik == mesh_z) bk = 0;
          else                   bk = ik;
          outfile.Printf(" %g", set.GetElement(bi, bj, bk));
          ++nvals;
          if (nvals == 5) {
            outfile.Printf("\n");
            nvals = 0;
          }
        }
      }
    }
  } else { // EXTENDED
    for (int ii = -1; ii <= mesh_x; ++ii) {
      bool zero_x = (ii < 0 || ii == mesh_x);
      for (int ij = -1; ij <= mesh_y; ++ij) {
        bool zero_y = (ij < 0 || ij == mesh_y);
        for (int ik = -1; ik <= mesh_z; ++ik) {
          if (zero_x || zero_y || ik < 0 || ik == mesh_z)
            outfile.Printf(" 0");
          else
            outfile.Printf(" %g", set.GetElement(ii, ij, ik));
          ++nvals;
          if (nvals == 5) {
            outfile.Printf("\n");
            nvals = 0;
          }
        }
      }
    }
  }
  if (nvals > 0) outfile.Printf("\n");
  return 0;
}

int DataIO_OpenDx::WriteGrid(DataSet const& setIn, CpptrajFile& outfile) const {
  DataSet_3D const& set = static_cast<DataSet_3D const&>( setIn );
  Vec3 oxyz = set.GridOrigin();
  if (gridWriteMode_ == BIN_CENTER)
    // Origin needs to be shifted to center of bin located at 0,0,0
    oxyz = set.BinCenter(0,0,0);
  // Print the OpenDX header
  WriteDxHeader(outfile, set.NX(), set.NY(), set.NZ(), set.NX(), set.NY(), set.NZ(),
                set.Ucell(), oxyz);
  // Now print out the data.
  size_t gridsize = set.Size();
  if (gridsize == 1)
    outfile.Printf("%g\n", set[0]);
  else if (gridsize == 2)
    outfile.Printf("%g %g\n", set[0], set[1]);
  else if (gridsize > 2) {
    // Data is already in row-major form (z-axis changes
    // fastest), so no need to do any kind of data adjustment
    for (size_t i = 0UL; i < gridsize - 2UL; i += 3UL)
      outfile.Printf("%g %g %g\n", set[i], set[i+1], set[i+2]);
    // Print out any points we may have missed
    switch (gridsize % 3) {
      case 2: outfile.Printf("%g %g\n", set[gridsize-2], set[gridsize-1]); break;
      case 1: outfile.Printf("%g\n", set[gridsize-1]); break;
    }
  }
  return 0;
}
