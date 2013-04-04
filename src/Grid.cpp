#include <cmath>   // floor
#include <cstring> // memset
#include <cstdio>  // sscanf
#include <cstdlib> // atof
#include <vector>
#include "Constants.h"
#include "CpptrajStdio.h"
#include "Grid.h"
#include "BufferedLine.h"
#include "PDBfile.h"

#define MIN(X, Y) ( ( (X) < (Y) ) ? (X) : (Y) )
#define MAX(X, Y) ( ( (X) < (Y) ) ? (Y) : (X) )

// CONSTRUCTOR
Grid::Grid() :
  increment_(1.0),
  mode_(ORIGIN),
  dx_(0),
  dy_(0),
  dz_(0),
  sx_(0),
  sy_(0),
  sz_(0),
  gridsize_(0),
  nx_(0),
  ny_(0),
  nz_(0),
  grid_(0),
  callingRoutine_("GRID")
{}

// DESTRUCTOR
Grid::~Grid() {
  if (grid_!=0) delete[] grid_;
}

// Copy Constructor
Grid::Grid(const Grid& rhs) :
  increment_(rhs.increment_),
  mode_(rhs.mode_),
  dx_(rhs.dx_),
  dy_(rhs.dy_),
  dz_(rhs.dz_),
  sx_(rhs.sx_),
  sy_(rhs.sy_),
  sz_(rhs.sz_),
  gridsize_(rhs.gridsize_),
  nx_(rhs.nx_),
  ny_(rhs.ny_),
  nz_(rhs.nz_),
  grid_(0),
  callingRoutine_(rhs.callingRoutine_)
{
  if (gridsize_ > 0)
    grid_ = new float[ gridsize_ ];
  memset(grid_, 0, gridsize_ * sizeof(float));
}

// ASSIGNMENT
Grid& Grid::operator=(const Grid& rhs) {
  if (this == &rhs) return *this;
  // Deallocate
  if (grid_!=0) delete[] grid_;
  grid_ = 0;
  increment_ = rhs.increment_;
  mode_ = rhs.mode_;
  dx_ = rhs.dx_;
  dy_ = rhs.dy_;
  dz_ = rhs.dz_;
  sx_ = rhs.sx_;
  sy_ = rhs.sy_;
  sz_ = rhs.sz_;
  gridsize_ = rhs.gridsize_;
  nx_ = rhs.nx_;
  ny_ = rhs.ny_;
  nz_ = rhs.nz_;
  callingRoutine_ = rhs.callingRoutine_;
  if (gridsize_ > 0) {
    grid_ = new float[ gridsize_ ];
    memcpy( grid_, rhs.grid_, gridsize_ * sizeof(float) );
  }
  return *this;
}

const char* Grid::HelpText = 
  "{nx dx ny dy nz dz [box|origin|center <mask>] [negative]} | {readdx <file>}";
 
// Grid::GridInit()
/** Initialize grid from argument list. Expected args are:
  * <nx> <dx> <ny> <dy> <nz> <dz> [box|origin|center <mask>] [negative]
  * Where n is the number of bin points and d is the bin spacing in
  * a given direction, box/origin/center <mask> specifies where the
  * grid should be centered, and [negative] builds negative instead of
  * positive density.
  * \return 1 on error, 0 on success.
  */
int Grid::GridInit(const char* callingRoutineIn, ArgList& argIn) {
  if (callingRoutineIn!=0)
    callingRoutine_.assign(callingRoutineIn);
  std::string dxfilename = argIn.GetStringKey("readdx");
  if (!dxfilename.empty()) 
    return InitFromFile( dxfilename, "DX" );
  // Get nx, dx, ny, dy, nz, dz
  nx_ = argIn.getNextInteger(-1);
  dx_ = argIn.getNextDouble(-1);
  ny_ = argIn.getNextInteger(-1);
  dy_ = argIn.getNextDouble(-1);
  nz_ = argIn.getNextInteger(-1);
  dz_ = argIn.getNextDouble(-1);
  if (nx_ < 0 || ny_ < 0 || nz_ < 0 ||
      dx_ < 0 || dy_ < 0 || dz_ < 0)
  {
    mprinterr("Error: %s: Invalid grid size/spacing.\n", callingRoutine_.c_str());
    mprinterr("       nx=%i ny=%i nz=%i | dx=%.3f dy=%.3f dz=%.3f\n",
              nx_, ny_, nz_, dx_, dy_, dz_);
    return 1;
  }
  if (nx_ % 2 == 1) {
    mprintf("Warning: %s: number of grid points must be even.\n",callingRoutine_.c_str());
    ++nx_;
    mprintf("         Incrementing NX by 1 to %i\n", nx_);
  }
  if (ny_ % 2 == 1) {
    mprintf("Warning: %s: number of grid points must be even.\n",callingRoutine_.c_str());
    ++ny_;
    mprintf("         Incrementing NY by 1 to %i\n", ny_);
  }
  if (nz_ % 2 == 1) {
    mprintf("Warning: %s: number of grid points must be even.\n",callingRoutine_.c_str());
    ++nz_;
    mprintf("         Incrementing NY by 1 to %i\n", nz_);
  }
  // Box/origin
  mode_ = ORIGIN;
  if (argIn.hasKey("box"))
    mode_ = BOX;
  else if (argIn.hasKey("origin"))
    mode_ = ORIGIN;
  else if (argIn.Contains("center")) {
    std::string maskexpr = argIn.GetStringKey("center");
    if (maskexpr.empty()) {
      mprinterr("Error: 'center' requires <mask>\n");
      return 1;
    }
    centerMask_.SetMaskString( maskexpr );
    mode_ = CENTER;
  }
  // Negative
  if (argIn.hasKey("negative"))
    increment_ = -1.0;
  // Allocate memory and calc half grid
  if (Allocate()) return 1;

  return 0;
}

// Grid::GridInitSizeRes()
/** Initialize grid from a given size and resolution */
int Grid::GridInitSizeRes(const char* callingRoutineIn, double size[3],
                          double res[3], std::string const& mode) {
  if (callingRoutineIn!=0)
    callingRoutine_.assign(callingRoutineIn);
  // Determine the mode
  if (mode == "center")
    mode_ = CENTER;
  else if (mode == "box")
    mode_ = BOX;
  else
    mode_ = ORIGIN;
  // Set our grid spacing
  dx_ = res[0]; dy_ = res[1]; dz_ = res[2];
  // Determine our bin count, and make sure it's even
  nx_ = (int) (size[0] / dx_);
  ny_ = (int) (size[1] / dy_);
  nz_ = (int) (size[2] / dz_);
  if (nx_ % 2 == 1) nx_++;
  if (ny_ % 2 == 1) ny_++;
  if (nz_ % 2 == 1) nz_++;
  // Now allocate the grid
  if (Allocate())
    return 1;
  return 0;
}

/** Allocate the grid and calculate the half grid. */
int Grid::Allocate() {
  // Calculate half grid
  sx_ = (double)nx_ * dx_/2.0;
  sy_ = (double)ny_ * dy_/2.0;
  sz_ = (double)nz_ * dz_/2.0;
  // Allocate memory
  gridsize_ = nx_ * ny_ * nz_;
  if (gridsize_ <= 0) {
    mprinterr("Error: %s: Grid size <= 0 (%i)\n",callingRoutine_.c_str(), gridsize_);
    return 1;
  }
  if (grid_!=0) delete[] grid_;
  grid_ = new float[ gridsize_ ];
  memset(grid_, 0, gridsize_ * sizeof(float));
  return 0;
}

// Grid::InitFromFile()
/** Initializes a grid from a density file instead of input arguments. The
 *  filetype must be a recognized file type (case-sensitive). So far only
 *  the following file types are recognized for input densities:
 *      DX
 */
int Grid::InitFromFile(std::string const& filename, std::string const& filetype)
{
  if (filetype == "DX") {
    BufferedLine infile;
    infile.OpenRead(filename);
    // Set some defaults
    mode_ = ORIGIN;
    increment_ = 1.0;
    // Skip comments
    std::string line = infile.GetLine();
    while (!line.empty() && line[0] == '#') line = infile.GetLine();
    if (line.empty()) {
      mprinterr("Error: Unexpected EOF in DX file %s\n", filename.c_str());
      return 1;
    }
    // object 1 class gridpositions counts nx ny nz
    if (sscanf(line.c_str(), "object 1 class gridpositions counts %d %d %d",
               &nx_, &ny_, &nz_) != 3)
    {
      mprinterr("Error: Reading grid counts from DX file %s\n", filename.c_str());
      return 1;
    }
    // origin xmin ymin zmin (unused), make sure it says origin
    line = infile.GetLine();
    if (line.compare(0, 6, "origin") != 0) 
      mprintf("Warning: DX file %s, expected 'origin ...', got [%s]\n",
              filename.c_str(), line.c_str());
    // 3x 'delta hx hy hz'
    double dxyz[3];
    for (int i = 0; i < 3; i++) {
      line = infile.GetLine();
      if (sscanf(line.c_str(), "delta %lg %lg %lg", dxyz, dxyz+1, dxyz+2) != 3) {
        mprinterr("Error: Reading delta line from DX file %s\n", filename.c_str());
        return 1;
      }
      // Check that only 1 of the 3 values is non-zero
      if (dxyz[i] != (dxyz[0] + dxyz[1] + dxyz[2])) {
        mprinterr("Error: Rotated basis in DX file %s not yet supported.\n", filename.c_str());
        return 1;
      }
      switch (i) {
        case 0: dx_ = dxyz[0]; break; 
        case 1: dy_ = dxyz[1]; break; 
        case 2: dz_ = dxyz[2]; break;
      }
    } 
    // object 2 class gridconnections counts nx ny nz
    int nxyz[3];
    line = infile.GetLine();
    if (sscanf(line.c_str(), "object 2 class gridconnections counts %d %d %d",
               nxyz, nxyz+1, nxyz+2) != 3)
    {
      mprintf("Error: Reading grid connections from DX file %s\n", filename.c_str());
      return 1; 
    }
    // Sanity check for conflicting grid dimensions
    if (nxyz[0] != nx_ || nxyz[1] != ny_ || nxyz[2] != nz_) {
      mprinterr("Error: Conflicting grid dimensions in input DX density file %s.\n",
                filename.c_str());
      mprinterr("Error: Grid positions: %d %d %d\n", nx_, ny_, nz_);
      mprinterr("Error: Grid connections: %d %d %d\n", nxyz[0], nxyz[1], nxyz[2]);
      return 1;
    }
    // object 3 class array type <type> rank <r> times <i>
    // This line describes whether data will be in binary or ascii format.
    line = infile.GetLine();
    if (line.compare(0, 8, "object 3") != 0) {
      mprinterr("Error: DX file %s; expected 'object 3 ...', got [%s]\n",
                filename.c_str(), line.c_str());
      return 1;
    }
    if (line.find("binary") != std::string::npos) {
      mprinterr("Error: DX file %s; binary DX files not yet supported.\n", filename.c_str());
      return 1;
    }
    // Allocate the grid and calc half grid
    if (Allocate()) return 1;
    // Read in data
    int ndata = 0;
    infile.SetupBuffer();
    while (ndata < gridsize_) {
      if (infile.Line() == 0) {
        mprinterr("Error: Unexpected EOF hit in %s\n", filename.c_str());
        infile.CloseFile();
        return 1;
      }
      int nTokens = infile.TokenizeLine(" \t");
      for (int j = 0; j < nTokens; j++) {
        if (ndata >= gridsize_) {
          mprinterr("Error: Too many grid points found!\n");
          infile.CloseFile();
          return 1;
        }
        grid_[ndata++] = (float)atof(infile.NextToken());
      }
    }
    // Finished
    infile.CloseFile();
  } else {
    mprinterr("Input density file type [%s] is not recognized\n",
              filetype.c_str());
    return 1;
  }
  return 0;
}

// Grid::GridInfo()
void Grid::GridInfo() {
  mprintf("    %s: Grid at", callingRoutine_.c_str());
  if (mode_ == BOX)
    mprintf(" box center");
  else if (mode_ == ORIGIN)
    mprintf(" origin");
  else if (mode_ == CENTER)
    mprintf(" center of atoms in mask [%s]\n", centerMask_.MaskString());
  mprintf(" will be calculated as");
  if (increment_ > 0)
    mprintf(" positive density\n");
  else
    mprintf(" negative density\n");
  mprintf("\tGrid points : %5i %5i %5i\n", nx_, ny_, nz_);
  mprintf("\tGrid spacing: %5.3f %5.3f %5.3f\n", dx_, dy_, dz_);
}

// Grid::GridSetup()
int Grid::GridSetup(Topology const& currentParm) {
  // Check box
  if (mode_ == BOX) {
    if (currentParm.BoxType()!=Box::ORTHO) {
      mprintf("Warning: %s: Code to shift to the box center is not yet\n",callingRoutine_.c_str());
      mprintf("Warning:\timplemented for non-orthorhomibic unit cells.\n");
      mprintf("Warning:\tShifting to the origin instead.\n");
      mode_ = ORIGIN;
    }
  } else if (mode_ == CENTER) {
    if ( currentParm.SetupIntegerMask( centerMask_ ) ) return 1;
    centerMask_.MaskInfo();
    if ( centerMask_.None() ) {
      mprinterr("Error: No atoms selected for grid center mask [%s]\n", centerMask_.MaskString());
      return 1;
    }
  }

  return 0;
}

// Grid::GridPoint()
/** Main grid routine. Check if position specified by coordinates
  * corresponds to a valid bin and if so increment the bin.
  * \return Bin index if binned, -1 if out of bounds.
  */
int Grid::GridPoint(Vec3 const& vIn) {
  double xx = vIn[0] + sx_;
  int i = (int) (xx / dx_) - 1;
  if (i >= 0 && i < nx_) {
    double yy = vIn[1] + sy_;
    int j = (int) (yy / dy_) - 1;
    if (j >= 0 && j < ny_) {
      double zz = vIn[2] + sz_;
      int k = (int) (zz / dz_) - 1;
      if (k >= 0 && k < nz_) {
        int idx = i*ny_*nz_ + j*nz_ + k;
        grid_[idx] += increment_;
        return idx;
      }
    }
  }
  return -1;
}

// Grid::GridFrame()
/** Grid coordinates in currentFrame according to mask. */
void Grid::GridFrame(Frame& currentFrame, AtomMask const& mask) {
  Vec3 center;
  if (mode_ == BOX) 
    center = currentFrame.BoxCrd().Center();
  else if (mode_ == CENTER)
    center = currentFrame.VGeometricCenter( centerMask_ );
  else // ORIGIN
    center.Zero();
  for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom)
    GridPoint( Vec3(currentFrame.XYZ( *atom )) - center );
}

// Grid::BinPoint
/** Grid routine used for backwards compatibility with ptraj Dipole
  * routine (Action_Dipole).
  */
int Grid::BinPoint(double xx, double yy, double zz) {
  int i = (int) (xx / dx_);
  if (i > 0 && i < nx_) {
    int j = (int) (yy / dy_);
    if (j > 0 && j < ny_) {
      int k = (int) (zz / dz_);
      if (k > 0 && k < nz_) {
        int idx = i*ny_*nz_ + j*nz_ + k;
        grid_[idx] += increment_;
        return idx;
      }
    }
  }
  return -1;
}

// Grid::PrintXplor()
void Grid::PrintXplor(std::string const& name, const char* title, 
                      std::string remark)
{
  CpptrajFile outfile;
  if (outfile.OpenWrite( name )) {
    mprinterr("Error: Could not open Xplor output file.\n");
    return;
  }
  // Title
  outfile.Printf("%s\n", title);
  // Remarks - only 1
  outfile.Printf("%8i\n%s\n",1,remark.c_str());
  // Header
  outfile.Printf("%8i%8i%8i",   nx_, -nx_/2 + 1, nx_/2 );
  outfile.Printf("%8i%8i%8i",   ny_, -ny_/2 + 1, ny_/2 );
  outfile.Printf("%8i%8i%8i\n", nz_, -nz_/2 + 1, nz_/2 );
  outfile.Printf("%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n",
                 (double)nx_ * dx_, (double)ny_ * dy_, (double)nz_ * dz_,
                 90.0, 90.0, 90.0);
  outfile.Printf("ZYX\n");
  // Print grid bins
  int NZ2 = nz_ / 2;
  for (int k = 0; k < nz_; ++k) {
    outfile.Printf("%8i\n", k - NZ2 + 1);
    for (int j = 0; j < ny_; ++j) {
      int col = 1;
      for (int i = 0; i < nx_; ++i) {
        outfile.Printf("%12.5f", GridVal(i, j, k));
        if ( (col % 6)==0 )
          outfile.Printf("\n");
        ++col;
      }
      if ( (col-1) % 6 != 0 )
        outfile.Printf("\n");
    }
  }
  outfile.CloseFile();
}

// Grid::PrintPDB()
void Grid::PrintPDB(std::string const& filename, double cut, double normIn) 
{
  double norm = normIn;
  // Calculate normalization if necessary
  if (norm <= 0) {
    for (int i = 0; i < gridsize_; ++i)
       if ((double)grid_[i] > norm)
         norm = (double)grid_[i];
    if (norm == 0) {
      mprinterr("Error: %s: Grid max is 0. No density for PDB write.\n",
                callingRoutine_.c_str());
      return;
    }
    mprintf("\t%s: Normalizing grid by %f\n", norm);
  }
  norm = 1.0 / norm;
  // Write PDB
  PDBfile pdbout;
  if (pdbout.OpenWrite(filename)) {
    mprinterr("Error: %s: Cannot open PDB output.\n", callingRoutine_.c_str()); 
    return;
  }
  mprintf("\tWriting PDB of grid points > %.3f of grid max.\n", cut);
  int res = 1;
  for (int k = 0; k < nz_; ++k) {
    for (int j = 0; j < ny_; ++j) {
      for (int i = 0; i < nx_; ++i) {
        double gridval = GridVal(i, j, k) * norm;
        if (gridval > cut)
          pdbout.WriteATOM(res++, Xcrd(i), Ycrd(j), Zcrd(k), "GRID", gridval);
      }
    }
  }
  // Write grid boundaries
  for (int k = 0; k <= nz_; k += nz_)
    for (int j = 0; j <= ny_; j += ny_)
      for (int i = 0; i <= nx_; i += nx_)
        pdbout.WriteHET(res, Xcrd(i), Ycrd(j), Zcrd(k));
}

// Grid::PrintDX
/** Use the default lower-left corner (first bin) as the origin when printing this
  * DX file
  */
void Grid::PrintDX(std::string const& filename) {
  PrintDX(filename, Xbin(0), Ybin(0), Zbin(0));
}

// Grid::PrintDX
/** This will print the grid in OpenDX format, commonly used by VMD, PBSA,
 * 3D-RISM, etc.
 */
void Grid::PrintDX(std::string const& filename, double xorig, double yorig, double zorig)
{
  CpptrajFile outfile;
  if (outfile.OpenWrite(filename)) {
    mprinterr("Error: Could not open OpenDX output file.\n");
    return;
  }
  // Print the OpenDX header
  outfile.Printf("object 1 class gridpositions counts %d %d %d\n",
                 nx_, ny_, nz_);
  outfile.Printf("origin %lg %lg %lg\n", xorig, yorig, zorig);
  outfile.Printf("delta %lg 0 0\n", dx_);
  outfile.Printf("delta 0 %lg 0\n", dy_);
  outfile.Printf("delta 0 0 %lg\n", dz_);
  outfile.Printf("object 2 class gridconnections counts %d %d %d\n",
                 nx_, ny_, nz_);
  outfile.Printf(
    "object 3 class array type double rank 0 items %d data follows\n",
    gridsize_);
  
  // Now print out the data. It is already in row-major form (z-axis changes
  // fastest), so no need to do any kind of data adjustment
  for (int i = 0; i < gridsize_ - 2; i += 3)
    outfile.Printf("%g %g %g\n", grid_[i], grid_[i+1], grid_[i+2]);
  // Print out any points we may have missed
  switch (gridsize_ % 3) {
    case 2: outfile.Printf("%g %g\n", grid_[gridsize_-2], grid_[gridsize_-1]); break;
    case 1: outfile.Printf("%g\n", grid_[gridsize_-1]); break;
  }
  
  // Print tail
  if (mode_ == CENTER)
    outfile.Printf("\nobject \"density (%s) [A^-3]\" class field\n",
                   centerMask_.MaskString());
  else
    outfile.Printf("\nobject \"density [A^-3]\" class field\n");

  outfile.CloseFile();
}

// Grid::GridPrint()
void Grid::PrintEntireGrid() {
  for (int i = 0; i < gridsize_; ++i)
    mprintf("\t%f\n", grid_[i]);
}

// Grid::ExtractPeaks()
/** Extract peaks from the current grid and return another Grid instance. This
  * works by taking every grid point and analyzing all grid points adjacent to
  * it (including diagonals). If any of those grid points have a higher value
  * (meaning there is a direction towards "increased" density) then that value
  * is _not_ a maximum. Any density peaks less than the minimum filter are
  * discarded
  * \return Grid instance with all non-peak grid points zeroed-out
  */
Grid& Grid::ExtractPeaks(double min_filter) {
   // Start out with a copy of myself and an integer vector for the background
  Grid& peaks = *this;                     // _nparray

  /* Zero out all points that are less than min_filter and any points whose
   * "adjacent" points in any dimension have a greater density (meaning it's not
   * a local max)
   */
  for (int i = 0; i < nx_; i++)
  for (int j = 0; j < ny_; j++)
  for (int k = 0; k < nz_; k++) {
    float val = GridVal(i, j, k);
    if (val < min_filter) {
      peaks.SetGridVal(i, j, k, 0.0f);
      continue;
    }
    for (int ii = MAX(0, i-1); ii <= MIN(nx_, i+1); ii++)
    for (int jj = MAX(0, j-1); jj <= MIN(ny_, j+1); jj++)
    for (int kk = MAX(0, k-1); kk <= MIN(nz_, k+1); kk++) {
      if (ii == i && jj == j && kk == k) continue;
      if (GridVal(ii, jj, kk) > val)
        peaks.SetGridVal(i, j, k, 0.0f);
    }
  }

  return peaks;
}
