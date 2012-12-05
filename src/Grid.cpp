#include <cstring> // memset
#include "Grid.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"
#include "PDBfile.h"

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

void Grid::Help() {
  mprintf(" nx dx ny dy nz dz [box|origin|center <mask>] [negative]");
}
 
// Grid::GridInit()
/** Initialize grid from argument list. */
int Grid::GridInit(const char* callingRoutineIn, ArgList& argIn) {
  if (callingRoutineIn!=NULL)
    callingRoutine_.assign(callingRoutineIn);

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
  // Calculate half grid
  sx_ = (double)nx_ * dx_/2.0;
  sy_ = (double)ny_ * dy_/2.0;
  sz_ = (double)nz_ * dz_/2.0;
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

  // Allocate memory
  gridsize_ = nx_ * ny_ * nz_;
  if (gridsize_ <= 0) {
    mprinterr("Error: %s: Grid size <= 0 (%i)\n",gridsize_, callingRoutine_.c_str());
    return 1;
  }
  if (grid_!=0) delete[] grid_;
  grid_ = new float[ gridsize_ ];
  memset(grid_, 0, gridsize_ * sizeof(float));

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
int Grid::GridSetup(Topology* currentParm) {
  // Check box
  if (mode_ == BOX) {
    if (currentParm->BoxType()!=Box::ORTHO) {
      mprintf("Warning: %s: Code to shift to the box center is not yet\n",callingRoutine_.c_str());
      mprintf("Warning:\timplemented for non-orthorhomibic unit cells.\n");
      mprintf("Warning:\tShifting to the origin instead.\n");
      mode_ = ORIGIN;
    }
  } else if (mode_ == CENTER) {
    if ( currentParm->SetupIntegerMask( centerMask_ ) ) return 1;
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

void Grid::GridFrame(Frame& currentFrame, AtomMask const& mask) {
  Vec3 center;
  switch (mode_) {
    case BOX: center = currentFrame.BoxCrd().Center(); break;
    case CENTER: center = currentFrame.VGeometricCenter( centerMask_ ); break;
    case ORIGIN: center.Zero(); break;
  }
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

// Grid::GridPrint()
void Grid::PrintEntireGrid() {
  for (int i = 0; i < gridsize_; ++i)
    mprintf("\t%f\n", grid_[i]);
}

