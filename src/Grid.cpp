#include <cstring> // memset
#include "Grid.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"

// CONSTRUCTOR
Grid::Grid() :
  increment_(1),
  box_(false),
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

// Grid::GridInit()
/** Initialize grid from argument list. Expected call:
  * <nx> <dx> <ny> <dy> <nz> <dz> [box|origin] [negative]  
  */
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
    mprinterr("       nx=%i ny=%i nz=%i | dx=%.3lf dy=%.3lf dz=%.3lf\n",
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
  if (argIn.hasKey("box"))
    box_ = true;
  else if (argIn.hasKey("origin"))
    box_ = false;
  // Negative
  if (argIn.hasKey("negative"))
    increment_ = -1;

  return 0;
}

// Grid::GridInfo()
void Grid::GridInfo() {
  mprintf("    %s: Grid at", callingRoutine_.c_str());
  if (box_)
    mprintf(" box center");
  else
    mprintf(" origin");
  mprintf(" will be calculated as");
  if (increment_==1)
    mprintf(" positive density\n");
  else
    mprintf(" negative density\n");
  mprintf("\tGrid points : %5i %5i %5i\n", nx_, ny_, nz_);
  mprintf("\tGrid spacing: %5.3lf %5.3lf %5.3lf\n", dx_, dy_, dz_);
}

// Grid::GridAllocate()
int Grid::GridAllocate() {
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

// Grid::GridSetup()
int Grid::GridSetup(Topology* currentParm) {
  // Check box
  if (box_) {
    if (currentParm->BoxType()!=Box::ORTHO) {
      mprintf("Warning: %s: Code to shift to the box center is not yet\n",callingRoutine_.c_str());
      mprintf("         implemented for non-orthorhomibic unit cells.\n");
      mprintf("         Shifting to the origin instead.\n");
      box_ = false;
    }
  }

  return 0;
}

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


// Grid::GridPrint()
void Grid::PrintEntireGrid() {
  for (int i = 0; i < gridsize_; ++i)
    mprintf("\t%lf\n", grid_[i]);
}

