#include <cstring> // memset
#include "Grid.h"
#include "CpptrajStdio.h"

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
//frames_(0),
  grid_(0),
//dipolex_(0),
//dipoley_(0),
//dipolez_(0),
  callingRoutine_("GRID")
{}

// DESTRUCTOR
Grid::~Grid() {
  if (grid_!=0) delete[] grid_;
  //if (dipolex_!=0) delete[] dipolex_;
  //if (dipoley_!=0) delete[] dipoley_;
  //if (dipolez_!=0) delete[] dipolez_;
}

// Grid::GridInit()
/** Initialize grid from argument list. Expected call:
  * <CMD> <filename> <nx> <dx> <ny> <dy> <nz> <dz> [box|origin] [negative]  
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
  // Mask expression
/*  char* maskexpr = argIn.getNextMask();
  if (maskexpr==NULL) {
    mprinterr("Error: %s: No mask specified.\n",callingRoutine_.c_str());
    return 1;
  }
  gridmask_.SetMaskString(maskexpr);*/
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
  //mprintf(" will be dumped to file %s as", filename_.c_str());
  mprintf(" will be calculated as");
  if (increment_==1)
    mprintf("positive density\n");
  else
    mprintf("negative density\n");
  //mprintf("\tMask expression: [%s]\n",gridmask_.MaskString());
  mprintf("\tGrid points : %5i %5i %5i\n", nx_, ny_, nz_);
  mprintf("\tGrid spacing: %5.3lf %5.3lf %5.3lf\n", dx_, dy_, dz_);
}

// Grid::GridAllocate()
int Grid::GridAllocate() {
  //gridsize_ =  (size_t)nx_;
  //gridsize_ *= (size_t)ny_;
  //gridsize_ *= (size_t)nz_;
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
  // Setup mask
/*  if (currentParm->SetupIntegerMask( gridmask_ ))
    return 1;
  mprintf("\t[%s] %i atoms selected.\n", gridmask_.MaskString(), gridmask_.Nselected());
  if (gridmask_.None()) {
    mprinterr("Error: %s: No atoms selected for parm %s\n", callingRoutine_.c_str(), 
              currentParm->c_str());
    return 1;
  }*/

  return 0;
}

// Grid::GridPrint()
// invert: iarg3
// smooth: darg3
// madura: darg2
//void Grid::GridPrint(bool invert, double smooth, double madura) 
void Grid::GridPrintHeader(CpptrajFile& outfile)
{
  //CpptrajFile outfile;
  //if (outfile.OpenWrite( filename_ )) return;
  outfile.Printf("This line is ignored\n%8i\nrdparm generated grid density\n", 1);
  outfile.Printf("%8i%8i%8i", nx_, -nx_/2 + 1, nx_/2 );
  outfile.Printf("%8i%8i%8i", ny_, -ny_/2 + 1, ny_/2 );
  outfile.Printf("%8i%8i%8i", nz_, -nz_/2 + 1, nz_/2 );
  outfile.Printf("%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n",
                 (double)nx_ * dx_, (double)ny_ * dy_, (double)nz_ * dz_,
                 90.0, 90.0, 90.0);
  outfile.Printf("ZYX\n");

/*  double gridMax = 0;
  for (int k = 0; k < nz_; ++k) {
    outfile.Printf("%8i\n", k - nz_/2 + 1);
    for (int j = 0; j < ny_; ++j) {
      int col = 1;
      for (int i = 0; i < nx_; ++i) {
        int idx = i*ny_*nz_ + j*nz_ + k;
        // ----- SMOOTHING -----
        if (smooth > 0.0) {
          double yy = grid_[idx] - smooth;
          double xx = yy*yy / (0.2 * smooth * smooth);
          xx = exp( -xx );
          if (invert) {
            if (grid_[idx] > smooth) // NOTE: Comparison OK? Needs cast?
              grid_[idx] = -5.0;
            else
              grid_[idx] -= grid_[idx] * xx;
            // COMMENTED OUT IN ORIGINAL PTRAJ CODE
            //if (gridInfo->grid[index] < action->darg3) {
            //  gridInfo->grid[index] = 0.0;
            //}
            //
            if (grid_[idx] >= 0)
              grid_[idx] = smooth - grid_[idx];
          } else {
            if (grid_[idx] < smooth)
              grid_[idx] = 0;
            else
              grid_[idx] -= grid_[idx] * xx;
            if (grid_[idx] < smooth)
              grid_[idx] = 0;
          }
        }

        // do the madura negative option to expose low density
        if ( madura > 0.0 && grid_[idx] > 0.0 && grid_[idx] < madura )
          outfile.Printf("%12.5f", -grid_[idx]);
        else
          outfile.Printf("%12.5f", grid_[idx]);

        if (col && (col%6 == 0))
          outfile.Printf("\n");
        ++col;

        if ( grid_[idx] > gridMax )
          gridMax = grid_[idx];
      } // END i loop over x
      if ( (col-1) % 6 != 0 ) // Unless a newline was just written...
        outfile.Printf("\n");
    } // END j loop over y
  } // END k loop over z
*/
}
