#include "GridAction.h"
#include "CpptrajStdio.h"

// GridAction::HelpText
const char* GridAction::HelpText =
  "{nx dx ny dy nz dz [box|origin|center <mask>] [negative]}";

// CheckEven()
static inline void CheckEven(int& N, std::string const& call) {
  if (N % 2 == 1) {
    mprintf("Warning: %s: number of grid points must be even.\n", call.c_str());
    ++N;
    mprintf("Warning: Incrementing NX by 1 to %u\n", N);
  }
}

// GridAction::GridInit()
DataSet_GridFlt* GridAction::GridInit(const char* callingRoutineIn, ArgList& argIn, 
                                      DataSetList& DSL) 
{
  DataSet_GridFlt* Grid = 0;
  std::string callingRoutine;
  if (callingRoutineIn!=0)
    callingRoutine.assign(callingRoutineIn);
  else
    callingRoutine.assign("Grid");
  // Check for existing dataset
  std::string dsname = argIn.GetStringKey("data");
  if (!dsname.empty()) {
    Grid = (DataSet_GridFlt*)DSL.FindSetOfType( dsname, DataSet::GRID_FLT );
    if (Grid == 0) {
      mprinterr("Error: %s: Could not find grid data set with name %s\n",
                callingRoutine.c_str(), dsname.c_str());
      return 0;
    }
    return Grid;
  }
  // Get nx, dx, ny, dy, nz, dz
  int nx = argIn.getNextInteger(-1);
  double dx = argIn.getNextDouble(-1);
  int ny = argIn.getNextInteger(-1);
  double dy = argIn.getNextDouble(-1);
  int nz = argIn.getNextInteger(-1);
  double dz = argIn.getNextDouble(-1);
  if (nx < 1 || ny < 1 || nz < 1 ||
      dx < 0 || dy < 0 || dz < 0)
  {
    mprinterr("Error: %s: Invalid grid size/spacing.\n", callingRoutine.c_str());
    mprinterr("       nx=%i ny=%i nz=%i | dx=%.3f dy=%.3f dz=%.3f\n",
              nx, ny, nz, dx, dy, dz);
    return 0;
  }
  // Check that # of grid points are even in each direction.
  CheckEven( nx, callingRoutine );
  CheckEven( ny, callingRoutine );
  CheckEven( nz, callingRoutine );
  // Set spacing and origin
  // FIXME: For now only allow actual origin to be consistent with previous grid.
  Grid->setOriginAndSpacing(nx,ny,nz,0.0,0.0,0.0,dx, dy, dz);
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
      return 0;
    }
    centerMask_.SetMaskString( maskexpr );
    mode_ = CENTER;
  }
  // Negative
  if (argIn.hasKey("negative"))
    increment_ = -1.0;
  // Allocate DataSet
  Grid = (DataSet_GridFlt*)DSL.AddSet( DataSet::GRID_FLT, argIn.GetStringKey("name"), "GRID" );
  if (Grid->Allocate3D(nx, ny, nz)) {
    DataSetList::const_iterator last = DSL.end();
    --last;
    DSL.erase( last );
    Grid = 0;
  }
  return Grid;
}
  
void GridAction::GridInfo(DataSet_GridFlt const& grid) {
  mprintf("\tGrid at");
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
  mprintf("\tGrid points : %5i %5i %5i\n", grid.NX(), grid.NY(), grid.NZ());
  mprintf("\tGrid spacing: %5.3f %5.3f %5.3f\n", grid.DX(), grid.DY(), grid.DZ());
}

int GridAction::GridSetup(Topology const& currentParm) {
  // Check box
  if (mode_ == BOX) {
    if (currentParm.BoxType()!=Box::ORTHO) {
      mprintf("Warning: Code to shift to the box center is not yet\n");
      mprintf("Warning: implemented for non-orthorhomibic unit cells.\n");
      mprintf("Warning: Shifting to the origin instead.\n");
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

