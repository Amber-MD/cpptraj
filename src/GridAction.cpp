#include "GridAction.h"
#include "CpptrajStdio.h"

const char* GridAction::HelpText =
  "{nx dx ny dy nz dz [box|origin|center <mask>] [negative]} | {readdx <file>}";

static inline void CheckEven(int& N, std::string const& call) {
  if (N % 2 == 1) {
    mprintf("Warning: %s: number of grid points must be even.\n", call.c_str());
    ++N;
    mprintf("Warning: Incrementing NX by 1 to %u\n", N);
  }
}

DataSet_3D* GridAction::GridInit(const char* callingRoutineIn, ArgList& argIn, DataSetList& DSL) {
  DataSet_3D* Grid = 0;
  std::string callingRoutine;
  if (callingRoutineIn!=0)
    callingRoutine.assign(callingRoutineIn);
  else
    callingRoutine.assign("Grid");
  // Check for existing dataset
  std::string dsname = argIn.GetStringKey("data");
  if (!dsname.empty()) {
    Grid = (DataSet_3D*)DSL.FindSetOfType( dsname, DataSet::GRID_FLT );
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
  // Allocate DataSet
  Grid = (DataSet_3D*)DSL.AddSet( DataSet::GRID_FLT, argIn.GetStringKey("name"), "GRID" );
  if (Grid->Allocate3D(nx, ny, nz)) {
    DataSetList::const_iterator last = DSL.end();
    --last;
    DSL.erase( last );
    Grid = 0;
  }
  return Grid;
}
  


