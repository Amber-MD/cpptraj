#include <cmath> // log
#include "Action_GridFreeEnergy.h"
#include "CpptrajStdio.h" // mprintf, mprinterr
#include "StringRoutines.h" // integerToString

// CONSTRUCTOR
Action_GridFreeEnergy::Action_GridFreeEnergy() :
  maxVoxelOccupancyCount_(600) // NOTE: See header for comments.
{}

void Action_GridFreeEnergy::Help() {
  mprintf("\t<filename>");
  Grid::Help();
  mprintf(" <mask>\n");
}

// Action_GridFreeEnergy::init()
Action::RetType Action_GridFreeEnergy::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get output filename
  filename_ = actionArgs.GetStringNext();
  if (filename_.empty()) {
    mprinterr("Error: GridFreeEnergy: no filename specified.\n");
    return Action::ERR;
  }
  // Get grid options (<nx> <dx> <ny> <dy> <nz> <dz> [box|origin] [negative])
  if (grid_.GridInit( "GridFreeEnergy", actionArgs ))
    return Action::ERR;

  // Get mask
  std::string maskexpr = actionArgs.GetMaskNext();
  if (maskexpr.empty()) {
    mprinterr("Error: GRID: No mask specified.\n");
    return Action::ERR;
  }
  mask_.SetMaskString(maskexpr);

  // Info
  grid_.GridInfo();
  mprintf("\tGrid will be printed to file %s\n",filename_.c_str());
  mprintf("\tMask expression: [%s]\n",mask_.MaskString());

  // Allocate grid
  //if (GridAllocate()) return 1;

  return Action::OK;
}

// Action_GridFreeEnergy::setup()
Action::RetType Action_GridFreeEnergy::Setup(Topology* currentParm, Topology** parmAddress) {
  // Setup grid, checks box info.
  if (grid_.GridSetup( *currentParm )) return Action::ERR;

  // Setup mask
  if (currentParm->SetupIntegerMask( mask_ ))
    return Action::ERR;
  mprintf("\t[%s] %i atoms selected.\n", mask_.MaskString(), mask_.Nselected());
  if (mask_.None()) {
    mprinterr("Error: GridFreeEnergy: No atoms selected for parm %s\n", currentParm->c_str());
    return Action::ERR;
  }

  return Action::OK;
}

// Action_GridFreeEnergy::action()
Action::RetType Action_GridFreeEnergy::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  grid_.GridFrame( *currentFrame, mask_ );
  return Action::OK;
}

// Action_GridFreeEnergy::print()
void Action_GridFreeEnergy::Print() {
  CpptrajFile outfile;
  /* How times does this occupancy count value arise?
   *    i.e. if  
   *                 voxelOccupancyCount[50] = 10 
   * 
   *         then there are 10 voxels with an occupancy count 50
   */
  std::vector<int> voxelOccupancyCount;
  // Largest occupancy count
  int currentLargestVoxelOccupancyCount;
  /* Value of most frequent voxel occupancy count
   *
   * i.e. if 
   *                 voxelOccupancyCount[50] = 10 
   *                 voxelOccupancyCount[51] = 100
   *                 voxelOccupancyCount[52] = 5
   *      then
   *                 mostFrequentVoxelOccupancy would be 51
   */
  int mostFrequentVoxelOccupancy;

  // Zero the voxelOccupancyCount array
  voxelOccupancyCount.assign(maxVoxelOccupancyCount_, 0);
  // Determine frequency for bin populations
  for (Grid::iterator gval = grid_.begin(); gval != grid_.end(); ++gval) {
        int bincount = (int) *gval;
        // Sanity check: if bin count >= size of voxelOccupancyCount, 
        //               increase size.
        if (bincount >= (int)voxelOccupancyCount.size())
          voxelOccupancyCount.resize( bincount+1, 0 );
        // Add one to the counter for the fequency of this bin population
        ++voxelOccupancyCount[ bincount ];
  }
  // Walk the voxelOccupancyCount[] array and determine which
  // occupancy is the most frequent (excluding 0).
  currentLargestVoxelOccupancyCount = voxelOccupancyCount[1];
  mostFrequentVoxelOccupancy = 1;
  //mprintf("CDBG: voxelOccupancyCount[1] = %i\n", voxelOccupancyCount[1]);
  for (int i = 2; i < (int)voxelOccupancyCount.size(); ++i) {
    //mprintf("CDBG: voxelOccupancyCount[%i] = %i\n", i, voxelOccupancyCount[i]);
    // Determine which voxel has the higest occupancy count
    // i.e. which is the most frequent value.
    if (voxelOccupancyCount[i] > currentLargestVoxelOccupancyCount) {
      currentLargestVoxelOccupancyCount = voxelOccupancyCount[i];
      mostFrequentVoxelOccupancy = i;
    }
  }
  mprintf("CDBG: Most frequent occupancy is %i (%i occurences)\n", mostFrequentVoxelOccupancy,
          currentLargestVoxelOccupancyCount);
  // The assumption here is that mostFrequentVoxelOccupancy corresponds to the
  // the occupancy count of bulk solvent
  int normalisationFactor = mostFrequentVoxelOccupancy;

  for (Grid::iterator gval = grid_.begin(); gval != grid_.end(); ++gval) {
    double gridval = (double)(*gval);
    gridval = -1.0 * log( (gridval / normalisationFactor) + 0.00000001 );
    *gval = (float)gridval;
  }

  grid_.PrintXplor( filename_, "", "REMARKS Change in Free energy from bulk solvent with bin normalisation of " + integerToString(normalisationFactor) );
}

