#include <cmath>
#include "Action_Grid.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Grid::Action_Grid() :
  max_(0.80),
  madura_(0),
  smooth_(0),
  invert_(false),
  dxform_(false)
{}

void Action_Grid::Help() {
  mprintf("\t<filename> %s <mask>\n", Grid::HelpText);
  mprintf("\t[max <fraction>] [smoothdensity <value>] [invert] [madura <madura>]\n");
  mprintf("\t[pdb <pdbout>] [opendx]\n");
  mprintf("\tBin atoms in <mask> into a 3D grid.\n");
  mprintf("\t<fraction>: Percent of max to write.\n");
  mprintf("\t<madura>  : Grid values lower than <madura> become flipped in sign, exposes low density.\n");
  mprintf("\t<value>   : Used to smooth density.\n");
  mprintf("\t[opendx]  : Write the density file in OpenDX format.\n");
}

// Action_Grid::init()
Action::RetType Action_Grid::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get output filename
  filename_ = actionArgs.GetStringNext();
  if (filename_.empty()) {
    mprinterr("Error: GRID: no filename specified.\n");
    return Action::ERR;
  }
  // Get grid options
  if (grid_.GridInit( "GRID", actionArgs ))
    return Action::ERR;

  // Get extra options
  max_ = actionArgs.getKeyDouble("max", 0.80);
  madura_ = actionArgs.getKeyDouble("madura", 0);
  smooth_ = actionArgs.getKeyDouble("smoothdensity", 0);
  invert_ = actionArgs.hasKey("invert");
  dxform_ = actionArgs.hasKey("opendx");
  pdbname_ = actionArgs.GetStringKey("pdb"); 

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
  if (pdbname_.empty())
    mprintf("\tPseudo-PDB will be printed to STDOUT.\n");
  else
    mprintf("\tPseudo-PDB will be printed to %s\n", pdbname_.c_str());
  // TODO: print extra options

  // Allocate grid
  //if (GridAllocate()) return 1;

  return Action::OK;
}

// Action_Grid::setup()
Action::RetType Action_Grid::Setup(Topology* currentParm, Topology** parmAddress) {
  // Setup grid, checks box info.
  if (grid_.GridSetup( *currentParm )) return Action::ERR;

  // Setup mask
  if (currentParm->SetupIntegerMask( mask_ ))
    return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprinterr("Error: GRID: No atoms selected for parm %s\n", currentParm->c_str());
    return Action::ERR;
  }

  return Action::OK;
}

// Action_Grid::action()
// TODO: Combine all below
Action::RetType Action_Grid::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  grid_.GridFrame( *currentFrame, mask_ );
  return Action::OK;
}

// Action_Grid::print()
void Action_Grid::Print() {
  // DEBUG
  //mprintf("CDBG: Printing grid.\n");
  //PrintEntireGrid();
  // END DEBUG

  // Perform normalization and find max
  double gridMax = 0;
  for (Grid::iterator gval = grid_.begin(); gval != grid_.end(); ++gval) {
    double gridval = (double)(*gval);
    // ----- SMOOTHING -----
    if (smooth_ > 0.0) {
      double yy = gridval - smooth_;
      double xx = yy*yy / (0.2 * smooth_ * smooth_);
      xx = exp( -xx );
      if (invert_) {
        if (gridval > smooth_) // NOTE: Comparison OK? Needs cast?
          gridval = -5.0;
        else
          gridval -= gridval * xx;
        /* COMMENTED OUT IN ORIGINAL PTRAJ CODE
        if (gridInfo->grid[index] < action->darg3) {
          gridInfo->grid[index] = 0.0;
        }
        */
        if (gridval >= 0)
          gridval = smooth_ - gridval;
      } else {
        if (gridval < smooth_)
          gridval = 0;
        else
          gridval -= gridval * xx;
        if (gridval < smooth_)
          gridval = 0;
      }
    }

    // do the madura negative option to expose low density
    if ( madura_ > 0.0 && gridval > 0.0 && gridval < madura_ )
      *gval = (float)-gridval;
    else
      *gval = (float) gridval;

    if ( gridval > gridMax )
      gridMax = gridval;
  } 

  // Write Xplor file
  if (dxform_)
    grid_.PrintDX(filename_);
  else
    grid_.PrintXplor( filename_, "This line is ignored", 
                      "rdparm generated grid density" );

  // PDBfile output
  mprintf("    GRID: grid max is %.3lf\n", gridMax);
  grid_.PrintPDB( pdbname_, max_, gridMax );
}

