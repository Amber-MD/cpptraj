#include <cmath>
#include "Action_Grid.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Grid::Action_Grid() :
  max_(0.80),
  madura_(0),
  smooth_(0),
  invert_(false)
{}

void Action_Grid::Help() {

}

// Action_Grid::init()
/** Usage: grid <filename> nx dx ny dy nz dz [origin] [negative] 
  *             [max <fraction>] [smoothdensity <value>] [invert] [madura <madura>]
  *             <mask>
  *
  * <fraction>: Percent of max to write.
  * <madura>  : Grid values lower than <madura> become flipped in sign, exposes low density.
  * <value>   : Used to smooth density.
  */
int Action_Grid::init() {
  // Get output filename
  filename_ = actionArgs.GetStringNext();
  if (filename_.empty()) {
    mprinterr("Error: GRID: no filename specified.\n");
    return 1;
  }
  // Get grid options
  if (grid_.GridInit( "GRID", actionArgs ))
    return 1;

  // Get extra options
  max_ = actionArgs.getKeyDouble("max", 0.80);
  madura_ = actionArgs.getKeyDouble("madura", 0);
  smooth_ = actionArgs.getKeyDouble("smoothdensity", 0);
  invert_ = actionArgs.hasKey("invert");
  pdbname_ = actionArgs.GetStringKey("pdb"); 

  // Get mask
  ArgList::ConstArg maskexpr = actionArgs.getNextMask();
  if (maskexpr==NULL) {
    mprinterr("Error: GRID: No mask specified.\n");
    return 1;
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

  return 0;
}

// Action_Grid::setup()
int Action_Grid::setup() {
  // Setup grid, checks box info.
  if (grid_.GridSetup( currentParm )) return 1;

  // Setup mask
  if (currentParm->SetupIntegerMask( mask_ ))
    return 1;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprinterr("Error: GRID: No atoms selected for parm %s\n", currentParm->c_str());
    return 1;
  }

  return 0;
}

// Action_Grid::action()
int Action_Grid::action() {
  if (grid_.GridBox()) {
    Vec3 boxCenter( currentFrame->BoxX() / 2.0,
                    currentFrame->BoxY() / 2.0,
                    currentFrame->BoxZ() / 2.0 );
    for (AtomMask::const_iterator atom = mask_.begin();
                                  atom != mask_.end(); ++atom)
    {
      Vec3 XYZ = currentFrame->XYZ( *atom );
      XYZ -= boxCenter;
      //mprintf("BATM %6i ", *atom + 1);
      grid_.GridPoint( XYZ[0], XYZ[1], XYZ[2] );
    }
  } else {
    for (AtomMask::const_iterator atom = mask_.begin();
                                  atom != mask_.end(); ++atom)
    {
      const double* XYZ = currentFrame->XYZ( *atom );
      //mprintf("ATOM %6i ", *atom + 1);
      grid_.GridPoint( XYZ[0], XYZ[1], XYZ[2] );
    }
  }

  return 0;
}

// Action_Grid::print()
void Action_Grid::print() {
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
  grid_.PrintXplor( filename_, "This line is ignored", 
                    "rdparm generated grid density" );

  // PDBfile output
  mprintf("    GRID: grid max is %.3lf\n", gridMax);
  grid_.PrintPDB( pdbname_, max_, gridMax );
}

