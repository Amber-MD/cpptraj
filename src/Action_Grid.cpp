#include <cmath>
#include "Action_Grid.h"
#include "CpptrajStdio.h"
#include "PDBfile.h"

// CONSTRUCTOR
Action_Grid::Action_Grid() :
  max_(0.80),
  madura_(0),
  smooth_(0),
  invert_(false)
{}

// Action_Grid::init()
int Action_Grid::init() {
  // Get output filename
  filename_ = actionArgs.GetStringNext();
  if (filename_.empty()) {
    mprinterr("Error: GRID: no filename specified.\n");
    return 1;
  }
  // Get grid options
  if (GridInit( "GRID", actionArgs ))
    return 1;

  // Get extra options
  max_ = actionArgs.getKeyDouble("max", 0.80);
  madura_ = actionArgs.getKeyDouble("madura", 0);
  invert_ = actionArgs.getKeyDouble("smoothdensity", 0);
  invert_ = actionArgs.hasKey("invert");
  pdbname_ = actionArgs.GetStringKey("pdb"); 

  // Get mask
  char* maskexpr = actionArgs.getNextMask();
  if (maskexpr==NULL) {
    mprinterr("Error: GRID: No mask specified.\n");
    return 1;
  }
  mask_.SetMaskString(maskexpr);

  // Info
  GridInfo();
  mprintf("\tGrid will be printed to file %s\n",filename_.c_str());
  mprintf("\tMask expression: [%s]\n",mask_.MaskString());
  if (pdbname_.empty())
    mprintf("\tPseudo-PDB will be printed to STDOUT.\n");
  else
    mprintf("\tPseudo-PDB will be printed to %s\n", pdbname_.c_str());
  // TODO: print extra options

  // Allocate grid
  if (GridAllocate()) return 1;

  return 0;
}

// Action_Grid::setup()
int Action_Grid::setup() {
  if (GridSetup( currentParm )) return 1;

  // Setup mask
  if (currentParm->SetupIntegerMask( mask_ ))
    return 1;
  mprintf("\t[%s] %i atoms selected.\n", mask_.MaskString(), mask_.Nselected());
  if (mask_.None()) {
    mprinterr("Error: GRID: No atoms selected for parm %s\n", currentParm->c_str());
    return 1;
  }

  return 0;
}

// Action_Grid::action()
int Action_Grid::action() {
  double XYZ[3], boxcrd[3];
  if (GridBox()) {
    currentFrame->BoxXYZ( boxcrd );
    boxcrd[0] /= 2.0;
    boxcrd[1] /= 2.0;
    boxcrd[2] /= 2.0;
    for (AtomMask::const_iterator atom = mask_.begin();
                                  atom != mask_.end(); ++atom)
    {
      currentFrame->GetAtomXYZ( XYZ, *atom );
      double xx = XYZ[0] - boxcrd[0];
      double yy = XYZ[1] - boxcrd[1];
      double zz = XYZ[2] - boxcrd[2];
      //mprintf("BATM %6i ", *atom + 1);
      GridPoint( xx, yy, zz );
    }
  } else {
    for (AtomMask::const_iterator atom = mask_.begin();
                                  atom != mask_.end(); ++atom)
    {
      currentFrame->GetAtomXYZ( XYZ, *atom );
      //mprintf("ATOM %6i ", *atom + 1);
      GridPoint( XYZ[0], XYZ[1], XYZ[2] );
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
  CpptrajFile outfile;
  if (outfile.OpenWrite( filename_ )) return;
  GridPrintHeader(outfile);

  double gridMax = 0;
  int NZ2 = NZ()/2;
  for (int k = 0; k < NZ(); ++k) {
    outfile.Printf("%8i\n", k - NZ2 + 1);
    for (int j = 0; j < NY(); ++j) {
      int col = 1;
      for (int i = 0; i < NX(); ++i) {
        double gridval = GridVal(i, j, k);
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
          outfile.Printf("%12.5f", -gridval);
        else
          outfile.Printf("%12.5f", gridval);

        if (col && (col%6 == 0))
          outfile.Printf("\n");
        ++col;

        if ( gridval > gridMax )
          gridMax = gridval;
      } // END i loop over x
      if ( (col-1) % 6 != 0 ) // Unless a newline was just written...
        outfile.Printf("\n");
    } // END j loop over y
  } // END k loop over z

  // PDBfile output
  // TODO: Grid pop in Bfactor/occ column
  PDBfile pdbout;
  if (pdbout.OpenPDB(pdbname_)) {
    mprinterr("Error: GRID: Cannot open PDB output.\n");
    return;
  }
  mprintf("    GRID: grid max is %.3lf\n", gridMax);
  mprintf("          dumping a pseudo-pdb representing all points > %.3f\n",
           0.80 * gridMax);
  int res = 1;
  for (int k = 0; k < NZ(); ++k) {
    for (int j = 0; j < NY(); ++j) {
      for (int i = 0; i < NX(); ++i) {
        double gridval = GridVal(i, j, k);
        if (gridval > (max_ * gridMax)) {
          pdbout.WriteHET(res, Xcrd(i), Ycrd(j), Zcrd(k));
          ++res;
        }
      }
    }
  }

}

