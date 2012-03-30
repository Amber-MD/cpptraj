// DISTRMSD
#include "Action_DistRmsd.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DistRmsd::DistRmsd() {
  drmsd=NULL;
}

// DistRmsd::init()
/** Called once before traj processing. Set up reference info.
  * Expected call: 
  * drmsd <name> <mask> [<refmask>] [out filename] 
  *       [ first | ref <filename> | refindex <#> | 
  *         reftraj <filename> [parm <parmname> | parmindex <#>] ] 
  */
int DistRmsd::init( ) {
  // Check for keywords
  char *rmsdFile = actionArgs.getKeyString("out",NULL);

  // Get the RMS mask string for target 
  char *mask0 = actionArgs.getNextMask();
  TgtMask.SetMaskString(mask0);

  // Initialize reference. If no reference mask given, mask0 will be used.
  if (RefInit(true, false, mask0, actionArgs, FL, PFL, NULL))
    return 1;

  // Set up the RMSD data set
  drmsd = DSL->Add(DataSet::DOUBLE, actionArgs.getNextString(),"DRMSD");
  if (drmsd==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(rmsdFile,drmsd);

  mprintf("    DISTRMSD: (%s), reference is",TgtMask.MaskString());
  RefInfo();
  mprintf("\n");

  return 0;
}

// DistRmsd::setup()
/** Called every time the trajectory changes. Set up TgtMask for the new 
  * parmtop and allocate space for selected atoms from the Frame.
  */
int DistRmsd::setup() {

  if ( currentParm->SetupIntegerMask(TgtMask) ) return 1;
  if ( TgtMask.None() ) {
    mprintf("    Error: DistRmsd::setup: No atoms in mask.\n");
    return 1;
  }
  // Allocate space for selected atoms in the frame. This will also put the
  // correct masses in based on the mask.
  SelectedTgt.SetupFrameFromMask(TgtMask, currentParm->Mass());

  // Reference setup
  if (RefSetup( currentParm )) return 1;

  // Check that num atoms in frame mask from this parm match ref parm mask
  if ( RefNselected() != TgtMask.Nselected() ) {
    mprintf( "    Error: Number of atoms in RMS mask (%i) does not \n",TgtMask.Nselected());
    mprintf( "           equal number of atoms in Ref mask (%i).\n",RefNselected());
    return 1;
  }

  return 0;
}

// DistRmsd::action()
/** Called every time a frame is read in. Calc distance RMSD.
  * If first is true, set the first frame read in as reference.
  */
int DistRmsd::action() {
  // Perform any needed reference actions
  RefAction(currentFrame, NULL);

  // Set selected frame atoms. Masses have already been set.
  SelectedTgt.SetCoordinates(*currentFrame, TgtMask);

  // DEBUG
/*  mprintf("  DEBUG: RMSD: First atom coord in SelectedTgt is : "); 
  SelectedTgt->printAtomCoord(0);
  mprintf("  DEBUG: RMSD: First atom coord in SelectedRef is : ");
  SelectedRef->printAtomCoord(0);
*/

  double DR = SelectedTgt.DISTRMSD( &SelectedRef_ );

  drmsd->Add(frameNum, &DR);

  return 0;
}

