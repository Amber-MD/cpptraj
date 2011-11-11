// Molsurf
#include <cstdlib> // Using malloc since interfacing with C code
#include <cstring> 
#include "Action_Molsurf.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Molsurf::Molsurf() {
  //fprintf(stderr,"Angle Con\n");
  sasa=NULL;
  atom = NULL;
  probe_rad = 1.4;
} 

// DESTRUCTOR
Molsurf::~Molsurf() { 
  if (atom!=NULL) free(atom);
}

// Molsurf::init()
/** Expected call: molsurf [<name>] [<mask1>] [out filename] [probe <probe_rad>]
  * Dataset name will be the last arg checked for. Check order is:
  *    1) Keywords
  *    2) Masks
  *    3) Dataset name
  */
int Molsurf::init() {
  char *mask1;
  char *molsurfFile;

  // Get keywords
  molsurfFile = actionArgs.getKeyString("out",NULL);
  probe_rad = actionArgs.getKeyDouble("probe",1.4);

  // Get Masks
  mask1 = actionArgs.getNextMask();
  Mask1.SetMaskString(mask1);

  // Dataset to store angles
  sasa = DSL->Add(DOUBLE, actionArgs.getNextString(),"MSURF");
  if (sasa==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(molsurfFile,sasa);

  mprintf("    MOLSURF: [%s] Probe Radius=%-8.3lf\n",Mask1.maskString,probe_rad);

  return 0;
}

// Molsurf::setup()
/** Set mask up for this parmtop.
  * currentParm is set in Action::Setup
  */
int Molsurf::setup() {

  if ( Mask1.SetupMask(currentParm,activeReference,debug) ) return 1;
  if (Mask1.None()) {
    mprintf("    Error: Molsurf::setup: Mask contains 0 atoms.\n");
    return 1;
  }

  mprintf("    MOLSURF: Calculating surface area for %i atoms.\n",Mask1.Nselected);
  // NOTE: If Mask is * dont include any solvent?

  // The ATOM structure is how molsurf organizes atomic data. Allocate
  // here and fill in parm info. Coords will be filled in during action. 
  atom = (ATOM*) malloc(Mask1.Nselected * sizeof(ATOM));
  if (atom==NULL) {
    mprinterr("Error: Molsurf::Setup Could not allocate memory for ATOMs.\n");
    return 1;
  }
  // Set up parm info for atoms in mask
  double *Radii = currentParm->GB_radii();
  for (int maskidx = 0; maskidx < Mask1.Nselected; maskidx++) {
    int parmatom = Mask1.Selected[maskidx];
    int nres = currentParm->atomToResidue(parmatom);
    atom[maskidx].anum = parmatom + 1; // anum is for debug output only, atoms start from 1
    strcpy(atom[maskidx].anam,currentParm->names[parmatom]);
    atom[maskidx].rnum = nres + 1; // again for debug output only, residues start from 1
    strcpy(atom[maskidx].rnam,currentParm->resnames[nres]);
    atom[maskidx].pos[0] = 0;
    atom[maskidx].pos[1] = 0;
    atom[maskidx].pos[2] = 0;
    atom[maskidx].q = currentParm->charge[parmatom];
    atom[maskidx].rad = Radii[parmatom];
  }

  return 0;  
}

// Molsurf::action()
int Molsurf::action() {
  double molsurf_sasa;

  // Set up coordinates for atoms in mask
  int i3 = 0;
  for (int maskidx = 0; maskidx < Mask1.Nselected; maskidx++) {
    int i3 = Mask1.Selected[maskidx] * 3;
    atom[maskidx].pos[0] = currentFrame->X[i3  ];
    atom[maskidx].pos[1] = currentFrame->X[i3+1];
    atom[maskidx].pos[2] = currentFrame->X[i3+2];
  }
  
  molsurf_sasa = molsurf( probe_rad, atom, Mask1.Nselected); 

  sasa->Add(frameNum, &molsurf_sasa);

  //fprintf(outfile,"%10i %10.4lf\n",frameNum,molsurf_sasa);
  
  return 0;
} 


