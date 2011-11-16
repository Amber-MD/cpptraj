// CheckStructure
#include <cmath>
#include "Action_CheckStructure.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
CheckStructure::CheckStructure() {
  //fprintf(stderr,"CheckStructure Con\n");
  noimage=false;
  imageType=0;
  bondoffset=1.0;
  nonbondcut2=0.64; // 0.8^2
} 

// DESTRUCTOR
CheckStructure::~CheckStructure() {
  //fprintf(stderr,"CheckStructure Destructor.\n");
  outfile.CloseFile();
}

// CheckStructure::init()
/** Expected call: check[structure] [<mask1>] [reportfile <report>] [noimage] 
  *                     [offset <offset>] [cut <cut>]
  */
// Dataset name will be the last arg checked for. Check order is:
//    1) Keywords
//    2) Masks
//    3) Dataset name
int CheckStructure::init( ) {
  char *mask1;
  char *reportFile;
  double nonbondcut;

  // Get Keywords
  noimage = actionArgs.hasKey("noimage");
  reportFile = actionArgs.getKeyString("reportfile",NULL);
  bondoffset = actionArgs.getKeyDouble("offset",1.0);
  nonbondcut = actionArgs.getKeyDouble("cut",0.8);

  // Get Masks
  mask1 = actionArgs.getNextMask();
  //fprintf(stdout,"    Mask 1: %s\n",mask1);
  Mask1.SetMaskString(mask1);

  mprintf("    CHECKSTRUCTURE: Checking atoms in mask [%s]",Mask1.maskString);
  if (noimage) 
    mprintf(", imaging off");
  if (reportFile!=NULL)
    mprintf(", output to %s",reportFile);
  mprintf(".\n");
  mprintf("                    Warnings will be printed for bond length > eq + %.2lf\n",
          bondoffset);
  mprintf("                    and non-bond distance < %.2lf\n",nonbondcut);
  nonbondcut2 = nonbondcut * nonbondcut;

  if (outfile.SetupFile(reportFile, WRITE, DATAFILE, STANDARD, debug))
    return 1;
  outfile.OpenFile();

  return 0;
}

// CheckStructure::setup()
/** Determine what atoms each mask pertains to for the current parm file.
  * Also determine whether imaging should be performed. Check if parm
  * has bonds.
  */
int CheckStructure::setup() {
  double req = 0;
  double rk = 0;

  if ( Mask1.SetupMask(currentParm,activeReference,debug) ) return 1;

  if (Mask1.None()) {
    mprintf("    Error: CheckStructure::setup: Mask has no atoms.\n");
    return 1;
  }

  // Check imaging - check box based on prmtop box
  imageType = 0;
  if (!noimage) {
    imageType = (int)currentParm->boxType;
    if (currentParm->boxType==NOBOX && debug>0) {
      mprintf("    Warning: No box info in %s, disabling imaging.\n",currentParm->parmName);
    }
  }

  // Check bonds
  req_array.resize( ((Mask1.Nselected * Mask1.Nselected)-Mask1.Nselected) / 2, 0);
  if ( currentParm->NbondsWithH + currentParm->NbondsWithoutH <= 0 ) {
    mprintf("    Warning: No bond info in %s, will not check bonds.\n",currentParm->parmName);
  } else {
    // Pre-fetch all necessary bond eq constants
    int nbon = 0;
    for (int maskidx1 = 0; maskidx1 < Mask1.Nselected - 1; maskidx1++) {
      int atom1 = Mask1.Selected[maskidx1];
      for (int maskidx2 = maskidx1 + 1; maskidx2 < Mask1.Nselected; maskidx2++) {
        int atom2 = Mask1.Selected[maskidx2];
        if (currentParm->GetBondParam(&rk, &req, atom1, atom2)) {
          double bondmax = req + bondoffset;
          bondmax *= bondmax;
          req_array[nbon] = bondmax;
        }
        nbon++;
      }
    }
  }      

  // Print imaging info for this parm
  mprintf("    CHECKSTRUCTURE: %s (%i atoms)",Mask1.maskString, Mask1.Nselected);
  if (imageType > 0)
    mprintf(", imaging on");
  else
    mprintf(", imaging off");
  mprintf(".\n");
        
  return 0;  
}

// CheckStructure::action()
int CheckStructure::action() {
  double ucell[9], recip[9], D2, D, bondmax;
  int nbon = 0;

  if (imageType>0) currentFrame->BoxToRecip(ucell,recip);

  int lastidx = Mask1.Nselected - 1;
  for (int maskidx1 = 0; maskidx1 < lastidx; maskidx1++) {
    int atom1 = Mask1.Selected[maskidx1];
    for (int maskidx2 = maskidx1 + 1; maskidx2 < Mask1.Nselected; maskidx2++) {
      int atom2 = Mask1.Selected[maskidx2];
      // Get distance^2
      D2 = currentFrame->DIST2(atom1,atom2,imageType,ucell,recip);
      // Check for long bond length
      bondmax = req_array[nbon];
      // req > 0 means atoms bonded, check if distance > pre-stored req (bondmax)
      if (bondmax>0) {
        // Check if distance2 > against (req+bondoffset)^2
        if (D2 > bondmax) {
          D = sqrt(D2);
          outfile.IO->Printf(
                  "%i\t Warning: Unusual bond length %i@%s to %i@%s (%.2lf)\n",
                  frameNum+OUTPUTFRAMESHIFT,
                  atom1+1, currentParm->names[atom1], atom2+1, currentParm->names[atom2],
                  D);
        }
      // req == 0 means atoms not bonded, check overlap
      } else {
        if (D2 < nonbondcut2) {
          D = sqrt(D2);
          outfile.IO->Printf(
                  "\t Warning: Atoms %i@%s and %i@%s are close (%.2lf\n",
                  frameNum+OUTPUTFRAMESHIFT,
                  atom1+1, currentParm->names[atom1], atom2+1, currentParm->names[atom2],D);
        }
      }
      nbon++;
    } // END second loop over mask atoms
  } // END first loop over mask atoms

  return 0;
}

