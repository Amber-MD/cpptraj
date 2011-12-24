// Closest
// Find closest waters to atoms in mask.
#include <cmath>
#include <algorithm> // sort
#include <cfloat> // DBL_MAX
#ifdef _OPENMP
#  include "omp.h"
#endif
#include "Action_Closest.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Closest::Closest() {
  //fprintf(stderr,"Closest Con\n");
  noimage=false;
  firstAtom=false;
  imageType=0;
  oldParm=NULL;
  newParm=NULL;
  outFile=NULL;
  outList=NULL;
  framedata=NULL;
  moldata=NULL;
  distdata=NULL;
  atomdata=NULL;
  Nclosest=0;
  prefix=NULL;
} 

// DESTRUCTOR
Closest::~Closest() {
  //fprintf(stderr,"Closest Destructor.\n");
  //this->ClearMaskList(); 
  if (newParm!=NULL) delete newParm;
  if (outList!=NULL) delete outList;
}

// Closest::init()
/** Expected call: closest <# to keep> <mask> [noimage] [first/oxygen] 
  *                [closestout <filename> [outprefix <parmprefix>]
  * Dataset name will be the last arg checked for. Check order is:
  *    1) Keywords
  *    2) Masks
  *    3) Dataset name
  */
int Closest::init( ) {
  char *mask1;

  // Get Keywords
  closestWaters = actionArgs.getNextInteger(-1);
  if (closestWaters < 0) {
    mprinterr("Error: Closest::init(): Invalid # solvent molecules to keep (%i).\n",
              closestWaters);
    return 1;
  }
  if ( actionArgs.hasKey("oxygen") || actionArgs.hasKey("first") )
    firstAtom=true;
  noimage = actionArgs.hasKey("noimage");
  prefix = actionArgs.getKeyString("outprefix",NULL);
  // Setup output file and sets if requested.
  // Will keep track of Frame, Mol#, Distance, and first solvent atom
  mask1 = actionArgs.getKeyString("closestout",NULL);
  if (mask1 != NULL) {
    // Set up datasets
    outList = new DataSetList();
    framedata = outList->Add(INT,(char*)"Frame\0","Frame");
    moldata   = outList->Add(INT,(char*)"Mol\0","Mol");
    distdata  = outList->Add(DOUBLE,(char*)"Dist\0","Dist");
    atomdata  = outList->Add(INT,(char*)"FirstAtm\0","FirstAtm");
    if (framedata==NULL || moldata==NULL || distdata==NULL || atomdata==NULL) {
      mprinterr("Error: Closest::init(): Could not setup data sets for output file %s\n",mask1);
      return 1;
    }
    // Add sets to datafile in list.
    outFile = DFL->Add(mask1, framedata);
    outFile = DFL->Add(mask1, moldata);
    outFile = DFL->Add(mask1, distdata);
    outFile = DFL->Add(mask1, atomdata);
    if (outFile==NULL) {
      mprintf("Error: Closest::init(): Could not setup output file %s\n",mask1);
      return 1;
    }
  }

  // Get Masks
  mask1 = actionArgs.getNextMask();
  if (mask1==NULL) {
    mprinterr("Error: Closest::init(): No mask specified.\n");
    return 1;
  }
  soluteMask.SetMaskString(mask1);

  mprintf("    CLOSEST: Finding closest %i solvent molecules to",closestWaters);
  //if (mask1==NULL)
  //  fprintf(stdout,"             all solute atoms\n");
  //else
    mprintf(" atoms in mask %s\n",soluteMask.MaskString());
  if (noimage) 
    mprintf("             Imaging will be turned off.\n");
  if (firstAtom)
    mprintf("             Only first atom of solvent molecule used for distance calc.\n");
  if (outFile!=NULL)
    mprintf("             Closest molecules will be saved to %s\n",outFile->Filename());
  if (prefix!=NULL)
    mprintf("             Stripped topology file will be written with prefix %s\n",prefix);

  return 0;
}

// Closest::setup()
/** Like the strip action, closest will modify the current parm keeping info
  * for atoms in mask plus the closestWaters solvent molecules. Set up the
  * vector of MolDist objects, one for every solvent molecule in the original
  * parm file. Atom masks for each solvent molecule will be set up.
  */
int Closest::setup() {
  int solventMol, NsolventAtoms;
  MolDist solvent;

  // If there are no solvent molecules this action is not valid.
  if (currentParm->solventMolecules==0) {
    mprintf("Warning: Closest::setup: Parm %s does not contain solvent.\n",currentParm->parmName);
    return 1;
  }
  // If # solvent to keep >= solvent in this parm the action is not valid.
  if (closestWaters >= currentParm->solventMolecules) {
    mprintf("Warning: Closest::setup: # solvent to keep (%i) >= # solvent molecules in\n",
            closestWaters);
    mprintf("                           %s (%i).\n",currentParm->parmName,
            currentParm->solventMolecules);
    return 1;
  } 
  // Check that all solvent molecules contain same # atoms. Solvent 
  // molecules must be identical for the command to work properly; 
  // the prmtop strip occurs only once so the solvent params become fixed.
  NsolventAtoms = currentParm->solventMoleculeStop[0] - currentParm->solventMoleculeStart[0];
  for (solventMol = 1; solventMol < currentParm->solventMolecules; solventMol++) {
    if ( NsolventAtoms != 
         (currentParm->solventMoleculeStop[solventMol] - 
          currentParm->solventMoleculeStart[solventMol]) ) 
    {
      mprintf("Warning: Closest::setup: Solvent molecules in %s are not of uniform size.\n",
              currentParm->parmName);
      mprintf("         First solvent mol=%i atoms, %i solvent mol=%i atoms.\n",
              NsolventAtoms, solventMol,
              currentParm->solventMoleculeStop[solventMol] - 
              currentParm->solventMoleculeStart[solventMol]);
      return 1;
    }
  }

  // Setup solute atom mask
  // NOTE: Should ensure that no solvent atoms are selected!
  if ( currentParm->SetupIntegerMask(soluteMask, activeReference) ) return 1;
  if (soluteMask.None()) {
    mprintf("Warning: Closest::setup: Mask %s contains no atoms.\n",soluteMask.MaskString());
    return 1;
  }

/*  if (currentParm->mass==NULL && useMass) {
    fprintf(stdout,"    Warning: Closest::setup: Mass for this parm is NULL.\n");
    fprintf(stdout,"             Geometric center of mass will be used.\n");
    useMass=false;
  }*/

  // Figure out imaging - check box based on prmtop box
  // NOTE: Should box be figured out from read-in coords?
  imageType = 0;
  if (!noimage) {
    imageType = (int)currentParm->boxType;
    if (currentParm->boxType==NOBOX) {
        mprintf("Warning: Closest::setup: ");
        mprintf(" Imaging specified but no box information in prmtop %s\n",currentParm->parmName);
        mprintf("             No imaging can occur..\n");
    }
  }

  // Figure out what the the total size of the selected solute atoms plus
  // the number of kept solvent atoms is in order to set up the stripped
  // parm.
  // solventMoleculeStop is really the first atom of the next molecule.
  // Want to keep the solute in addition to the closest solvent.
  stripMask = soluteMask;
  //mprintf("DEBUG:\t");
  //stripMask.PrintMaskAtoms("stripMask");
  NsolventAtoms = currentParm->solventMoleculeStop[closestWaters-1] - 
                  currentParm->solventMoleculeStart[0];
  mprintf("\tKeeping %i solvent atoms.\n",NsolventAtoms);
  // Put solvent atom #s at end of array for creating stripped parm
  stripMask.AddAtomRange(currentParm->solventMoleculeStart[0], 
                         currentParm->solventMoleculeStop[closestWaters-1]);
  //mprintf("DEBUG:\t");
  //stripMask.PrintMaskAtoms("stripMaskWsolvent");

  // Store old parm
  oldParm = currentParm;
 
  // Get atom masks for all solvent molecules in old parm.
  SolventMols.clear();
  solvent.D=0.0;
  solvent.mol=0;
  SolventMols.resize(oldParm->solventMolecules, solvent);
  int firstSolventMol = oldParm->FirstSolventMol();
  for (solventMol=0; solventMol < oldParm->solventMolecules; solventMol++) {
    SolventMols[solventMol].mol = firstSolventMol + solventMol;
    // Setup solvent molecule mask
    SolventMols[solventMol].mask.AddAtomRange(oldParm->solventMoleculeStart[solventMol],
                                               oldParm->solventMoleculeStop[solventMol]  );
    //SolventMols[solventMol].mask.PrintMaskAtoms("solvent");
  }

  // Create stripped Parm
  if (newParm!=NULL) delete newParm;
  newParm = currentParm->modifyStateByMask(stripMask.Selected, prefix);
  if (newParm==NULL) {
    mprinterr("Error: Closest::setup: Could not create new parmtop.\n");
    return 1;
  }
  newParm->Summary();

  // Allocate space for new frame
  newFrame.SetupFrame(newParm->natom, newParm->mass);

  // If prefix given then output stripped parm
  if (prefix!=NULL) {
    mprintf("\tWriting out amber topology file %s\n",newParm->parmName);
    if ( newParm->WriteAmberParm(newParm->parmName) ) {
      mprinterr("Error: CLOSEST: Could not write out stripped parm file %s\n",
              newParm->parmName);
    }
  }

  // Set parm
  currentParm = newParm;

  return 0;  
}

// Closest::action()
/** Find the minimum distance between atoms in soluteMask and each solvent Mask.
  */
int Closest::action() {
  int solventMol, solventAtom, maskPosition, atom, maxSolventMolecules;
  double Dist, maxD, ucell[9], recip[9];

  if (imageType>0) {
    currentFrame->BoxToRecip(ucell, recip);
    // Calculate max possible imaged distance
    maxD = currentFrame->box[0] + currentFrame->box[1] + currentFrame->box[2];
    maxD = maxD * maxD;
  } else {
    // If not imaging, set max distance to an arbitrarily large number
    maxD = DBL_MAX;
  }

  // Loop over all solvent molecules in original frame
  maxSolventMolecules = oldParm->solventMolecules;
  // DEBUG
  //mprintf("Closest: Begin parallel loop for %i\n",frameNum);
  // DEBUG
#ifdef _OPENMP
#pragma omp parallel private(solventMol,atom,Dist,solventAtom)
{
  //mprintf("OPENMP: %i threads\n",omp_get_num_threads());
#pragma omp for
#endif
  for (solventMol=0; solventMol < maxSolventMolecules; solventMol++) {
    //mprintf("[%i] Calculating distance for molecule %i\n",omp_get_thread_num(),solventMol);
    // Set the initial minimum distance for this solvent mol to be the
    // max possible distance.
    SolventMols[solventMol].D = maxD;

    // DEBUG - show solvent mask
    //fprintf(stdout,"      Solvent %i %i %i\n", MaskList[solventMol]->Selected[0]+1,
    //        MaskList[solventMol]->Selected[1]+1,MaskList[solventMol]->Selected[2]+1);

    // Calculate distance between each atom in soluteMask and atoms in solvent Mask
    for (atom = 0; atom < soluteMask.Nselected; atom++) {
      Dist = currentFrame->DIST2(soluteMask.Selected[atom], 
                                 SolventMols[solventMol].mask.Selected[0], 
                                 imageType, ucell, recip);
      if (Dist < SolventMols[solventMol].D) SolventMols[solventMol].D=Dist;
      //fprintf(stdout,"D atom %i %i = %lf image %i\n",soluteMask.Selected[atom],
      //        MaskList[solventMol]->Selected[0],minD,imageType);
      // Check the rest of the solvent atoms if specified
      if (!firstAtom) {
        for (solventAtom=1; solventAtom<SolventMols[solventMol].mask.Nselected; solventAtom++) 
        {
          Dist = currentFrame->DIST2(soluteMask.Selected[atom], 
                               SolventMols[solventMol].mask.Selected[solventAtom],
                               imageType, ucell, recip);
          if (Dist < SolventMols[solventMol].D) SolventMols[solventMol].D=Dist;
        }
      }
    }

    // DEBUG - Print distances
    //mprintf("DEBUG:\tMol %8i minD= %lf\n",solventMol, SolventMols[solventMol].D);
  } // END for loop over solventMol
#ifdef _OPENMP
} // END pragma omp parallel
#endif
  // DEBUG
  //mprintf("Closest: End parallel loop for %i, got %i Distances.\n",frameNum,(int)SolventMols.size());
  // DEBUG

  // Sort distances
  sort( SolventMols.begin(), SolventMols.end(), moldist_cmp() );
  // DEBUG
  //mprintf("Closest: Distances sorted for %i\n",frameNum);
  // DEBUG

  // Add first closestWaters solvent atoms to stripMask
  solventMol = 0;
  maskPosition = soluteMask.Nselected;
  for ( std::vector<MolDist>::iterator solvent = SolventMols.begin();
                                       solvent != SolventMols.end();
                                       solvent++ ) 
  {
    //mprintf("DEBUG:\tmol %i ",(*solvent).mol);
    //(*solvent).mask.PrintMaskAtoms("Mask");
    for (solventAtom = 0; solventAtom < (*solvent).mask.Nselected; solventAtom++)
      stripMask.Selected[maskPosition++] = (*solvent).mask.Selected[solventAtom];
    // Record which water molecules are closest if requested
    if (outFile!=NULL) {
      framedata->Add(Nclosest, &frameNum);
      moldata->Add(Nclosest, &((*solvent).mol));
      Dist = sqrt( (*solvent).D );
      distdata->Add(Nclosest, &Dist);
      atomdata->Add(Nclosest, &((*solvent).mask.Selected[0]));
      Nclosest++;
    }
    // DEBUG - print first closestWaters distances
    //mprintf("DEBUG: Mol %i   D2= %lf   Atom0= %i\n",(*it).mol, (*it).D, (*it).mask->Selected[0]);
    solventMol++;
    if (solventMol==closestWaters) break;
  }

  // Modify and set frame
  //mprintf("DEBUG:\t");
  //stripMask.PrintMaskAtoms("action_stripMask");
  newFrame.SetFrameFromMask(currentFrame, &stripMask);
  currentFrame = &newFrame;

  return 0;
} 

// Closest::print()
/** Set up the closest output file for writing. Since the datasets used in
  * closest are local and not part of the master data set list, call Sync
  * for them here. Also set so that the X column (which is just an index)
  * is not written.
  */
void Closest::print() {
  if (outFile==NULL) return;
  // Sync up the dataset list here since it is not part of the master 
  // dataset list.
  outList->Sync();
  // Set specific datafile options
  outFile->SetNoXcol();
}

