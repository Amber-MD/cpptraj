// Closest
// Find closest waters to atoms in mask.
#include <cstdlib>
#include <cmath>
#include <cstring> //strlen
#include <cstdio> // sprintf
#ifdef _OPENMP
#  include "omp.h"
#endif
#include "Action_Closest.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Closest::Closest() {
  //fprintf(stderr,"Closest Con\n");
  noimage=false;
  tempMask=NULL;
  //useMass=true;
  firstAtom=false;
  imageType=0;
  oldParm=NULL;
  newParm=NULL;
  newFrame=NULL;
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
  this->ClearMaskList(); 
  if (newParm!=NULL) delete newParm;
  if (newFrame!=NULL) delete newFrame;
  if (tempMask!=NULL) delete tempMask;
  if (outList!=NULL) delete outList;
}

/*
 * Closest::ClearMaskList()
 */
void Closest::ClearMaskList() {
  if (!MaskList.empty()) {
    // Clear away old array
    for (std::vector<AtomMask*>::iterator it=MaskList.begin(); it!=MaskList.end(); it++) 
      delete *it;
  }
  MaskList.clear();
}

/*
 * Closest::init()
 * Expected call: closest <# to keep> <mask> [noimage] [first/oxygen]
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int Closest::init( ) {
  char *mask1;

  // Get Keywords
  closestWaters = A->getNextInteger(-1);
  if (closestWaters < 0) {
    mprintf("Error: Closest::init(): Invalid # solvent molecules to keep (%i).\n",
            closestWaters);
    return 1;
  }
  if ( A->hasKey("oxygen") || A->hasKey("first") )
    firstAtom=true;
  noimage = A->hasKey("noimage");
  prefix = A->getKeyString("outprefix",NULL);
  // Setup output file and sets if requested.
  // Will keep track of Frame, Mol#, Distance, and first solvent atom
  mask1 = A->getKeyString("closestout",NULL);
  if (mask1 != NULL) {
    // Set up datasets
    outList = new DataSetList();
    framedata = outList->Add(INT,(char*)"Frame\0","Frame");
    moldata   = outList->Add(INT,(char*)"Mol\0","Mol");
    distdata  = outList->Add(DOUBLE,(char*)"Dist\0","Dist");
    atomdata  = outList->Add(INT,(char*)"FirstAtm\0","FirstAtm");
    if (framedata==NULL || moldata==NULL || distdata==NULL || atomdata==NULL) {
      mprintf("Error: Closest::init(): Could not setup data sets for output file %s\n",mask1);
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
  mask1 = A->getNextMask();
  if (mask1==NULL) {
    mprintf("Error: Closest::init(): No mask specified.\n");
    return 1;
  }
  Mask1.SetMaskString(mask1);

  mprintf("    CLOSEST: Finding closest %i solvent molecules to\n",closestWaters);
  //if (mask1==NULL)
  //  fprintf(stdout,"             all solute atoms\n");
  //else
    mprintf("             atoms in mask %s\n",Mask1.maskString);
  if (noimage) 
    mprintf("             non-imaged");
  else
    mprintf("             imaged");
  if (firstAtom)
    mprintf(", using first atom in solvent molecule for distance calc");
  mprintf(".\n");
  if (outFile!=NULL)
    mprintf("             Closest molecules will be saved to %s\n",outFile->filename);
  if (prefix!=NULL)
    mprintf("             Stripped topology file will be written with prefix %s\n",prefix);

  return 0;
}

/*
 * Closest::setup()
 * Like the strip action, closest will modify the current parm keeping info
 * for atoms in mask plus the closestWaters solvent molecules.
 */
int Closest::setup() {
  int solventMol, solventAtom, NsolventAtoms;
  AtomMask *SolventMask;

  // If there are no solvent molecules this action is not valid.
  if (P->solventMolecules==0) {
    mprintf("    Error: Closest::setup: Parm %s does not contain solvent.\n",P->parmName);
    return 1;
  }
  // If # solvent to keep >= solvent in this parm the action is not valid.
  if (closestWaters >= P->solventMolecules) {
    mprintf("    Error: Closest::setup: # solvent to keep (%i) >= # solvent molecules in\n",
            closestWaters);
    mprintf("                           %s (%i).\n",P->parmName,P->solventMolecules);
    return 1;
  } 
  // Check that all solvent molecules contain same # atoms. Solvent 
  // molecules must be identical for the command to work properly; 
  // the prmtop strip occurs only once so the solvent params become fixed.
  NsolventAtoms = P->solventMoleculeStop[0] - P->solventMoleculeStart[0];
  for (solventMol = 1; solventMol < P->solventMolecules; solventMol++) {
    if ( NsolventAtoms != 
         (P->solventMoleculeStop[solventMol] - P->solventMoleculeStart[solventMol]) ) {
      mprintf("    Error: Closest::setup: Solvent molecules in %s are not of uniform size.\n",
              P->parmName);
      mprintf("           First solvent mol=%i atoms, %i solvent mol=%i atoms.\n",
              NsolventAtoms, solventMol,
              P->solventMoleculeStop[solventMol] - P->solventMoleculeStart[solventMol]);
      return 1;
    }
  }

  if ( Mask1.SetupMask(P,debug) ) return 1;
  if (Mask1.None()) {
    mprintf("    Error: Closest::setup: Mask %s contains no atoms.\n",Mask1.maskString);
    return 1;
  }

/*  if (P->mass==NULL && useMass) {
    fprintf(stdout,"    Warning: Closest::setup: Mass for this parm is NULL.\n");
    fprintf(stdout,"             Geometric center of mass will be used.\n");
    useMass=false;
  }*/

  // Figure out imaging - check box based on prmtop box
  // NOTE: Should box be figured out from read-in coords?
  imageType = 0;
  if (!noimage) {
    imageType = P->BoxType;
    if (P->BoxType==0) {
        mprintf("    Warning: Closest::setup: ");
        mprintf(" Imaging specified but no box information in prmtop %s\n",P->parmName);
        mprintf("             No imaging can occur..\n");
    }
  }

  // Figure out the total size of selected solvent molecules in atoms.
  // solventMoleculeStop is really the first atom of the next molecule.
  // Use temporary mask to keep original mask unmodified.
  if (tempMask!=NULL) delete tempMask;
  tempMask = Mask1.Copy();
  NsolventAtoms = P->solventMoleculeStop[closestWaters-1] - P->solventMoleculeStart[0];
  tempMask->Selected = (int*) realloc( tempMask->Selected, 
                                       (tempMask->Nselected + NsolventAtoms) * sizeof(int));
  mprintf("    CLOSEST: Keeping %i solvent atoms.\n",NsolventAtoms);
  // Put solvent atom #s at end of temporary array for creating stripped parm
  for ( solventMol=0; solventMol < closestWaters; solventMol++) {
    for ( solventAtom = P->solventMoleculeStart[solventMol];
          solventAtom < P->solventMoleculeStop[solventMol];
          solventAtom++ ) {
      //fprintf(stdout,"      Saving solvent atom %i\n",solventAtom+1);
      tempMask->Selected[tempMask->Nselected++] = solventAtom;
    }
  }

  // Store old parm
  oldParm = P;

  // Create integer masks for all solvent molecules in old parm.
  this->ClearMaskList();
  for (solventMol=0; solventMol < oldParm->solventMolecules; solventMol++) {
    // Setup solvent molecule mask
    SolventMask = new AtomMask();
    NsolventAtoms = oldParm->solventMoleculeStop[solventMol] -
                    oldParm->solventMoleculeStart[solventMol];
    SolventMask->Selected = (int*) malloc(NsolventAtoms * sizeof(int));
    SolventMask->Nselected = 0;
    for (solventAtom = oldParm->solventMoleculeStart[solventMol];
         solventAtom < oldParm->solventMoleculeStop[solventMol];
         solventAtom++) {
      SolventMask->Selected[SolventMask->Nselected++] = solventAtom;
    }
    MaskList.push_back(SolventMask);
  }

  // Create stripped Parm
  if (newParm!=NULL) delete newParm;
  newParm = P->modifyStateByMask(tempMask->Selected, tempMask->Nselected);
  if (newParm==NULL) {
    mprintf("    Error: Closest::setup: Could not create new parmtop.\n");
    return 1;
  }

  // Allocate space for new frame
  if (newFrame!=NULL) delete newFrame;
  newFrame = new Frame(newParm->natom, newParm->mass);

  // If prefix given then output stripped parm
  if (prefix!=NULL && newParm->parmName==NULL) {
    newParm->parmName=(char*)malloc((strlen(oldParm->parmName)+strlen(prefix)+2)*sizeof(char));
    sprintf(newParm->parmName,"%s.%s",prefix,oldParm->parmName);
    mprintf("             Writing out amber topology file %s\n",newParm->parmName);
    if ( newParm->WriteAmberParm() ) {
      mprintf("      Error: CLOSEST: Could not write out stripped parm file %s\n",
              newParm->parmName);
    }

  // Otherwise Set stripped parm name only, default prefix closest 
  } else if ( newParm->parmName==NULL ) {
    newParm->parmName=(char*)malloc((strlen(oldParm->parmName)+9)*sizeof(char));
    sprintf(newParm->parmName,"closest.%s",oldParm->parmName);
  }

  // Set parm
  P = newParm;

  return 0;  
}

/*
 * Closest::action()
 * Find the minimum distance between atoms in Mask1 and each solvent Mask.
 */
int Closest::action() {
  std::list<MolDist> Distances;
  MolDist moldist;
  int solventMol, solventAtom, maskPosition, atom;
  double minD, Dist, maxD, ucell[9], recip[9];

  if (imageType>0) F->BoxToRecip(ucell, recip);
  Distances.clear();
  // Calculate max distance
  maxD = F->box[0] + F->box[1] + F->box[2];
  maxD = maxD * maxD;
  // Loop over all solvent molecules in original frame
#ifdef _OPENMP
#pragma omp parallel private(solventMol,moldist,minD,atom,Dist,solventAtom)
{
#pragma omp for
#endif
  for (solventMol=0; solventMol < oldParm->solventMolecules; solventMol++) {
    moldist.mol = oldParm->firstSolvMol + solventMol;

    // DEBUG - show solvent mask
    //fprintf(stdout,"      Solvent %i %i %i\n", MaskList[solventMol]->Selected[0]+1,
    //        MaskList[solventMol]->Selected[1]+1,MaskList[solventMol]->Selected[2]+1);

    // Set initial target minimimum to be larger than any possible 
    // imaged distance.
    minD = maxD;

    // Calculate distance between each atom in Mask1 and atoms in solvent Mask
    for (atom = 0; atom < Mask1.Nselected; atom++) {
      Dist = F->DIST2(Mask1.Selected[atom], MaskList[solventMol]->Selected[0], imageType, ucell, recip);
      if (Dist < minD) minD=Dist;
      //fprintf(stdout,"D atom %i %i = %lf image %i\n",Mask1.Selected[atom],
      //        MaskList[solventMol]->Selected[0],minD,imageType);
      if (!firstAtom) {
        for (solventAtom=1; solventAtom<MaskList[solventMol]->Nselected; solventAtom++) {
          Dist = F->DIST2(Mask1.Selected[atom], MaskList[solventMol]->Selected[solventAtom],imageType,
                          ucell, recip);
          if (Dist < minD) minD=Dist;
        }
      }
    }
    moldist.D = minD; 

    // DEBUG - Print distances
    //fprintf(stdout,"Mol %8i minD= %lf\n",solventMol, moldist.D);
 
    moldist.mask = MaskList[solventMol];
    Distances.push_back(moldist);
  } // END for loop over solventMol
#ifdef _OPENMP
} // END pragma omp parallel
#endif

  // Sort distances
  Distances.sort( moldist_cmp() );

  // Add first closestWaters solvent atoms to tempMask
  solventMol = 0;
  maskPosition = Mask1.Nselected;
  for ( std::list<MolDist>::iterator it=Distances.begin();
        it != Distances.end();
        it++ ) {
    for (solventAtom = 0; solventAtom < (*it).mask->Nselected; solventAtom++)
      tempMask->Selected[maskPosition++] = (*it).mask->Selected[solventAtom];
    // Record which water molecules are closest if requested
    if (outFile!=NULL) {
      framedata->Add(Nclosest, &currentFrame);
      moldata->Add(Nclosest, &((*it).mol));
      Dist = sqrt( (*it).D );
      distdata->Add(Nclosest, &Dist);
      atomdata->Add(Nclosest, &((*it).mask->Selected[0]));
      Nclosest++;
    }
    // DEBUG - print first closestWaters distances
    //mprintf("DEBUG: Mol %i   D2= %lf   Atom0= %i\n",(*it).mol, (*it).D, (*it).mask->Selected[0]);
    solventMol++;
    if (solventMol==closestWaters) break;
  }

  // Modify and set frame
  newFrame->SetFrameFromMask(F, tempMask);
  F = newFrame;

  return 0;
} 

/*
 * Closest::print()
 * Write info on what solvent molecules were kept if specified.
 * NOTE: Get this into the master file list!
 */
void Closest::print() {
  if (outFile==NULL) return;
  // Sync up the dataset list here since it is not part of the master 
  // dataset list.
  outList->Sync();
  // Set specific datafile options
  outFile->SetNoXcol();
}

