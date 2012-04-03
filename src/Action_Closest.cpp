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
#include "ParmFile.h"

// CONSTRUCTOR
Closest::Closest() {
  //fprintf(stderr,"Closest Con\n");
  useImage=true;
  firstAtom=false;
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
  useImage = !(actionArgs.hasKey("noimage"));
  prefix = actionArgs.getKeyString("outprefix",NULL);
  // Setup output file and sets if requested.
  // Will keep track of Frame, Mol#, Distance, and first solvent atom
  mask1 = actionArgs.getKeyString("closestout",NULL);
  if (mask1 != NULL) {
    // Set up datasets
    outList = new DataSetList();
    framedata = outList->Add(DataSet::INT,(char*)"Frame\0","Frame");
    moldata   = outList->Add(DataSet::INT,(char*)"Mol\0","Mol");
    distdata  = outList->Add(DataSet::DOUBLE,(char*)"Dist\0","Dist");
    atomdata  = outList->Add(DataSet::INT,(char*)"FirstAtm\0","FirstAtm");
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
  if (!useImage) 
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
  MolDist solvent;

  // If there are no solvent molecules this action is not valid.
  if (currentParm->Nsolvent()==0) {
    mprintf("Warning: Closest::setup: Parm %s does not contain solvent.\n",currentParm->c_str());
    return 1;
  }
  // If # solvent to keep >= solvent in this parm the action is not valid.
  // TODO: Just use max # waters?
  if (closestWaters >= currentParm->Nsolvent()) {
    mprintf("Warning: Closest::setup: # solvent to keep (%i) >= # solvent molecules in\n",
            closestWaters);
    mprintf("                           %s (%i).\n",currentParm->c_str(),
            currentParm->Nsolvent());
    return 1;
  } 
  // Check that all solvent molecules contain same # atoms. Solvent 
  // molecules must be identical for the command to work properly; 
  // the prmtop strip occurs only once so the solvent params become fixed.
  Topology::mol_iterator solvmol = currentParm->SolventStart();
  int NsolventAtoms = (*solvmol).NumAtoms();
  ++solvmol;
  for (; solvmol != currentParm->SolventEnd(); solvmol++) {
    if ( NsolventAtoms != (*solvmol).NumAtoms() ) {
      mprinterr("Error: Closest::setup: Solvent molecules in %s are not of uniform size.\n",
                currentParm->c_str());
      mprinterr("       First solvent mol = %i atoms, this solvent mol = %i atoms.\n",
                NsolventAtoms, (*solvmol).NumAtoms());
      return 1;
    }
  }

  // Setup solute atom mask
  // NOTE: Should ensure that no solvent atoms are selected!
  if ( currentParm->SetupIntegerMask(soluteMask) ) return 1;
  if (soluteMask.None()) {
    mprintf("Warning: Closest::setup: Mask %s contains no atoms.\n",soluteMask.MaskString());
    return 1;
  }

  // Figure out what the the total size of the selected solute atoms plus
  // the number of kept solvent atoms is in order to set up the stripped
  // parm.
  // solventMoleculeStop is really the first atom of the next molecule.
  // Want to keep the solute in addition to the closest solvent.
  stripMask = soluteMask;
  //mprintf("DEBUG:\t");
  //stripMask.PrintMaskAtoms("stripMask");
  // Since we have ensured all solvent atoms are of uniform size this is just
  // the desired number of kept waters * number solvent atoms in each mol.
  NsolventAtoms *= closestWaters;
  mprintf("\tKeeping %i solvent atoms.\n",NsolventAtoms);
  // Add space for kept solvent atom #s at end of mask.
  stripMask.AddAtomRange( stripMask.Nselected(), stripMask.Nselected() + NsolventAtoms );
  //stripMask.AddAtomRange(currentParm->solventMoleculeStart[0], 
  //                       currentParm->solventMoleculeStop[closestWaters-1]);
  //mprintf("DEBUG:\t");
  //stripMask.PrintMaskAtoms("stripMaskWsolvent");

  // Store old parm
  oldParm = currentParm;
 
  // Get atom masks for all solvent molecules in old parm.
  SolventMols.clear();
  solvent.D=0.0;
  solvent.mol=0;
  SolventMols.resize(oldParm->Nsolvent(), solvent);
  int solventMolNum = oldParm->FirstSolventMol() + 1;
  std::vector<MolDist>::iterator mdist = SolventMols.begin();
  for (Topology::mol_iterator solvmol = oldParm->SolventStart(); 
                              solvmol!= oldParm->SolventEnd(); solvmol++) 
  {
    (*mdist).mol = solventMolNum++;
    // NOTE: EndAtom - 1?
    (*mdist).mask.AddAtomRange( (*solvmol).BeginAtom(), (*solvmol).EndAtom() );
    //SolventMols[solventMol].mask.PrintMaskAtoms("solvent");
    ++mdist;
  }

  // Create stripped Parm
  if (newParm!=NULL) delete newParm;
  newParm = currentParm->modifyStateByMask(stripMask);
  if (newParm==NULL) {
    mprinterr("Error: Closest::setup: Could not create new parmtop.\n");
    return 1;
  }
  newParm->Summary();

  // Allocate space for new frame
  newFrame.SetupFrame(newParm->Natom(), newParm->Mass());

  // If prefix given then output stripped parm
  if (prefix!=NULL) {
    std::string newfilename(prefix);
    newfilename += ".";
    newfilename += oldParm->OriginalFilename();
    mprintf("\tWriting out amber topology file %s to %s\n",newParm->c_str(),newfilename.c_str());
    ParmFile pfile;
    pfile.SetDebug( debug );
    if ( pfile.Write(*newParm, newfilename, ParmFile::AMBERPARM ) ) {
    //if ( newParm->WriteTopology(newParm->parmName) ) {
      mprinterr("Error: CLOSEST: Could not write out stripped parm file %s\n",
              newParm->c_str());
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
  int solventMol, maskPosition, maxSolventMolecules;
  double Dist, maxD, ucell[9], recip[9];
  AtomMask::const_iterator solute_atom, solvent_atom;

  if (imageType != NOIMAGE) {
    currentFrame->BoxToRecip(ucell, recip);
    // Calculate max possible imaged distance
    maxD = currentFrame->MaxImagedDistance();
  } else {
    // If not imaging, set max distance to an arbitrarily large number
    maxD = DBL_MAX;
  }

  // Loop over all solvent molecules in original frame
  maxSolventMolecules = oldParm->Nsolvent();
  // DEBUG
  //mprintf("Closest: Begin parallel loop for %i\n",frameNum);
  // DEBUG
#ifdef _OPENMP
#pragma omp parallel private(solventMol,solute_atom,Dist,solvent_atom)
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
    solvent_atom = SolventMols[solventMol].mask.begin();
    for (solute_atom = soluteMask.begin(); solute_atom != soluteMask.end(); solute_atom++)
    {
      Dist = currentFrame->DIST2(*solute_atom, *solvent_atom, imageType, ucell, recip);
      if (Dist < SolventMols[solventMol].D) 
        SolventMols[solventMol].D = Dist;
      //fprintf(stdout,"D atom %i %i = %lf image %i\n",soluteMask.Selected[atom],
      //        MaskList[solventMol]->Selected[0],minD,imageType);
      // Check the rest of the solvent atoms if specified
      if (!firstAtom) {
        ++solvent_atom;
        for (; solvent_atom != SolventMols[solventMol].mask.end(); solvent_atom++) 
        {
          Dist = currentFrame->DIST2(*solute_atom, *solvent_atom, imageType, ucell, recip);
          if (Dist < SolventMols[solventMol].D) 
            SolventMols[solventMol].D = Dist;
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
  maskPosition = soluteMask.Nselected();
  for ( std::vector<MolDist>::iterator solvent = SolventMols.begin();
                                       solvent != SolventMols.end();
                                       solvent++ ) 
  {
    //mprintf("DEBUG:\tmol %i ",(*solvent).mol);
    //(*solvent).mask.PrintMaskAtoms("Mask");
    stripMask.AddMaskAtPosition( (*solvent).mask, maskPosition );

    // Record which water molecules are closest if requested
    if (outFile!=NULL) {
      framedata->Add(Nclosest, &frameNum);
      moldata->Add(Nclosest, &((*solvent).mol));
      Dist = sqrt( (*solvent).D );
      distdata->Add(Nclosest, &Dist);
      solvent_atom = (*solvent).mask.begin();
      int solvent_first_atom = *solvent_atom; 
      atomdata->Add(Nclosest, &solvent_first_atom);
      Nclosest++;
    }
    // DEBUG - print first closestWaters distances
    //mprintf("DEBUG: Mol %i   D2= %lf   Atom0= %i\n",(*it).mol, (*it).D, (*it).mask->Selected[0]);
    ++solventMol;
    if (solventMol==closestWaters) break;
  }

  // Modify and set frame
  //mprintf("DEBUG:\t");
  //stripMask.PrintMaskAtoms("action_stripMask");
  newFrame.SetFrame(*currentFrame, stripMask);
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
  outFile->ProcessArgs("noxcol");
}

