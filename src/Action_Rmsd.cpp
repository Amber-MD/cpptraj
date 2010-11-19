#include <cstdlib>
#include "Action_Rmsd.h"
using namespace std;

// CONSTRUCTOR
Rmsd::Rmsd() {
  rmsd=NULL;
  PerResRMSD=NULL;
  nres=0;
  ResRange=NULL;
  RefRange=NULL;
  perresout=NULL;
  nofit=false;
  first=false;
  useMass=false;
  perres=false;
  RefFrame=NULL;
  RefParm=NULL;
  SelectedRef=NULL;
  SelectedFrame=NULL;
  perresmask=NULL;
  perrescenter=false;
  perresinvert=false;
  //currentType=RMSD;
}

// DESTRUCTOR
Rmsd::~Rmsd() {
  AtomMask *ResMask;
  // If first, ref Frame was allocd (not assigned from reference Frame List)
  if (first && RefFrame!=NULL)
    delete RefFrame;
  if (SelectedRef!=NULL) delete SelectedRef;
  if (SelectedFrame!=NULL) delete SelectedFrame;
  if (PerResRMSD!=NULL) delete PerResRMSD;
  while ( !PerResMask.empty() ) {
    ResMask = PerResMask.back();
    delete ResMask;
    PerResMask.pop_back();
  }
  while ( !PerRefMask.empty() ) {
    ResMask = PerRefMask.back();
    delete ResMask;
    PerRefMask.pop_back();
  }
  if (ResRange!=NULL) delete ResRange;
  if (RefRange!=NULL) delete RefRange;
}

/*
 * Rmsd::SetRefMask()
 * Setup reference mask based on maskRef. Requires RefParm to be set. Should 
 * only be called once.
 * If reference, this is called from init. If first, this is called from setup.
 */
int Rmsd::SetRefMask() {
  if ( RefMask.SetupMask(RefParm,debug) ) return 1;
  if (RefMask.None()) {
    fprintf(stdout,"    Error: Rmsd::SetRefMask: No atoms in reference mask.\n");
    return 1;
  }
  // Allocate frame for selected reference atoms
  SelectedRef = new Frame(RefMask.Nselected, NULL);

  return 0;
}

/*
 * Rmsd::init()
 * Called once before traj processing. Set up reference info.
 * Expected call: 
 * rmsd <name> <mask> [first | ref <filename> | refindex <#>] [<refmask>] [out filename] [nofit]
 *      [perres] [perresout <filename>] [range <res range>] [perresmask <addtl mask>]
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int Rmsd::init( ) {
  char *referenceName, *mask0, *maskRef;
  char *rmsdFile, *range, *refrange;
  int refindex, referenceKeyword;

  // Check for keywords
  referenceKeyword=A->hasKey("reference"); // For compatibility with ptraj
  referenceName=A->getKeyString("ref",NULL);
  refindex=A->getKeyInt("refindex",-1);
  nofit = A->hasKey("nofit");
  first = A->hasKey("first");
  useMass = A->hasKey("mass");
  perres = A->hasKey("perres");
  rmsdFile = A->getKeyString("out",NULL);
  // Per-res keywords
  perresout = A->getKeyString("perresout",NULL);
  perresinvert = A->hasKey("perresinvert");
  range = A->getKeyString("range",NULL);
  ResRange = A->NextArgToRange(range);
  refrange = A->getKeyString("refrange",NULL);
  if (refrange!=NULL) 
    RefRange = A->NextArgToRange(refrange);
  perresmask = A->getKeyString("perresmask",(char*)"");
  perrescenter = A->hasKey("perrescenter");

  // Get the RMS mask string for frames
  mask0 = A->getNextMask();
  FrameMask.SetMaskString(mask0);
  // Get RMS mask string for reference
  maskRef = A->getNextMask();
  // If no reference mask specified, make same as RMS mask
  if (maskRef==NULL) maskRef=mask0; 
  RefMask.SetMaskString(maskRef);

  // Set up the RMSD data set
  rmsd = DSL->Add(DOUBLE, A->getNextString());
  if (rmsd==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(rmsdFile,rmsd);

  if (!first && referenceName==NULL && refindex==-1 && referenceKeyword==0) {
    fprintf(stdout,
            "    Warning: Rmsd::init: No reference structure given. Defaulting to first.\n");
    first=true;
  }

  if (!first) {
    // Attempt to get reference index by name
    if (referenceName!=NULL)
      refindex=FL->GetFrameIndex(referenceName);

    // For compatibility with ptraj, if 'reference' specified use first reference structure
    if (referenceKeyword) refindex=0;

    // Get reference frame by index
    RefFrame=FL->GetFrame(refindex);
    if (RefFrame==NULL) {
      fprintf(stdout,"    Error: Rmsd::init: Could not get reference index %i\n",refindex);
      return 1;
    }
    // Set reference parm
    RefParm=FL->GetFrameParm(refindex);
    // Setup reference mask here since reference frame/parm are allocated
    if ( SetRefMask() ) return 1;
    //RefFrame->printAtomCoord(0);
    //fprintf(stderr,"  NATOMS IN REF IS %i\n",RefFrame->P->natom); // DEBUG
  }

  //rmsd->Info();
  fprintf(stdout,"    RMSD: (%s), reference is ",FrameMask.maskString);
  if (RefFrame==NULL)
    fprintf(stdout,"first frame");
  else if (referenceName!=NULL)
    fprintf(stdout,"%s",referenceName);
  else
    fprintf(stdout,"reference index %i",refindex);
  fprintf(stdout," (%s)",RefMask.maskString);
  if (nofit)
    fprintf(stdout,", no fitting");
  else
    fprintf(stdout,", with fitting");
  fprintf(stdout,".\n");
  if (useMass) 
    fprintf(stdout,"          Center of mass will be used.\n");
  else
    fprintf(stdout,"          Geometric center will be used.\n");
  if (perres) {
    fprintf(stdout,"          No-fit RMSD will also be calculated for ");
    if (range==NULL) 
      fprintf(stdout,"each solute residue");
    else
      fprintf(stdout,"residues %s",range);
    if (refrange!=NULL)
      fprintf(stdout," (reference residues %s)",refrange);
    fprintf(stdout," using mask [:X%s].\n",perresmask);
    fprintf(stdout,"          WARNING: Currently residues are set up based on the first trajectory read in!\n");
    if (perresout==NULL) {
      fprintf(stdout,"Error: perres specified but no output filename given (perresout).\n");
      perres=false;
      return 1;
    }
    fprintf(stdout,"          Per-residue output file is %s\n",perresout);
    if (perrescenter)
      fprintf(stdout,"          perrescenter: Each residue will be centered prior to RMS calc.\n");
    if (perresinvert)
      fprintf(stdout,"          perresinvert: Frames will be written in rows instead of columns.\n");
  }

  return 0;
}

/*
 * Rmsd::setup()
 * Called every time the trajectory changes. Set up FrameMask for the new parmtop
 * and allocate space for selected atoms from the Frame.
 */
int Rmsd::setup() {
  char resArg[1024];
  char refArg[1024];
  char resName[5];
  AtomMask *ResMask;
  list<int>::iterator it;
  //int RefNselected;
  int refRes;;

  if ( FrameMask.SetupMask(P,debug) ) return 1;
  if ( FrameMask.None() ) {
    fprintf(stdout,"    Error: Rmsd::setup: No atoms in mask.\n");
    return 1;
  }
  // Allocate space for selected atoms in the frame
  // NOTE: Should masses be passed in to trigger the Alloc even though the vals will
  //       be wrong initially until set by SetFrameFromMask in action?
  if (SelectedFrame!=NULL) delete SelectedFrame;
  SelectedFrame = new Frame(FrameMask.Nselected, NULL);
  
  // first: If RefParm not set, set it here and set the reference mask.
  //        Should only occur once.
  if (first && RefParm==NULL) {
    RefParm=P;
    if ( SetRefMask() ) return 1;
  } 

  // Check that num atoms in frame mask from this parm match ref parm mask
  if ( RefMask.Nselected != FrameMask.Nselected ) {
    fprintf(stdout, "    Error: Number of atoms in RMS mask (%i) does not \n",FrameMask.Nselected);
    fprintf(stdout, "           equal number of atoms in Ref mask (%i).\n",RefMask.Nselected);
    return 1;
  }

  // Check Mass
  if (useMass && P->mass==NULL) {
    fprintf(stdout,"    Warning: Rmsd::setup: usemass: Parmtop %s does not contain mass info.\n",
            P->parmName);
    fprintf(stdout,"             Geometric center will be used instead.\n");
    useMass=false;
  }

  // ---===  PER RESIDUE RMSD OPTION ===---
  // If perres was specified, need a data set for each residue
  // Currently perres will only work for the first parmtop used
  // NOTE THAT ALL RESIDUES FROM INPUT SHOULD BE SHIFTED BY -1
  if (PerResRMSD==NULL && perres) {
    // If no range specified do all solute residues.
    if (ResRange==NULL) {
      if (P->finalSoluteRes>0)
        nres = P->finalSoluteRes;
      else
        nres=P->nres;
      // Has to match user input where resnums start from 1
      // Since ArgToRange returns up to and including, want 1-nres
      sprintf(resArg,"%i-%i",1,nres);
      ResRange = A->NextArgToRange(resArg);
    } 
    fprintf(stdout,"      RMSD: PerRes: Setting up for %i residues.\n",(int)ResRange->size());
    PerResRMSD=new DataSetList();
    resName[4]='\0';
    for (it=ResRange->begin(); it!=ResRange->end(); it++) {
      // Get corresponding reference resnum - if none specified use current res
      if (RefRange==NULL) 
        refRes = (*it);
      else {
        refRes = RefRange->front();
        RefRange->pop_front();
      }
      //res = *it - 1; // res is the internal resnumber, *it the user resnumber
      // Setup Dataset Name to be name of this residue 
      P->ResName(resName,(*it)-1);
      
      // Setup mask strings - masks are based off user resnums
      sprintf(resArg,":%i%s",*it,perresmask);
      sprintf(refArg,":%i%s",refRes,perresmask);
      //fprintf(stdout,"DEBUG: RMSD: PerRes: Mask %s RefMask %s\n",resArg,refArg);

      // Check that mask can be set up for reference
      ResMask=new AtomMask();
      ResMask->SetMaskString(refArg);
      if ( ResMask->SetupMask(RefParm,0) ) {
        fprintf(stdout,"Warning: RMSD: PerRes: Could not set up reference for residue %i\n",*it);
        delete ResMask;
        continue;
      }
      if (ResMask->None()) {
        fprintf(stdout,"Warning: RMSD: PerRes: No atoms in reference for residue %i\n",*it);
        delete ResMask;
        continue;
      }
      //RefNselected = RefMask->Nselected;
      PerRefMask.push_back(ResMask);

      // Set up mask for this parm - If unable make sure reference mask is popped as well
      ResMask=new AtomMask();
      ResMask->SetMaskString(resArg);
      if ( ResMask->SetupMask(P,0) ) { // NOTE: Allow debug value in here?
        fprintf(stdout,"Warning: RMSD: PerRes: Could not set up mask for residue %i\n",*it);
        delete ResMask;
        ResMask = PerRefMask.back();
        delete ResMask;
        PerRefMask.pop_back();
        continue;
      }
      if (ResMask->None()) {
        fprintf(stdout,"Warning: RMSD: PerRes: No atoms in mask for residue %i\n",*it);
        delete ResMask;
        ResMask = PerRefMask.back();
        delete ResMask;
        PerRefMask.pop_back();
        continue;
      }
      // Check that number of atoms selected in parm is same as reference
      if ( (PerRefMask.back())->Nselected != ResMask->Nselected) {
        fprintf(stdout,"Warning: RMSD: PerRes: # atoms in mask for residue %i (%i) not equal\n",
                *it,ResMask->Nselected);
        fprintf(stdout,"                       to # atoms in reference mask (%i).\n",
                (PerRefMask.back())->Nselected);
        delete ResMask;
        ResMask = PerRefMask.back();
        delete ResMask;
        PerRefMask.pop_back();
        continue;
      }
      PerResMask.push_back(ResMask);

      // DEBUG
      //fprintf(stdout,"PERRES_RMS: Mask %s RefMask %s\n",(PerResMask.back())->maskString,
      //        (PerRefMask.back())->maskString);

      // Setup Dataset for this residue
      sprintf(resArg,"%4s%i",resName,*it);
      // TEST - add all datasets to the same output file
      // NOTE - eventually give this its own output file and make a print routine
      DFL->Add(perresout,PerResRMSD->Add(DOUBLE, resArg));
    }
    // Set output file to be inverted if requested
    if (perresinvert) {
      DataFile *Current = DFL->GetDataFile(perresout);
      if (Current!=NULL)
        Current->SetInverted();
    }
  }

  return 0;
}

/*
 * Rmsd::action()
 * Called every time a frame is read in. Calc RMSD. If RefFrame is NULL
 * at this point we want to set the first frame read in as reference.
 */
int Rmsd::action() {
  double R, U[9], Trans[6];
  vector<AtomMask*>::iterator it;
  int res;

  // first: If Ref is NULL, allocate this frame as reference
  //        Should only occur once.
  // NOTE: For MPI this will currently result in different references between threads.
  if (first && RefFrame==NULL) 
    RefFrame = F->Copy();

  // Set selected reference atoms - always done since RMS fit modifies SelectedRef 
  SelectedRef->SetFrameFromMask(RefFrame, &RefMask);

  // Set selected frame atoms
  SelectedFrame->SetFrameFromMask(F, &FrameMask);

  // DEBUG
/*  fprintf(stdout,"  DEBUG: RMSD: First atom coord in SelectedFrame is : "); 
  SelectedFrame->printAtomCoord(0);
  fprintf(stdout,"  DEBUG: RMSD: First atom coord in SelectedRef is : ");
  SelectedRef->printAtomCoord(0);
*/

  if (nofit) {
    R = SelectedFrame->RMSD(SelectedRef, useMass);
  } else {
    R = SelectedFrame->RMSD(SelectedRef, U, Trans, useMass);
    F->Translate(Trans);
    F->Rotate(U);
    F->Translate(Trans+3);
  }

  rmsd->Add(currentFrame, &R);

  // Per Residue RMSD - Set reference and selected frame for each mask in PerResMask
  if (perres) {
    res=0;
    for (it = PerResMask.begin(); it!=PerResMask.end(); it++) {
      SelectedRef->SetFrameFromMask(RefFrame, PerRefMask[res]);
      SelectedFrame->SetFrameFromMask(F, (*it));
      if (perrescenter) 
        SelectedFrame->ShiftToCenter(SelectedRef); 
      R = SelectedFrame->RMSD(SelectedRef,useMass);
      //fprintf(stdout,"DEBUG:           Res %i nofit RMSD = %lf\n",res,R);
      // NOTE: Should check for error on AddData?
      PerResRMSD->AddData(currentFrame, &R, res);
      res++;
    }
  }
  

  return 0;
}

