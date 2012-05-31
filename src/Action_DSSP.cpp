#include "Action_DSSP.h"
#include "CpptrajStdio.h"
/// Hbond energy calc prefactor
// From ptraj actions.c:transformSecStruct 0.42*0.20*332
const double Action_DSSP::DSSP_fac = 27.888;

// CONSTRUCTOR
Action_DSSP::Action_DSSP() {
  outfilename_=NULL;
  dssp_=NULL; 
  Nres_=0;
  Nframe_=0;
  SSline_=0;
  printString_=false;
  SSdata_=0;
  dsspData_=0;
}

// DESTRUCTOR
Action_DSSP::~Action_DSSP() {
//  debugout.CloseFile(); // DEBUG
  if (SSline_!=0) delete[] SSline_;
  if (SSdata_!=0) delete SSdata_;
  if (dsspData_!=0) delete dsspData_;
}

const char Action_DSSP::SSchar[7]={ '0', 'b', 'B', 'G', 'H', 'I', 'T' };
const char Action_DSSP::SSname[7][6]={"None", "Para", "Anti", "3-10", "Alpha", "Pi", "Turn"};

// Action_DSSP::init()
/** Expected call: secstruct [out <filename>] [<mask>] [sumout <filename>]
  *                          [ptrajformat]
  * If sumout is not specified the filename specified by out is used with .sum suffix. 
  */
// For now dont allow NULL(stdout) filename for output
int Action_DSSP::init() {
  // DEBUG
//  debugout.SetupFile((char*)"dsspdebug.out",WRITE,UNKNOWN_FORMAT,STANDARD,0);
//  debugout.OpenFile();

  // Get keywords
  outfilename_ = actionArgs.getKeyString("out",NULL);
  char* temp = actionArgs.getKeyString("sumout",NULL);
  if (temp!=NULL) {
    sumOut.assign(temp);
  } else if (outfilename_!=NULL) {
    sumOut.assign(outfilename_);
    sumOut += ".sum";
  } 
  if (actionArgs.hasKey("ptrajformat")) printString_=true;
  // Get masks
  char* mask = actionArgs.getNextMask();
  Mask_.SetMaskString(mask);

  // Set up the DSSP data set
  if (printString_) {
    dssp_ = DSL->Add(DataSet::STRING, actionArgs.getNextString(),"DSSP");
    if (dssp_==NULL) return 1;
    DFL->Add(outfilename_, dssp_);
  } else
    SSdata_ = new DataSetList;

  mprintf( "    SECSTRUCT: Calculating secondary structure using mask [%s]\n",Mask_.MaskString());
  if (outfilename_!=NULL) 
    mprintf("               Dumping results to %s\n", outfilename_);
  if (!sumOut.empty())
    mprintf("               Sum results to %s\n",sumOut.c_str());
  if (printString_) 
    mprintf("               SS data for each residue will be stored as a string.\n");
  else
    mprintf("               SS data for each residue will be stored as integers.\n");

  return 0;
}

// Action_DSSP::setup()
/** Set up secondary structure calculation for all residues selected by the
  * mask expression. A residue is selected if at least one of the following 
  * atoms named "C   ", "O   ", "N   ", or "H   " (i.e. standard atom protein 
  * BB atom names) is selected. The residue is only initialized if it has not
  * been previously selected and set up by a prior call to setup.
  */
// NOTE: Currently relatively memory-intensive. Eventually set up so that SecStruct and
// CO_HN_Hbond members exist only for selected residues? Use Map?
int Action_DSSP::setup() {
  Residue RES;

  // Set up mask for this parm
  if ( currentParm->SetupIntegerMask( Mask_ ) ) return 1;
  if ( Mask_.None() ) {
    mprintf("Warning: DSSP: Mask has no atoms.\n");
    return 1;
  }

  // Initially mark all residues already set up as not selected and 
  // reset all atom numbers
  for (int res = 0; res < (int)SecStruct_.size(); res++) {
    SecStruct_[res].isSelected = false;
    SecStruct_[res].C = -1;
    SecStruct_[res].H = -1;
    SecStruct_[res].N = -1;
    SecStruct_[res].O = -1;
  }

  // Set up SecStruct for each solute residue
  Nres_ = currentParm->FinalSoluteRes();
  if (debug>0) mprintf("\tDSSP: Setting up for %i residues.\n",Nres_);

  // Set up for each residue of the current Parm if not already set-up.
  RES.sstype=SECSTRUCT_NULL;
  RES.isSelected=false;
  RES.C=-1;
  RES.O=-1;
  RES.N=-1;
  RES.H=-1;
  RES.CO_HN_Hbond.assign( Nres_, 0 );
  RES.SSprob[0]=0;
  RES.SSprob[1]=0;
  RES.SSprob[2]=0;
  RES.SSprob[3]=0;
  RES.SSprob[4]=0;
  RES.SSprob[5]=0;
  RES.SSprob[6]=0;
  RES.resDataSet = NULL;
  // Only resize SecStruct if current # residues > previous # residues
  if (Nres_ > (int) SecStruct_.size())
    SecStruct_.resize(Nres_, RES);

  // Go through all atoms in mask. Determine which residues have their C,
  // O, N, or H atoms selected. Store the actual coordinate index 
  // (i.e. atom# * 3) for use with COORDDIST routine.
  for (AtomMask::const_iterator atom = Mask_.begin(); atom!=Mask_.end(); ++atom) {
    int atom_res = (*currentParm)[*atom].ResNum();
    // If residue is out of bounds skip it
    if ( atom_res >= Nres_ ) continue;
    //fprintf(stdout,"DEBUG: Atom %i Res %i [%s]\n",*atom,atom_res,P->names[*atom]);
    SecStruct_[atom_res].isSelected = true;
    if (      (*currentParm)[*atom].Name() == "C   ")
      SecStruct_[atom_res].C = (*atom) * 3;
    else if ( (*currentParm)[*atom].Name() == "O   ")
      SecStruct_[atom_res].O = (*atom) * 3;
    else if ( (*currentParm)[*atom].Name() == "N   ")
      SecStruct_[atom_res].N = (*atom) * 3;
    else if ( (*currentParm)[*atom].Name() == "H   ")
      SecStruct_[atom_res].H = (*atom) * 3;
  }

  // For each residue selected in the mask, check if residue is already 
  // set up in SecStruct. If so update the atom indices, otherwise set it up.
  int selected = 0;
  for (int res = 0; res < Nres_; ++res) {
    if (!SecStruct_[res].isSelected) continue;
    // Residue needs at least C=O or N-H, Check?
    // Set up dataset if necessary 
    if (!printString_ && SecStruct_[res].resDataSet==NULL) {
      // Setup dataset name for this residue
      SecStruct_[res].resDataSet = SSdata_->AddMultiN(DataSet::INT,"",
                                                    currentParm->ResidueName(res),res+1);
      if (SecStruct_[res].resDataSet!=NULL) DFL->Add(outfilename_, SecStruct_[res].resDataSet);
    }
    ++selected;
  }

  // Count number of selected residues
  mprintf("\tMask [%s] corresponds to %i residues.\n",Mask_.MaskString(),selected);

  // Set up output buffer to hold string
  if (printString_) {
    if (SSline_!=0) delete[] SSline_;
    SSline_ = new char[ (2*selected) + 1 ];
    SSline_[(2*selected)]='\0';
  }

  // DEBUG - Print atom nums for each residue set up
//  for (res=0; res < Nres_; res++) {
//    if (SecStruct_[res].isSelected)
//      fprintf(stdout,"DEBUG: %i C=%i O=%i N=%i H=%i\n",res,SecStruct_[res].C/3,
//              SecStruct_[res].O/3, SecStruct_[res].N/3, SecStruct_[res].H/3);
//  }

  return 0;
}

// Action_DSSP::isBonded()
/** Return 1 if residue 1 CO bonded to residue 2 NH.
  * Ensure residue numbers are valid and residues are selected.
  */
int Action_DSSP::isBonded(int res1, int res2) {
  if (res1<0) return 0;
  if (res2<0) return 0;
  if (res1>=Nres_) return 0;
  if (res2>=Nres_) return 0;
  if (!SecStruct_[res1].isSelected) return 0;
  if (!SecStruct_[res2].isSelected) return 0;
  if ( SecStruct_[res1].CO_HN_Hbond[res2] ) return 1;
  return 0;
}

// Action_DSSP::SSassign()
/** Assign all residues from res1 to res2-1 the secondary structure sstype
  * only if not already assigned.
  * Assumes given residue range is valid.
  */
void Action_DSSP::SSassign(int res1, int res2, SStype typeIn) {
  for (int res=res1; res<res2; res++) {
    if (res==Nres_) break;
    if (!SecStruct_[res].isSelected) continue;
    if (SecStruct_[res].sstype==SECSTRUCT_NULL)
      SecStruct_[res].sstype=typeIn;
  }
}   
 
// Action_DSSP::action()
/** Determine secondary structure by hydrogen bonding pattern. */    
int Action_DSSP::action() {
  int resi, resj;
  int C, O, H, N;
  double rON, rCH, rOH, rCN, E;

  // Determine C=0 to H-N hydrogen bonds for each residue to each other residue
#ifdef _OPENMP
#pragma omp parallel private(resi,resj,C,O,H,N,rON, rCH, rOH, rCN, E)
{
#pragma omp for
#endif
  for (resi=0; resi < Nres_; resi++) {
    if (!SecStruct_[resi].isSelected) continue;
    // Reset previous SS assignment
    SecStruct_[resi].sstype=SECSTRUCT_NULL;
    SecStruct_[resi].CO_HN_Hbond.assign( Nres_, 0 );    

    if (SecStruct_[resi].C==-1 || SecStruct_[resi].O==-1) continue;
    C = SecStruct_[resi].C;
    O = SecStruct_[resi].O;
    for (resj=0; resj < Nres_; resj++) {
      if (!SecStruct_[resj].isSelected) continue;
// DEBUG
//      debugout.IO->Printf("\n%i Res%i-Res%i:",frameNum,resi,resj);
//      debugout.IO->Printf(" C=%i O=%i | N=%i H=%i:",C,O,SecStruct_[resj].N,SecStruct_[resj].H);
// DEBUG 
      if (resi==resj) continue;
      
      // NOTE: Should check all atoms here?
      if (SecStruct_[resj].H==-1 || SecStruct_[resj].N==-1) continue;
      N = SecStruct_[resj].N;
      H = SecStruct_[resj].H;

      rON = currentFrame->COORDDIST(O, N);
      rCH = currentFrame->COORDDIST(C, H);
      rOH = currentFrame->COORDDIST(O, H);
      rCN = currentFrame->COORDDIST(C, N);

      E = DSSP_fac * (1/rON + 1/rCH - 1/rOH - 1/rCN);
      if (E < -0.5)
        SecStruct_[resi].CO_HN_Hbond[resj] = 1;
//      if ( SecStruct_[resi].CO_HN_Hbond[resj] ) debugout.IO->Printf(" HBONDED!"); // DEBUG
    }
  }
#ifdef _OPENMP
} // END pragma omp parallel
#endif

  // Determine Secondary Structure based on Hbonding pattern.
  // In case of structural overlap, priority is given to the structure first in this list: 
  //   H, B, (E), G, I, T  (s. p. 2595 in the Kabsch & Sander paper)
  for (resi=0; resi < Nres_; resi++) {
    if (!SecStruct_[resi].isSelected) continue;
    // Alpha helices
    if ( isBonded( resi - 1, resi+3 ) && isBonded( resi, resi + 4) ) {
      SSassign(resi,resi+4,SECSTRUCT_ALPHA);
      continue;
    }

    // Beta sheets - only needed if SS not already assigned
    if ( SecStruct_[resi].sstype == SECSTRUCT_NULL ) {
      for (resj=0; resj < Nres_; resj++) {
        if (!SecStruct_[resj].isSelected) continue;
        // Only consider residues spaced more than 2 apart
        int abs_resi_resj = resi - resj;
        if (abs_resi_resj<0) abs_resi_resj = -abs_resi_resj;
        if (abs_resi_resj > 2) {
          // Parallel
          if ( (isBonded(resi-1, resj) && isBonded(resj, resi+1)) ||
               (isBonded(resj-1, resi) && isBonded(resi, resj+1)) ) {
            SecStruct_[resi].sstype = SECSTRUCT_PARA;
            break;
          // Anti-parallel
          } else if ( (isBonded(resi-1, resj+1) && isBonded(resj-1, resi+1)) ||
                      (isBonded(resi,   resj  ) && isBonded(resj,   resi  )) ) {
            SecStruct_[resi].sstype = SECSTRUCT_ANTI;
            break;
          }
        }
      }
      if (SecStruct_[resi].sstype!=SECSTRUCT_NULL) continue; 
    }

    // 3-10 helix
    if ( isBonded( resi - 1, resi+2 ) && isBonded( resi, resi + 3) ) {
      SSassign(resi,resi+3,SECSTRUCT_3_10);
      continue;
    } 

    // Pi helix
    if ( isBonded( resi - 1, resi+4 ) && isBonded( resi, resi + 5) ) {
      SSassign(resi,resi+5,SECSTRUCT_PI);
      continue;
    }
  } // End Initial SS assignment over all residues

  // Assign Turn structure
  for (resi=0; resi < Nres_; resi++) {
    if (!SecStruct_[resi].isSelected) continue;

    for (resj=5; resj > 2; resj--) {
//      fprintf(stdout,"DEBUG: %i Res %i and %i+%i are",frameNum,resi,resi,resj);
//      if ( isBonded( resi, resi+resj) )
//        fprintf(stdout," BONDED!\n");
//      else
//        fprintf(stdout," not bonded.\n");
      if ( isBonded( resi, resi+resj) ) {
        SSassign(resi+1, resi+resj, SECSTRUCT_TURN);
        break;
      }
    }
  }

  // Store data for each residue 
  //fprintf(stdout,"%10i ",frameNum);
  // String data set
  if (printString_) {
    resj = 0;
    for (resi=0; resi < Nres_; resi++) {
      if (!SecStruct_[resi].isSelected) continue;
      SecStruct_[resi].SSprob[SecStruct_[resi].sstype]++;
      SSline_[resj++] = SSchar[SecStruct_[resi].sstype];
      SSline_[resj++] = ' ';
    }
    dssp_->Add(frameNum, SSline_);
  // Integer data sets
  } else {
    for (resi=0; resi < Nres_; resi++) {
      if (!SecStruct_[resi].isSelected) continue;
      //fprintf(stdout,"%c",SSchar[SecStruct_[resi].sstype]);
      SecStruct_[resi].SSprob[SecStruct_[resi].sstype]++;
      SecStruct_[resi].resDataSet->Add(frameNum, &(SecStruct_[resi].sstype));
    }
  }
  //fprintf(stdout,"\n");
  ++Nframe_;

  return 0;
}

// Action_DSSP::print()
/** Calculate the average of each secondary structure type across all residues.
  * Prepare for output via the master data file list.
  */
void Action_DSSP::print() {
  DataFile *dsspFile;
  int resi, ss;
  double avg;

  // If not printing a string, sync the integer dataset here since
  // it is not part of the master dataset list.
  if (!printString_ && SSdata_!=0) SSdata_->Sync();

  if (sumOut.empty()) return;
  // Set up dataset list to store averages
  dsspData_ = new DataSetList(); 
  // Set up a dataset for each SS type
  for (ss=1; ss<7; ss++) 
    dsspFile = DFL->Add((char*)sumOut.c_str(), 
                        dsspData_->Add(DataSet::DOUBLE, (char*)SSname[ss], "SS") );
  // Change the X label to Residue
  dsspFile->ProcessArgs("xlabel Residue");
  // Dont print empty frames
  // NOTE: Obsolete?
  dsspFile->ProcessArgs("noemptyframes");

  // Calc the avg structure of each type for each selected residue 
  for (resi=0; resi < Nres_; resi++) {
    if (!SecStruct_[resi].isSelected) continue;
    for (ss=1; ss<7; ss++) {
      avg = (double)SecStruct_[resi].SSprob[ss];
      avg /= (double)Nframe_;
      DataSet *tempDS = dsspData_->GetDataSetN(ss-1);
      tempDS->Add(resi, &avg);
    }
  }
}
