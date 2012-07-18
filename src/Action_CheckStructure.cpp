// CheckStructure
#include <cmath>
#include <algorithm> // sort
#include "Action_CheckStructure.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
CheckStructure::CheckStructure() : 
  bondoffset_(1.0),
  nonbondcut2_(0.64), // 0.8^2
  isSeparate_(false)
{ 
  //fprintf(stderr,"CheckStructure Con\n");
  useImage_ = true;
}

// DESTRUCTOR
CheckStructure::~CheckStructure() {
  //fprintf(stderr,"CheckStructure Destructor.\n");
  outfile_.CloseFile();
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
  // Get Keywords
  useImage_ = !(actionArgs.hasKey("noimage"));
  char *reportFile = actionArgs.getKeyString("reportfile",NULL);
  bondoffset_ = actionArgs.getKeyDouble("offset",1.0);
  double nonbondcut = actionArgs.getKeyDouble("cut",0.8);

  // Get Masks
  char *mask1 = actionArgs.getNextMask();
  //fprintf(stdout,"    Mask 1: %s\n",mask1);
  Mask1_.SetMaskString(mask1);

  mprintf("    CHECKSTRUCTURE: Checking atoms in mask [%s]",Mask1_.MaskString());
  if (!useImage_) 
    mprintf(", imaging off");
  if (reportFile!=NULL)
    mprintf(", output to %s",reportFile);
  mprintf(".\n");
  mprintf("                    Warnings will be printed for bond length > eq + %.2lf\n",
          bondoffset_);
  mprintf("                    and non-bond distance < %.2lf\n",nonbondcut);
  nonbondcut2_ = nonbondcut * nonbondcut;

  if (outfile_.SetupWrite(reportFile, debug))
    return 1;
  outfile_.OpenFile();

  return 0;
}

// CheckStructure::SeparateInit()
void CheckStructure::SeparateInit(double bondoffsetIn, double nonbondcutIn, int debugIn) 
{
  double nonbondcut;
  isSeparate_ = true;
  debug = debugIn;
  useImage_ = false;
  if (bondoffsetIn < 0)
    bondoffset_ = 1.0;
  else
    bondoffset_ = bondoffsetIn;
  if (nonbondcutIn < 0)
    nonbondcut = 0.8;
  else 
    nonbondcut = nonbondcutIn;
  nonbondcut2_ = nonbondcut * nonbondcut;
  Mask1_.SetMaskString(NULL);
}

/** Set up bond arrays in a sorted list for easy access during loop
  * over all pairs of atoms. Only use bonds for which both atoms are in
  * the mask. It is expected that BndLst has 2 atom indices (i.e.
  * atom# * 3) followed by parameter index that starts from 1.
  */
void CheckStructure::SetupBondlist(std::vector<int> const& BndLst) {
  bond_list bnd;
  for (std::vector<int>::const_iterator bondatom = BndLst.begin();
                                        bondatom != BndLst.end();
                                        bondatom++) 
  {
    bnd.atom1 = *bondatom / 3;
    ++bondatom;
    if ( !Mask1_.AtomInCharMask(bnd.atom1) ) { ++bondatom; continue; }
    bnd.atom2 = *bondatom / 3;
    ++bondatom;
    if ( !Mask1_.AtomInCharMask(bnd.atom2) ) { continue; }
    bnd.param = *bondatom - 1; // Amber indices start from 1
    bondL_.push_back(bnd); 
  }
}

// CheckStructure::setup()
/** Determine what atoms each mask pertains to for the current parm file.
  * Also determine whether imaging should be performed. Check if parm
  * has bonds. If so, set up a list of the bonds between atoms in mask,
  * along with the expected bond lengths. Store bond lengths as 
  * (req + bondoffset)^2 for easy comparison with calculated distance^2.
  */
int CheckStructure::setup() {
  double req = 0;
  double rk = 0;
  unsigned int totalbonds;

  // Initially set up character mask to easily determine whether
  // both atoms of a bond are in the mask.
  if ( currentParm->SetupCharMask(Mask1_) ) return 1;
  if (Mask1_.None()) {
    mprintf("    Error: CheckStructure::setup: Mask has no atoms.\n");
    return 1;
  }

  // Check bonds
  bondL_.clear();
  SetupBondlist( currentParm->Bonds() );
  SetupBondlist( currentParm->BondsH() );
  if (bondL_.empty()) {
    mprintf("    Warning: No bond info in %s, will not check bonds.\n",currentParm->c_str());
    totalbonds = 0;
  } else {  
    // Since in the loop atom1 always < atom2, enforce this with parameters.
    // Sort by atom1, then by atom2
    sort( bondL_.begin(), bondL_.end(), bond_list_cmp() );
    // Fill in (req + offset)^2 values
    for (std::vector<bond_list>::iterator it = bondL_.begin(); it!=bondL_.end(); it++) {
      if ( currentParm->GetBondParamIdx((*it).param, rk, req) ) {
        // If no parameters exist for bond. Get cutoff from atom names.
        req = currentParm->GetBondedCutoff( (*it).atom1, (*it).atom2 );
      }
      req = (req + bondoffset_);
      req *= req;
      (*it).req = req;
    }
    // DEBUG
    if (debug>0) {
      mprintf("DEBUG:\tSorted bond list:\n");
      for (std::vector<bond_list>::iterator it = bondL_.begin(); it!=bondL_.end(); it++)
        mprintf("\t%8i %8i %8i %6.2lf\n",(*it).atom1,(*it).atom2,(*it).param, (*it).req);
    }
    totalbonds = bondL_.size();
  }
  // Insert a placeholder. Since atoms -1 -1 dont exist the bond will never
  // actually be accessed. If bond info is present this serves to indicate 
  // to the pair loop there are no more bonds to check. If no bond info is
  // present this serves to skip the bond check entirely.
  bond_list bnd; 
  bnd.atom1 = -1;
  bnd.atom2 = -1;
  bondL_.push_back(bnd);
      
  // Reset to integer mask.
  if ( currentParm->SetupIntegerMask( Mask1_ ) ) return 1;
  // Print imaging info for this parm
  if (!isSeparate_) {
    mprintf("    CHECKSTRUCTURE: %s (%i atoms, %u bonds)",Mask1_.MaskString(), Mask1_.Nselected(),
            totalbonds);
    if (imageType_ != Frame::NOIMAGE)
      mprintf(", imaging on");
    else
      mprintf(", imaging off");
    mprintf(".\n");
  }
        
  return 0;  
}

// CheckStructure::action()
int CheckStructure::action() {
  double ucell[9], recip[9], D2, D, bondmax;
  std::vector<bond_list>::iterator currentBond = bondL_.begin();

  if (imageType_==Frame::NONORTHO) currentFrame->BoxToRecip(ucell,recip);

  int lastidx = Mask1_.Nselected() - 1;
  for (int maskidx1 = 0; maskidx1 < lastidx; maskidx1++) {
    int atom1 = Mask1_[maskidx1];
    for (int maskidx2 = maskidx1 + 1; maskidx2 < Mask1_.Nselected(); maskidx2++) {
      int atom2 = Mask1_[maskidx2];
      // Get distance^2
      D2 = currentFrame->DIST2(atom1,atom2,imageType_,ucell,recip);
      // Are these atoms bonded?
      if ( (atom1==(*currentBond).atom1) && (atom2==(*currentBond).atom2) ) {
        // req has been precalced to (req + bondoffset)^2
        bondmax = (*currentBond).req;
        // Check for long bond length; distance2 > (req+bondoffset)^2
        if (D2 > bondmax) {
          D = sqrt(D2);
          outfile_.Printf(
                  "%i\t Warning: Unusual bond length %i@%s to %i@%s (%.2lf)\n",
                  frameNum+OUTPUTFRAMESHIFT,
                  atom1+1, (*currentParm)[atom1].c_str(), 
                  atom2+1, (*currentParm)[atom2].c_str(), D);
        }
        // Next bond
        currentBond++;
      // Atoms not bonded, check overlap
      } else {
        if (D2 < nonbondcut2_) {
          D = sqrt(D2);
          outfile_.Printf(
                  "%i\t Warning: Atoms %i@%s and %i@%s are close (%.2lf)\n",
                  frameNum+OUTPUTFRAMESHIFT,
                  atom1+1, (*currentParm)[atom1].c_str(), 
                  atom2+1, (*currentParm)[atom2].c_str(), D);
        }
      }
    } // END second loop over mask atoms
  } // END first loop over mask atoms

  return 0;
}

// CheckStructure::SeparateAction()
/// Used when you want to check all input coordinates.
int CheckStructure::SeparateAction(Frame *frameIn) {
  double D2, bondmax;
  std::vector<bond_list>::iterator currentBond = bondL_.begin();
  int Nproblems = 0;

  //if (imageType>0) frameIn->BoxToRecip(ucell,recip);

  int frame_natom = frameIn->Natom();
  int lastatom = frame_natom - 1;
  int idx1 = 0;
  for (int atom1 = 0; atom1 < lastatom; atom1++) {
    int idx2 = idx1 + 3;
    for (int atom2 = atom1 + 1; atom2 < frame_natom; atom2++) {
      // Get distance^2
      D2 = frameIn->COORDDIST2(idx1,idx2);
      // Are these atoms bonded?
      if ( (atom1==(*currentBond).atom1) && (atom2==(*currentBond).atom2) ) {
        // req has been precalced to (req + bondoffset)^2
        bondmax = (*currentBond).req;
        // Check for long bond length; distance2 > (req+bondoffset)^2
        if (D2 > bondmax) {
          Nproblems++; 
          if (debug>0)
            mprintf("\t\tWarning: Unusual bond length %i@%s to %i@%s (%.2lf)\n",
                    atom1+1, (*currentParm)[atom1].c_str(), 
                    atom2+1, (*currentParm)[atom2].c_str(), sqrt(D2));
        }
        // Next bond
        currentBond++;
      // Atoms not bonded, check overlap
      } else {
        if (D2 < nonbondcut2_) {
          Nproblems++;
          if (debug>0)
            mprintf("\t\tWarning: Atoms %i@%s and %i@%s are close (%.2lf)\n",
                    atom1+1, (*currentParm)[atom1].c_str(), 
                    atom2+1, (*currentParm)[atom1].c_str(), sqrt(D2));
        }
      }
      idx2+=3;
    } // END second loop over mask atoms
    idx1+=3;
  } // END first loop over mask atoms

  return Nproblems;
}

