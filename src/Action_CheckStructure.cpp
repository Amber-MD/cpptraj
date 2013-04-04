// Action_CheckStructure
#include <cmath>
#include <algorithm> // sort
#include "Action_CheckStructure.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_CheckStructure::Action_CheckStructure() : 
  bondoffset_(1.15),
  nonbondcut2_(0.64), // 0.8^2
  bondcheck_(true),
  silent_(false),
  CurrentParm_(0),
  debug_(0)
{} 

void Action_CheckStructure::Help() {
  mprintf("\t[<mask1>] [reportfile <report>] [noimage]\n");
  mprintf("\t[offset <offset>] [cut <cut>] [nobondcheck]\n");
  mprintf("\tCheck frames for atomic overlaps and unusual bond lengths\n");
}

// DESTRUCTOR
Action_CheckStructure::~Action_CheckStructure() {
  //fprintf(stderr,"CheckStructure Destructor.\n");
  outfile_.CloseFile();
}

// Action_CheckStructure::init()
Action::RetType Action_CheckStructure::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  // Get Keywords
  InitImaging( !(actionArgs.hasKey("noimage")) );
  std::string reportFile = actionArgs.GetStringKey("reportfile");
  bondoffset_ = actionArgs.getKeyDouble("offset",1.15);
  double nonbondcut = actionArgs.getKeyDouble("cut",0.8);
  bondcheck_ = !actionArgs.hasKey("nobondcheck");
  // Hidden option, for use when used by other actions
  silent_ = actionArgs.hasKey("silent");

  // Get Masks
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );

  mprintf("    CHECKSTRUCTURE: Checking atoms in mask [%s]",Mask1_.MaskString());
  if (!UseImage()) 
    mprintf(", imaging off");
  if (!reportFile.empty())
    mprintf(", output to %s",reportFile.c_str());
  mprintf(".\n");
  if (!bondcheck_) {
    mprintf("\tChecking inter-atom distances only.\n");
    mprintf("\tWarnings will be printed for non-bond distances < %.2f Ang.\n", nonbondcut);
  } else {
    mprintf("\tChecking inter-atom and bond distances.\n");
    mprintf("\tWarnings will be printed for bond lengths > eq + %.2f Ang\n",
            bondoffset_);
    mprintf("\tand non-bond distances < %.2f Ang.\n",nonbondcut);
  }
  nonbondcut2_ = nonbondcut * nonbondcut;

  if (outfile_.SetupWrite(reportFile, debug_))
    return Action::ERR;
  outfile_.OpenFile();

  return Action::OK;
}

/** Set up bond arrays in a sorted list for easy access during loop
  * over all pairs of atoms. Only use bonds for which both atoms are in
  * the mask. It is expected that BndLst has 2 atom indices (i.e.
  * atom# * 3) followed by parameter index that starts from 1.
  */
void Action_CheckStructure::SetupBondlist(std::vector<int> const& BndLst) {
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
  if (bondL_.empty())
    mprintf("Warning: No bond info in parm %s, will not check bonds.\n",CurrentParm_->c_str());
}

// Action_CheckStructure::setup()
/** Determine what atoms each mask pertains to for the current parm file.
  * Also determine whether imaging should be performed. Check if parm
  * has bonds. If so, set up a list of the bonds between atoms in mask,
  * along with the expected bond lengths. Store bond lengths as 
  * (req + bondoffset)^2 for easy comparison with calculated distance^2.
  */
Action::RetType Action_CheckStructure::Setup(Topology* currentParm, Topology** parmAddress) {
  double req = 0;
  double rk = 0;
  unsigned int totalbonds;

  CurrentParm_ = currentParm;      
  // Initially set up character mask to easily determine whether
  // both atoms of a bond are in the mask.
  if ( currentParm->SetupCharMask(Mask1_) ) return Action::ERR;
  if (Mask1_.None()) {
    mprintf("    Error: CheckStructure::setup: Mask has no atoms.\n");
    return Action::ERR;
  }

  SetupImaging( currentParm->BoxType() );

  // Check bonds
  bondL_.clear();
    totalbonds = 0;
  if (bondcheck_) {
    SetupBondlist( currentParm->Bonds() );
    SetupBondlist( currentParm->BondsH() );
  }
  if (!bondL_.empty()) {
    // Since in the loop atom1 always < atom2, enforce this with parameters.
    // Sort by atom1, then by atom2
    std::sort( bondL_.begin(), bondL_.end(), bond_list_cmp() );
    // Fill in (req + offset)^2 values
    for (std::vector<bond_list>::iterator it = bondL_.begin(); it!=bondL_.end(); it++) {
      if ( currentParm->GetBondParamIdx((*it).param, rk, req) ) {
        // If no parameters exist for bond. Get cutoff from atom names.
        req = currentParm->GetBondedCutoff( (*it).atom1, (*it).atom2 );
      }
      req += bondoffset_;
      req *= req;
      (*it).req = req;
    }
    // DEBUG
    if (debug_>0) {
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
  if ( currentParm->SetupIntegerMask( Mask1_ ) ) return Action::ERR;
  // Print imaging info for this parm
  mprintf("    CHECKSTRUCTURE: %s (%i atoms, %u bonds)",Mask1_.MaskString(), Mask1_.Nselected(),
          totalbonds);
  if (ImagingEnabled())
    mprintf(", imaging on");
  else
    mprintf(", imaging off");
  mprintf(".\n");
  return Action::OK;  
}

// Action_CheckStructure::CheckFrame()
int Action_CheckStructure::CheckFrame(int frameNum, Frame const& currentFrame) {
  Matrix_3x3 ucell, recip;
  int Nproblems = 0;
  // Start bond list 
  std::vector<bond_list>::iterator currentBond = bondL_.begin();
  // Get imaging info for non-orthogonal box
  if (ImageType()==NONORTHO) currentFrame.BoxCrd().ToRecip(ucell, recip);
  // Begin loop
  int lastidx = Mask1_.Nselected() - 1;
  for (int maskidx1 = 0; maskidx1 < lastidx; maskidx1++) {
    int atom1 = Mask1_[maskidx1];
    for (int maskidx2 = maskidx1 + 1; maskidx2 < Mask1_.Nselected(); maskidx2++) {
      int atom2 = Mask1_[maskidx2];
      // Get distance^2
      double D2 = DIST2(currentFrame.XYZ(atom1), currentFrame.XYZ(atom2),
                        ImageType(), currentFrame.BoxCrd(), ucell, recip);
      if ( (atom1==(*currentBond).atom1) && (atom2==(*currentBond).atom2) ) {
        // Atoms bonded, check bond length.
        // req has been precalced to (req + bondoffset)^2
        double bondmax = (*currentBond).req;
        // Check for long bond length; distance2 > (req+bondoffset)^2
        if (D2 > bondmax) {
          ++Nproblems;
          outfile_.Printf(
                  "%i\t Warning: Unusual bond length %i:%s to %i:%s (%.2lf)\n", frameNum,
                  atom1+1, CurrentParm_->TruncResAtomName(atom1).c_str(), 
                  atom2+1, CurrentParm_->TruncResAtomName(atom2).c_str(), sqrt(D2));
        }
        // Next bond
        currentBond++;
      } else {
        // Atoms not bonded, check overlap
        if (D2 < nonbondcut2_) {
          ++Nproblems;
          outfile_.Printf(
                  "%i\t Warning: Atoms %i:%s and %i:%s are close (%.2lf)\n", frameNum,
                  atom1+1, CurrentParm_->TruncResAtomName(atom1).c_str(), 
                  atom2+1, CurrentParm_->TruncResAtomName(atom2).c_str(), sqrt(D2));
        }
      }
    } // END second loop over mask atoms
  } // END first loop over mask atoms
  return Nproblems;
}

// Action_CheckStructure::action()
Action::RetType Action_CheckStructure::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  CheckFrame(frameNum+1, *currentFrame);
  return Action::OK;
}
