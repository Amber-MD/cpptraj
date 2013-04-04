// Action_Hbond 
#include <cmath> // sqrt
#include <algorithm> // sort
#include "Action_Hbond.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // DigitWidth
#include "DistRoutines.h"
#include "TorsionRoutines.h"
#include "Constants.h" // RADDEG, DEGRAD

// CONSTRUCTOR
Action_Hbond::Action_Hbond() :
  debug_(0),
  Nframes_(0),
  hasDonorMask_(false),
  hasDonorHmask_(false),
  hasAcceptorMask_(false),
  hasSolventDonor_(false),
  hasSolventAcceptor_(false),
  calcSolvent_(false),
  acut_(0),
  dcut2_(0),
  CurrentParm_(0),
  series_(false),
  NumHbonds_(0),
  NumSolvent_(0),
  NumBridge_(0),
  BridgeID_(0),
  masterDSL_(0)
{}

void Action_Hbond::Help() {
  mprintf("\t[out <filename>] <mask> [angle <cut>] [dist <cut>] [series]\n");
  mprintf("\t[donormask <mask> [donorhmask <mask>]] [acceptormask <mask>]\n");
  mprintf("\t[avgout <filename>]\n");
  mprintf("\t[solventdonor <mask>] [solventacceptor <mask>]\n");
  mprintf("\t[solvout <filename>] [bridgeout <filename>]\n");
  mprintf("\tSearch for hydrogen bonds using atoms in the region specified by mask.\n");
  mprintf("\tIf just <mask> specified donors and acceptors will be automatically searched for.\n");
  mprintf("\tIf donormask is specified but not acceptormask, acceptors will be\n");
  mprintf("\tautomatically searched for in <mask>.\n");
  mprintf("\tIf acceptormask is specified but not donormask, donors will be automatically\n");
  mprintf("\tsearched for in <mask>.\n");
  mprintf("\tIf both donormask and acceptor mask are specified no automatic searching will occur.\n");
  mprintf("\tIf donorhmask is specified atoms in that mask will be paired with atoms in\n");
  mprintf("\tdonormask instead of automatically searching for hydrogen atoms.\n");
}

// Action_Hbond::Init()
Action::RetType Action_Hbond::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  // Get keywords
  DataFile* DF = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  series_ = actionArgs.hasKey("series");
  avgout_ = actionArgs.GetStringKey("avgout");
  solvout_ = actionArgs.GetStringKey("solvout");
  bridgeout_ = actionArgs.GetStringKey("bridgeout");
  acut_ = actionArgs.getKeyDouble("angle",135.0);
  // Convert angle cutoff to radians
  acut_ *= DEGRAD;
  double dcut = actionArgs.getKeyDouble("dist",3.0);
  dcut2_ = dcut * dcut;
  // Get donor mask
  std::string mask = actionArgs.GetStringKey("donormask");
  if (!mask.empty()) {
    DonorMask_.SetMaskString(mask);
    hasDonorMask_=true;
    // Get donorH mask (if specified)
    mask = actionArgs.GetStringKey("donorhmask");
    if (!mask.empty()) {
      DonorHmask_.SetMaskString(mask);
      hasDonorHmask_=true;
    }
  }
  // Get acceptor mask
  mask = actionArgs.GetStringKey("acceptormask");
  if (!mask.empty()) {
    AcceptorMask_.SetMaskString(mask);
    hasAcceptorMask_=true;
  }
  // Get solvent donor mask
  mask = actionArgs.GetStringKey("solventdonor");
  if (!mask.empty()) {
    SolventDonorMask_.SetMaskString(mask);
    hasSolventDonor_ = true;
    calcSolvent_ = true;
  }
  // Get solvent acceptor mask
  mask = actionArgs.GetStringKey("solventacceptor");
  if (!mask.empty()) {
    SolventAcceptorMask_.SetMaskString(mask);
    hasSolventAcceptor_ = true;
    calcSolvent_ = true;
  }
  // Get generic mask
  Mask_.SetMaskString(actionArgs.GetMaskNext());

  // If calculating solvent and avgout filename is specified but 
  // solvout/bridgeout is not, set solvout = avgout and 
  // bridgeout = solvout.
  if ( calcSolvent_ ) {
    if (solvout_.empty() && !avgout_.empty())
      solvout_ = avgout_;
    if (bridgeout_.empty() && !solvout_.empty())
      bridgeout_ = solvout_;
  }

  // Setup datasets
  hbsetname_ = actionArgs.GetStringNext();
  if (hbsetname_.empty())
    hbsetname_ = DSL->GenerateDefaultName("HB");
  NumHbonds_ = DSL->AddSetAspect(DataSet::INT, hbsetname_, "UU");
  if (NumHbonds_==0) return Action::ERR;
  if (DF != 0) DF->AddSet( NumHbonds_ );
  if (calcSolvent_) {
    NumSolvent_ = DSL->AddSetAspect(DataSet::INT, hbsetname_, "UV");
    if (NumSolvent_ == 0) return Action::ERR;
    if (DF != 0) DF->AddSet( NumSolvent_ );
    NumBridge_ = DSL->AddSetAspect(DataSet::INT, hbsetname_, "Bridge");
    if (NumBridge_ == 0) return Action::ERR;
    if (DF != 0) DF->AddSet( NumBridge_ );
    BridgeID_ = DSL->AddSetAspect(DataSet::STRING, hbsetname_, "ID");
    if (BridgeID_ == 0) return Action::ERR;
    if (DF != 0) DF->AddSet( BridgeID_ );
  } 

  mprintf( "  HBOND: ");
  if (!hasDonorMask_ && !hasAcceptorMask_)
    mprintf("Searching for Hbond donors/acceptors in region specified by %s\n",
            Mask_.MaskString());
  else if (hasDonorMask_ && !hasAcceptorMask_)
    mprintf("Donor mask is %s, acceptors will be searched for in region specified by %s\n",
            DonorMask_.MaskString(), Mask_.MaskString());
  else if (hasAcceptorMask_ && !hasDonorMask_)
    mprintf("Acceptor mask is %s, donors will be searched for in a region specified by %s\n",
            AcceptorMask_.MaskString(), Mask_.MaskString());
  else
    mprintf("Donor mask is %s, Acceptor mask is %s\n",
            DonorMask_.MaskString(), AcceptorMask_.MaskString());
  if (hasDonorHmask_)
    mprintf("\tSeparate donor H mask is %s\n", DonorHmask_.MaskString() );
  if (hasSolventDonor_)
    mprintf("\tWill search for hbonds between solute and solvent donors in [%s]\n",
            SolventDonorMask_.MaskString());
  if (hasSolventAcceptor_)
    mprintf("\tWill search for hbonds between solute and solvent acceptors in [%s]\n",
            SolventAcceptorMask_.MaskString());
  mprintf("\tDistance cutoff = %.3lf, Angle Cutoff = %.3lf\n",dcut,acut_*RADDEG);
  if (DF != 0) 
    mprintf("\tDumping # Hbond v time results to %s\n", DF->DataFilename().base());
  if (!avgout_.empty())
    mprintf("\tDumping Hbond avgs to %s\n",avgout_.c_str());
  if (calcSolvent_ && !solvout_.empty())
    mprintf("\tDumping solute-solvent hbond avgs to %s\n", solvout_.c_str());
  if (calcSolvent_ && !bridgeout_.empty())
    mprintf("\tDumping solvent bridging info to %s\n", bridgeout_.c_str());
  if (series_)
    mprintf("\tTime series data for each hbond will be saved for analysis.\n");
  masterDSL_ = DSL;
  return Action::OK;
}

// Action_Hbond::SearchAcceptor()
/** Search for hbond acceptors X in the region specified by amask.
  * If Auto is true select acceptors based on the rule that "Hydrogen 
  * bonds are FON"
  */
void Action_Hbond::SearchAcceptor(HBlistType& alist, AtomMask& amask, bool Auto) {
  bool isAcceptor;
  // Set up acceptors: F, O, N
  // NOTE: Attempt to determine electronegative carbons?
  for (AtomMask::const_iterator atom = amask.begin();
                                atom != amask.end(); ++atom)
  {
    isAcceptor=true;
    // If auto searching, only consider acceptor atoms as F, O, N
    if (Auto) {
      isAcceptor=false;
      if ( (*CurrentParm_)[*atom].Element() == Atom::FLUORINE ||
           (*CurrentParm_)[*atom].Element() == Atom::OXYGEN ||
           (*CurrentParm_)[*atom].Element() == Atom::NITROGEN    )
       isAcceptor=true;
    }
    if (isAcceptor)
      alist.push_back(*atom);
  }
}

// Action_Hbond::SearchDonor()
/** Search for hydrogen bond donors X-H in the region specified by dmask.
  * If Auto is true select donors based on the rule that "Hydrogen bonds 
  * are FON". If useHmask is true pair each atom in dmask with atoms
  * in DonorHmask.
  */
void Action_Hbond::SearchDonor(HBlistType& dlist, AtomMask& dmask, bool Auto,
                               bool useHmask) 
{
  bool isDonor;
  // Set up donors: F-H, O-H, N-H
  AtomMask::const_iterator donorhatom = DonorHmask_.begin();
  for (AtomMask::const_iterator donoratom = dmask.begin();
                                donoratom != dmask.end(); ++donoratom)
  {
    // If this is already an H atom continue
    if ( (*CurrentParm_)[*donoratom].Element() == Atom::HYDROGEN ) continue;
    isDonor = true;
    // If auto searching, only consider donor atoms as F, O, N
    if (Auto) {
      isDonor=false;
      if ( (*CurrentParm_)[*donoratom].Element() == Atom::FLUORINE ||
           (*CurrentParm_)[*donoratom].Element() == Atom::OXYGEN ||
           (*CurrentParm_)[*donoratom].Element() == Atom::NITROGEN   )
        isDonor=true;
    }
    if (isDonor) {
      // If no bonds to this atom assume it is an ion. Only do this if !Auto
      if (!Auto && (*CurrentParm_)[*donoratom].Nbonds() == 0) {
        dlist.push_back(*donoratom);
        dlist.push_back(*donoratom);
      } else {
        if (!useHmask) {
          // Get list of hydrogen atoms bonded to this atom
          for (Atom::bond_iterator batom = (*CurrentParm_)[*donoratom].bondbegin();
                                   batom != (*CurrentParm_)[*donoratom].bondend();
                                   batom++)
          {
            if ( (*CurrentParm_)[*batom].Element() == Atom::HYDROGEN ) {
              //mprintf("BOND H: %i@%s -- %i@%s\n",*donoratom+1,(*CurrentParm_)[*donoratom].c_str(),
              //        *batom+1,(*CurrentParm_)[*batom].c_str());
              dlist.push_back(*donoratom);
              dlist.push_back(*batom);
            }
          }
        } else {
          // Use next atom in donor h atom mask
          if (donorhatom == DonorHmask_.end())
            mprintf("Warning: Donor %i:%s: Ran out of atoms in donor H mask.\n",
                    *donoratom + 1, (*CurrentParm_)[*donoratom].c_str());
          else {
            dlist.push_back(*donoratom);
            dlist.push_back(*(donorhatom++));
          }
        }
      }
    } // END atom is potential donor
  } // END loop over selected atoms
}

// Action_Hbond::Setup()
/** Search for hbond donors and acceptors. */
Action::RetType Action_Hbond::Setup(Topology* currentParm, Topology** parmAddress) {
  CurrentParm_ = currentParm;
  // Set up mask
  if (!hasDonorMask_ || !hasAcceptorMask_) {
    if ( currentParm->SetupIntegerMask( Mask_ ) ) return Action::ERR;
    if ( Mask_.None() ) {
      mprintf("Warning: Hbond::setup: Mask has no atoms.\n");
      return Action::ERR;
    }
  }
  // Set up donor mask
  if (hasDonorMask_) {
    if ( currentParm->SetupIntegerMask( DonorMask_ ) ) return Action::ERR;
    if (DonorMask_.None()) {
      mprintf("Warning: Hbond: DonorMask has no atoms.\n");
      return Action::ERR;
    }
    if ( hasDonorHmask_ ) {
      if ( currentParm->SetupIntegerMask( DonorHmask_ ) ) return Action::ERR;
      if ( DonorHmask_.None() ) {
        mprintf("Warning: Hbond: Donor H mask has no atoms.\n");
        return Action::ERR;
      }
      if ( DonorHmask_.Nselected() != DonorMask_.Nselected() ) {
        mprinterr("Error: There is not a 1 to 1 correspondance between donor and donorH masks.\n");
        mprinterr("Error: donor (%i atoms), donorH (%i atoms).\n",DonorMask_.Nselected(),
                  DonorHmask_.Nselected());
        return Action::ERR;
      }
    }
  }
  // Set up acceptor mask
  if (hasAcceptorMask_) {
    if ( currentParm->SetupIntegerMask( AcceptorMask_ ) ) return Action::ERR;
    if (AcceptorMask_.None()) {
      mprintf("Warning: Hbond: AcceptorMask has no atoms.\n");
      return Action::ERR;
    }
  }
  // Set up solvent donor/acceptor masks
  if (hasSolventDonor_) {
    if (currentParm->SetupIntegerMask( SolventDonorMask_ )) return Action::ERR;
    if (SolventDonorMask_.None()) {
      mprintf("Warning: Hbond: SolventDonorMask has no atoms.\n");
      return Action::ERR;
    }
  }
  if (hasSolventAcceptor_) {
    if (currentParm->SetupIntegerMask( SolventAcceptorMask_ )) return Action::ERR;
    if (SolventAcceptorMask_.None()) {
      mprintf("Warning: Hbond: SolventAcceptorMask has no atoms.\n");
      return Action::ERR;
    }
  }

  // OK TO CLEAR?
  Acceptor_.clear();
  Donor_.clear();
  // SOLUTE: Four cases:
  // 1) DonorMask and AcceptorMask null: donors and acceptors automatically searched for.
  if (!hasDonorMask_ && !hasAcceptorMask_) {
    SearchAcceptor(Acceptor_, Mask_,true);
    SearchDonor(Donor_, Mask_, true, false);
  
  // 2) DonorMask only: acceptors automatically searched for in Mask
  } else if (hasDonorMask_ && !hasAcceptorMask_) {
    SearchAcceptor(Acceptor_, Mask_,true);
    SearchDonor(Donor_, DonorMask_, false, hasDonorHmask_);

  // 3) AcceptorMask only: donors automatically searched for in Mask
  } else if (!hasDonorMask_ && hasAcceptorMask_) {
    SearchAcceptor(Acceptor_, AcceptorMask_, false);
    SearchDonor(Donor_, Mask_, true, false);

  // 4) Both DonorMask and AcceptorMask: No automatic search.
  } else {
    SearchAcceptor(Acceptor_, AcceptorMask_, false);
    SearchDonor(Donor_, DonorMask_, false, hasDonorHmask_);
  }

  // Print acceptor/donor information
  mprintf("\tSet up %zu acceptors:\n", Acceptor_.size() );
  if (debug_>0) {
    for (HBlistType::iterator accept = Acceptor_.begin(); accept!=Acceptor_.end(); accept++)
      mprintf("        %8i: %4s\n",*accept+1,(*currentParm)[*accept].c_str());
  }
  mprintf("\tSet up %zu donors:\n", Donor_.size()/2 );
  if (debug_>0) {
    for (HBlistType::iterator donor = Donor_.begin(); donor!=Donor_.end(); donor++) {
      int atom = (*donor);
      ++donor;
      int a2   = (*donor);
      mprintf("        %8i:%4s - %8i:%4s\n",atom+1,(*currentParm)[atom].c_str(),
              a2+1,(*currentParm)[a2].c_str()); 
    } 
  }
  if ( Acceptor_.empty() && Donor_.empty() ) {
    mprinterr("Error: No HBond donors or acceptors.\n");
    return Action::ERR;
  }

  // SOLVENT:
  if (hasSolventAcceptor_) {
    SolventAcceptor_.clear();
    SearchAcceptor(SolventAcceptor_, SolventAcceptorMask_, false);
    mprintf("\tSet up %zu solvent acceptors\n", SolventAcceptor_.size() );
  }
  if (hasSolventDonor_) {
    SolventDonor_.clear();
    SearchDonor(SolventDonor_, SolventDonorMask_, false, false);
    mprintf("\tSet up %zu solvent donors\n", SolventDonor_.size()/2 );
  }

  return Action::OK;
}

// Action_Hbond::AtomsAreHbonded()
/** Used to determine if solute atoms are bonded to solvent atoms. */
int Action_Hbond::AtomsAreHbonded(Frame const& currentFrame, int frameNum, 
                                  int a_atom, int d_atom, int h_atom, 
                                  int hbidx, bool solutedonor) 
{
  std::string hblegend;
  HbondType HB;
  double angle;

  if (a_atom == d_atom) return 0;
  double dist2 = DIST2_NoImage(currentFrame.XYZ(a_atom), currentFrame.XYZ(d_atom));
  if (dist2 > dcut2_) return 0;
  /*mprintf("DEBUG: Donor %i@%s -- acceptor %i@%s = %lf",
         d_atom+1, (*currentParm)[d_atom].c_str(),
         a_atom+1, (*currentParm)[a_atom].c_str(), sqrt(dist2));*/
  // For ions, donor atom will be same as h atom so no angle needed.
  if (d_atom != h_atom) {
    angle = CalcAngle( currentFrame.XYZ(a_atom), 
                       currentFrame.XYZ(h_atom),
                       currentFrame.XYZ(d_atom) );
    if (angle < acut_) return 0;
  }
  double dist = sqrt(dist2);
  //mprintf( "A-D HBOND[%6i]: %6i@%-4s ... %6i@%-4s-%6i@%-4s Dst=%6.2lf Ang=%6.2lf\n", hbidx, 
  //        a_atom, (*currentParm)[a_atom].c_str(),
  //        h_atom, (*currentParm)[h_atom].c_str(), 
  //        d_atom, (*currentParm)[d_atom].c_str(), dist, angle*RADDEG);
  // Find hbond in map
  HBmapType::iterator entry = SolventMap_.find( hbidx );
  if (entry == SolventMap_.end() ) {
    // New Hbond
    if (solutedonor) {
      HB.A = -1; // Do not care about which solvent acceptor
      HB.D = d_atom;
      HB.H = h_atom;
      hblegend = CurrentParm_->TruncResAtomName(h_atom) + "-V";
    } else {
      HB.A = a_atom;
      HB.D = -1; // Do not care about solvent donor heavy atom
      HB.H = -1; // Do not care about solvent donor H atom
      hblegend = CurrentParm_->TruncResAtomName(a_atom) + "-V";
    }
    HB.Frames = 1;
    HB.dist = dist;
    HB.angle = angle;
    if (series_) {
      HB.data_ = (DataSet_integer*) masterDSL_->AddSetIdxAspect( DataSet::INT, hbsetname_, 
                                                          hbidx, "solventhb" );
      //mprinterr("Created Solvent HB data frame %i idx %i %p\n",frameNum,hbidx,HB.data_);
      HB.data_->Resize( masterDSL_->MaxFrames() );
      HB.data_->SetLegend( hblegend );
      (*HB.data_)[ frameNum ] = 1;
    }
    SolventMap_.insert( entry, std::pair<int,HbondType>(hbidx, HB) );
  } else {
    (*entry).second.Frames++;
    (*entry).second.dist += dist;
    (*entry).second.angle += angle;
    if (series_) {
      //mprinterr("Adding Solvent HB data frame %i idx %i %p\n",frameNum,hbidx,(*entry).second.data_);
      (*(*entry).second.data_)[ frameNum ] = 1;
    }
  }     
  return 1;
}

// Action_Hbond::DoAction()
/** Calculate distance between all donors and acceptors. Store Hbond info.
  */    
Action::RetType Action_Hbond::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  // accept ... H-D
  int D, H;
  double dist, dist2, angle;
  HBmapType::iterator it;
  HbondType HB;

  // SOLUTE-SOLUTE HBONDS
  int hbidx = 0; 
  int numHB=0;
  for (HBlistType::iterator donor = Donor_.begin(); donor!=Donor_.end(); ++donor) {
    D = (*donor);
    ++donor;
    H = (*donor);
    for (HBlistType::iterator accept = Acceptor_.begin(); 
                              accept != Acceptor_.end(); ++accept, ++hbidx) 
    {
      if (*accept == D) continue;
      dist2 = DIST2_NoImage(currentFrame->XYZ(*accept), currentFrame->XYZ(D));
      if (dist2 > dcut2_) continue;
      angle = CalcAngle( currentFrame->XYZ(*accept), 
                         currentFrame->XYZ(H),
                         currentFrame->XYZ(D)       );
      if (angle < acut_) continue;
      ++numHB;
      dist = sqrt(dist2);
      //mprintf("HBOND[%i] A=%i ... H=%i D=%i  Dist=%f Angle=%f\n", hbidx, // DEBUG 
      //        *accept + 1, H + 1, D + 1, dist, angle*RADDEG);            // DEBUG
      // Find hbond in map
      it = HbondMap_.find( hbidx );
      if (it == HbondMap_.end() ) {
        // New Hbond
        HB.A = *accept;
        HB.D = D;
        HB.H = H;
        HB.Frames = 1;
        HB.dist = dist;
        HB.angle = angle;
        if (series_) {
          std::string hblegend = CurrentParm_->TruncResAtomName(*accept) + "-" +
                                 CurrentParm_->TruncResAtomName(D) + "-" +
                                 (*CurrentParm_)[H].Name().Truncated();
          HB.data_ = (DataSet_integer*) masterDSL_->AddSetIdxAspect( DataSet::INT, hbsetname_,
                                                              hbidx, "solutehb" );
          //mprinterr("Created solute Hbond dataset index %i\n", hbidx);
          HB.data_->Resize( masterDSL_->MaxFrames() );
          HB.data_->SetLegend( hblegend );
          (*HB.data_)[ frameNum ] = 1 ;
        }
        HbondMap_.insert( it, std::pair<int,HbondType>(hbidx, HB) );
      } else {
        (*it).second.Frames++;
        (*it).second.dist += dist;
        (*it).second.angle += angle;
        if (series_)
          (*(*it).second.data_)[ frameNum ] = 1;
      }
    }
  }
  NumHbonds_->Add(frameNum, &numHB);
  //mprintf("HBOND: Scanned %i hbonds.\n",hbidx);
  
  if (calcSolvent_) {
    // Contains info about which residue(s) a Hbonding solvent mol is
    // Hbonded to.
    std::map< int, std::set<int> > solvent2solute;
    int solventHbonds = 0;
    // SOLUTE DONOR-SOLVENT ACCEPTOR
    // Index by solute H atom. 
    if (hasSolventAcceptor_) {
      numHB = 0;
      for (HBlistType::iterator donor = Donor_.begin(); 
                                donor != Donor_.end(); ++donor) 
      {
        D = (*donor);
        ++donor;
        H = (*donor);
        for (HBlistType::iterator accept = SolventAcceptor_.begin(); 
                                  accept != SolventAcceptor_.end(); ++accept)
        { 
          if (AtomsAreHbonded( *currentFrame, frameNum, *accept, D, H, H, true )) {
            ++numHB;
            int soluteres = (*CurrentParm_)[D].ResNum();
            int solventmol = (*CurrentParm_)[*accept].ResNum();
            solvent2solute[solventmol].insert( soluteres );
            //mprintf("DBG:\t\tSolvent Res %i bonded to solute res %i\n",solventmol+1,soluteres+1);
          }
        }
      }
      //mprintf("DEBUG: # Solvent Acceptor to Solute Donor Hbonds is %i\n", numHB);
      solventHbonds += numHB;
    }
    // SOLVENT DONOR-SOLUTE ACCEPTOR
    // Index by solute acceptor atom
    if (hasSolventDonor_) {
      numHB = 0;
      for (HBlistType::iterator donor = SolventDonor_.begin();
                                donor != SolventDonor_.end(); ++donor)
      {
        D = (*donor);
        ++donor;
        H = (*donor);
        for (HBlistType::iterator accept = Acceptor_.begin();
                                  accept != Acceptor_.end(); ++accept)
        {
          if (AtomsAreHbonded( *currentFrame, frameNum, *accept, D, H, *accept, false )) {
            ++numHB;
            int soluteres = (*CurrentParm_)[*accept].ResNum();
            int solventmol = (*CurrentParm_)[D].ResNum();
            solvent2solute[solventmol].insert( soluteres );
            //mprintf("DBG:\t\tSolvent Res %i bonded to solute res %i\n",solventmol+1,soluteres+1);
          }
        }
      }
      //mprintf("DEBUG: # Solvent Donor to Solute Acceptor Hbonds is %i\n", numHB);
      solventHbonds += numHB;
    }
    NumSolvent_->Add(frameNum, &solventHbonds);

    // Determine number of bridging waters.
    numHB = 0;
    std::string bridgeID;
    for (std::map< int, std::set<int> >::iterator bridge = solvent2solute.begin();
                                                  bridge != solvent2solute.end();
                                                  ++bridge)
    {
      // If solvent molecule is bound to 2 or more different residues,
      // it is bridging. 
      if ( (*bridge).second.size() > 1) {
        ++numHB;
        bridgeID.append(integerToString( (*bridge).first+1 ) + "("); // Bridging Solvent res 
        for (std::set<int>::iterator res = (*bridge).second.begin();
                                     res != (*bridge).second.end(); ++res)
          bridgeID.append( integerToString( *res+1 ) + "+" ); // Solute res being bridged
        bridgeID.append("),");
        // Find bridge in map based on this combo of residues (bridge.second)
        BridgeType::iterator b_it = BridgeMap_.find( (*bridge).second );
        if (b_it == BridgeMap_.end() ) // New Bridge 
          BridgeMap_.insert( b_it, std::pair<std::set<int>,int>((*bridge).second, 1) );
        else                           // Increment bridge #frames
          (*b_it).second++;
      }
    }
    if (bridgeID.empty())
      bridgeID.assign("None");
    NumBridge_->Add(frameNum, &numHB);
    BridgeID_->Add(frameNum, (char*)bridgeID.c_str()); // FIXME: Fix cast
  }

  ++Nframes_;

  return Action::OK;
}

// Action_Hbond::print()
/** Print average occupancies over all frames for all detected Hbonds
  */
void Action_Hbond::Print() {
  std::vector<HbondType> HbondList; // For sorting
  std::string Aname, Hname, Dname;
  CpptrajFile outfile;

  if (CurrentParm_ == 0) return;
  // Calculate necessary column width for strings based on how many residues.
  // ResName+'_'+ResNum+'@'+AtomName | NUM = 4+1+R+1+4 = R+10
  int NUM = DigitWidth( CurrentParm_->Nres() ) + 10;

  // Solute Hbonds 
  if (!avgout_.empty()) { 
    if (outfile.OpenWrite( avgout_ )) return;
    // Place all detected Hbonds in a list and sort. Free memory as we go. 
    for (HBmapType::iterator it = HbondMap_.begin(); it!=HbondMap_.end(); ++it) 
      HbondList.push_back( (*it).second );
    HbondMap_.clear();
    sort( HbondList.begin(), HbondList.end(), hbond_cmp() );
    // Calculate averages and print
    //outfile.Printf("#Solute Hbonds:\n");
    outfile.Printf("%-*s %*s %*s %8s %12s %12s %12s\n", NUM, "#Acceptor", 
                   NUM, "DonorH", NUM, "Donor", "Frames", "Frac", "AvgDist", "AvgAng");
    for (std::vector<HbondType>::iterator hbond = HbondList.begin(); 
                                          hbond!=HbondList.end(); ++hbond ) 
    {
      double avg = (double) (*hbond).Frames;
      avg = avg / ((double) Nframes_);
      double dist = (double) (*hbond).dist;
      dist = dist / ((double) (*hbond).Frames);
      double angle = (double) (*hbond).angle;
      angle = angle / ((double) (*hbond).Frames);
      angle *= RADDEG;

      Aname = CurrentParm_->TruncResAtomName((*hbond).A);
      Hname = CurrentParm_->TruncResAtomName((*hbond).H);
      Dname = CurrentParm_->TruncResAtomName((*hbond).D);

      outfile.Printf("%-*s %*s %*s %8i %12.4lf %12.4lf %12.4lf\n",
                     NUM, Aname.c_str(), NUM, Hname.c_str(), NUM, Dname.c_str(),
                     (*hbond).Frames, avg, dist, angle);
    }
    outfile.CloseFile();
  }

  // Solute-solvent Hbonds 
  if (!solvout_.empty() && calcSolvent_) {
    if (solvout_ == avgout_) {
      if (outfile.OpenAppend( solvout_ )) return;
    } else {
      if (outfile.OpenWrite( solvout_)) return;
    }
    HbondList.clear();
    for (HBmapType::iterator it = SolventMap_.begin(); it != SolventMap_.end(); ++it)
      HbondList.push_back( (*it).second );
    SolventMap_.clear();
    sort( HbondList.begin(), HbondList.end(), hbond_cmp() );
    // Calc averages and print
    outfile.Printf("#Solute-Solvent Hbonds:\n");
    outfile.Printf("%-*s %*s %*s %8s %12s %12s %12s\n", NUM, "#Acceptor", 
                   NUM, "DonorH", NUM, "Donor", "Count", "Frac", "AvgDist", "AvgAng");
    for (std::vector<HbondType>::iterator hbond = HbondList.begin();
                                          hbond != HbondList.end(); ++hbond )
    {
      // Average has slightly diff meaning since for any given frame multiple
      // solvent can bond to the same solute.
      double avg = (double) (*hbond).Frames;
      avg = avg / ((double) Nframes_);
      double dist = (double) (*hbond).dist;
      dist = dist / ((double) (*hbond).Frames);
      double angle = (double) (*hbond).angle;
      angle = angle / ((double) (*hbond).Frames);
      angle *= RADDEG;

      if ((*hbond).A==-1) // Solvent acceptor
        Aname = "SolventAcc";
      else
        Aname = CurrentParm_->TruncResAtomName((*hbond).A);
      if ((*hbond).D==-1) { // Solvent donor
        Dname = "SolventDnr";
        Hname = "SolventH";
      } else {
        Dname = CurrentParm_->TruncResAtomName((*hbond).D);
        Hname = CurrentParm_->TruncResAtomName((*hbond).H);
      }

      outfile.Printf("%-*s %*s %*s %8i %12.4lf %12.4lf %12.4lf\n",
                     NUM, Aname.c_str(), NUM, Hname.c_str(), NUM, Dname.c_str(),
                     (*hbond).Frames, avg, dist, angle);
    }
    outfile.CloseFile();
  }

  // BRIDGING INFO
  if (!bridgeout_.empty() && calcSolvent_) {
    if (bridgeout_ == avgout_ || bridgeout_ == solvout_) {
      if (outfile.OpenAppend( bridgeout_ )) return;
    } else {
      if (outfile.OpenWrite( bridgeout_ )) return; 
    }
    outfile.Printf("#Bridging Solute Residues:\n");
    // Place bridging values in a vector for sorting
    std::vector<std::pair< std::set<int>, int> > bridgevector;
    for (BridgeType::iterator it = BridgeMap_.begin(); 
                              it != BridgeMap_.end(); ++it) 
      bridgevector.push_back( *it );
    std::sort( bridgevector.begin(), bridgevector.end(), bridge_cmp() );
    for (std::vector<std::pair< std::set<int>, int> >::iterator bv = bridgevector.begin();
                                                                bv != bridgevector.end(); ++bv)
    {
      outfile.Printf("Bridge Res");
      for (std::set<int>::iterator res = (*bv).first.begin();
                                   res != (*bv).first.end(); ++res)
        outfile.Printf(" %i:%s", *res+1, CurrentParm_->Res( *res ).c_str());
      outfile.Printf(", %i frames.\n", (*bv).second);
    } 
    outfile.CloseFile();
  } 
}
