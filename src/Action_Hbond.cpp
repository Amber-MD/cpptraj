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
  avgout_(0), solvout_(0), bridgeout_(0),
  UUseriesout_(0), UVseriesout_(0),
  useAtomNum_(false),
  hasDonorMask_(false),
  hasDonorHmask_(false),
  hasAcceptorMask_(false),
  hasSolventDonor_(false),
  hasSolventAcceptor_(false),
  calcSolvent_(false),
  noIntramol_(false),
  acut_(0),
  dcut2_(0),
  CurrentParm_(0),
  series_(false),
  seriesUpdated_(false),
  NumHbonds_(0),
  NumSolvent_(0),
  NumBridge_(0),
  BridgeID_(0),
  masterDSL_(0)
{}

void Action_Hbond::Help() const {
  mprintf("\t[<dsname>] [out <filename>] [<mask>] [angle <acut>] [dist <dcut>]\n"
          "\t[donormask <dmask> [donorhmask <dhmask>]] [acceptormask <amask>]\n"
          "\t[avgout <filename>] [printatomnum] [nointramol] [image]\n"
          "\t[solventdonor <sdmask>] [solventacceptor <samask>]\n"
          "\t[solvout <filename>] [bridgeout <filename>]\n"
          "\t[series [uuseries <filename>] [uvseries <filename>]]\n"
          "  Hydrogen bond is defined as A-HD, where A is acceptor heavy atom, H is\n"
          "  hydrogen, D is donor heavy atom. Hydrogen bond is formed when\n"
          "  A to D distance < dcut and A-H-D angle > acut; if acut < 0 it is ignored.\n"
          "  Search for hydrogen bonds using atoms in the region specified by mask.\n"
          "  If just <mask> specified donors and acceptors will be automatically searched for.\n"
          "  If donormask is specified but not acceptormask, acceptors will be\n"
          "  automatically searched for in <mask>.\n"
          "  If acceptormask is specified but not donormask, donors will be automatically\n"
          "  searched for in <mask>.\n"
          "  If both donormask and acceptor mask are specified no automatic searching will occur.\n"
          "  If donorhmask is specified atoms in that mask will be paired with atoms in\n"
          "  donormask instead of automatically searching for hydrogen atoms.\n");
}

// Action_Hbond::Init()
Action::RetType Action_Hbond::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  // Get keywords
  Image_.InitImaging( (actionArgs.hasKey("image")) );
  DataFile* DF = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  series_ = actionArgs.hasKey("series");
  if (series_) {
    UUseriesout_ = init.DFL().AddDataFile(actionArgs.GetStringKey("uuseries"), actionArgs);
    UVseriesout_ = init.DFL().AddDataFile(actionArgs.GetStringKey("uvseries"), actionArgs);
    init.DSL().SetDataSetsPending(true);
  }
  std::string avgname = actionArgs.GetStringKey("avgout");
  std::string solvname = actionArgs.GetStringKey("solvout");
  if (solvname.empty()) solvname = avgname;
  std::string bridgename = actionArgs.GetStringKey("bridgeout");
  if (bridgename.empty()) bridgename = solvname;
  
  useAtomNum_ = actionArgs.hasKey("printatomnum");
  acut_ = actionArgs.getKeyDouble("angle",135.0);
  noIntramol_ = actionArgs.hasKey("nointramol");
  // Convert angle cutoff to radians
  acut_ *= Constants::DEGRAD;
  double dcut = actionArgs.getKeyDouble("dist",3.0);
  dcut = actionArgs.getKeyDouble("distance", dcut); // for PTRAJ compat.
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

  // Setup datasets
  hbsetname_ = actionArgs.GetStringNext();
  if (hbsetname_.empty())
    hbsetname_ = init.DSL().GenerateDefaultName("HB");
  NumHbonds_ = init.DSL().AddSet(DataSet::INTEGER, MetaData(hbsetname_, "UU"));
  if (NumHbonds_==0) return Action::ERR;
  if (DF != 0) DF->AddDataSet( NumHbonds_ );
  avgout_ = init.DFL().AddCpptrajFile(avgname, "Avg. solute-solute HBonds");
  if (calcSolvent_) {
    NumSolvent_ = init.DSL().AddSet(DataSet::INTEGER, MetaData(hbsetname_, "UV"));
    if (NumSolvent_ == 0) return Action::ERR;
    if (DF != 0) DF->AddDataSet( NumSolvent_ );
    NumBridge_ = init.DSL().AddSet(DataSet::INTEGER, MetaData(hbsetname_, "Bridge"));
    if (NumBridge_ == 0) return Action::ERR;
    if (DF != 0) DF->AddDataSet( NumBridge_ );
    BridgeID_ = init.DSL().AddSet(DataSet::STRING, MetaData(hbsetname_, "ID"));
    if (BridgeID_ == 0) return Action::ERR;
    if (DF != 0) DF->AddDataSet( BridgeID_ );
    solvout_ = init.DFL().AddCpptrajFile(solvname,"Avg. solute-solvent HBonds");
    bridgeout_ = init.DFL().AddCpptrajFile(bridgename,"Solvent bridging info");
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
  if (noIntramol_)
    mprintf("\tOnly looking for intermolecular hydrogen bonds.\n");
  if (hasSolventDonor_)
    mprintf("\tWill search for hbonds between solute and solvent donors in [%s]\n",
            SolventDonorMask_.MaskString());
  if (hasSolventAcceptor_)
    mprintf("\tWill search for hbonds between solute and solvent acceptors in [%s]\n",
            SolventAcceptorMask_.MaskString());
  mprintf("\tDistance cutoff = %.3lf, Angle Cutoff = %.3lf\n",dcut,acut_*Constants::RADDEG);
  if (DF != 0) 
    mprintf("\tWriting # Hbond v time results to %s\n", DF->DataFilename().full());
  if (avgout_ != 0)
    mprintf("\tWriting Hbond avgs to %s\n",avgout_->Filename().full());
  if (calcSolvent_ && solvout_ != 0)
    mprintf("\tWriting solute-solvent hbond avgs to %s\n", solvout_->Filename().full());
  if (calcSolvent_ && bridgeout_ != 0)
    mprintf("\tWriting solvent bridging info to %s\n", bridgeout_->Filename().full());
  if (useAtomNum_)
    mprintf("\tAtom numbers will be written to output.\n");
  if (series_) {
    mprintf("\tTime series data for each hbond will be saved for analysis.\n");
    if (UUseriesout_ != 0) mprintf("\tWriting solute-solute time series to %s\n",
                                   UUseriesout_->DataFilename().full());
    if (UVseriesout_ != 0) mprintf("\tWriting solute-solvent time series to %s\n",
                                   UVseriesout_->DataFilename().full());
  }
  if (Image_.UseImage())
    mprintf("\tImaging enabled.\n");
  masterDSL_ = init.DslPtr();
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
    Atom const& AcceptorAtom = (*CurrentParm_)[*atom];
    isAcceptor=true;
    // If auto searching, only consider acceptor atoms as F, O, N
    if (Auto) {
      isAcceptor=false;
      // Ignore solvent when auto-searching 
      if (!CurrentParm_->Mol(AcceptorAtom.MolNum()).IsSolvent()) {
        if ( AcceptorAtom.Element() == Atom::FLUORINE ||
             AcceptorAtom.Element() == Atom::OXYGEN ||
             AcceptorAtom.Element() == Atom::NITROGEN    )
         isAcceptor=true;
      }
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
    Atom const& DonorAtom = (*CurrentParm_)[*donoratom];
    // If this is already an H atom continue
    if ( DonorAtom.Element() == Atom::HYDROGEN ) continue;
    isDonor = true;
    // If auto searching, only consider donor atoms as F, O, N
    if (Auto) {
      isDonor=false;
      // Ignore solvent when auto-searching
      if (!CurrentParm_->Mol(DonorAtom.MolNum()).IsSolvent()) {
        if ( DonorAtom.Element() == Atom::FLUORINE ||
             DonorAtom.Element() == Atom::OXYGEN ||
             DonorAtom.Element() == Atom::NITROGEN   )
          isDonor=true;
      }
    }
    if (isDonor) {
      // If no bonds to this atom assume it is an ion. Only do this if !Auto
      if (!Auto && DonorAtom.Nbonds() == 0) {
        dlist.push_back(*donoratom);
        dlist.push_back(*donoratom);
      } else {
        if (!useHmask) {
          // Get list of hydrogen atoms bonded to this atom
          for (Atom::bond_iterator batom = DonorAtom.bondbegin();
                                   batom != DonorAtom.bondend(); ++batom)
          {
            if ( (*CurrentParm_)[*batom].Element() == Atom::HYDROGEN ) {
              //mprintf("BOND H: %i@%s -- %i@%s\n",*donoratom+1,DonorAtom.c_str(),
              //        *batom+1,(*CurrentParm_)[*batom].c_str());
              dlist.push_back(*donoratom);
              dlist.push_back(*batom);
            }
          }
        } else {
          // Use next atom in donor h atom mask
          if (donorhatom == DonorHmask_.end())
            mprintf("Warning: Donor %i:%s: Ran out of atoms in donor H mask.\n",
                    *donoratom + 1, DonorAtom.c_str());
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
Action::RetType Action_Hbond::Setup(ActionSetup& setup) {
  CurrentParm_ = setup.TopAddress();
  Image_.SetupImaging( setup.CoordInfo().TrajBox().Type() );
  // Set up mask
  if (!hasDonorMask_ || !hasAcceptorMask_) {
    if ( setup.Top().SetupIntegerMask( Mask_ ) ) return Action::ERR;
    if ( Mask_.None() ) {
      mprintf("Warning: Mask has no atoms.\n");
      return Action::SKIP;
    }
  }
  // Set up donor mask
  if (hasDonorMask_) {
    if ( setup.Top().SetupIntegerMask( DonorMask_ ) ) return Action::ERR;
    if (DonorMask_.None()) {
      mprintf("Warning: DonorMask has no atoms.\n");
      return Action::SKIP;
    }
    if ( hasDonorHmask_ ) {
      if ( setup.Top().SetupIntegerMask( DonorHmask_ ) ) return Action::ERR;
      if ( DonorHmask_.None() ) {
        mprintf("Warning: Donor H mask has no atoms.\n");
        return Action::SKIP;
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
    if ( setup.Top().SetupIntegerMask( AcceptorMask_ ) ) return Action::ERR;
    if (AcceptorMask_.None()) {
      mprintf("Warning: AcceptorMask has no atoms.\n");
      return Action::SKIP;
    }
  }
  if (calcSolvent_) {
    // Set up solvent donor/acceptor masks
    if (hasSolventDonor_) {
      if (setup.Top().SetupIntegerMask( SolventDonorMask_ )) return Action::ERR;
      if (SolventDonorMask_.None()) {
        mprintf("Warning: SolventDonorMask has no atoms.\n");
        return Action::SKIP;
      }
    }
    if (hasSolventAcceptor_) {
      if (setup.Top().SetupIntegerMask( SolventAcceptorMask_ )) return Action::ERR;
      if (SolventAcceptorMask_.None()) {
        mprintf("Warning: SolventAcceptorMask has no atoms.\n");
        return Action::SKIP;
      }
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
    for (HBlistType::const_iterator accept = Acceptor_.begin(); accept!=Acceptor_.end(); accept++)
      mprintf("        %8i: %4s\n",*accept+1,setup.Top()[*accept].c_str());
  }
  mprintf("\tSet up %zu donors:\n", Donor_.size()/2 );
  if (debug_>0) {
    for (HBlistType::const_iterator donor = Donor_.begin(); donor!=Donor_.end(); donor++) {
      int atom = (*donor);
      ++donor;
      int a2   = (*donor);
      mprintf("        %8i:%4s - %8i:%4s\n",atom+1,setup.Top()[atom].c_str(),
              a2+1,setup.Top()[a2].c_str()); 
    } 
  }
  if ( Acceptor_.empty() && Donor_.empty() ) {
    mprintf("Warning: No HBond donors or acceptors.\n");
    return Action::SKIP;
  }

  // SOLVENT:
  if (calcSolvent_) {
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
  }
  if (Image_.ImagingEnabled())
    mprintf("\tImaging on.\n");
  else
    mprintf("\tImaging off.\n");
  // Estimate maximum memory usage by hbond data.
  // Max # of U-U pairs.
  size_t nPairs = (Acceptor_.size() * (Donor_.size()/2));
  // If calculating solvent every U acceptor can have V donor etc.
  if (calcSolvent_)
    nPairs += (Acceptor_.size() + Donor_.size()/2);
  mprintf("\tEstimated max potential memory usage: %s\n",
          MemoryUsage(nPairs, masterDSL_->MaxFrames()).c_str());

  return Action::OK;
}

std::string Action_Hbond::MemoryUsage(size_t nPairs, size_t nFrames) const {
  static const size_t HBmapTypeElt = 32 + sizeof(int) + 
                                     (2*sizeof(double) + sizeof(DataSet_integer*) + 4*sizeof(int));
  static const size_t BridgeTypeElt = 32 + sizeof(std::set<int>) + sizeof(int);
  size_t memTotal = nPairs * HBmapTypeElt;
  // If calculating series every hbond will have time series.
  // NOTE: This does not include memory used by DataSet.
  if (series_ && nFrames > 0) {
    size_t seriesSet = (nFrames * sizeof(int)) + sizeof(std::vector<int>);
    memTotal += (seriesSet * nPairs);
  }
  // Current memory used by bridging solvent
  for (BridgeType::const_iterator it = BridgeMap_.begin(); it != BridgeMap_.end(); ++it)
    memTotal += (it->first.size() * sizeof(int));
  memTotal += (BridgeMap_.size() * BridgeTypeElt);
  return ByteString( memTotal, BYTE_DECIMAL );
}

double Action_Hbond::ImagedAngle(const double* xyz_a, const double* xyz_h, const double* xyz_d) const
{
  double angle;
  Vec3 VH = Vec3(xyz_h);
  Vec3 H_A = MinImagedVec(VH, Vec3(xyz_a), ucell_, recip_);
  Vec3 H_D = Vec3(xyz_d) - VH;
  double rha = H_A.Magnitude2();
  double rhd = H_D.Magnitude2();
  if (rha > Constants::SMALL && rhd > Constants::SMALL) {
    angle = (H_A * H_D) / sqrt(rha * rhd);
    if      (angle >  1.0) angle =  1.0;
    else if (angle < -1.0) angle = -1.0;
    angle = acos(angle);
  } else
    angle = 0.0;
  return angle;
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
  double dist2 = DIST2(currentFrame.XYZ(a_atom), currentFrame.XYZ(d_atom),
                       Image_.ImageType(), currentFrame.BoxCrd(),
                       ucell_, recip_);
  if (dist2 > dcut2_) return 0;
  /*mprintf("DEBUG: Donor %i@%s -- acceptor %i@%s = %lf",
         d_atom+1, (*currentParm)[d_atom].c_str(),
         a_atom+1, (*currentParm)[a_atom].c_str(), sqrt(dist2));*/
  // For ions, donor atom will be same as h atom so no angle needed.
  if (d_atom != h_atom) {
    if (Image_.ImageType() == NOIMAGE) {
      angle = CalcAngle( currentFrame.XYZ(a_atom), 
                         currentFrame.XYZ(h_atom),
                         currentFrame.XYZ(d_atom) );
    } else {
      angle = ImagedAngle( currentFrame.XYZ(a_atom),
                           currentFrame.XYZ(h_atom),
                           currentFrame.XYZ(d_atom) );
    }
    if (angle < acut_) return 0;
  } else
    angle = 0.0;
  double dist = sqrt(dist2);
  //mprintf( "A-D HBOND[%6i]: %6i@%-4s ... %6i@%-4s-%6i@%-4s Dst=%6.2lf Ang=%6.2lf\n", hbidx, 
  //        a_atom, (*currentParm)[a_atom].c_str(),
  //        h_atom, (*currentParm)[h_atom].c_str(), 
  //        d_atom, (*currentParm)[d_atom].c_str(), dist, angle*Constants::RADDEG);
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
      HB.data_ = (DataSet_integer*)
                 masterDSL_->AddSet( DataSet::INTEGER,
                                     MetaData(hbsetname_, "solventhb", hbidx) ); 
      if (UVseriesout_ != 0) UVseriesout_->AddDataSet( HB.data_ );
      //mprinterr("Created Solvent HB data frame %i idx %i %p\n",frameNum,hbidx,HB.data_);
      HB.data_->SetLegend( hblegend );
      HB.data_->AddVal( frameNum, 1 );
    }
    SolventMap_.insert( entry, std::pair<int,HbondType>(hbidx, HB) );
  } else {
    entry->second.Frames++;
    entry->second.dist += dist;
    entry->second.angle += angle;
    if (series_) {
      //mprinterr("Adding Solvent HB data frame %i idx %i %p\n",frameNum,hbidx,(*entry).second.data_);
      entry->second.data_->AddVal( frameNum, 1 );
    }
  }     
  return 1;
}

/** This macro is used instead of a function to ensure it is inlined, which
  * actually results in a 1.05x speedup vs an attempted inline function call.
  * This is also faster (1.08x) than having an 'if' statement in the solute 
  * hbond loop.
  */
#define SoluteHbond(a_atom, d_atom, h_atom)  { \
  dist2 = DIST2(frm.Frm().XYZ(a_atom), frm.Frm().XYZ(d_atom), Image_.ImageType(), frm.Frm().BoxCrd(), ucell_, recip_); \
  if (dist2 > dcut2_) continue; \
  if (Image_.ImageType() == NOIMAGE) \
    angle = CalcAngle(frm.Frm().XYZ(a_atom), frm.Frm().XYZ(h_atom), frm.Frm().XYZ(d_atom)); \
  else \
    angle = ImagedAngle(frm.Frm().XYZ(a_atom), frm.Frm().XYZ(h_atom), frm.Frm().XYZ(d_atom)); \
  if (angle < acut_) continue; \
  ++numHB; \
  dist = sqrt(dist2); \
  it = HbondMap_.find( hbidx ); \
  if (it == HbondMap_.end() ) { \
    HB.A = a_atom; \
    HB.D = d_atom; \
    HB.H = h_atom; \
    HB.Frames = 1; \
    HB.dist = dist; \
    HB.angle = angle; \
    if (series_) { \
      std::string hblegend = CurrentParm_->TruncResAtomName(a_atom) + "-" + \
                             CurrentParm_->TruncResAtomName(d_atom) + "-" + \
                             (*CurrentParm_)[h_atom].Name().Truncated(); \
      HB.data_ = (DataSet_integer*) \
                 masterDSL_->AddSet( DataSet::INTEGER, \
                                     MetaData(hbsetname_, "solutehb", hbidx) ); \
      if (UUseriesout_ != 0) UUseriesout_->AddDataSet( HB.data_ ); \
      HB.data_->SetLegend( hblegend ); \
      HB.data_->AddVal( frameNum, 1 ); \
    } \
    HbondMap_.insert( it, std::pair<int,HbondType>(hbidx, HB) ); \
  } else { \
    it->second.Frames++; \
    it->second.dist += dist; \
    it->second.angle += angle; \
    if (series_) \
      it->second.data_->AddVal( frameNum, 1 ); \
  } \
} 

// Action_Hbond::DoAction()
/** Calculate distance between all donors and acceptors. Store Hbond info.
  */    
Action::RetType Action_Hbond::DoAction(int frameNum, ActionFrame& frm) {
  int D, H;
  // accept ... H-D
  if (Image_.ImagingEnabled())
    frm.Frm().BoxCrd().ToRecip(ucell_, recip_);
  // SOLUTE-SOLUTE HBONDS
# ifdef HBOND_OPENMP
  int dAt, hAt, aAt, didx, aidx, hbidx, mol1 = -1, numHB = 0;
  double dist, dist2, angle;
  HBmapType::iterator it;
  HbondType HB;

  for (didx = 0; didx < (int)Donor_.size(); didx += 2) {
    dAt = Donor_[didx  ];
    hAt = Donor_[didx+1];
    if (noIntramol_)
      mol1 = (*CurrentParm_)[dAt].MolNum();
    for (aidx = 0; aidx < (int)Acceptor_.size(); aidx++) {
      aAt = Acceptor_[aidx];
      if (aAt != dAt && mol1 != (*CurrentParm_)[aAt].MolNum()) {
        dist2 = DIST2(frm.Frm().XYZ(aAt), frm.Frm().XYZ(dAt), Image_.ImageType(),
                      frm.Frm().BoxCrd(), ucell_, recip_);
        if (dist2 > dcut2_) continue;
        if (Image_.ImageType() == NOIMAGE)
          angle = CalcAngle(frm.Frm().XYZ(aAt), frm.Frm().XYZ(hAt),
                            frm.Frm().XYZ(dAt));
        else
          angle = ImagedAngle(frm.Frm().XYZ(aAt), frm.Frm().XYZ(hAt),
                              frm.Frm().XYZ(dAt));
        if (angle < acut_) continue;
        ++numHB;
        hbidx = ((didx / 2) * Acceptor_.size()) + aidx;
        dist = sqrt(dist2);
        it = HbondMap_.find( hbidx );
        if (it == HbondMap_.end() ) {
          // New hbond
          HB.A = aAt; HB.D = dAt; HB.H = hAt;
          HB.Frames = 1; HB.dist = dist; HB.angle = angle;
          if (series_) {
            std::string hblegend = CurrentParm_->TruncResAtomName(aAt) + "-" +
                                   CurrentParm_->TruncResAtomName(dAt) + "-" +
                                   (*CurrentParm_)[hAt].Name().Truncated();
            HB.data_ = (DataSet_integer*)
                       masterDSL_->AddSet(DataSet::INTEGER,
                                          MetaData(hbsetname_,"solutehb",hbidx) );
            if (UUseriesout_ != 0) UUseriesout_->AddSet( HB.data_ );
            HB.data_->SetLegend( hblegend );
            HB.data_->AddVal( frameNum, 1 );
          }
          HbondMap_.insert( it, std::pair<int,HbondType>(hbidx, HB) );
        } else {
          // Update existing hbond
          it->second.Frames++;
          it->second.dist += dist;
          it->second.angle += angle;
          if (series_) it->second.data_->AddVal( frameNum, 1 );
        }
      }
    } // END acceptor loop
  } // END donor loop
  NumHbonds_->Add(frameNum, &numHB);
# else
  double dist, dist2, angle;
  HBmapType::iterator it;
  HbondType HB;

  int hbidx = 0; 
  int numHB=0;
  if (noIntramol_) {
    for (HBlistType::const_iterator donor = Donor_.begin(); donor!=Donor_.end(); ++donor) {
      D = (*donor);
      ++donor;
      H = (*donor);
      int mol1 = (*CurrentParm_)[D].MolNum();
      for (HBlistType::const_iterator accept = Acceptor_.begin();
                                      accept != Acceptor_.end(); ++accept, ++hbidx)
      {
        if (*accept == D || mol1 == (*CurrentParm_)[*accept].MolNum()) continue;
        // The following macro handles hbond determination.
        SoluteHbond(*accept, D, H);
      }
    }
  } else {
    for (HBlistType::const_iterator donor = Donor_.begin(); donor!=Donor_.end(); ++donor) {
      D = (*donor);
      ++donor;
      H = (*donor);
      for (HBlistType::const_iterator accept = Acceptor_.begin();
                                accept != Acceptor_.end(); ++accept, ++hbidx)
      {
        if (*accept == D) continue;
        // The following macro handles hbond determination.
        SoluteHbond(*accept, D, H);
      }
    }
  } // END noIntramol
  NumHbonds_->Add(frameNum, &numHB);
  //mprintf("HBOND: Scanned %i hbonds.\n",hbidx);
# endif 
  if (calcSolvent_) {
    // Contains info about which residue(s) a Hbonding solvent mol is
    // Hbonded to.
    std::map< int, std::set<int> > solvent2solute;
    int solventHbonds = 0;
    // SOLUTE DONOR-SOLVENT ACCEPTOR
    // Index by solute H atom. 
    if (hasSolventAcceptor_) {
      numHB = 0;
      for (HBlistType::const_iterator donor = Donor_.begin();
                                      donor != Donor_.end(); ++donor)
      {
        D = (*donor);
        ++donor;
        H = (*donor);
        for (HBlistType::const_iterator accept = SolventAcceptor_.begin();
                                        accept != SolventAcceptor_.end(); ++accept)
        { 
          if (AtomsAreHbonded( frm.Frm(), frameNum, *accept, D, H, H, true )) {
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
      for (HBlistType::const_iterator donor = SolventDonor_.begin();
                                      donor != SolventDonor_.end(); ++donor)
      {
        D = (*donor);
        ++donor;
        H = (*donor);
        for (HBlistType::const_iterator accept = Acceptor_.begin();
                                        accept != Acceptor_.end(); ++accept)
        {
          if (AtomsAreHbonded( frm.Frm(), frameNum, *accept, D, H, *accept, false )) {
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
    for (std::map< int, std::set<int> >::const_iterator bridge = solvent2solute.begin();
                                                        bridge != solvent2solute.end();
                                                      ++bridge)
    {
      // If solvent molecule is bound to 2 or more different residues,
      // it is bridging. 
      if ( bridge->second.size() > 1) {
        ++numHB;
        bridgeID.append(integerToString( bridge->first+1 ) + "("); // Bridging Solvent res 
        for (std::set<int>::const_iterator res = bridge->second.begin();
                                           res != bridge->second.end(); ++res)
          bridgeID.append( integerToString( *res+1 ) + "+" ); // Solute res being bridged
        bridgeID.append("),");
        // Find bridge in map based on this combo of residues (bridge.second)
        BridgeType::iterator b_it = BridgeMap_.find( bridge->second );
        if (b_it == BridgeMap_.end() ) // New Bridge 
          BridgeMap_.insert( b_it, std::pair<std::set<int>,int>(bridge->second, 1) );
        else                           // Increment bridge #frames
          b_it->second++;
      }
    }
    if (bridgeID.empty())
      bridgeID.assign("None");
    NumBridge_->Add(frameNum, &numHB);
    BridgeID_->Add(frameNum, bridgeID.c_str());
  }

  ++Nframes_;

  return Action::OK;
}

#ifdef MPI
// TODO Use in other hbond functions
static inline std::string CreateHBlegend(Topology const& topIn, int a_atom, int h_atom, int d_atom)
{
  if (a_atom == -1)
    return (topIn.TruncResAtomName(h_atom) + "-V");
  else if (d_atom == -1)
    return (topIn.TruncResAtomName(a_atom) + "-V");
  else
    return (topIn.TruncResAtomName(a_atom) + "-" +
            topIn.TruncResAtomName(d_atom) + "-" +
            topIn[h_atom].Name().Truncated());
}


void Action_Hbond::SyncMap(HBmapType& mapIn, const int* rank_frames, const int* rank_offsets,
                           const char* aspect, Parallel::Comm const& commIn)
const
{
  // Need to know how many hbonds on each thread.
  int num_hb = (int)mapIn.size();
  int* nhb_on_rank = 0;
  if (commIn.Master())
    nhb_on_rank = new int[ commIn.Size() ];
  commIn.GatherMaster( &num_hb, 1, MPI_INT, nhb_on_rank );
  std::vector<double> dArray;
  std::vector<int> iArray;
  if (commIn.Master()) {
    for (int rank = 1; rank < commIn.Size(); rank++) {
      if (nhb_on_rank[rank] > 0) {
        mprintf("DEBUG:\tReceiving %i hbonds from rank %i.\n", nhb_on_rank[rank], rank);
        dArray.resize( 2 * nhb_on_rank[rank] );
        iArray.resize( 5 * nhb_on_rank[rank] );
        commIn.Recv( &(dArray[0]), dArray.size(), MPI_DOUBLE, rank, 1300 );
        commIn.Recv( &(iArray[0]), iArray.size(), MPI_INT,    rank, 1301 );
        HbondType HB;
        int ii = 0, id = 0;
        for (int in = 0; in != nhb_on_rank[rank]; in++, ii += 5, id += 2) {
          HBmapType::iterator it = mapIn.find( iArray[ii] ); // hbidx
          if (it == mapIn.end() ) {
            // Hbond on rank that has not been found on master
            HB.dist  = dArray[id  ];
            HB.angle = dArray[id+1];
            HB.A      = iArray[ii+1];
            HB.H      = iArray[ii+2];
            HB.D      = iArray[ii+3];
            HB.Frames = iArray[ii+4];
            HB.data_ = 0;
            mprintf("\tNEW Hbond %i: %i-%i-%i D=%g A=%g %i frames", iArray[ii],
                    HB.A+1, HB.H+1, HB.D+1, HB.dist,
                    HB.angle, HB.Frames);
            if (series_) {
              HB.data_ = (DataSet_integer*)
                         masterDSL_->AddSet( DataSet::INTEGER,
                                             MetaData(hbsetname_, aspect, iArray[ii]) );
              // FIXME: This may be incorrect if CurrentParm_ has changed
              HB.data_->SetLegend( CreateHBlegend(*CurrentParm_, HB.A, HB.H, HB.D) );
              mprintf(" \"%s\"", HB.data_->legend());
            }
            mprintf("\n");
            mapIn.insert( it, std::pair<int,HbondType>(iArray[ii], HB) );
          } else {
            // Hbond on rank and master. Update on master.
            mprintf("\tAPPENDING Hbond %i: %i-%i-%i D=%g A=%g %i frames\n", iArray[ii],
                    it->second.A+1, it->second.H+1, it->second.D+1, dArray[id],
                    dArray[id+1], iArray[ii+4]);
            it->second.dist  += dArray[id  ];
            it->second.angle += dArray[id+1];
            it->second.Frames += iArray[ii+4];
            if (series_)
              HB.data_ = it->second.data_;
          }
          if (series_) {
            HB.data_->Resize( Nframes_ );
            int* d_beg = HB.data_->Ptr() + rank_offsets[ rank ];
            mprintf("\tResizing hbond series data to %i, starting frame %i, # frames %i\n",
                    Nframes_, rank_offsets[rank], rank_frames[rank]);
            commIn.Recv( d_beg, rank_frames[ rank ], MPI_INT, rank, 1304 + in );
            HB.data_->SetSynced();
          }
        } // END master loop over hbonds from rank
      }
    } // END master loop over ranks
    delete[] nhb_on_rank;
  } else {
    if (mapIn.size() > 0) {
      dArray.reserve( 2 * mapIn.size() );
      iArray.reserve( 5 * mapIn.size() );
      for (HBmapType::const_iterator hb = mapIn.begin(); hb != mapIn.end(); ++hb) {
        dArray.push_back( hb->second.dist );
        dArray.push_back( hb->second.angle );
        iArray.push_back( hb->first );
        iArray.push_back( hb->second.A );
        iArray.push_back( hb->second.H );
        iArray.push_back( hb->second.D );
        iArray.push_back( hb->second.Frames );
      }
      commIn.Send( &(dArray[0]), dArray.size(), MPI_DOUBLE, 0, 1300 );
      commIn.Send( &(iArray[0]), iArray.size(), MPI_INT,    0, 1301 );
      // Send series data to master
      if (series_) {
        int in = 0; // For tag
        for (HBmapType::const_iterator hb = mapIn.begin(); hb != mapIn.end(); ++hb, in++) {
          commIn.Send( hb->second.data_->Ptr(), hb->second.data_->Size(),
                                  MPI_INT, 0, 1304 + in );
          hb->second.data_->SetSynced();
        }
      }
    }
  }
}

/** PARALLEL NOTES:
  * The following tags are used for MPI send/receive:
  *   1300  : Array containing hbond double info on rank.
  *   1301  : Array containing hbond integer info on rank.
  *   1302  : Number of bridges to expect from rank.
  *   1303  : Array containing bridge integer info on rank.
  *   1304+X: Array of hbond X series info from rank.
  */
int Action_Hbond::SyncAction(Parallel::Comm const& commIn) {
  // Make sure all time series are updated at this point.
  UpdateSeries();
  // Get total number of frames.
  int* rank_frames = new int[ commIn.Size() ];
  commIn.GatherMaster( &Nframes_, 1, MPI_INT, rank_frames );
  mprintf("DEBUG: Master= %i frames\n", Nframes_);
  for (int rank = 1; rank < commIn.Size(); rank++) {
    mprintf("DEBUG: Rank%i= %i frames\n", rank, rank_frames[ rank ]);
    Nframes_ += rank_frames[ rank ];
  }
  mprintf("DEBUG: Total= %i frames.\n", Nframes_);
  // Convert rank frames to offsets.
  int* rank_offsets = new int[ commIn.Size() ];
  rank_offsets[0] = 0;
  for (int rank = 1; rank < commIn.Size(); rank++) {
    rank_offsets[rank] = rank_offsets[rank-1] + rank_frames[rank-1];
    mprintf("DEBUG:\t\tRank %i offset is %i\n", rank, rank_offsets[rank]);
  }
  //int total_frames = 0;
  //commIn.Reduce( &total_frames, &Nframes_, 1, MPI_INT, MPI_SUM );
  //if (commIn.Master())
  //  Nframes_ = total_frames;
  // Need to send hbond data from all ranks to master.
  SyncMap( HbondMap_, rank_frames, rank_offsets, "solutehb", commIn );
  if (calcSolvent_) {
    SyncMap( SolventMap_, rank_frames, rank_offsets, "solventhb", commIn );
    // iArray will contain for each bridge: Nres, res1, ..., resN, Frames
    std::vector<int> iArray;
    int iSize;
    if (commIn.Master()) {
      for (int rank = 1; rank < commIn.Size(); rank++)
      {
        // Receive size of iArray
        commIn.Recv( &iSize,           1, MPI_INT, rank, 1302 );
        mprintf("DEBUG: Receiving %i bridges from rank %i\n", iSize, rank);
        iArray.resize( iSize );
        commIn.Recv( &(iArray[0]), iSize, MPI_INT, rank, 1303 );
        unsigned int idx = 0;
        while (idx < iArray.size()) {
          std::set<int> residues;
          unsigned int i2 = idx + 1;
          for (int ir = 0; ir != iArray[idx]; ir++, i2++)
            residues.insert( iArray[i2] );
          BridgeType::iterator b_it = BridgeMap_.find( residues );
          if (b_it == BridgeMap_.end() ) // New Bridge 
            BridgeMap_.insert( b_it, std::pair<std::set<int>,int>(residues, iArray[i2]) );
          else                           // Increment bridge #frames
            b_it->second += iArray[i2];
          idx = i2 + 1;
        }
      }
    } else {
       // Construct bridge info array.
       for (BridgeType::const_iterator b = BridgeMap_.begin(); b != BridgeMap_.end(); ++b)
       {
         iArray.push_back( b->first.size() ); // # of bridging res
         for ( std::set<int>::const_iterator r = b->first.begin(); r != b->first.end(); ++r)
           iArray.push_back( *r ); // Bridging res
         iArray.push_back( b->second ); // # frames
      }
      // Since the size of each bridge can be different (i.e. differing #s of
      // residues may be bridged), first send size of the transport array.
      iSize = (int)iArray.size();
      commIn.Send( &iSize,           1, MPI_INT, 0, 1302 );
      commIn.Send( &(iArray[0]), iSize, MPI_INT, 0, 1303 );
    }
  }
  delete[] rank_frames;
  delete[] rank_offsets;
  return 0;
}
#endif

/** Calculate average distance and angle for hbond. */
void Action_Hbond::HbondTypeCalcAvg(HbondType& hb) {
  double dFrames = (double)hb.Frames;
  hb.dist /= dFrames;
  hb.angle /= dFrames;
  hb.angle *= Constants::RADDEG;
}

void Action_Hbond::UpdateSeries() {
  if (seriesUpdated_) return;
  if (series_ && Nframes_ > 0) {
    const int ZERO = 0;
    for (HBmapType::iterator hb = HbondMap_.begin(); hb != HbondMap_.end(); ++hb)
    {
      if ((int)hb->second.data_->Size() < Nframes_)
        hb->second.data_->Add( Nframes_-1, &ZERO );
    }
    for (HBmapType::iterator hb = SolventMap_.begin(); hb != SolventMap_.end(); ++hb)
    {
      if ((int)hb->second.data_->Size() < Nframes_)
        hb->second.data_->Add( Nframes_-1, &ZERO );
    }
  }
  // Should only be called once.
  seriesUpdated_ = true;
} 

// Action_Hbond::Print()
/** Print average occupancies over all frames for all detected Hbonds. */
void Action_Hbond::Print() {
  std::vector<HbondType> HbondList; // For sorting
  std::string Aname, Hname, Dname;

  // Final memory usage
  mprintf("    HBOND: Actual memory usage is %s\n",
          MemoryUsage(HbondMap_.size()+SolventMap_.size(), Nframes_).c_str());
  mprintf("\t%zu solute-solute hydrogen bonds.\n", HbondMap_.size());
  if (calcSolvent_) {
   mprintf("\t%zu solute-solvent hydrogen bonds.\n", SolventMap_.size());
   mprintf("\t%zu unique solute-solvent bridging interactions.\n", BridgeMap_.size());
  }

  // Ensure all series have been updated for all frames.
  UpdateSeries();

  if (CurrentParm_ == 0) return;
  // Calculate necessary column width for strings based on how many residues.
  // ResName+'_'+ResNum+'@'+AtomName | NUM = 4+1+R+1+4 = R+10
  int NUM = DigitWidth( CurrentParm_->Nres() ) + 10;
  // If useAtomNum_ +'_'+AtomNum += 1+A
  if (useAtomNum_) NUM += ( DigitWidth( CurrentParm_->Natom() ) + 1 );

  // Solute Hbonds 
  if (avgout_ != 0) { 
    // Place all detected Hbonds in a list and sort.
    for (HBmapType::const_iterator it = HbondMap_.begin(); it!=HbondMap_.end(); ++it) {
      HbondList.push_back( (*it).second );
      // Calculate average distance and angle for this hbond.
      HbondTypeCalcAvg( HbondList.back() );
    }
    HbondMap_.clear();
    // Sort and Print 
    sort( HbondList.begin(), HbondList.end(), hbond_cmp() );
    avgout_->Printf("%-*s %*s %*s %8s %12s %12s %12s\n", NUM, "#Acceptor", 
                   NUM, "DonorH", NUM, "Donor", "Frames", "Frac", "AvgDist", "AvgAng");
    for (std::vector<HbondType>::const_iterator hbond = HbondList.begin(); 
                                                hbond != HbondList.end(); ++hbond ) 
    {
      double avg = ((double)(*hbond).Frames) / ((double) Nframes_);
      Aname = CurrentParm_->TruncResAtomName((*hbond).A);
      Hname = CurrentParm_->TruncResAtomName((*hbond).H);
      Dname = CurrentParm_->TruncResAtomName((*hbond).D);
      if (useAtomNum_) {
        Aname.append("_" + integerToString((*hbond).A+1));
        Hname.append("_" + integerToString((*hbond).H+1));
        Dname.append("_" + integerToString((*hbond).D+1));
      }
      avgout_->Printf("%-*s %*s %*s %8i %12.4f %12.4f %12.4f\n",
                     NUM, Aname.c_str(), NUM, Hname.c_str(), NUM, Dname.c_str(),
                     (*hbond).Frames, avg, (*hbond).dist, (*hbond).angle);
    }
  }

  // Solute-solvent Hbonds 
  if (solvout_ != 0 && calcSolvent_) {
    HbondList.clear();
    for (HBmapType::const_iterator it = SolventMap_.begin(); it != SolventMap_.end(); ++it) {
      HbondList.push_back( (*it).second );
      // Calculate average distance and angle for this hbond.
      HbondTypeCalcAvg( HbondList.back() );
    }
    SolventMap_.clear();
    sort( HbondList.begin(), HbondList.end(), hbond_cmp() );
    // Calc averages and print
    solvout_->Printf("#Solute-Solvent Hbonds:\n");
    solvout_->Printf("%-*s %*s %*s %8s %12s %12s %12s\n", NUM, "#Acceptor", 
                   NUM, "DonorH", NUM, "Donor", "Count", "Frac", "AvgDist", "AvgAng");
    for (std::vector<HbondType>::iterator hbond = HbondList.begin();
                                          hbond != HbondList.end(); ++hbond )
    {
      // Average has slightly diff meaning since for any given frame multiple
      // solvent can bond to the same solute.
      double avg = ((double)(*hbond).Frames) / ((double) Nframes_);
      if ((*hbond).A==-1) // Solvent acceptor
        Aname = "SolventAcc";
      else {
        Aname = CurrentParm_->TruncResAtomName((*hbond).A);
        if (useAtomNum_) Aname.append("_" + integerToString((*hbond).A+1));
      }
      if ((*hbond).D==-1) { // Solvent donor
        Dname = "SolventDnr";
        Hname = "SolventH";
      } else {
        Dname = CurrentParm_->TruncResAtomName((*hbond).D);
        Hname = CurrentParm_->TruncResAtomName((*hbond).H);
        if (useAtomNum_) {
          Dname.append("_" + integerToString((*hbond).D+1));
          Hname.append("_" + integerToString((*hbond).H+1));
        }
      }
      solvout_->Printf("%-*s %*s %*s %8i %12.4f %12.4f %12.4f\n",
                     NUM, Aname.c_str(), NUM, Hname.c_str(), NUM, Dname.c_str(),
                     (*hbond).Frames, avg, (*hbond).dist, (*hbond).angle);
    }
  }

  // BRIDGING INFO
  if (bridgeout_ != 0 && calcSolvent_) {
    bridgeout_->Printf("#Bridging Solute Residues:\n");
    // Place bridging values in a vector for sorting
    std::vector<std::pair< std::set<int>, int> > bridgevector;
    for (BridgeType::const_iterator it = BridgeMap_.begin();
                                    it != BridgeMap_.end(); ++it)
      bridgevector.push_back( *it );
    std::sort( bridgevector.begin(), bridgevector.end(), bridge_cmp() );
    for (std::vector<std::pair< std::set<int>, int> >::const_iterator bv = bridgevector.begin();
                                                                      bv != bridgevector.end();
                                                                    ++bv)
    {
      bridgeout_->Printf("Bridge Res");
      for (std::set<int>::const_iterator res = bv->first.begin();
                                         res != bv->first.end(); ++res)
        bridgeout_->Printf(" %i:%s", *res+1, CurrentParm_->Res( *res ).c_str());
      bridgeout_->Printf(", %i frames.\n", bv->second);
    } 
  } 
}
