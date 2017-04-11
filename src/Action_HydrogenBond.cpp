#include "Action_HydrogenBond.h"
#include "CpptrajStdio.h"
#include "Constants.h"

// CONSTRUCTOR
Action_HydrogenBond::Action_HydrogenBond() :
  CurrentParm_(0),
  masterDSL_(0),
  NumHbonds_(0),
  NumSolvent_(0),
  NumBridge_(0),
  BridgeID_(0),
  UUseriesout_(0),
  UVseriesout_(0),
  avgout_(0),
  solvout_(0),
  bridgeout_(0),
  dcut2_(0.0),
  acut_(0.0),
  debug_(0),
  series_(0),
  useAtomNum_(false),
  noIntramol_(false),
  hasDonorMask_(false),
  hasDonorHmask_(false),
  hasAcceptorMask_(false),
  hasSolventDonor_(false),
  calcSolvent_(false),
  hasSolventAcceptor_(false)
{}

// void Action_HydrogenBond::Help()
void Action_HydrogenBond::Help() const {
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

// Action_HydrogenBond::Init()
Action::RetType Action_HydrogenBond::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  trajComm_ = init.TrajComm();
# endif
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
  mprintf("\tDistance cutoff = %.3f, Angle Cutoff = %.3f\n",dcut,acut_*Constants::RADDEG);
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

inline bool IsFON(Atom const& atm) {
  return (atm.Element() == Atom::FLUORINE ||
          atm.Element() == Atom::OXYGEN ||
          atm.Element() == Atom::NITROGEN);
}

// Action_HydrogenBond::Setup()
Action::RetType Action_HydrogenBond::Setup(ActionSetup& setup) {
  CurrentParm_ = setup.TopAddress();
  Image_.SetupImaging( setup.CoordInfo().TrajBox().Type() );

  // Set up generic mask
  if (!hasDonorMask_ || !hasAcceptorMask_) {
    if ( setup.Top().SetupIntegerMask( Mask_ ) ) return Action::ERR;
    if ( Mask_.None() ) {
      mprintf("Warning: Mask has no atoms.\n");
      return Action::SKIP;
    }
  }

  // ACCEPTOR MASK SETUP
  CharMask Amask;
  if (hasAcceptorMask_) {
    // Acceptor mask specified
    if ( setup.Top().SetupIntegerMask( AcceptorMask_ ) ) return Action::ERR;
    if (AcceptorMask_.None()) {
      mprintf("Warning: AcceptorMask has no atoms.\n");
      return Action::SKIP;
    }
    Acceptor_ = AcceptorMask_.Selected();
    // Set up the char mask
    Amask = CharMask(AcceptorMask_.ConvertToCharMask(), AcceptorMask_.Nselected());
  } else {
    // No specified acceptor mask; search generic mask.
    // Setup the char mask
    Amask = CharMask( setup.Top().Natom() );
    for (AtomMask::const_iterator at = Mask_.begin(); at != Mask_.end(); ++at) {
      if (IsFON( setup.Top()[*at] )) {
        Amask.Select(*at);
        Acceptor_.push_back( *at );
      }
    }
  }

  // DONOR/ACCEPTOR SETUP
  if (hasDonorMask_) {
    // Donor heavy atom mask specified
    if ( setup.Top().SetupIntegerMask( DonorMask_ ) ) return Action::ERR;
    if (DonorMask_.None()) {
      mprintf("Warning: DonorMask has no atoms.\n");
      return Action::SKIP;
    }
    if ( hasDonorHmask_ ) {
      // Donor hydrogen mask also specified
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
      for (int midx = 0; midx != DonorMask_.Nselected(); midx++)
        if ( Amask.AtomInCharMask(DonorMask_[midx]) )
          Both_.push_back( Site(DonorMask_[midx], DonorHmask_[midx]) );
        else
          Donor_.push_back( Site(DonorMask_[midx], DonorHmask_[midx]) );
    } else {
      // No donor hydrogen mask; use any hydrogens bonded to donor heavy atoms.
      for (AtomMask::const_iterator at = DonorMask_.begin(); at != DonorMask_.end(); ++at) {
        Iarray Hatoms;
        for (Atom::bond_iterator batom = setup.Top()[*at].bondbegin();
                                 batom != setup.Top()[*at].bondend(); ++batom)
          if (setup.Top()[*batom].Element() == Atom::HYDROGEN)
            Hatoms.push_back( *batom );
        if (!Hatoms.empty()) {
          if ( Amask.AtomInCharMask(*at) )
            Both_.push_back( Site(*at, Hatoms) );
          else
            Donor_.push_back( Site(*at, Hatoms) );
        }  
      }
    }
  } else {
    // No specified donor mask; search generic mask.
    for (AtomMask::const_iterator at = Mask_.begin(); at != Mask_.end(); ++at) {
      if (IsFON( setup.Top()[*at] )) {
        Iarray Hatoms;
        for (Atom::bond_iterator batom = setup.Top()[*at].bondbegin();
                                 batom != setup.Top()[*at].bondend(); ++batom)
          if (setup.Top()[*batom].Element() == Atom::HYDROGEN)
            Hatoms.push_back( *batom );
        if (!Hatoms.empty()) {
          if ( Amask.AtomInCharMask(*at) )
            Both_.push_back( Site(*at, Hatoms) );
          else
            Donor_.push_back( Site(*at, Hatoms) );
        }
      }
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

  mprintf("Acceptor atoms (%zu):\n", Acceptor_.size());
  for (Iarray::const_iterator at = Acceptor_.begin(); at != Acceptor_.end(); ++at)
    mprintf("\t%20s %8i\n", setup.Top().TruncResAtomName(*at).c_str(), *at+1);
  mprintf("Donor/acceptor sites (%zu):\n", Both_.size());
  for (Sarray::const_iterator si = Both_.begin(); si != Both_.end(); ++si) {
    mprintf("\t%20s %8i", setup.Top().TruncResAtomName(si->Idx()).c_str(), si->Idx()+1);
    for (Iarray::const_iterator at = si->Hbegin(); at != si->Hend(); ++at)
      mprintf(" %s", setup.Top()[*at].c_str());
    mprintf("\n");
  }
  mprintf("Donor sites (%zu):\n", Donor_.size());
  for (Sarray::const_iterator si = Donor_.begin(); si != Donor_.end(); ++si) {
    mprintf("\t%20s %8i", setup.Top().TruncResAtomName(si->Idx()).c_str(), si->Idx()+1);
    for (Iarray::const_iterator at = si->Hbegin(); at != si->Hend(); ++at)
      mprintf(" %s", setup.Top()[*at].c_str());
    mprintf("\n");
  }


  return Action::OK;
}

// Action_HydrogenBond::DoAction()
Action::RetType Action_HydrogenBond::DoAction(int frameNum, ActionFrame& frm) {

  return Action::OK;
}

// Action_HydrogenBond::Print()
void Action_HydrogenBond::Print() {
}
