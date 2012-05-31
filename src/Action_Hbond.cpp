// Action_Hbond 
#include <cmath> // sqrt
#include <algorithm> // sort 
#include "Action_Hbond.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG, DEGRAD

// CONSTRUCTOR
Action_Hbond::Action_Hbond() :
  Nframes_(0),
  avgout_(NULL),
  hasDonorMask_(false),
  hasAcceptorMask_(false),
  acut_(0),
  dcut2_(0),
  NumHbonds_(NULL),
  HBavg_(NULL)
{}

// DESTRUCTOR
Action_Hbond::~Action_Hbond() {
  if (HBavg_!=NULL) delete HBavg_;
}

// Action_Hbond::init()
/** Expected call: hbond [out <filename>] <mask> [angle <cut>] [dist <cut>] [avgout <filename>]
  *                      [donormask <mask>] [acceptormask <mask>]
  * Search for Hbonding atoms in region specified by mask. 
  * Arg. check order is:
  * - Keywords
  * - Masks
  * If just <mask> is specified donors and acceptors will be automatically
  * searched for.
  * If donormask is specified but not acceptormask, acceptors will be 
  * automatically searched for in <mask>.
  * If acceptormask is specified but not donormask, donors will be automatically
  * searched for in <mask>.
  * If both donormask and acceptor mask are specified no searching will occur.
  */
int Action_Hbond::init() {
  // Get keywords
  char* outfilename = actionArgs.getKeyString("out",NULL);
  avgout_ = actionArgs.getKeyString("avgout",NULL);
  acut_ = actionArgs.getKeyDouble("angle",135.0);
  // Convert angle cutoff to radians
  acut_ *= DEGRAD;
  double dcut = actionArgs.getKeyDouble("dist",3.0);
  dcut2_ = dcut * dcut;
  // Get donor mask
  char* mask = actionArgs.getKeyString("donormask",NULL);
  if (mask!=NULL) {
    DonorMask_.SetMaskString(mask);
    hasDonorMask_=true;
  }
  // Get acceptor mask
  mask = actionArgs.getKeyString("acceptormask",NULL);
  if (mask!=NULL) {
    AcceptorMask_.SetMaskString(mask);
    hasAcceptorMask_=true;
  }
  // Get generic mask
  mask = actionArgs.getNextMask();
  Mask_.SetMaskString(mask);

  // Setup datasets
  NumHbonds_ = DSL->Add(DataSet::INT, actionArgs.getNextString(),"NumHB");
  if (NumHbonds_==NULL) return 1;
  DFL->Add(outfilename,NumHbonds_);

  mprintf( "  HBOND: ");
  if (!hasDonorMask_ && !hasAcceptorMask_)
    mprintf("Searching for Hbond donors/acceptors in region specified by %s\n",Mask_.MaskString());
  else if (hasDonorMask_ && !hasAcceptorMask_)
    mprintf("Donor mask is %s, acceptors will be searched for in region specified by %s\n",
            DonorMask_.MaskString(), Mask_.MaskString());
  else if (hasAcceptorMask_ && !hasDonorMask_)
    mprintf("Acceptor mask is %s, donors will be searched for in a region specified by %s\n",
            AcceptorMask_.MaskString(), Mask_.MaskString());
  else
    mprintf("Donor mask is %s, Acceptor mask is %s\n",
            DonorMask_.MaskString(), AcceptorMask_.MaskString());
  mprintf( "         Distance cutoff = %.3lf, Angle Cutoff = %.3lf\n",dcut,acut_*RADDEG);
  if (outfilename!=NULL) 
    mprintf( "         Dumping # Hbond v time results to %s\n", outfilename);
  if (avgout_!=NULL)
    mprintf( "         Dumping Hbond avgs to %s\n",avgout_);

  return 0;
}

// Action_Hbond::SearchAcceptor()
/** Search for hbond acceptors X in the region specified by amask.
  * If Auto is true select acceptors based on the rule that "Hydrogen 
  * bonds are FON"
  */
void Action_Hbond::SearchAcceptor(AtomMask& amask, bool Auto) {
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
      if ( (*currentParm)[*atom].Element() == Atom::FLUORINE ||
           (*currentParm)[*atom].Element() == Atom::OXYGEN ||
           (*currentParm)[*atom].Element() == Atom::NITROGEN    )
       isAcceptor=true;
    }
    if (isAcceptor)
      Acceptor_.push_back(*atom);
  }
}

// Action_Hbond::SearchDonor()
/** Search for hydrogen bond donors X-H in the region specified by dmask.
  * If Auto is true select donors based on the rule that "Hydrogen bonds 
  * are FON"
  */
void Action_Hbond::SearchDonor(AtomMask& dmask, bool Auto) {
  bool isDonor;
  // Set up donors: F-H, O-H, N-H
  for (AtomMask::const_iterator donoratom = dmask.begin();
                                donoratom != dmask.end(); ++donoratom)
  {
    // If this is already an H atom continue
    if ( (*currentParm)[*donoratom].Element() == Atom::HYDROGEN ) continue;
    isDonor = true;
    // If auto searching, only consider donor atoms as F, O, N
    if (Auto) {
      isDonor=false;
      if ( (*currentParm)[*donoratom].Element() == Atom::FLUORINE ||
           (*currentParm)[*donoratom].Element() == Atom::OXYGEN ||
           (*currentParm)[*donoratom].Element() == Atom::NITROGEN   )
        isDonor=true;
    }
    if (isDonor) {
      // Get list of hydrogen atoms bonded to this atom
      for (Atom::bond_iterator batom = (*currentParm)[*donoratom].bondbegin();
                               batom != (*currentParm)[*donoratom].bondend();
                               batom++)
      {
        if ( (*currentParm)[*batom].Element() == Atom::HYDROGEN ) {
          //mprintf("BOND TO H: %i@%s -- %i@%s\n",*donoratom+1,(*currentParm)[*donoratom].c_str(),
          //        *batom+1,(*currentParm)[*batom].c_str());
          Donor_.push_back(*donoratom);
          Donor_.push_back(*batom);
        }
      }
    } // END atom is potential donor
  } // END loop over selected atoms
}

// Action_Hbond::setup()
/** Search for hbond donors and acceptors. */
int Action_Hbond::setup() {
  // Set up mask
  if (!hasDonorMask_ || !hasAcceptorMask_) {
    if ( currentParm->SetupIntegerMask( Mask_ ) ) return 1;
    if ( Mask_.None() ) {
      mprintf("Warning: Hbond::setup: Mask has no atoms.\n");
      return 1;
    }
  }
  // Set up donor mask
  if (hasDonorMask_) {
    if ( currentParm->SetupIntegerMask( DonorMask_ ) ) return 1;
    if (DonorMask_.None()) {
      mprintf("Warning: Hbond: DonorMask has no atoms.\n");
      return 1;
    }
  }
  // Set up acceptor mask
  if (hasAcceptorMask_) {
    if ( currentParm->SetupIntegerMask( AcceptorMask_ ) ) return 1;
    if (AcceptorMask_.None()) {
      mprintf("Warning: Hbond: AcceptorMask has no atoms.\n");
      return 1;
    }
  }

  // Four cases:
  // 1) DonorMask and AcceptorMask NULL: donors and acceptors automatically searched for.
  if (!hasDonorMask_ && !hasAcceptorMask_) {
    SearchAcceptor(Mask_,true);
    SearchDonor(Mask_,true);
  
  // 2) DonorMask only: acceptors automatically searched for in Mask
  } else if (hasDonorMask_ && !hasAcceptorMask_) {
    SearchAcceptor(Mask_,true);
    SearchDonor(DonorMask_, false);

  // 3) AcceptorMask only: donors automatically searched for in Mask
  } else if (!hasDonorMask_ && hasAcceptorMask_) {
    SearchAcceptor(AcceptorMask_, false);
    SearchDonor(Mask_,true);

  // 4) Both DonorMask and AcceptorMask: No automatic search.
  } else {
    SearchAcceptor(AcceptorMask_, false);
    SearchDonor(DonorMask_, false);
  }

  // Print acceptor/donor information
  mprintf("\tSet up %i acceptors:\n",(int)Acceptor_.size());
  if (debug>0) {
    for (HBlistType::iterator accept = Acceptor_.begin(); accept!=Acceptor_.end(); accept++)
      mprintf("        %8i: %4s\n",*accept+1,(*currentParm)[*accept].c_str());
  }
  mprintf("\tSet up %i donors:\n",((int)Donor_.size())/2);
  if (debug>0) {
    for (HBlistType::iterator donor = Donor_.begin(); donor!=Donor_.end(); donor++) {
      int atom = (*donor);
      ++donor;
      int a2   = (*donor);
      mprintf("        %8i:%4s - %8i:%4s\n",atom+1,(*currentParm)[atom].c_str(),
              a2+1,(*currentParm)[a2].c_str()); 
    } 
  }

  return 0;
}

// Action_Hbond::action()
/** Calculate distance between all donors and acceptors. Store Hbond info.
  */    
int Action_Hbond::action() {
  // accept ... H-D
  int D, H;
  double dist, dist2, angle;//, ucell[9], recip[9];
  std::map<int,HbondType>::iterator it;
  HbondType HB;

  int Nhb = 0; 
  int numHB=0;
  for (HBlistType::iterator donor = Donor_.begin(); donor!=Donor_.end(); ++donor) {
    D = (*donor);
    ++donor;
    H = (*donor);
    for (HBlistType::iterator accept = Acceptor_.begin(); accept!=Acceptor_.end(); ++accept, ++Nhb) 
    {
      if (*accept == D) continue;
      dist2 = currentFrame->DIST2(*accept, D);
      //dist2 = currentFrame->DIST2(*accept, D, (int)P->boxType, ucell, recip);
      if (dist2 > dcut2_) continue;
      angle = currentFrame->ANGLE(*accept, H, D);
      if (angle < acut_) continue;
//      mprintf( "HBOND[%i]: %i:%s ... %i:%s-%i:%s Dist=%lf Angle=%lf\n", 
//              Nhb, *accept, P->names[*accept],
//              H, P->names[H], D, P->names[D], dist, angle);
      ++numHB;
      dist = sqrt(dist2);
      // Find hbond in map
      it = HbondMap_.find( Nhb );
      if (it == HbondMap_.end() ) {
        // New Hbond
        HB.A=*accept;
        HB.D=D;
        HB.H=H;
        HB.Frames = 1;
        HB.dist=dist;
        HB.angle=angle;
        HbondMap_.insert( it, std::pair<int,HbondType>(Nhb, HB) );
      } else {
        (*it).second.Frames++;
        (*it).second.dist+=dist;
        (*it).second.angle+=angle;
      }
    }
  }
  NumHbonds_->Add(frameNum, &numHB);
//  mprintf("HBOND: Scanned %i hbonds.\n",Nhb);
  ++Nframes_;

  return 0;
}

// Action_Hbond::print()
/** Print average occupancies over all frames for all detected Hbonds
  */
void Action_Hbond::print() {
  std::vector<HbondType> HbondList;
  DataFile *hbavgFile;

  // If avgout is NULL no averaging.
  if (avgout_==NULL) return;

  // Set up data set list for all avg-related data.
  HBavg_ = new DataSetList(); 
  DFL->Add(avgout_, HBavg_->Add(DataSet::STRING, (char*)"Acceptor", "Acceptor"));
  DFL->Add(avgout_, HBavg_->Add(DataSet::STRING, (char*)"DonorH", "DonorH"));
  DFL->Add(avgout_, HBavg_->Add(DataSet::STRING, (char*)"Donor", "Donor"));
  DFL->Add(avgout_, HBavg_->Add(DataSet::INT, (char*)"Frames", "Frames"));
  DFL->Add(avgout_, HBavg_->Add(DataSet::DOUBLE, (char*)"Frac", "Frac"));
  DFL->Add(avgout_, HBavg_->Add(DataSet::DOUBLE, (char*)"AvgDist", "AvgDist"));
  hbavgFile = DFL->Add(avgout_, HBavg_->Add(DataSet::DOUBLE, (char*)"AvgAng", "AvgAng"));

  // Place all detected Hbonds in a list and sort 
  for (std::map<int,HbondType>::iterator it = HbondMap_.begin(); it!=HbondMap_.end(); ++it) 
    HbondList.push_back( (*it).second );
  sort( HbondList.begin(), HbondList.end(), hbond_cmp() );
  // Calculate averages
  int hbondnum=0;
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

    std::string Aname = currentParm->ResAtomName((*hbond).A);
    std::string Hname = currentParm->ResAtomName((*hbond).H);
    std::string Dname = currentParm->ResAtomName((*hbond).D);
    // TODO: DataSetList should accept string
    HBavg_->AddData(hbondnum, (char*)Aname.c_str(), 0);
    HBavg_->AddData(hbondnum, (char*)Hname.c_str(), 1);
    HBavg_->AddData(hbondnum, (char*)Dname.c_str(), 2);
    HBavg_->AddData(hbondnum, &((*hbond).Frames), 3);
    HBavg_->AddData(hbondnum, &avg, 4);
    HBavg_->AddData(hbondnum, &dist, 5);
    HBavg_->AddData(hbondnum, &angle, 6);
    ++hbondnum;
  }
  hbavgFile->ProcessArgs("noxcol");
}
