#include <cmath> // sqrt
#include <cfloat> // DBL_MAX
#include "Action_NativeContacts.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"

// CONSTRUCTOR
Action_NativeContacts::Action_NativeContacts() :
  distance_(7.0),
  debug_(0),
  first_(false),
  byResidue_(false),
  numnative_(0),
  nonnative_(0),
  mindist_(0),
  maxdist_(0),
  CurrentParm_(0)
{}

void Action_NativeContacts::Help() {
  mprintf("\t[<mask1> [<mask2>]]\n"
          "\t[noimage] [distance <cut>] [first] [out <filename>]\n"
          "\t[name <dsname>] [mindist] [maxdist] [byresidue]\n");
}

// Action_NativeContacts::SetupList()
Action_NativeContacts::Marray Action_NativeContacts::SetupList(Topology const& parmIn, 
                                        AtomMask const& mask) const
{
  Marray listOut;
  if (byResidue_) {
    AtomMask::const_iterator selected = mask.begin();
    for (Topology::res_iterator res = parmIn.ResStart(); res != parmIn.ResEnd(); ++res) {
      AtomMask resMask;
      for (int atnum = res->FirstAtom(); atnum != res->LastAtom(); ++atnum) {
        // If we are at the last selected atom we are done
        if (selected == mask.end()) break;
        if ( *selected == atnum ) {
          resMask.AddAtom( atnum );
          ++selected;
        }
      }
      if (!resMask.None())
        listOut.push_back( resMask );
      // If we are at the last selected atom we are done
      if (selected == mask.end()) break;
    }
  } else { // FIXME: Just use masks? If so get rid of this AtomMask constructor
    for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom)
      listOut.push_back( AtomMask( *atom ) );
  }
  return listOut;
}

// DEBUG
static void DebugContactList(std::vector<AtomMask> const& clist, Topology const& parmIn) {
  for (std::vector<AtomMask>::const_iterator c = clist.begin();
                                             c != clist.end(); ++c)
  {
    mprintf("\tPotential Contact %u:", c - clist.begin());
    for (AtomMask::const_iterator m = c->begin(); m != c->end(); ++m)
      mprintf(" %s", parmIn.AtomMaskName(*m).c_str());
    mprintf("\n");
  }
}

/** Determine the minimum and maximum distances given two masks.
  * \return the minimum distance between the masks.
  */
double Action_NativeContacts::GetMinMax(AtomMask const& m1, AtomMask const& m2,
                                        Frame const& frameIn, 
                                        double& minDist2, double& maxDist2) const
{
  double min = DBL_MAX;
  for (AtomMask::const_iterator a1 = m1.begin(); a1 != m1.end(); ++a1) {
    for (AtomMask::const_iterator a2 = m2.begin(); a2 != m2.end(); ++a2) {
      double Dist2 = DIST2(frameIn.XYZ(*a1), frameIn.XYZ(*a2), image_.ImageType(),
                           frameIn.BoxCrd(), ucell_, recip_);
      if (Dist2 < min) {
        min = Dist2;
        minDist2 = Dist2;
      }
      if (Dist2 > maxDist2) maxDist2 = Dist2;
    }
  }
  return min;
}

/** Based on select atoms in Mask1 (and optionally Mask2), set up
  * contact lists by atom or residue.
  */
int Action_NativeContacts::SetupContactLists(Topology const& parmIn, Frame const& frameIn) {
  if ( parmIn.SetupIntegerMask( Mask1_, frameIn ) ) return 1;
  Mask1_.MaskInfo();
  // Setup first contact list
  contactList1_ = SetupList(parmIn, Mask1_);
  if (contactList1_.empty()) {
    mprinterr("Warning: Nothing selected for '%s'\n", Mask1_.MaskString());
    return 1;
  }
  if (debug_ > 0) DebugContactList( contactList1_, parmIn );
  // Setup second contact list if necessary
  if ( Mask2_.MaskStringSet() ) {
    if (parmIn.SetupIntegerMask( Mask2_, frameIn ) ) return 1;
    Mask2_.MaskInfo();
    // Warn if masks overlap
    int nOverlap = Mask1_.NumAtomsInCommon( Mask2_ );
    if (nOverlap > 0)
      mprintf("Warning: Masks '%s' and '%s' overlap by %i atoms.\n", 
              Mask1_.MaskString(), Mask2_.MaskString(), nOverlap);
    contactList2_ = SetupList(parmIn, Mask2_);
    if (contactList2_.empty()) {
      mprinterr("Warning: Nothing selected for '%s'\n", Mask2_.MaskString());
      return 1;
    }
    if (debug_ > 0) DebugContactList( contactList2_, parmIn );
  }

  // Loop over contact lists. Determine what satisfies the cutoff.
  double maxDist2 = 0.0, minDist2 = 0.0;
  nativeContacts_.clear();
  if ( Mask2_.MaskStringSet() ) {
    for (unsigned int c1 = 0; c1 < contactList1_.size(); ++c1)
      for (unsigned int c2 = 0; c2 < contactList2_.size(); ++c2)
      {
        double Dist2 = GetMinMax(contactList1_[c1], contactList2_[c2], frameIn, minDist2, maxDist2);
        if (Dist2 < distance_)
          nativeContacts_.insert( contactType(c1, c2) );
      }
  } else {
    for (unsigned int c1 = 0; c1 < contactList1_.size(); ++c1)
      for (unsigned int c2 = c1 + 1; c2 < contactList1_.size(); ++c2)
      {
        double Dist2 = GetMinMax(contactList1_[c1], contactList1_[c2], frameIn, minDist2, maxDist2);
        //mprintf("\tContact %u to %u mindist = %f\n", c1, c2, sqrt(Dist2));
        if (Dist2 < distance_)
          nativeContacts_.insert( contactType(c1, c2) );
      }
  }
  //mprintf("\tMinimum observed distance= %f, maximum observed distance= %f\n",
  //        sqrt(minDist2), sqrt(maxDist2));
  // DEBUG - Print contacts
  mprintf("\tSetup %zu native contacts:\n", nativeContacts_.size());
  for (contactListType::const_iterator contact = nativeContacts_.begin();
                                       contact != nativeContacts_.end(); ++contact)
  {
    int a1 = contactList1_[contact->first][0];
    int a2;
    if (Mask2_.MaskStringSet())
      a2 = contactList2_[contact->second][0];
    else
      a2 = contactList1_[contact->second][0];
    if (byResidue_) {
      int r1 = parmIn[a1].ResNum();
      int r2 = parmIn[a2].ResNum();
      mprintf("\t\tResidue %s to %s\n",
              parmIn.TruncResNameNum(r1).c_str(),
              parmIn.TruncResNameNum(r2).c_str());
    } else
      mprintf("\t\tAtom %i[%s] to %i[%s]\n", a1+1,
              parmIn.AtomMaskName(a1).c_str(),
              parmIn.AtomMaskName(a2).c_str());
  }
    
  return 0;
}
// -----------------------------------------------------------------------------
// Action_NativeContacts::Init()
Action::RetType Action_NativeContacts::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  // Get Keywords
  image_.InitImaging( !(actionArgs.hasKey("noimage")) );
  double dist = actionArgs.getKeyDouble("distance", 7.0);
  byResidue_ = actionArgs.hasKey("byresidue");
  // Square the cutoff
  distance_ = dist * dist;
  first_ = actionArgs.hasKey("first");
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  // Get reference
  ReferenceFrame REF = FL->GetFrameFromArgs( actionArgs );
  if (!first_ && REF.error()) return Action::ERR;
  // Create data sets
  std::string name = actionArgs.GetStringKey("name");
  if (name.empty())
    name = DSL->GenerateDefaultName("Contacts");
  numnative_ = DSL->AddSetAspect(DataSet::INTEGER, name, "Native");
  nonnative_ = DSL->AddSetAspect(DataSet::INTEGER, name, "NonNative");
  if (outfile != 0) {
    outfile->AddSet(numnative_);
    outfile->AddSet(nonnative_);
  }
  if (numnative_ == 0 || nonnative_ == 0) return Action::ERR;
  if (actionArgs.hasKey("mindist")) {
    mindist_ = DSL->AddSetAspect(DataSet::DOUBLE, name, "MinDist");
    if (mindist_ == 0) return Action::ERR;
    if (outfile != 0) outfile->AddSet(mindist_);
  }
  if (actionArgs.hasKey("maxdist")) {
    maxdist_ = DSL->AddSetAspect(DataSet::DOUBLE, name, "MaxDist");
    if (maxdist_ == 0) return Action::ERR;
    if (outfile != 0) outfile->AddSet(maxdist_);
  }
  // Get Masks
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );
  std::string mask2 = actionArgs.GetMaskNext();
  if (!mask2.empty())
    Mask2_.SetMaskString( mask2 );
  // Set up reference if necessary.
  if (!first_ && REF.empty()) {
    mprintf("\tNo reference structure specified. Defaulting to first.\n");
    first_ = true;
  }
  if (!first_) {
    // Set up imaging info for ref parm
    image_.SetupImaging( REF.Parm()->BoxType() );
    if (image_.ImageType() == NONORTHO)
      REF.Coord()->BoxCrd().ToRecip(ucell_, recip_);
    if (SetupContactLists( *(REF.Parm()), *(REF.Coord()) )) return Action::ERR;
  }
  mprintf("    NATIVECONTACTS: Mask1='%s'", Mask1_.MaskString());
  if (Mask2_.MaskStringSet())
    mprintf(" Mask2='%s'", Mask2_.MaskString());
  mprintf(", contacts set up based on");
  if (first_)
    mprintf(" first frame.\n");
  else
    mprintf("'%s'.\n", REF.FrameName().base());
  if (byResidue_)
    mprintf("\tContacts will be determined by residue\n");
  mprintf("\tDistance cutoff is %g Ang,", sqrt(distance_));
  if (image_.ImageType() == NOIMAGE)
    mprintf(" imaging is off.\n");
  else
    mprintf(" imaging is on.\n");
  mprintf("\tData set name: %s\n", name.c_str());
  if (maxdist_ != 0)
    mprintf("\tSaving maximum observed distances to '%s'\n", maxdist_->Legend().c_str());
  if (mindist_ != 0)
    mprintf("\tSaving minimum observed distances to '%s'\n", mindist_->Legend().c_str());
  if (outfile != 0)
    mprintf("\tOutput to '%s'\n", outfile->DataFilename().full());

  return Action::OK;
}

// Action_NativeContacts::Setup()
Action::RetType Action_NativeContacts::Setup(Topology* currentParm, Topology** parmAddress) {
  // Set up imaging info for this parm
  image_.SetupImaging( currentParm->BoxType() );
  if (image_.ImagingEnabled())
    mprintf("\tImaging enabled.\n");
  else
    mprintf("\tImaging disabled.\n");
  CurrentParm_ = currentParm;
  return Action::OK;
}

// Action_NativeContacts::DoAction()
Action::RetType Action_NativeContacts::DoAction(int frameNum, Frame* currentFrame,
                                                Frame** frameAddress)
{
  if (image_.ImageType() == NONORTHO) currentFrame->BoxCrd().ToRecip(ucell_, recip_);
  if (first_) {
    if (SetupContactLists( *CurrentParm_, *currentFrame )) return Action::ERR;
    first_ = false;
  }
  // Loop over all potential contacts
  int Nnative = 0;
  int NnonNative = 0;
  double maxDist2 = 0.0;
  double minDist2 = 0.0;
  if ( Mask2_.MaskStringSet() ) {
    for (unsigned int c1 = 0; c1 < contactList1_.size(); ++c1)
      for (unsigned int c2 = 0; c2 < contactList2_.size(); ++c2)
      {
        double Dist2 = GetMinMax(contactList1_[c1], contactList2_[c2],
                                 *currentFrame, minDist2, maxDist2);
        if (Dist2 < distance_) { // Potential contact
          if (nativeContacts_.find( contactType(c1, c2) ) != nativeContacts_.end())
            ++Nnative;    // Native contact
          else
            ++NnonNative; // Non-native contact
        }
      }
  } else {
    for (unsigned int c1 = 0; c1 < contactList1_.size(); ++c1)
      for (unsigned int c2 = c1 + 1; c2 < contactList1_.size(); ++c2)
      {
        double Dist2 = GetMinMax(contactList1_[c1], contactList1_[c2],
                                 *currentFrame, minDist2, maxDist2);
        //mprintf("DEBUG: Dist for %s <-> %s: %f", 
        //        CurrentParm_->AtomMaskName(contactList1_[c1][0]).c_str(),
        //        CurrentParm_->AtomMaskName(contactList1_[c2][0]).c_str(), sqrt(Dist2));
        if (Dist2 < distance_) { // Potential contact
          if (nativeContacts_.find( contactType(c1, c2) ) != nativeContacts_.end()) {
            ++Nnative;    // Native contact
            //mprintf(" NATIVE");
          } else {
            ++NnonNative; // Non-native contact
            //mprintf(" NON-NATIVE");
          }
        }
        //mprintf("\n");
      }
  }
  numnative_->Add(frameNum, &Nnative);
  nonnative_->Add(frameNum, &NnonNative);
  if (mindist_ != 0) {
    minDist2 = sqrt(minDist2);
    mindist_->Add(frameNum, &minDist2);
  }
  if (maxdist_ != 0) {
    maxDist2 = sqrt(maxDist2);
    maxdist_->Add(frameNum, &maxDist2);
  }
  return Action::OK;
}
