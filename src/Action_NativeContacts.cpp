#include <cmath> // sqrt
#include <cfloat> // DBL_MAX
#include <cstdlib> // abs, intel 11 compilers choke on std::abs
#include <algorithm> // std::max, std::sort
#include "Action_NativeContacts.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "PDBfile.h"

// CONSTRUCTOR
Action_NativeContacts::Action_NativeContacts() :
  distance_(7.0),
  pdbcut_(0.0),
  debug_(0),
  matrix_min_(0),
  resoffset_(1),
  nframes_(0),
  first_(false),
  byResidue_(false),
  includeSolvent_(false),
  series_(false),
  usepdbcut_(false),
  seriesUpdated_(false),
  saveNonNative_(false),
  cfile_(0), pfile_(0), nfile_(0), rfile_(0),
  seriesout_(0),
  seriesNNout_(0),
  numnative_(0),
  nonnative_(0),
  mindist_(0),
  maxdist_(0),
  nativeMap_(0),
  nonnatMap_(0),
  CurrentParm_(0),
  refParm_(0),
  masterDSL_(0)
{}
// TODO: mapout, avg contacts over traj, 1=native, -1=nonnative
void Action_NativeContacts::Help() const {
  mprintf("\t[<mask1> [<mask2>]] [writecontacts <outfile>] [resout <resfile>]\n"
          "\t[noimage] [distance <cut>] [out <filename>] [includesolvent]\n"
          "\t[ first | %s ]\n"
          "\t[resoffset <n>] [contactpdb <file>] [pdbcut <cut>] [mindist] [maxdist]\n"
          "\t[name <dsname>] [byresidue] [map [mapout <mapfile>]] [series [seriesout <file>]]\n"
          "\t[savenonnative [seriesnnout <file>] [nncontactpdb <file>]]\n"
          "  Calculate number of contacts in <mask1>, or between <mask1> and <mask2>\n"
          "  if both are specified. Native contacts are determined based on the given\n"
          "  reference structure (or first frame if not specified) and the specified\n"
          "  distance cut-off (7.0 Ang. default). If [byresidue] is specified contacts\n"
          "  between two residues spaced <resoffset> residues apart are ignored, and\n"
          "  the map (if specified) is written per-residue.\n", DataSetList::RefArgs);
}

/** Set up atom/residue indices corresponding to atoms selected in mask.
  * This is done to make creating an atom/residue contact map easier.
  */
Action_NativeContacts::Iarray Action_NativeContacts::SetupContactIndices(
                                AtomMask const& mask, Topology const& parmIn)
{
  Iarray contactIdx;
  for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom)
    if (byResidue_)
      contactIdx.push_back( parmIn[*atom].ResNum() );
    else
      contactIdx.push_back( *atom );
  return contactIdx;
}

// DEBUG
static void DebugContactList(AtomMask const& mask, Topology const& parmIn)
{
  for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom)
    mprintf("\tPotential Contact %u: %s\n", atom - mask.begin(),
            parmIn.AtomMaskName(*atom).c_str());
}

/** Remove any selected solvent atoms from mask. */
static void removeSelectedSolvent( Topology const& parmIn, AtomMask& mask ) {
  AtomMask newMask = mask;
  newMask.ClearSelected();
  for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom) {
    int molnum = parmIn[*atom].MolNum();
    if (!parmIn.Mol(molnum).IsSolvent())
      newMask.AddSelectedAtom( *atom );
  }
  mask = newMask;
}

/** Based on selected atoms in Mask1 (and optionally Mask2), set up
  * potential contact lists. Also set up atom/residue indices corresponding
  * to each potential contact.
  */
int Action_NativeContacts::SetupContactLists(Topology const& parmIn, Frame const& frameIn)
{
  // Setup first contact list
  if ( parmIn.SetupIntegerMask( Mask1_, frameIn ) ) return 1;
  if (!includeSolvent_) removeSelectedSolvent( parmIn, Mask1_ );
  Mask1_.MaskInfo();
  if (Mask1_.None())
  {
    mprinterr("Warning: Nothing selected for '%s'\n", Mask1_.MaskString());
    return 1;
  }
  if (debug_ > 0) DebugContactList( Mask1_, parmIn );
  contactIdx1_ = SetupContactIndices( Mask1_, parmIn );
  // Setup second contact list if necessary
  if ( Mask2_.MaskStringSet() ) {
    if (parmIn.SetupIntegerMask( Mask2_, frameIn ) ) return 1;
    if (!includeSolvent_) removeSelectedSolvent( parmIn, Mask2_ );
    Mask2_.MaskInfo();
    if (Mask2_.None())
    {
      mprinterr("Warning: Nothing selected for '%s'\n", Mask2_.MaskString());
      return 1;
    }
    // Warn if masks overlap
    int nOverlap = Mask1_.NumAtomsInCommon( Mask2_ );
    if (nOverlap > 0) {
      mprintf("Warning: Masks '%s' and '%s' overlap by %i atoms.\n"
              "Warning: Some contacts may be double-counted.\n", 
              Mask1_.MaskString(), Mask2_.MaskString(), nOverlap);
      if (mindist_ != 0)
        mprintf("Warning: Minimum distance will always be 0.0\n");
    }
    if (debug_ > 0) DebugContactList( Mask2_, parmIn );
    contactIdx2_ = SetupContactIndices( Mask2_, parmIn );
  }
  return 0;
}

/** This macro is used by DetermineNativeContacts to set up a new contact
  * if it is valid.
  */
#define SetNativeContact() { \
        if (ValidContact(*c1, *c2, parmIn)) { \
          double Dist2 = DIST2(fIn.XYZ(*c1), fIn.XYZ(*c2), image_.ImageType(), \
                               fIn.BoxCrd(), ucell_, recip_); \
          minDist2 = std::min( Dist2, minDist2 ); \
          maxDist2 = std::max( Dist2, maxDist2 ); \
          if (Dist2 < distance_) { \
            std::string legend(parmIn.AtomMaskName(*c1) + "_" + parmIn.AtomMaskName(*c2)); \
            int r1 = parmIn[*c1].ResNum(); \
            int r2 = parmIn[*c2].ResNum(); \
            ret = nativeContacts_.insert( Mpair(Cpair(*c1,*c2), contactType(legend,r1,r2)) ); \
            if (ret.second && series_) {\
              MetaData md(numnative_->Meta().Name(), "NC", nativeContacts_.size()); \
              md.SetLegend( legend ); \
              DataSet* ds = masterDSL_->AddSet(DataSet::INTEGER, md); \
              ret.first->second.SetData( ds ); \
              if (seriesout_ != 0) seriesout_->AddDataSet( ds ); \
            } \
          } \
        } \
}

// Action_NativeContacts::DetermineNativeContacts()
/** Determine potential contacts for given Topology and Frame, then determine 
  * which pairs of contacts satisfy the cutoff and set those as native contacts.
  * Should only be called once.
  */
#ifdef MPI
int Action_NativeContacts::DetermineNativeContacts(Topology const& parmIn, Frame const& frameIn)
#else
int Action_NativeContacts::DetermineNativeContacts(Topology const& parmIn, Frame const& fIn)
#endif
{
# ifdef MPI
  Frame fIn = frameIn;
  if (trajComm_.Size() > 1) {
    // Ensure all threads use same reference
    if (trajComm_.Master())
      for (int rank = 1; rank < trajComm_.Size(); rank++)
        fIn.SendFrame( rank, trajComm_ );
    else
      fIn.RecvFrame( 0, trajComm_ );
  }
# endif
  if (pfile_ != 0 || nfile_ != 0) {
    refFrame_ = fIn; // Save frame for later PDB output.
    refParm_ = &parmIn;  // Save parm for later PDB output.
  }
  if ( SetupContactLists(parmIn, fIn) ) return 1;
  // If specified, set up contacts maps; base size on atom masks.
  if (nativeMap_ != 0) {
    int matrix_max;
    if (Mask2_.MaskStringSet()) {
      matrix_min_ = std::min( Mask1_[0], Mask2_[0] );
      matrix_max = std::max( Mask1_.back(), Mask2_.back() );
    } else {
      matrix_min_ = Mask1_[0];
      matrix_max = Mask1_.back();
    }
    std::string label("Atom");
    if (byResidue_) {
      matrix_min_ = parmIn[matrix_min_].ResNum();
      matrix_max = parmIn[matrix_max].ResNum();
      label.assign("Residue");
    }
    int matrix_cols = matrix_max - matrix_min_ + 1;
    if (nativeMap_->AllocateHalf( matrix_cols )) return 1;
    if (nonnatMap_->AllocateHalf( matrix_cols )) return 1;
    Dimension matrix_dim( matrix_min_+1, 1, label );
    nativeMap_->SetDim(Dimension::X, matrix_dim);
    nativeMap_->SetDim(Dimension::Y, matrix_dim);
    nonnatMap_->SetDim(Dimension::X, matrix_dim);
    nonnatMap_->SetDim(Dimension::Y, matrix_dim);
  }
  double maxDist2 = 0.0;
  double minDist2 = DBL_MAX;
  nativeContacts_.clear();
  std::pair<contactListType::iterator, bool> ret; 
  if ( Mask2_.MaskStringSet() ) {
    for (AtomMask::const_iterator c1 = Mask1_.begin(); c1 != Mask1_.end(); ++c1)
      for (AtomMask::const_iterator c2 = Mask2_.begin(); c2 != Mask2_.end(); ++c2)
      {
        SetNativeContact();
      }
  } else {
    for (AtomMask::const_iterator c1 = Mask1_.begin(); c1 != Mask1_.end(); ++c1)
      for (AtomMask::const_iterator c2 = c1 + 1; c2 != Mask1_.end(); ++c2)
      {
        SetNativeContact();
      }
  }
  //mprintf("\tMinimum observed distance= %f, maximum observed distance= %f\n",
  //        sqrt(minDist2), sqrt(maxDist2));
  // Print contacts
  mprintf("\tSetup %zu native contacts:\n", nativeContacts_.size());
  for (contactListType::const_iterator contact = nativeContacts_.begin();
                                       contact != nativeContacts_.end(); ++contact)
  {
    int a1 = contact->first.first;
    int a2 = contact->first.second;
    mprintf("\t\tAtom '%s' to '%s'\n", parmIn.AtomMaskName(a1).c_str(),
            parmIn.AtomMaskName(a2).c_str());
  }
  return 0;  
}
// -----------------------------------------------------------------------------
// Action_NativeContacts::Init()
Action::RetType Action_NativeContacts::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  trajComm_ = init.TrajComm();
# endif
  masterDSL_ = init.DslPtr();
  debug_ = debugIn;
  // Get Keywords
  image_.InitImaging( !(actionArgs.hasKey("noimage")) );
  double dist = actionArgs.getKeyDouble("distance", 7.0);
  byResidue_ = actionArgs.hasKey("byresidue");
  resoffset_ = actionArgs.getKeyInt("resoffset", 0) + 1;
  if (resoffset_ < 1) {
    mprinterr("Error: Residue offset must be >= 0\n");
    return Action::ERR;
  }
  includeSolvent_ = actionArgs.hasKey("includesolvent");
  series_ = actionArgs.hasKey("series");
  saveNonNative_ = actionArgs.hasKey("savenonnative");
# ifdef MPI
  if (saveNonNative_) {
    mprinterr("Error: Saving non-native contact data not yet supported for MPI\n");
    return Action::ERR;
  }
# endif
  distance_ = dist * dist; // Square the cutoff
  first_ = actionArgs.hasKey("first");
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  if (series_) {
    seriesout_ = init.DFL().AddDataFile(actionArgs.GetStringKey("seriesout"), actionArgs);
    init.DSL().SetDataSetsPending( true );
    if (saveNonNative_)
      seriesNNout_ = init.DFL().AddDataFile(actionArgs.GetStringKey("seriesnnout"), actionArgs);
  }
  cfile_ = init.DFL().AddCpptrajFile(actionArgs.GetStringKey("writecontacts"), "Native Contacts",
                               DataFileList::TEXT, true);
  pfile_ = init.DFL().AddCpptrajFile(actionArgs.GetStringKey("contactpdb"), "Contact PDB",
                               DataFileList::PDB);
  if (saveNonNative_)
    nfile_ = init.DFL().AddCpptrajFile(actionArgs.GetStringKey("nncontactpdb"),
                                       "Non-native Contact PDB", DataFileList::PDB);
  rfile_ = init.DFL().AddCpptrajFile(actionArgs.GetStringKey("resout"), "Contact Res Pairs",
                               DataFileList::TEXT, true);
  if (cfile_ == 0 || rfile_ == 0) return Action::ERR;
  pdbcut_ = (float)actionArgs.getKeyDouble("pdbcut", -1.0);
  usepdbcut_ = (pdbcut_ > -1.0);
  // Get reference
  ReferenceFrame REF = init.DSL().GetReferenceFrame( actionArgs );
  if (!first_) {
    if (REF.error()) return Action::ERR;
    if (REF.empty()) {
      mprintf("Warning: No reference structure specified. Defaulting to first.\n");
      first_ = true;
    }
  } else {
    if (!REF.empty()) {
      mprinterr("Error: Must only specify 'first' or a reference structure, not both.\n");
      return Action::ERR;
    }
  }
  // Create data sets
  std::string name = actionArgs.GetStringKey("name");
  if (name.empty())
    name = init.DSL().GenerateDefaultName("Contacts");
  numnative_ = init.DSL().AddSet(DataSet::INTEGER, MetaData(name, "native"));
  nonnative_ = init.DSL().AddSet(DataSet::INTEGER, MetaData(name, "nonnative"));
  if (outfile != 0) {
    outfile->AddDataSet(numnative_);
    outfile->AddDataSet(nonnative_);
  }
  if (numnative_ == 0 || nonnative_ == 0) return Action::ERR;
  if (actionArgs.hasKey("mindist")) {
    mindist_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(name, "mindist"));
    if (mindist_ == 0) return Action::ERR;
    if (outfile != 0) outfile->AddDataSet(mindist_);
  }
  if (actionArgs.hasKey("maxdist")) {
    maxdist_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(name, "maxdist"));
    if (maxdist_ == 0) return Action::ERR;
    if (outfile != 0) outfile->AddDataSet(maxdist_);
  }
  DataFile *natmapfile = 0, *nonmapfile = 0;
  if (actionArgs.hasKey("map")) {
    nativeMap_ = (DataSet_MatrixDbl*)init.DSL().AddSet(DataSet::MATRIX_DBL, MetaData(name, "nativemap"));
    if (nativeMap_ == 0) return Action::ERR;
    nonnatMap_ = (DataSet_MatrixDbl*)init.DSL().AddSet(DataSet::MATRIX_DBL, MetaData(name, "nonnatmap"));
    if (nonnatMap_ == 0) return Action::ERR;
    FileName mapFilename;
    mapFilename.SetFileName( actionArgs.GetStringKey("mapout") );
    if (!mapFilename.empty()) {
      natmapfile = init.DFL().AddDataFile(mapFilename.PrependFileName("native."));
      if (natmapfile != 0) natmapfile->AddDataSet(nativeMap_);
      nonmapfile = init.DFL().AddDataFile(mapFilename.PrependFileName("nonnative."));
      if (nonmapfile != 0) nonmapfile->AddDataSet(nonnatMap_);
    }
  }
  // Get Masks
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );
  std::string mask2 = actionArgs.GetMaskNext();
  if (!mask2.empty())
    Mask2_.SetMaskString( mask2 );
  mprintf("    NATIVECONTACTS: Mask1='%s'", Mask1_.MaskString());
  if (Mask2_.MaskStringSet())
    mprintf(" Mask2='%s'", Mask2_.MaskString());
  mprintf(", contacts set up based on");
  if (first_)
    mprintf(" first frame.\n");
  else
    mprintf("'%s'.\n", REF.refName());
  if (byResidue_) {
    mprintf("\tContacts will be ignored for residues spaced < %i apart.\n", resoffset_);
    if (nativeMap_ != 0)
      mprintf("\tMap will be printed by residue.\n");
  }
  if (saveNonNative_)
    mprintf("\tSaving non-native contacts as well (may use a lot of memory).\n");
  mprintf("\tDistance cutoff is %g Angstroms,", sqrt(distance_));
  if (!image_.UseImage())
    mprintf(" imaging is off.\n");
  else
    mprintf(" imaging is on.\n");
  if (includeSolvent_)
    mprintf("\tMask selection will including solvent.\n");
  else
    mprintf("\tMask selection will not include solvent.\n");
  mprintf("\tData set name: %s\n", name.c_str());
  if (maxdist_ != 0)
    mprintf("\tSaving maximum observed distances in set '%s'\n", maxdist_->legend());
  if (mindist_ != 0)
    mprintf("\tSaving minimum observed distances in set '%s'\n", mindist_->legend());
  if (outfile != 0)
    mprintf("\tOutput to '%s'\n", outfile->DataFilename().full());
  mprintf("\tContact stats will be written to '%s'\n", cfile_->Filename().full());
  mprintf("\tContact res pairs will be written to '%s'\n", rfile_->Filename().full());
  if (pfile_ != 0) {
    mprintf("\tContact PDB will be written to '%s'\n", pfile_->Filename().full());
    if (usepdbcut_) mprintf("\tOnly atoms with values > %g will be written to PDB.\n", pdbcut_);
  }
  if (nfile_ != 0) {
    mprintf("\tNon-native contact PDB will be written to '%s'\n", nfile_->Filename().full());
    if (usepdbcut_) mprintf("\tOnly atoms with values > %g will be written to PDB.\n", pdbcut_);
  }
  if (nativeMap_ != 0) {
    mprintf("\tNative contacts map will be saved as set '%s'\n"
            "\tNon-native contacts map will be saved as set '%s'\n",
            nativeMap_->legend(), nonnatMap_->legend());
    if (natmapfile!=0) mprintf("\tNative map output to '%s'\n",natmapfile->DataFilename().full());
    if (nonmapfile!=0) mprintf("\tNative map output to '%s'\n",nonmapfile->DataFilename().full());
  }
  if (series_) {
    mprintf("\tSaving native contact time series %s[NC].\n", name.c_str());
    if (seriesout_ != 0) mprintf("\tWriting native contact time series to %s\n",
                                 seriesout_->DataFilename().full());
    if (saveNonNative_) {
      mprintf("\tSaving non-native contact time series %s[NN]\n", name.c_str());
      if (seriesNNout_ != 0) mprintf("\tWriting non-native contact time series to %s\n",
                                     seriesNNout_->DataFilename().full());
    }
  }
  // Set up reference if necessary.
  if (!first_) {
    // Set up imaging info for ref parm
    image_.SetupImaging( REF.CoordsInfo().TrajBox().Type() );
    if (image_.ImageType() == NONORTHO)
      REF.Coord().BoxCrd().ToRecip(ucell_, recip_);
    if (DetermineNativeContacts( REF.Parm(), REF.Coord() )) return Action::ERR;
  }
  return Action::OK;
}

// Action_NativeContacts::Setup()
Action::RetType Action_NativeContacts::Setup(ActionSetup& setup) {
  // Setup potential contact lists for this topology
  if (SetupContactLists( setup.Top(), Frame()))
    return Action::SKIP;
  mprintf("\t%zu potential contact sites for '%s'\n", Mask1_.Nselected(), Mask1_.MaskString());
  if (Mask2_.MaskStringSet())
    mprintf("\t%zu potential contact sites for '%s'\n", Mask2_.Nselected(), Mask2_.MaskString());
  // Set up imaging info for this parm
  image_.SetupImaging( setup.CoordInfo().TrajBox().Type() );
  if (image_.ImagingEnabled())
    mprintf("\tImaging enabled.\n");
  else
    mprintf("\tImaging disabled.\n");
  CurrentParm_ = setup.TopAddress();
  return Action::OK;
}

/// \return true if valid contact; when by residue atoms cannot be in same residue
bool Action_NativeContacts::ValidContact(int a1, int a2, Topology const& parmIn) const {
  if (byResidue_) {
    if ( abs(parmIn[a1].ResNum() - parmIn[a2].ResNum()) < resoffset_ )
      return false;
  }
  return true;
}

/** This macro is used by DoAction to check if a contact is valid, formed,
  * if it is native, and if so update it.
  */
#define UpdateNativeContact(M1_, M2_, CI1_, CI2_) { \
        if (ValidContact(M1_[c1], M2_[c2], *CurrentParm_)) { \
          double Dist2 = DIST2(frm.Frm().XYZ(M1_[c1]), frm.Frm().XYZ(M2_[c2]), \
                               image_.ImageType(), frm.Frm().BoxCrd(), ucell_, recip_); \
          minDist2 = std::min( Dist2, minDist2 ); \
          maxDist2 = std::max( Dist2, maxDist2 ); \
          if (Dist2 < distance_) { \
            Cpair contact_pair(M1_[c1], M2_[c2]); \
            contactListType::iterator it = nativeContacts_.find( contact_pair ); \
            if (it != nativeContacts_.end()) \
            { \
              ++Nnative; \
              it->second.Increment(frameNum, sqrt(Dist2), Dist2); \
              if (nativeMap_ != 0) nativeMap_->Element(CI1_[c1] - matrix_min_, \
                                                       CI2_[c2] - matrix_min_) += 1; \
            } else { \
              ++NnonNative; \
              if (nonnatMap_ != 0) nonnatMap_->Element(CI1_[c1] - matrix_min_, \
                                                       CI2_[c2] - matrix_min_) += 1; \
              if (saveNonNative_) { \
                it = nonNativeContacts_.find( contact_pair ); \
                if (it == nonNativeContacts_.end()) { \
                  int at1 = contact_pair.first; \
                  int at2 = contact_pair.second; \
                  std::string legend(CurrentParm_->AtomMaskName(at1) + "_" + \
                                     CurrentParm_->AtomMaskName(at2)); \
                  int r1 = (*CurrentParm_)[at1].ResNum(); \
                  int r2 = (*CurrentParm_)[at2].ResNum(); \
                  std::pair<contactListType::iterator, bool> ret; \
                  ret = nonNativeContacts_.insert( Mpair(contact_pair, \
                                                         contactType(legend,r1,r2)) ); \
                  it = ret.first; \
                  if (series_) { \
                    MetaData md(numnative_->Meta().Name(), "NN", nonNativeContacts_.size()); \
                    md.SetLegend( legend ); \
                    DataSet* ds = masterDSL_->AddSet(DataSet::INTEGER, md); \
                    it->second.SetData( ds ); \
                    if (seriesNNout_ != 0) seriesNNout_->AddDataSet( ds ); \
                  } \
                } \
                it->second.Increment(frameNum, sqrt(Dist2), Dist2); \
              } \
            } \
          } \
        } \
}

// Action_NativeContacts::DoAction()
Action::RetType Action_NativeContacts::DoAction(int frameNum, ActionFrame& frm) {
  if (image_.ImageType() == NONORTHO) frm.Frm().BoxCrd().ToRecip(ucell_, recip_);
  if (first_) {
    mprintf("\tUsing first frame to determine native contacts.\n");
    if (DetermineNativeContacts( *CurrentParm_, frm.Frm() )) return Action::ERR;
    first_ = false;
  }
  nframes_++;
  // Loop over all potential contacts
  int Nnative = 0;
  int NnonNative = 0;
  double maxDist2 = 0.0;
  double minDist2 = DBL_MAX;
  if ( Mask2_.MaskStringSet() ) {
    for (int c1 = 0; c1 != Mask1_.Nselected(); c1++)
      for (int c2 = 0; c2 != Mask2_.Nselected(); c2++)
      {
        UpdateNativeContact(Mask1_, Mask2_, contactIdx1_, contactIdx2_);
      }
  } else {
    for (int c1 = 0; c1 != Mask1_.Nselected(); c1++)
      for (int c2 = c1 + 1; c2 != Mask1_.Nselected(); c2++)
      {
        UpdateNativeContact(Mask1_, Mask1_, contactIdx1_, contactIdx1_);
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

void Action_NativeContacts::UpdateSeries() {
  if (seriesUpdated_) return;
  if (series_ && nframes_ > 0) {
    const int ZERO = 0;
    // Ensure all series have been updated for all frames.
    for (contactListType::iterator it = nativeContacts_.begin();
                                   it != nativeContacts_.end(); ++it)
      if (it->second.Data().Size() < nframes_)
        it->second.Data().Add( nframes_ - 1, &ZERO );
    for (contactListType::iterator it = nonNativeContacts_.begin();
                                   it != nonNativeContacts_.end(); ++it)
      if (it->second.Data().Size() < nframes_)
        it->second.Data().Add( nframes_ - 1, &ZERO );
  }
  // Should only be called once.
  seriesUpdated_ = true;
}

#ifdef MPI
int Action_NativeContacts::SyncAction() {
  // Make sure all time series are updated at this point.
  UpdateSeries();
  // Get total number of frames.
  int N = (int)nframes_;
  int total_frames;
  trajComm_.Reduce( &total_frames, &N, 1, MPI_INT, MPI_SUM );
  if (trajComm_.Master())
    nframes_ = (unsigned int)total_frames;
  // Should have the same number of contacts on each thread since reference is shared.
  for (contactListType::iterator it = nativeContacts_.begin(); it != nativeContacts_.end(); ++it)
  {
    double total[3], buf[3];
    buf[0] = it->second.Avg();
    buf[1] = it->second.Stdev();
    buf[2] = (double)it->second.Nframes();
    trajComm_.Reduce( total, buf, 3, MPI_DOUBLE, MPI_SUM );
    if (trajComm_.Master())
      it->second.SetValues( total[0], total[1], (int)total[2] );
  }
  return 0;
}
#endif

// Action_NativeContacts::WriteContacts()
void Action_NativeContacts::WriteContacts(contactListType& ContactsIn) {
  if (ContactsIn.empty()) return;
  // Map of residue pairs to total contact values.
  typedef std::map<Cpair, resContact> resContactMap;
  resContactMap ResContacts;
  std::pair<resContactMap::iterator, bool> ret;
  // Normalize native contacts. Place them into an array where they will
  // be sorted. Sum up total contact over residue pairs.
  std::vector<contactType> sortedList;
  for (contactListType::iterator it = ContactsIn.begin();
                                 it != ContactsIn.end(); ++it)
  {
    it->second.Finalize();
    sortedList.push_back( it->second );
    ret = ResContacts.insert( Rpair(Cpair(it->second.Res1(),it->second.Res2()),
                                    resContact(it->second.Nframes())) );
    if (!ret.second) // residue pair exists, update it.
      ret.first->second.Increment( it->second.Nframes() );
  }
  std::sort( sortedList.begin(), sortedList.end() );
  // Place residue pairs into an array to be sorted.
  std::vector<Rpair> ResList;
  for (resContactMap::const_iterator it = ResContacts.begin(); it != ResContacts.end(); ++it)
    ResList.push_back( *it );
  std::sort( ResList.begin(), ResList.end(), res_cmp() );
  // Print out total fraction frames for residue pairs.
  rfile_->Printf("%-8s %8s %10s %10s\n", "#Res1", "#Res2", "TotalFrac", "Contacts");
  //for (resContactMap::const_iterator it = ResContacts.begin(); it != ResContacts.end(); ++it)
  for (std::vector<Rpair>::const_iterator it = ResList.begin();
                                          it != ResList.end(); ++it)
    rfile_->Printf("%-8i %8i %10g %10i\n", it->first.first+1, it->first.second+1,
                  (double)it->second.Nframes()/(double)nframes_,
                  it->second.Ncontacts());
  // Print out sorted atom contacts.
  cfile_->Printf("%-8s %20s %8s %8s %8s %8s\n", "#", "Contact", "Nframes", "Frac.", "Avg", "Stdev");
  unsigned int num = 1;
  for (std::vector<contactType>::const_iterator NC = sortedList.begin();
                                                NC != sortedList.end(); ++NC, ++num)
  { 
    double fracPresent = (double)NC->Nframes() / (double)nframes_;
    cfile_->Printf("%8u %20s %8i %8.3g %8.3g %8.3g\n", num, NC->id(),
                   NC->Nframes(), fracPresent, NC->Avg(), NC->Stdev());
  }
}

// Action_NativeContacts::WriteContactPDB()
void Action_NativeContacts::WriteContactPDB( contactListType& ContactsIn, CpptrajFile* fileIn)
{
    std::vector<double> atomContactFrac(refParm_->Natom(), 0.0);
    double norm = 1.0 / ((double)nframes_ * 2.0);
    for (contactListType::const_iterator it = ContactsIn.begin();
                                         it != ContactsIn.end(); ++it)
    {
      int a1 = it->first.first;
      int a2 = it->first.second;
      contactType const& NC = it->second;
      double fracShared = (double)NC.Nframes() * norm;
      atomContactFrac[a1] += fracShared;
      atomContactFrac[a2] += fracShared;
    }
    // Normalize so the strongest contact value is 100.00
    norm = (double)*std::max_element(atomContactFrac.begin(), atomContactFrac.end());
    norm = 100.00 / norm;
    PDBfile& contactPDB = static_cast<PDBfile&>( *fileIn );
    mprintf("\tWriting contacts PDB to '%s'\n", fileIn->Filename().full());
    contactPDB.WriteTITLE( numnative_->Meta().Name() + " " + Mask1_.MaskExpression() + " " +
                           Mask2_.MaskExpression() );
    int cidx = 0;
    for (int aidx = 0; aidx != refParm_->Natom(); aidx++, cidx += 3) {
      float bfac = (float)(atomContactFrac[aidx] * norm);
      if (!usepdbcut_ || (bfac > pdbcut_)) {
        int resnum = (*refParm_)[aidx].ResNum();
        contactPDB.WriteCoord(PDBfile::ATOM, aidx+1, (*refParm_)[aidx].Name(),
                              refParm_->Res(resnum).Name(), resnum+1,
                              refFrame_[cidx], refFrame_[cidx+1], refFrame_[cidx+2],
                              1.0, bfac, (*refParm_)[aidx].ElementName(), 0);
      }
    }
}

// Action_NativeContacts::Print()
void Action_NativeContacts::Print() {
  if (nativeMap_ != 0) {
    // Normalize maps by number of frames.
    double norm = 1.0 / (double)nframes_;
    for (DataSet_MatrixDbl::iterator m = nativeMap_->begin(); m != nativeMap_->end(); ++m)
      *m *= norm;
    for (DataSet_MatrixDbl::iterator m = nonnatMap_->begin(); m != nonnatMap_->end(); ++m)
      *m *= norm;
  }
  // Ensure all series have been updated for all frames.
  UpdateSeries();
  if (!cfile_->IsStream()) {
    mprintf("    CONTACTS: %s: Writing native contacts to file '%s'\n",
            numnative_->Meta().Name().c_str(), cfile_->Filename().full());
    cfile_->Printf("# Contacts: %s\n", numnative_->Meta().Name().c_str());
    cfile_->Printf("# Native contacts determined from mask '%s'", Mask1_.MaskString());
    if (Mask2_.MaskStringSet())
      cfile_->Printf(" and mask '%s'", Mask2_.MaskString());
    cfile_->Printf("\n");
  } else
    mprintf("    CONTACTS: %s\n", numnative_->Meta().Name().c_str());
  WriteContacts( nativeContacts_ );
  if (saveNonNative_) {
    if (!cfile_->IsStream()) {
      mprintf("              %s: Writing non-native contacts to file '%s'\n",
              numnative_->Meta().Name().c_str(), cfile_->Filename().full());
      cfile_->Printf("# Non-native Contacts: %s\n", numnative_->Meta().Name().c_str());
    } else
      mprintf("      ------- Non-native %s -------\n", numnative_->Meta().Name().c_str());
    WriteContacts( nonNativeContacts_ );
  }
  // Break down contacts by atom, write to PDB.
  if (pfile_ != 0)
    WriteContactPDB( nativeContacts_, pfile_ );
  if (nfile_ != 0)
    WriteContactPDB( nonNativeContacts_, nfile_ );
}
// -----------------------------------------------------------------------------
void Action_NativeContacts::contactType::Finalize() {
  if (nframes_ > 0) {
    dist_ /= (double)nframes_;
    dist2_ /= (double)nframes_;
    dist2_ -= (dist_ * dist_);
    if (dist2_ > 0)
      dist2_ = sqrt(dist2_);
    else
      dist2_ = 0.0;
  }
}
