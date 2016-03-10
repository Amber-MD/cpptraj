#include <cmath>
#include <cstdio> // sscanf
#include <algorithm> // sort
#include "Action_NMRrst.h"
#include "DataSet_double.h"
#include "DataSet_float.h"
#include "StringRoutines.h"
#include "CpptrajStdio.h"
#include "AtomMap.h"

// CONSTRUCTOR
Action_NMRrst::Action_NMRrst() : Action(HIDDEN),
   findOutput_(0), specOutput_(0), masterDSL_(0), numNoePairs_(0), max_cut_(6.0),
   strong_cut_(2.9), medium_cut_(3.5), weak_cut_(5.0),
   resOffset_(0), debug_(0), nframes_(0), useMass_(false),
   findNOEs_(false), series_(false) {} 

void Action_NMRrst::Help() const {
  mprintf("\t[<name>] file <rstfile> [name <dataname>] [geom] [noimage] [resoffset <r>]\n"
          "  Calculate distances based on entries in the given NMR restraint file.\n");
}

/// \return true if first character is a 'skippable' one.
static inline bool SkipChar( const char* ptr ) {
  return (ptr != 0 && (*ptr == '#' || *ptr == '!' || *ptr == '\n' || *ptr == '\r'));
}

static inline std::string MaskExpression(int resnum, std::string& aname) {
  return std::string( ":" + integerToString(resnum) + "@" + aname);
}

static void TranslateAmbiguous(std::string& aname) {
  if (aname == "QA") // Gly a-methylene
    aname.assign("HA=");
}
// -----------------------------------------------------------------------------
// Action_NMRrst::Init()
Action::RetType Action_NMRrst::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  if (init.TrajComm().Size() > 1) {
    mprinterr("Error: 'nmrrst' action does not work with > 1 thread (%i threads currently).\n",
              init.TrajComm().Size());
    return Action::ERR;
  }
# endif
  debug_ = debugIn;
  // Get Keywords
  Image_.InitImaging( !(actionArgs.hasKey("noimage")) );
  useMass_ = !(actionArgs.hasKey("geom"));
  findNOEs_ = actionArgs.hasKey("findnoes");
  findOutput_ = init.DFL().AddCpptrajFile(actionArgs.GetStringKey("findout"), "Found NOEs",
                                    DataFileList::TEXT, true);
  specOutput_ = init.DFL().AddCpptrajFile(actionArgs.GetStringKey("specout"), "Specified NOEs",
                                    DataFileList::TEXT, true);
  if (findOutput_ == 0 || specOutput_ == 0) return Action::ERR;
  resOffset_ = actionArgs.getKeyInt("resoffset", 0);
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  max_cut_ = actionArgs.getKeyDouble("cut", 6.0);
  strong_cut_ = actionArgs.getKeyDouble("strongcut", 2.9);
  medium_cut_ = actionArgs.getKeyDouble("mediumcut", 3.5);
  weak_cut_ = actionArgs.getKeyDouble("weakcut", 5.0);
  series_ = actionArgs.hasKey("series");
  std::string rstfilename = actionArgs.GetStringKey("file");
  setname_ = actionArgs.GetStringKey("name");
  if (setname_.empty())
    setname_ = init.DSL().GenerateDefaultName("NMR");
  nframes_ = 0;

  // Atom Mask
  Mask_.SetMaskString( actionArgs.GetMaskNext() );

  // Pairs specified on command line.
  std::string pair1 = actionArgs.GetStringKey("pair");
  while (!pair1.empty()) {
    std::string pair2 = actionArgs.GetStringNext();
    if (pair2.empty()) {
      mprinterr("Error: Only one mask specified for pair (%s)\n", pair1.c_str());
      return Action::ERR;
    }
    Pairs_.push_back( MaskPairType(AtomMask(pair1), AtomMask(pair2)) );
    pair1 = actionArgs.GetStringKey("pair");
  }

  // Check that something will be done
  if (!findNOEs_ && rstfilename.empty() && Pairs_.empty()) {
    mprinterr("Error: Must specify restraint file, 'pair', and/or 'findnoes'.\n");
    return Action::ERR;
  }

  // Read in NMR restraints.
  if (!rstfilename.empty()) {
    if (ReadNmrRestraints( rstfilename )) return Action::ERR;
  }

  // Set up distances.
  int num_noe = 1;
  for (noeDataArray::iterator noe = NOEs_.begin(); noe != NOEs_.end(); ++noe, ++num_noe) {
    // Translate any ambiguous atom names
    TranslateAmbiguous( noe->aName1_ ); 
    TranslateAmbiguous( noe->aName2_ );
    // Create mask expressions from resnum/atom name
    noe->dMask1_.SetMaskString( MaskExpression( noe->resNum1_, noe->aName1_ ) ); 
    noe->dMask2_.SetMaskString( MaskExpression( noe->resNum2_, noe->aName2_ ) ); 
    // Dataset to store distances
    AssociatedData_NOE noeData(noe->bound_, noe->boundh_, noe->rexp_);
    MetaData md(setname_, "NOE", num_noe);
    md.SetLegend( noe->dMask1_.MaskExpression() + " and " + noe->dMask2_.MaskExpression());
    md.SetScalarMode( MetaData::M_DISTANCE );
    md.SetScalarType( MetaData::NOE );
    noe->dist_ = init.DSL().AddSet(DataSet::DOUBLE, md);
    if (noe->dist_==0) return Action::ERR;
    noe->dist_->AssociateData( &noeData );
    // Add dataset to data file
    if (outfile != 0) outfile->AddDataSet( noe->dist_ );
  }

  masterDSL_ = init.DslPtr();
 
  mprintf("Warning: *** THIS ACTION IS EXPERIMENTAL. ***\n");
  mprintf("    NMRRST: %zu NOEs from NMR restraint file.\n", NOEs_.size());
  mprintf("\tShifting residue numbers in restraint file by %i\n", resOffset_);
  // DEBUG - print NOEs
  for (noeDataArray::const_iterator noe = NOEs_.begin(); noe != NOEs_.end(); ++noe)
    mprintf("\t'%s'  %f < %f < %f\n", noe->dist_->legend(),
            noe->bound_, noe->rexp_, noe->boundh_);
  if (findNOEs_) {
    mprintf("\tSearching for potential NOEs. Max cutoff is %g Ang.\n", max_cut_);
    mprintf("\tNOE distance criteria (Ang.): S= %g, M= %g, W= %g\n",
            strong_cut_, medium_cut_, weak_cut_);
    if (series_)
      mprintf("\tDistance data for NOEs less than cutoff will be saved as '%s[foundNOE]'.\n",
              setname_.c_str());
    mprintf("\tFound NOEs will be written to '%s'\n", findOutput_->Filename().full());
  }
  if (!Pairs_.empty()) {
    mprintf("\tSpecified NOE pairs:\n");
    for (MaskPairArray::const_iterator mp = Pairs_.begin(); mp != Pairs_.end(); ++mp)
      mprintf("\t\t[%s] to [%s]\n", mp->first.MaskString(), mp->second.MaskString());
    mprintf("\tSpecified NOE data will be written to '%s'\n", specOutput_->Filename().full());
  }
  if (!Image_.UseImage()) 
    mprintf("\tNon-imaged");
  else
    mprintf("\tImaged");
  if (useMass_) 
    mprintf(", center of mass.\n");
  else
    mprintf(", geometric center.\n");

  return Action::OK;
}
// -----------------------------------------------------------------------------
int Action_NMRrst::ReadNmrRestraints( std::string const& rstfilename )
{
  BufferedLine infile;
  if (infile.OpenFileRead( rstfilename )) return 1;
  // Try to determine what kind of file.
  const char* ptr = infile.Line();
  // Try to skip past any blank lines and comments
  while ( SkipChar( ptr ) )
    ptr = infile.Line();
  if (ptr == 0) {
    mprinterr("Error: Unexpected end of restraint file.\n");
    return 1; 
  }
  std::string inputLine( ptr );
  infile.CloseFile();
  // Re-open file
  if (infile.OpenFileRead( rstfilename )) return 1;
  int err = 0;
  if ( inputLine.compare(0, 7, "*HEADER")==0 ||
       inputLine.compare(0, 6, "*TITLE")==0 ||
       inputLine.compare(0, 6, "assign")==0 )
    // XPLOR
    err = ReadXplor( infile );
  else
    // Assume DIANA/Amber
    err = ReadAmber( infile );
  infile.CloseFile();
  if (err != 0) {
    mprinterr("Error: Could not parse restraint file.\n");
    return 1;
  }
  return 0;
}

// -----------------------------------------------------------------------------
/** Read Amber/DIANA style restraints. */
int Action_NMRrst::ReadAmber( BufferedLine& infile ) {
  noeDataType NOE;
  const char* ptr = infile.Line();
  if (ptr == 0) {
    mprinterr("Error: Unexpected end of Amber restraint file.\n");
    return 1;
  }
  char rname1[16], rname2[16], aname1[16], aname2[16];
  double l_bound, u_bound;
  while (ptr != 0) {
    if (!SkipChar(ptr)) {
      int cols = sscanf( ptr, "%d %s %s %d %s %s %lf %lf",
                         &NOE.resNum1_, rname1, aname1,
                         &NOE.resNum2_, rname2, aname2,
                         &l_bound, &u_bound );
      if (cols == 7) { // 7-column, upper-bound only
        NOE.bound_ = 0.0;
        NOE.boundh_ = l_bound;
      } else if (cols == 8) { // 8-column, lower/upper bounds
        NOE.bound_ = l_bound;
        NOE.boundh_ = u_bound;
      } else {
        mprinterr("Error: Expected only 7 or 8 columns in Amber restraint file, got %i.\n", cols);
        return 1;
      }
      NOE.rexp_ = -1.0;
      NOE.dist_ = 0;
      NOE.resNum1_ += resOffset_;
      NOE.resNum2_ += resOffset_;
      if (NOE.resNum1_ < 1 || NOE.resNum2_ < 1) {
        mprinterr("Error: One or both residue numbers are out of bounds (%i, %i)\n"
                  "Error: Line: %s", NOE.resNum1_, NOE.resNum2_, ptr);
      } else {
        NOE.aName1_.assign( aname1 );
        NOE.aName2_.assign( aname2 );
        NOEs_.push_back( NOE );
      }
    }
    ptr = infile.Line();
  }
  return 0;
}
// -----------------------------------------------------------------------------
/// Convert next Xplor-style selection 'resid X name A' resnum/atom name 
static inline int GetAssignSelection(std::string& aName, ArgList& line, int offset)
{
  int resnum = line.getKeyInt("resid",0) + offset;
  if (resnum < 1) return -1;
  aName = line.GetStringKey("name");
  return resnum;
}

/** Read Xplor-style restraint file. */
int Action_NMRrst::ReadXplor( BufferedLine& infile ) {
  noeDataType NOE;
  const char* ptr = infile.Line();
  if (ptr == 0) {
    mprinterr("Error: Unexpected end of XPLOR restraint file.\n");
    return 1;
  }
  while ( ptr != 0 ) {
    if (ptr[0] == 'a' && ptr[1] == 's' && ptr[2] == 's' &&
        ptr[3] == 'i' && ptr[4] == 'g' && ptr[5] == 'n'   )
    {
      // 'assign' statement
      ArgList line(ptr, " ()");
      if (line.empty()) {
        mprinterr("Error: Could not parse XPLOR 'assign' line:\n\t%s",ptr);
      } else {
        line.MarkArg(0); // Mark 'assign'
        // Get 2 Masks
        NOE.resNum1_ = GetAssignSelection( NOE.aName1_, line, resOffset_ );
        NOE.resNum2_ = GetAssignSelection( NOE.aName2_, line, resOffset_ );
        if (NOE.resNum1_ < 1 || NOE.resNum2_ < 1) {
          mprinterr("Error: Could not get masks from line:\n\t%s", ptr);
          mprinterr("Error: Check if residue number + offset is out of bounds.\n");
        } else {
          // Check for noe bounds
          NOE.rexp_ = line.getNextDouble(-1.0);
          if ( NOE.rexp_ < 0.0 ) {
            // No more on this line, assume jcoupling
            ptr = infile.Line();
            line.SetList(ptr, " ()");
            // Get 2 more masks and jcoupling values 
          } else {
            // NOE
            NOE.boundh_ = NOE.rexp_ + line.getNextDouble(0.0);
            NOE.bound_ = NOE.rexp_ - line.getNextDouble(0.0);
            NOE.dist_ = 0;
            NOEs_.push_back( NOE );
          }
        }
      }
    }
    ptr = infile.Line();
  }
  return 0;
}
// -----------------------------------------------------------------------------
// Action_NMRrst::Setup()
/** Determine what atoms each mask pertains to for the current parm file.
  */
Action::RetType Action_NMRrst::Setup(ActionSetup& setup) {
  // ---------------------------------------------
  // Set up NOEs from file.
  for (noeDataArray::iterator noe = NOEs_.begin(); noe != NOEs_.end(); ++noe) {
    if (setup.Top().SetupIntegerMask( noe->dMask1_ )) return Action::ERR;
    if (setup.Top().SetupIntegerMask( noe->dMask2_ )) return Action::ERR;
    if (noe->dMask1_.None() || noe->dMask2_.None()) {
      mprintf("Warning: One or both masks for NOE '%s' have no atoms (%i and %i).\n",
              noe->dist_->legend(), noe->dMask1_.Nselected(),
              noe->dMask2_.Nselected());
      noe->active_ = false; 
    } else
      noe->active_ = true;
  }
  // ---------------------------------------------
  // Set up potential NOE sites.
  if (findNOEs_) {
    if (setup.Top().SetupCharMask( Mask_ )) return Action::ERR;
    Mask_.MaskInfo();
    if (Mask_.None()) return Action::SKIP;
    SiteArray potentialSites; // .clear();
    AtomMap resMap;
    resMap.SetDebug( debug_ );
    std::vector<bool> selected;
    Range soluteRes = setup.Top().SoluteResidues();
    for (Range::const_iterator res = soluteRes.begin(); res != soluteRes.end(); ++res)
    {
      int res_first_atom = setup.Top().Res(*res).FirstAtom();
      selected.assign( setup.Top().Res(*res).NumAtoms(), false );
      // Find symmetric atom groups.
      AtomMap::AtomIndexArray symmGroups;
      if (resMap.SymmetricAtoms(setup.Top(), symmGroups, *res)) return Action::ERR;
      // DEBUG
      if (debug_ > 0) {
        mprintf("DEBUG: Residue %i: symmetric atom groups:\n", *res + 1);
        for (AtomMap::AtomIndexArray::const_iterator grp = symmGroups.begin();
                                                     grp != symmGroups.end(); ++grp)
        {
          mprintf("\t\t");
          for (AtomMap::Iarray::const_iterator at = grp->begin();
                                               at != grp->end(); ++at)
            mprintf(" %s", setup.Top().TruncAtomNameNum( *at ).c_str());
          mprintf("\n");
        }
      }
      // Each symmetric hydrogen atom group is a site.
      for (AtomMap::AtomIndexArray::const_iterator grp = symmGroups.begin();
                                                   grp != symmGroups.end(); ++grp)
      { // NOTE: If first atom is H all should be H.
        if ( setup.Top()[ grp->front() ].Element() == Atom::HYDROGEN )
        {
          Iarray symmAtomGroup;
          for (Iarray::const_iterator at = grp->begin();
                                      at != grp->end(); ++at)
            if (Mask_.AtomInCharMask( *at ))
              symmAtomGroup.push_back( *at );
          if (!symmAtomGroup.empty()) {
            potentialSites.push_back( Site(*res, symmAtomGroup) );
            // Mark symmetric atoms as selected.
            for (AtomMap::Iarray::const_iterator at = grp->begin();
                                                 at != grp->end(); ++at)
              selected[ *at - res_first_atom ] = true;
          }
        }
      }
      // All other non-selected hydrogens bonded to same heavy atom are sites.
      for (int ratom = res_first_atom; ratom != setup.Top().Res(*res).LastAtom(); ++ratom)
      {
        if ( setup.Top()[ratom].Element() != Atom::HYDROGEN ) {
          Iarray heavyAtomGroup;
          for (Atom::bond_iterator ba = setup.Top()[ratom].bondbegin();
                                   ba != setup.Top()[ratom].bondend(); ++ba)
            if ( Mask_.AtomInCharMask(*ba) && 
                 *ba >= res_first_atom && 
                 *ba < setup.Top().Res(*res).LastAtom() )
            {
              if ( !selected[ *ba - res_first_atom ] &&
                   setup.Top()[ *ba ].Element() == Atom::HYDROGEN )
                heavyAtomGroup.push_back( *ba );
            }
          if (!heavyAtomGroup.empty())
            potentialSites.push_back( Site(*res, heavyAtomGroup) );
        }
      }
    }
    mprintf("\t%zu potential NOE sites:\n", potentialSites.size());
    for (SiteArray::const_iterator site = potentialSites.begin();
                                   site != potentialSites.end(); ++site)
    {
      mprintf("  %u\tRes %i:", site - potentialSites.begin(), site->ResNum()+1);
      for (unsigned int idx = 0; idx != site->Nindices(); ++idx)
        mprintf(" %s", setup.Top().TruncAtomNameNum( site->Idx(idx) ).c_str());
      mprintf("\n");
    }
    if (noeArray_.empty()) {
      size_t siteArraySize = 0;
      // Set up all potential NOE pairs. Keep track of size.
      for (SiteArray::const_iterator site1 = potentialSites.begin();
                                     site1 != potentialSites.end(); ++site1)
      {
        for (SiteArray::const_iterator site2 = site1 + 1;
                                       site2 != potentialSites.end(); ++site2)
        {
          if (site1->ResNum() != site2->ResNum()) {
            std::string legend = site1->SiteLegend(setup.Top()) + "--" +
                                 site2->SiteLegend(setup.Top());
            DataSet* ds = 0;
            if (series_) {
              ds = masterDSL_->AddSet(DataSet::FLOAT,
                                      MetaData(setname_, "foundNOE", noeArray_.size()));
              if (ds == 0) return Action::ERR;
              // Construct a data set name.
              ds->SetLegend(legend);
            }
            noeArray_.push_back( NOEtype(*site1, *site2, ds, legend) );
            siteArraySize += (2 * sizeof(int) * site1->Nindices()) +
                             (2 * sizeof(int) * site2->Nindices());
          }
        }
      }
      numNoePairs_ = noeArray_.size();
      size_t siteSize = sizeof(int) + (2 * sizeof(Iarray)) + sizeof(Site);
      size_t noeSize = (2 * siteSize) + sizeof(DataSet*) + sizeof(double) 
                       + sizeof(NOEtype);
      if (series_) noeSize += sizeof(std::vector<float>);
      size_t noeArraySize = (noeSize * numNoePairs_) + siteArraySize;
      if (series_)
        noeArraySize += (setup.Nframes() * numNoePairs_ * sizeof(float));
      mprintf("\t%zu potential NOE pairs. Estimated memory usage is %s\n",
              numNoePairs_, ByteString(noeArraySize, BYTE_DECIMAL).c_str());
    } else if (numNoePairs_ != potentialSites.size()) {
      mprinterr("Warning: Found NOE matrix has already been set up for %zu potential\n"
                "Warning:   NOEs, but %zu NOEs currently found.\n", numNoePairs_,
                potentialSites.size());
      return Action::SKIP;
    }
  }
  // ---------------------------------------------
  // Set up NOEs specified on the command line
  if (!Pairs_.empty()) {
    if (!specifiedNOEs_.empty()) {
      mprintf("Warning: Specifying NOEs currently only works with first topology used.\n");
      return Action::SKIP;
    }
    for (MaskPairArray::iterator mp = Pairs_.begin(); mp != Pairs_.end(); mp++) {
      if (setup.Top().SetupIntegerMask( mp->first )) return Action::ERR;
      int res1 = CheckSameResidue(setup.Top(), mp->first);
      if (res1 < 0) continue;
      if (setup.Top().SetupIntegerMask( mp->second )) return Action::ERR;
      int res2 = CheckSameResidue(setup.Top(), mp->second);
      if (res2 < 0) continue;
      Site site1( res1, mp->first.Selected() );
      Site site2( res2, mp->second.Selected() );
      std::string legend = site1.SiteLegend(setup.Top()) + "--" +
                           site2.SiteLegend(setup.Top());
      DataSet* ds = 0;
      if (series_) {
        ds = masterDSL_->AddSet(DataSet::FLOAT,
                                MetaData(setname_, "specNOE", specifiedNOEs_.size()));
        if (ds == 0) return Action::ERR;
        ds->SetLegend(legend);
      }
      specifiedNOEs_.push_back( NOEtype(site1, site2, ds, legend) );
    }
  } 
  // Set up imaging info for this parm
  Image_.SetupImaging( setup.CoordInfo().TrajBox().Type() );
  if (Image_.ImagingEnabled())
    mprintf("\tImaged.\n");
  else
    mprintf("\tImaging off.\n");

  return Action::OK;  
}

/** Check that all atoms in mask belong to same residue. */
int Action_NMRrst::CheckSameResidue(Topology const& top, AtomMask const& mask) const {
  if (mask.None()) return -1;
  int resnum = top[mask[0]].ResNum();
  for (AtomMask::const_iterator at = mask.begin(); at != mask.end(); ++at) {
    int r = top[*at].ResNum();
    if (r != resnum) {
      mprintf("Warning: Mask atom %i %s not in same residue as %i %s\n",
              *at + 1, top.AtomMaskName(*at).c_str(),
              mask[0] + 1, top.AtomMaskName(mask[0]).c_str());
    }
  }
  return resnum;
}

/** Process an NOEtypeArray */
void Action_NMRrst::ProcessNoeArray(NOEtypeArray& Narray, Frame const& frameIn, int frameNum) const
{
  for (NOEtypeArray::iterator my_noe = Narray.begin();
                              my_noe != Narray.end(); ++my_noe)
  {
    unsigned int shortest_idx1 = 0, shortest_idx2 = 0;
    double shortest_dist2 = -1.0;
    for (unsigned int idx1 = 0; idx1 != my_noe->Site1().Nindices(); ++idx1)
    {
      for (unsigned int idx2 = 0; idx2 != my_noe->Site2().Nindices(); ++idx2)
      {
        double dist2 = DIST2(frameIn.XYZ(my_noe->Site1().Idx(idx1)),
                             frameIn.XYZ(my_noe->Site2().Idx(idx2)),
                             Image_.ImageType(), frameIn.BoxCrd(),
                             ucell_, recip_);
        if (shortest_dist2 < 0.0 || dist2 < shortest_dist2) {
          shortest_dist2 = dist2;
          shortest_idx1 = idx1;
          shortest_idx2 = idx2;
        }
      }
    }
    // NOTE: Saving d^2
    my_noe->UpdateNOE(frameNum, shortest_dist2, shortest_idx1, shortest_idx2);
  }
}

// Action_NMRrst::DoAction()
Action::RetType Action_NMRrst::DoAction(int frameNum, ActionFrame& frm) {
  double Dist;
  Vec3 a1, a2;

  if (Image_.ImageType() == NONORTHO)
    frm.Frm().BoxCrd().ToRecip(ucell_, recip_);
  // NOEs from file.
  for (noeDataArray::iterator noe = NOEs_.begin(); noe != NOEs_.end(); ++noe) {
    if ( noe->active_ ) {
      if (useMass_) {
        a1 = frm.Frm().VCenterOfMass( noe->dMask1_ );
        a2 = frm.Frm().VCenterOfMass( noe->dMask2_ );
      } else {
        a1 = frm.Frm().VGeometricCenter( noe->dMask1_ );
        a2 = frm.Frm().VGeometricCenter( noe->dMask2_ );
      }

      switch ( Image_.ImageType() ) {
        case NONORTHO: Dist = DIST2_ImageNonOrtho(a1, a2, ucell_, recip_); break;
        case ORTHO: Dist = DIST2_ImageOrtho(a1, a2, frm.Frm().BoxCrd()); break;
        case NOIMAGE: Dist = DIST2_NoImage(a1, a2); break;
      }
      Dist = sqrt(Dist);
      noe->dist_->Add(frameNum, &Dist);
    }
  }
  // Find NOEs
  ProcessNoeArray( noeArray_, frm.Frm(), frameNum );
  // Specified NOEs
  ProcessNoeArray( specifiedNOEs_, frm.Frm(), frameNum );

  ++nframes_;
  return Action::OK;
}

/** Analyze NOEtypeArray */
void Action_NMRrst::AnalyzeNoeArray(NOEtypeArray& Narray, CpptrajFile* out) const {
    // Calculate <r^-6>^(-1/6) 
    for (NOEtypeArray::iterator my_noe = Narray.begin();
                                my_noe != Narray.end(); ++my_noe)
    {
      double r6Avg = pow( my_noe->R6_Avg() / (double)nframes_, -1.0/6.0 );
      my_noe->SetR6Avg( r6Avg );
    }
    // Sort by R6 avg.
    std::sort( Narray.begin(), Narray.end() );
    // Remove sets greater than the cutoff.
    NOEtypeArray::iterator bad_noe = Narray.begin();
    while (bad_noe != Narray.end() && bad_noe->R6_Avg() < max_cut_)
      ++bad_noe;
    size_t newSize = bad_noe - Narray.begin();
    for (NOEtypeArray::iterator my_noe = bad_noe;
                                my_noe != Narray.end(); ++my_noe)
    {
      if (debug_ > 0)
        mprintf("\tRemoving: %s (%g Ang)\n", my_noe->PrintNOE().c_str(), my_noe->R6_Avg());
      if (my_noe->Data() != 0)
        masterDSL_->RemoveSet( my_noe->Data() );
    }
    Narray.resize( newSize );
    // Print Final found NOEs. Calculate distances from d^2 if necessary.
    std::vector<unsigned int> Bins(4, 0); // Strong, weak, medium, none
    double Cutoffs[3] = {strong_cut_, medium_cut_, weak_cut_};
    unsigned int current_cut = 0;
    const char* Labels[4] = {"STRONG", "MEDIUM", "WEAK", "NONE"};
    while (current_cut < 3 && Narray.front().R6_Avg() > Cutoffs[current_cut])
      ++current_cut;
    out->Printf("#Format: <r1>:{ @<a1X>(c1X) ... } -- <r2>:{ @<a2X>(c2X) ... }"
                   " <Avg Dist. (Ang)> <Label>\n"
                   "# r1, r2: Residue Numbers\n"
                   "# a1X, a2X: Group atom numbers\n"
                   "# c1X, c2X: Number of times atom was part of shortest distance.\n");
    out->Printf("#Final NOEs (%zu):\n#   %s\n", Narray.size(), Labels[current_cut]);
    for (NOEtypeArray::iterator my_noe = Narray.begin();
                                my_noe != Narray.end(); ++my_noe)
    {
      unsigned int old_cut = current_cut;
      while (current_cut < 3 && my_noe->R6_Avg() > Cutoffs[current_cut]) 
        current_cut++;
      if (old_cut != current_cut)
        out->Printf("#   %s\n", Labels[current_cut]);
      out->Printf("\t %s %g \"%s\"\n", my_noe->PrintNOE().c_str(),
                     my_noe->R6_Avg(), my_noe->legend());
      if (my_noe->Data() != 0) {
        DataSet_float& data = static_cast<DataSet_float&>( *(my_noe->Data()) );
        for (unsigned int i = 0; i != data.Size(); i++)
          data[i] = sqrt(data[i]);
      }
      // Bin
      if      (my_noe->R6_Avg() < strong_cut_) ++Bins[0];
      else if (my_noe->R6_Avg() < medium_cut_) ++Bins[1];
      else if (my_noe->R6_Avg() < weak_cut_)   ++Bins[2];
      else                                     ++Bins[3];
    }
    out->Printf("#Totals: %u strong, %u medium, %u weak, %u none.\n",
            Bins[0], Bins[1], Bins[2], Bins[3]);
}

// Action_NMRrst::Print()
void Action_NMRrst::Print() {
  if (findNOEs_ || !specifiedNOEs_.empty()) {
    mprintf("    NMRRST:\n");
    if (nframes_ < 1) {
      mprintf("Warning: No frames processed.\n");
      return;
    }
    if (findNOEs_)
      AnalyzeNoeArray(noeArray_, findOutput_);
    if (!specifiedNOEs_.empty()) {
      max_cut_  = 1000.00; // FIXME: Make sure specified NOEs do not get deleted.
      AnalyzeNoeArray(specifiedNOEs_, specOutput_);
    }
  }
}

// -----------------------------------------------------------------------------
std::string Action_NMRrst::Site::SiteLegend(Topology const& top) const {
  std::string legend( top.TruncResNameNum(resNum_) + "(" );
  for (Iarray::const_iterator atom = indices_.begin(); atom != indices_.end(); ++atom) {
    if (atom != indices_.begin()) legend.append(",");
    legend.append( top[*atom].Name().Truncated() );
  }
  legend.append(")");
  return legend;
}

std::string Action_NMRrst::NOEtype::PrintNOE() const {
  std::string noe( integerToString(Site1().ResNum()+1) + ":{" );
//  if (Site1().Nindices() > 1) {
    for (unsigned int idx = 0; idx != Site1().Nindices(); idx++)
      noe.append(" @" + integerToString(Site1().Idx(idx)+1) +
                  "(" + integerToString(Site1().Count(idx)) + ")");
//  } else
//    noe.append(" @" + integerToString(Site1().Idx(0)+1));
  noe.append(" } -- " + integerToString(Site2().ResNum()+1) + ":{");
//  if (Site2().Nindices() > 1) {
    for (unsigned int idx = 0; idx != Site2().Nindices(); idx++)
      noe.append(" @" + integerToString(Site2().Idx(idx)+1) +
                  "(" + integerToString(Site2().Count(idx)) + ")");
//  } else
//    noe.append(" @" + integerToString(Site2().Idx(0)+1));
  noe.append(" }");
  return noe;
}
