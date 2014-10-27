#include <cmath>
#include <cstdio> // sscanf
#include "Action_NMRrst.h"
#include "DataSet_double.h"
#include "StringRoutines.h"
#include "CpptrajStdio.h"
#include "AtomMap.h"

// CONSTRUCTOR
Action_NMRrst::Action_NMRrst() :
   useMass_(false), findNOEs_(false), resOffset_(0), debug_(0) {} 

void Action_NMRrst::Help() {
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

// Action_NMRrst::Init()
Action::RetType Action_NMRrst::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  // Get Keywords
  Image_.InitImaging( !(actionArgs.hasKey("noimage")) );
  useMass_ = !(actionArgs.hasKey("geom"));
  findNOEs_ = actionArgs.hasKey("findnoes");
  resOffset_ = actionArgs.getKeyInt("resoffset", 0);
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  std::string rstfilename = actionArgs.GetStringKey("file");
  std::string setname = actionArgs.GetStringKey("name");
  if (setname.empty())
    setname = DSL->GenerateDefaultName("NMR");
  if (!findNOEs_ && rstfilename.empty()) {
    mprinterr("Error: Must specify restraint file and/or 'findnoes'.\n");
    return Action::ERR;
  }

  // Read in NMR restraints.
  if (!rstfilename.empty()) {
    if (ReadNmrRestraints( rstfilename )) return Action::ERR;
  }

  // Set up distances.
  int num_noe = 1;
  for (noeArray::iterator noe = NOEs_.begin(); noe != NOEs_.end(); ++noe, ++num_noe) {
     // Translate any ambiguous atom names
     TranslateAmbiguous( noe->aName1_ ); 
     TranslateAmbiguous( noe->aName2_ );
     // Create mask expressions from resnum/atom name
     noe->dMask1_.SetMaskString( MaskExpression( noe->resNum1_, noe->aName1_ ) ); 
     noe->dMask2_.SetMaskString( MaskExpression( noe->resNum2_, noe->aName2_ ) ); 
     // Dataset to store distances
    noe->dist_ = DSL->AddSetIdxAspect(DataSet::DOUBLE, setname, num_noe, "NOE");
    if (noe->dist_==0) return Action::ERR;
    noe->dist_->SetScalar( DataSet::M_DISTANCE, DataSet::NOE );
    ((DataSet_double*)noe->dist_)->SetNOE(noe->bound_, noe->boundh_, noe->rexp_);
    noe->dist_->SetLegend( noe->dMask1_.MaskExpression() + " and " + 
                           noe->dMask2_.MaskExpression());
    // Add dataset to data file
    if (outfile != 0) outfile->AddSet( noe->dist_ );
  }
 
  mprintf("    NMRRST: %zu NOEs from NMR restraint file.\n", NOEs_.size());
  mprintf("\tShifting residue numbers in restraint file by %i\n", resOffset_);
  // DEBUG - print NOEs
  for (noeArray::const_iterator noe = NOEs_.begin(); noe != NOEs_.end(); ++noe)
    mprintf("\t'%s'  %f < %f < %f\n", noe->dist_->Legend().c_str(),
            noe->bound_, noe->rexp_, noe->boundh_);
  if (findNOEs_)
    mprintf("\tSearching for potential NOEs.\n");
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
Action::RetType Action_NMRrst::Setup(Topology* currentParm, Topology** parmAddress) {
  // Set up NOEs from file.
  for (noeArray::iterator noe = NOEs_.begin(); noe != NOEs_.end(); ++noe) {
    if (currentParm->SetupIntegerMask( noe->dMask1_ )) return Action::ERR;
    if (currentParm->SetupIntegerMask( noe->dMask2_ )) return Action::ERR;
    //mprintf("\t%s (%i atoms) to %s (%i atoms)",Mask1_.MaskString(), Mask1_.Nselected(),
    //        Mask2_.MaskString(),Mask2_.Nselected());
    if (noe->dMask1_.None() || noe->dMask2_.None()) {
      mprintf("Warning: One or both masks for NOE '%s' have no atoms (%i and %i).\n",
              noe->dist_->Legend().c_str(), noe->dMask1_.Nselected(),
              noe->dMask2_.Nselected());
      noe->active_ = false; 
    } else
      noe->active_ = true;
  }
  // Set up potential NOE sites.
  if (findNOEs_) {
    potentialSites_.clear();
    AtomMap resMap;
    resMap.SetDebug( debug_ );
    std::vector<bool> selected;
    Range soluteRes = currentParm->SoluteResidues();
    for (Range::const_iterator res = soluteRes.begin(); res != soluteRes.end(); ++res)
    {
      int res_first_atom = currentParm->Res(*res).FirstAtom();
      selected.assign( currentParm->Res(*res).NumAtoms(), false );
      // Find symmetric atom groups.
      AtomMap::AtomIndexArray symmGroups;
      if (resMap.SymmetricAtoms(*currentParm, symmGroups, *res)) return Action::ERR;
      // DEBUG
      mprintf("DEBUG: Residue %i: symmetric atom groups:\n", *res + 1);
      for (AtomMap::AtomIndexArray::const_iterator grp = symmGroups.begin();
                                                   grp != symmGroups.end(); ++grp)
      {
        mprintf("\t\t");
        for (AtomMap::Iarray::const_iterator at = grp->begin();
                                             at != grp->end(); ++at)
          mprintf(" %s", currentParm->TruncAtomNameNum( *at ).c_str());
        mprintf("\n");
      }
      // Each symmetric hydrogen atom group is a site.
      for (AtomMap::AtomIndexArray::const_iterator grp = symmGroups.begin();
                                                   grp != symmGroups.end(); ++grp)
      { // NOTE: If first atom is H all should be H.
        if ( (*currentParm)[ grp->front() ].Element() == Atom::HYDROGEN )
        {
          potentialSites_.push_back( Site(*res, *grp) );
          // Mark symmetric atoms as selected.
          for (AtomMap::Iarray::const_iterator at = grp->begin();
                                               at != grp->end(); ++at)
            selected[ *at - res_first_atom ] = true;
        }
      }
      // All other non-selected hydrogens bonded to same heavy atom are sites.
      for (int ratom = res_first_atom; ratom != currentParm->Res(*res).LastAtom(); ++ratom)
      {
        if ( (*currentParm)[ratom].Element() != Atom::HYDROGEN ) {
          Iarray heavyAtomGroup;
          for (Atom::bond_iterator ba = (*currentParm)[ratom].bondbegin();
                                   ba != (*currentParm)[ratom].bondend(); ++ba)
            if ( *ba >= res_first_atom &&
                 !selected[ *ba - res_first_atom ] &&
                 (*currentParm)[ *ba ].Element() == Atom::HYDROGEN )
                heavyAtomGroup.push_back( *ba );
          if (!heavyAtomGroup.empty())
            potentialSites_.push_back( Site(*res, heavyAtomGroup) );
        }
      }
    }
    mprintf("DEBUG: Potential NOE sites:\n");
    for (SiteArray::const_iterator site = potentialSites_.begin();
                                   site != potentialSites_.end(); ++site)
    {
      mprintf("\tRes %i:", site->ResNum()+1);
      for (Site::Idx_it it = site->IdxBegin(); it != site->IdxEnd(); ++it)
        mprintf(" %s", currentParm->TruncAtomNameNum( *it ).c_str());
      mprintf("\n");
    }
  }
  // Set up imaging info for this parm
  Image_.SetupImaging( currentParm->BoxType() );
  if (Image_.ImagingEnabled())
    mprintf("\tImaged.\n");
  else
    mprintf("\tImaging off.\n");

  return Action::OK;  
}

// Action_NMRrst::DoAction()
Action::RetType Action_NMRrst::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  double Dist;
  Matrix_3x3 ucell, recip;
  Vec3 a1, a2;

  for (noeArray::iterator noe = NOEs_.begin(); noe != NOEs_.end(); ++noe) {
    if ( (*noe).active_ ) {
      if (useMass_) {
        a1 = currentFrame->VCenterOfMass( (*noe).dMask1_ );
        a2 = currentFrame->VCenterOfMass( (*noe).dMask2_ );
      } else {
        a1 = currentFrame->VGeometricCenter( (*noe).dMask1_ );
        a2 = currentFrame->VGeometricCenter( (*noe).dMask2_ );
      }

      switch ( Image_.ImageType() ) {
        case NONORTHO:
          currentFrame->BoxCrd().ToRecip(ucell, recip);
          Dist = DIST2_ImageNonOrtho(a1, a2, ucell, recip);
          break;
        case ORTHO:
          Dist = DIST2_ImageOrtho(a1, a2, currentFrame->BoxCrd());
          break;
        case NOIMAGE:
          Dist = DIST2_NoImage(a1, a2);
          break;
      }
      Dist = sqrt(Dist);
      (*noe).dist_->Add(frameNum, &Dist);
    }
  }
  return Action::OK;
} 
