#include <cmath> // sqrt
#include <algorithm> // std::fill, std::min, std::max
#include <set> // SET
#include "Action_DSSP.h"
#include "CpptrajStdio.h"
#include "DataSet.h"
#include "DistRoutines.h"
#ifdef DSSPDEBUG
#include "Constants.h"
#endif
#ifdef _OPENMP
# include <omp.h>
#endif

/** From the original Kabsch & Sander 1983 paper. Obtained via
  *   q1 = 0.42e
  *   q2 = 0.20e
  *   f = 332 (approximate conversion factor to get energies in kcal/mol)
  *   fac = q1*q2*f
  */
const double Action_DSSP::DSSP_fac_ = 27.888;

/** DSSP default cutoff for determining hbonds, from 1983 Kabsch & Sander paper. */
const double Action_DSSP::DSSP_cut_ = -0.5;

/** DSSP 1 character SS assignment. */
const char  Action_DSSP::DSSP_char_[] = { ' ', 'E', 'B', 'G', 'H', 'I', 'T', 'S' };

/** Used for output to STRING data set. */
const char* Action_DSSP::SSchar_[]    = { "0", "E", "B", "G", "H", "I", "T", "S" };

/** Full SS names. */
const char* Action_DSSP::DSSP_name_[]={"None", "Extended", "Bridge", "3-10", "Alpha", "Pi", "Turn", "Bend"};


// ----- SSres -----------------------------------------------------------------
/// CONSTRUCTOR
Action_DSSP::SSres::SSres() :
  resDataSet_(0),
  chirality_(0),
  bend_(0),
  sstype_(NONE),
  num_(-1),
  C_(-1),
  O_(-1),
  N_(-1),
  H_(-1),
  CA_(-1),
  prevIdx_(-1),
  nextIdx_(-1),
  resChar_(' '),
  isSelected_(false)
{
  std::fill(SScount_, SScount_ + NSSTYPE_, 0);
  std::fill(Bcount_, Bcount_ + NBRIDGETYPE_, 0);
  std::fill(turnChar_, turnChar_ + NTURNTYPE_, ' ');
}

/// COPY CONSTRUCTOR
Action_DSSP::SSres::SSres(SSres const& rhs) :
  resDataSet_(rhs.resDataSet_),
  chirality_(rhs.chirality_),
  bend_(rhs.bend_),
  sstype_(rhs.sstype_),
  num_(rhs.num_),
  C_(rhs.C_),
  O_(rhs.O_),
  N_(rhs.N_),
  H_(rhs.H_),
  CA_(rhs.CA_),
  prevIdx_(rhs.prevIdx_),
  nextIdx_(rhs.nextIdx_),
  bridges_(rhs.bridges_),
  resChar_(rhs.resChar_),
  isSelected_(rhs.isSelected_)
{
  std::copy(rhs.SScount_, rhs.SScount_ + NSSTYPE_, SScount_);
  std::copy(rhs.Bcount_, rhs.Bcount_ + NBRIDGETYPE_, Bcount_);
  std::copy(rhs.turnChar_, rhs.turnChar_ + NTURNTYPE_, turnChar_);
}

/// ASSIGNMENT
Action_DSSP::SSres& Action_DSSP::SSres::operator=(SSres const& rhs) {
  if (this == &rhs) return *this;
  resDataSet_ = rhs.resDataSet_;
  chirality_ = rhs.chirality_;
  bend_ = rhs.bend_;
  sstype_ = rhs.sstype_;
  num_ = rhs.num_;
  C_ = rhs.C_;
  O_ = rhs.O_;
  N_ = rhs.N_;
  H_ = rhs.H_;
  CA_ = rhs.CA_;
  prevIdx_ = rhs.prevIdx_;
  nextIdx_ = rhs.nextIdx_;
  bridges_ = rhs.bridges_;
  resChar_ = rhs.resChar_;
  isSelected_ = rhs.isSelected_;
  std::copy(rhs.SScount_, rhs.SScount_ + NSSTYPE_, SScount_);
  std::copy(rhs.Bcount_, rhs.Bcount_ + NBRIDGETYPE_, Bcount_);
  std::copy(rhs.turnChar_, rhs.turnChar_ + NTURNTYPE_, turnChar_);
  return *this;
}

/** Reset atom indices and mark as not selected. */
void Action_DSSP::SSres::Deselect() {
  isSelected_ = false;
  C_ = -1;
  H_ = -1;
  N_ = -1;
  O_ = -1;
  CA_ = -1;
  prevIdx_ = -1;
  nextIdx_ = -1;
}

/** Reset assignment status. */
void Action_DSSP::SSres::Unassign() {
  sstype_ = NONE;
  bridges_.clear();
  std::fill(turnChar_, turnChar_ + NTURNTYPE_, ' ');
}

/** Accumulate SS data. */
void Action_DSSP::SSres::AccumulateData(int frameNum, bool useString, bool betaDetail)
{
  SScount_[sstype_]++;
  int idata = (int)sstype_;
  const char* sdata = SSchar_[sstype_];
  if (sstype_ == EXTENDED || sstype_ == BRIDGE) {
    // Record which types are present and in what amounts. Only record each
    // unique type once.
    int npara = 0;
    int nanti = 0;
    for (BridgeArray::const_iterator it = bridges_.begin(); it != bridges_.end(); ++it)
    {
      if (it->Btype() == PARALLEL)
        npara++;
      else if (it->Btype() == ANTIPARALLEL)
        nanti++;
    }
    if (npara > 0)
      Bcount_[PARALLEL]++;
    if (nanti > 0)
      Bcount_[ANTIPARALLEL]++;
    if (betaDetail) {
      // Determine parallel or antiparallel by the greater number of bridges,
      // tie goes to antiparallel. If both are zero, do nothing to match
      // previous behavior.
      if (npara > nanti) {
        idata = (int)EXTENDED;
        sdata = "b";
      } else if (nanti > 0) {
        idata = (int)BRIDGE;
        sdata = "B";
      }
    }
    // TODO bulge?
  } else
    Bcount_[NO_BRIDGE]++;
  if (useString)
    resDataSet_->Add(frameNum, sdata);
  else
    resDataSet_->Add(frameNum, &idata);
}

/** Set turn beginning. */
void Action_DSSP::SSres::SetTurnBegin(TurnType typeIn) {
  if (turnChar_[typeIn] == '<')
    turnChar_[typeIn] = 'X';
  else
    turnChar_[typeIn] = '>';
}

/** Set turn end. */
void Action_DSSP::SSres::SetTurnEnd(TurnType typeIn) {
  if (turnChar_[typeIn] == '>')
    turnChar_[typeIn] = 'X';
  else
    turnChar_[typeIn] = '<';
}

/** Set as a turn of given type. */
void Action_DSSP::SSres::SetTurn(TurnType typeIn) {
  // Do not overwrite an existing end character
  if (turnChar_[typeIn] == ' ') {
    if      (typeIn == T3) turnChar_[typeIn] = '3';
    else if (typeIn == T4) turnChar_[typeIn] = '4';
    else if (typeIn == T5) turnChar_[typeIn] = '5';
  }
}

/** \return Priority of given s.s. type; higher # == higher priority. */
int Action_DSSP::SSres::ssPriority(SStype typeIn) {
  switch (typeIn) {
    case ALPHA    : return 8;
    case BRIDGE   : return 7;
    case EXTENDED : return 6;
    case H3_10    : return 5;
    case HPI      : return 4;
    case TURN     : return 3;
    case BEND     : return 2;
    case NONE     : return 1;
  }
  return 0;
}

/** \return Priority of the current assigned s.s. type. */
int Action_DSSP::SSres::SSpriority() const {
  return ssPriority(sstype_);
}

/** Set residue with given s.s. type. */
void Action_DSSP::SSres::SetSS(SStype typeIn) {
  // TODO check if the priority check is necessary
  //if (ssPriority(typeIn) > ssPriority(sstype_))
    sstype_ = typeIn;
}

/** \return true if this residue is the start of a turn. */
bool Action_DSSP::SSres::HasTurnStart(TurnType typeIn) const {
  if (turnChar_[typeIn] == '>' ||
      turnChar_[typeIn] == 'X')
    return true;
  return false;
}

/** Set this residue as bridging to the specified index with specified type. */
void Action_DSSP::SSres::SetBridge(int idx, BridgeType btypeIn) {
  bridges_.push_back( Bridge(idx, btypeIn) );
  if (bridges_.size() > 2) {
    mprintf("Warning: Residue %i has more than 2 bridges (currently %zu).\n", Num()+1, bridges_.size());
  }
}

/** Determine the dominant bridge type among all bridges. Tie goes to antiparallel. */
Action_DSSP::BridgeType Action_DSSP::SSres::DominantBridgeType() const {
  if (bridges_.empty()) return NO_BRIDGE;
  int npara = 0;
  int nanti = 0;
  for (BridgeArray::const_iterator it = bridges_.begin(); it != bridges_.end(); ++it)
  {
    if (it->Btype() == PARALLEL)
      npara++;
    else if (it->Btype() == ANTIPARALLEL)
      nanti++;
  }
  if (npara > nanti)
    return PARALLEL;
  else
    return ANTIPARALLEL;
}

/** \return True if this residue has at least 1 bridge. */
bool Action_DSSP::SSres::HasBridge() const {
  return !(bridges_.empty());
}

/** \return True if this residue is bridged with the specified index. */
bool Action_DSSP::SSres::IsBridgedWith(int idx2) const {
  for (BridgeArray::const_iterator it = bridges_.begin(); it != bridges_.end(); ++it)
    if (it->Idx() == idx2) return true;
  return false;
}

/*char Action_DSSP::SSres::StrandChar() const {
  // TODO ever b2?
  return ssChar_[B1];
}*/

/** Print a shorthand for the current SS assignment to STDOUT. */
void Action_DSSP::SSres::PrintSSchar() const {
  static const char btypeChar[] = { ' ', 'p', 'A' };
  mprintf("\t%6i %c %c %c %c", num_+1, resChar_, turnChar_[T3], turnChar_[T4], turnChar_[T5]);
  // For backwards compat, always print at least 2.
  unsigned int maxIdx = 2;
  if (bridges_.size() > maxIdx)
    maxIdx = bridges_.size();
  for (unsigned int idx = 0; idx < maxIdx; idx++) {
    if (idx >= bridges_.size())
      mprintf(" %c(%6i)", ' ', 0);
    else
      mprintf(" %c(%6i)", btypeChar[bridges_[idx].Btype()], bridges_[idx].Idx()+1);
  }
  //for (BridgeArray::const_iterator it = bridges_.begin(); it != bridges_.end(); ++it)
  //  mprintf(" %c(%6i)", btypeChar[it->Btype()], it->Idx()+1);
  mprintf(" %6i %6i %6i %c\n", Bcount_[NO_BRIDGE], Bcount_[PARALLEL], Bcount_[ANTIPARALLEL],
          DSSP_char_[sstype_]);
}

#ifdef MPI
/** Sync residue SS and beta counts to master. */
void Action_DSSP::SSres::SyncToMaster(Parallel::Comm const& commIn) {
  int SSprob[NSSTYPE_ + NBRIDGETYPE_];
  int Buffer[NSSTYPE_ + NBRIDGETYPE_];
  std::copy( SScount_, SScount_ + NSSTYPE_,    SSprob );
  std::copy( Bcount_,  Bcount_ + NBRIDGETYPE_, SSprob + NSSTYPE_ );
  commIn.ReduceMaster( Buffer, SSprob, NSSTYPE_ + NBRIDGETYPE_, MPI_INT, MPI_SUM );
  if (commIn.Master()) {
    std::copy( Buffer, Buffer + NSSTYPE_, SScount_ );
    std::copy( Buffer + NSSTYPE_, Buffer + NSSTYPE_ + NBRIDGETYPE_, Bcount_ );
  }
}
#endif

// ----- Action_DSSP -----------------------------------------------------------

Action_DSSP::Action_DSSP() :
  debug_(0),
  BB_N_("N"),
  BB_H_("H"),
  BB_C_("C"),
  BB_O_("O"),
  BB_CA_("CA"),
  SG_("SG"),
  outfile_(0),
  dsspFile_(0),
  assignout_(0),
  printString_(false),
  betaDetail_(false)
{}

// Action_DSSP::Help()
void Action_DSSP::Help() const {
  mprintf("\t[<name>] [out <filename>] [<mask>] [sumout <filename>]\n"
          "\t[assignout <filename>] [totalout <filename>] [ptrajformat]\n"
          "\t[betadetail]\n"
          "\t[namen <N name>] [nameh <H name>] [nameca <CA name>]\n"
          "\t[namec <C name>] [nameo <O name>] [namesg <sulfur name>]\n"
          "  Calculate secondary structure (SS) content for residues in <mask>.\n"
          "  The 'out' file will contain SS vs frame.\n"
          "  The 'sumout' file will contain total SS content for each residue, by SS type.\n"
          "   If sumout not specified, the filename specified by out is used with .sum suffix.\n"
          "  The 'assignout' file will contain the SS assignment foe each residue based\n"
          "   on the majority SS type.\n"
          "  The 'totalout' file will contain overall SS content vs frame, by SS type.\n"
          "  The 'ptrajformat' keyword will use characters instead of #s in the 'out' file.\n"
          "  The 'betadetail' keyword will print parallel/anti-parallel beta instead of\n"
          "   extended/bridge.\n"
          "  The backbone N, H, CA, C, and O atom names can be specifed with 'nameX' keywords.\n"
          "  The disulfide sulfur atom name can be specified with the 'nameSG' keyword.\n"
         );
          
}

// Action_DSSP::Init()
Action::RetType Action_DSSP::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  Nframes_ = 0;
  // Get keywords
  outfile_ = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  std::string temp = actionArgs.GetStringKey("sumout");
  if (temp.empty() && outfile_ != 0) 
    temp = outfile_->DataFilename().Full() + ".sum";
  dsspFile_ = init.DFL().AddDataFile( temp );
  DataFile* totalout = init.DFL().AddDataFile( actionArgs.GetStringKey("totalout"), actionArgs );
  assignout_ = init.DFL().AddCpptrajFile(actionArgs.GetStringKey("assignout"), "SS assignment");
  printString_ = actionArgs.hasKey("ptrajformat");
  betaDetail_ = actionArgs.hasKey("betadetail");
  temp = actionArgs.GetStringKey("namen");
  if (!temp.empty()) BB_N_ = temp;
  temp = actionArgs.GetStringKey("nameh");
  if (!temp.empty()) BB_H_ = temp;
  temp = actionArgs.GetStringKey("namec");
  if (!temp.empty()) BB_C_ = temp;
  temp = actionArgs.GetStringKey("nameo");
  if (!temp.empty()) BB_O_ = temp;
  temp = actionArgs.GetStringKey("nameca");
  if (!temp.empty()) BB_CA_ = temp;
  temp = actionArgs.GetStringKey("namesg");
  if (!temp.empty()) SG_ = temp;
  // Get masks
  if (Mask_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;

  // Set up the DSSP data set
  dsetname_ = actionArgs.GetStringNext();
  // Set default name if none specified
  if (dsetname_.empty())
    dsetname_ = init.DSL().GenerateDefaultName("DSSP");
  // Set up Z labels
  if (outfile_ != 0) {
    if (betaDetail_)
      outfile_->ProcessArgs("zlabels None,Para,Anti,3-10,Alpha,Pi,Turn,Bend");
    else
      outfile_->ProcessArgs("zlabels None,Ext,Bridge,3-10,Alpha,Pi,Turn,Bend");
  }
  // Create data sets for total fraction SS vs time.
  for (int i = 0; i < NSSTYPE_; i++) {
    const char* aspect = DSSP_name_[i];
    if (betaDetail_ && (SStype)i == EXTENDED)
      aspect = "Para";
    else if (betaDetail_ && (SStype)i == BRIDGE)
      aspect = "Anti";
    SSname_.push_back( std::string(aspect) );
    totalDS_[i] = init.DSL().AddSet(DataSet::FLOAT, MetaData(dsetname_, aspect));
    if (totalDS_[i] == 0) {
      mprinterr("Error: Could not create DSSP total frac v time data set.\n");
      return Action::ERR;
    }
    // For now dont add 'None' so colors match up.
    if (totalout != 0 && i > 0) totalout->AddDataSet( totalDS_[i] );
  }
# ifdef _OPENMP
  // Each thread needs temp. space to store found hbonds every frame
  // to avoid memory clashes when adding/updating in map.
# pragma omp parallel
  {
# pragma omp master
  {
  CO_NH_bondsArray_.resize( omp_get_num_threads() );
  }
  }
# endif

  mprintf( "    SECSTRUCT: Calculating secondary structure using mask [%s]\n",Mask_.MaskString());
# ifdef _OPENMP
  mprintf("\tParallelizing calculation with %zu threads.\n", CO_NH_bondsArray_.size());
# endif
  if (outfile_ != 0) 
    mprintf("\tDumping results to %s\n", outfile_->DataFilename().full());
  if (dsspFile_ != 0)
    mprintf("\tSum results to %s\n", dsspFile_->DataFilename().full());
  if (betaDetail_)
    mprintf("\tWill print parallel/anti-parallel beta in place of extended/bridge\n");
  if (printString_) { 
    mprintf("\tSS data for each residue will be stored as a string.\n");
    for (int i = 0; i < NSSTYPE_; i++)
      mprintf("\t\t%s = %s\n", SSchar_[i], SSname_[i].c_str());
  } else {
    mprintf("\tSS data for each residue will be stored as integers.\n");
    for (int i = 0; i < NSSTYPE_; i++)
      mprintf("\t\t%i = %s\n", i, SSname_[i].c_str());
  }
  if (assignout_ != 0)
    mprintf("\tOverall assigned SS will be written to %s\n", assignout_->Filename().full());
  mprintf("\tBackbone Atom Names: N=[%s]  H=[%s]  C=[%s]  O=[%s]  CA=[%s]\n",
          *BB_N_, *BB_H_, *BB_C_, *BB_O_, *BB_CA_ );
  mprintf("\tDisulfide sulfur atom name: %s\n", *SG_);
  mprintf("# Citation: Kabsch, W.; Sander, C.; \"Dictionary of Protein Secondary Structure:\n"
          "#           Pattern Recognition of Hydrogen-Bonded and Geometrical Features.\"\n"
          "#           Biopolymers (1983), V.22, pp.2577-2637.\n" );
  init.DSL().SetDataSetsPending(true);
  Init_ = init;
  return Action::OK;
}

/** Print the atom name and "mask name" to stdout. */
static inline void PrintAtom(Topology const& top, NameType const& name, int idx) {
  if (idx > -1)
    mprintf(" '%s'=%-12s", *name, top.AtomMaskName(idx/3).c_str());
  else
    mprintf(" '%s'=%-12s", *name, "NONE");
}

// TODO use Num()?
static inline int AbsResDelta(int idx1, int idx2) {
  int resDelta =  idx1 - idx2;
  if (resDelta < 0) resDelta = -resDelta;
  return resDelta;
}

/** Print a warning about multiple selected residues.
  * \return the index that minimizes the absolute delta between residues.
  */ 
static inline int MultiResWarning(int resNum, int currNum, int newNum, const char* typeStr)
{
  int currentDelta = AbsResDelta(resNum, currNum);
  int newDelta = AbsResDelta(resNum, newNum);
  mprintf("Warning: Multiple %s residues for res %i (current %s is %i; potential %s is %i).\n",
          typeStr, resNum+1, typeStr, currNum+1, typeStr, newNum+1);
  mprintf("Warning: Old delta= %i  New delta= %i\n", currentDelta, newDelta);
  int selectedNum, rejectedNum;
  if (currentDelta < newDelta) {
    selectedNum = currNum;
    rejectedNum = newNum;
  } else {
    selectedNum = newNum;
    rejectedNum = currNum;
  }
  mprintf("Warning: Residues %i and %i will not be considered consecutive.\n",
          resNum+1, rejectedNum+1);
  return selectedNum;
}

/** Find the residue number of the atom with target name bonded to specified atom/residue.
  * \param topIn Current topology
  * \param at0 Specified atom
  * \param res0 Specified residue
  * \param tgtName Name of the desired bonded atom.
  * \param typeStr Indicate if we're looking for previous/next
  * \return Residue index of atom with target name bonded to specified atom.
  */
static inline int FindResNum(Topology const& topIn, int at0, int res0, NameType const& tgtName, const char* typeStr)
{
  int tgtres = -1;
  // Find tgtName in atoms bonded to this res.
  for (Atom::bond_iterator ib = topIn[ at0 ].bondbegin(); ib != topIn[ at0 ].bondend(); ++ib)
  {
    if ( topIn[*ib].Name() == tgtName && topIn[*ib].ResNum() != res0 ) {
      int bres = topIn[*ib].ResNum();
      if (tgtres == -1)
        tgtres = bres;
      else
        tgtres = MultiResWarning(res0, tgtres, bres, typeStr);
    }
  }
  return tgtres;
}

// Action_DSSP::Setup()
/** Set up secondary structure calculation for all residues selected by the
  * mask. A residue is selected if at least one atom in the residue is
  * selected by the mask. The coordinate indices (i.e. atom # * 3) for
  * the C, O, N, H, and CA atoms are set up if those atoms are present.
  * The residue is only initialized if it has not been previously selected
  * and set up by a prior call to Setup().
  */
Action::RetType Action_DSSP::Setup(ActionSetup& setup)
{
  // Set up mask for this parm
  if ( setup.Top().SetupCharMask( Mask_ ) ) return Action::ERR;
  if ( Mask_.None() ) {
    mprintf("Warning: DSSP: Mask has no atoms.\n");
    return Action::SKIP;
  }

  // Deselect any existing residues
  for (SSarrayType::iterator it = Residues_.begin(); it != Residues_.end(); ++it)
    it->Deselect();
  // Set up for each solute residue 
  Range soluteRes = setup.Top().SoluteResidues();
  mprintf("\tSetting up for %i solute residues.\n", soluteRes.Size());
  if ((unsigned int)soluteRes.Size() > Residues_.size())
    Residues_.resize( soluteRes.Size() );
  MetaData md(dsetname_, "res");
  DataSet::DataType dt;
  if (printString_)
    dt = DataSet::STRING;
  else
    dt = DataSet::INTEGER;
  unsigned int nResSelected = 0;
  SSarrayType::iterator Res = Residues_.begin();
  for (Range::const_iterator ridx = soluteRes.begin(); ridx != soluteRes.end(); ++ridx, ++Res)
  {
    Residue const& thisRes = setup.Top().Res( *ridx );
    if (Res->Num() != -1) {
      // Residue has previously been set up. Check that indices match.
      if (Res->Num() != *ridx) {
        mprinterr("Error: Solute residue index %i does not match previously setup\n"
                  "Error: index %i\n", *ridx, Res->Num());
        return Action::ERR;
      }
    } else {
      // Set up Residue. TODO also molecule index?
      Res->SetNum( *ridx );
      Res->SetResChar( thisRes.SingleCharName() );
    }
    // Search for atoms in residue.
    bool hasAtoms = false;
    int prevresnum = -1;
    int nextresnum = -1;
    int nAtSelected = 0;
    for (int at = thisRes.FirstAtom(); at != thisRes.LastAtom(); at++)
    {
      if (Mask_.AtomInCharMask( at )) {
        nAtSelected++;
        if      ( setup.Top()[at].Name() == BB_C_ )  { Res->SetC( at*3 ); hasAtoms = true; }
        else if ( setup.Top()[at].Name() == BB_O_ )  { Res->SetO( at*3 ); hasAtoms = true; }
        else if ( setup.Top()[at].Name() == BB_N_ )  { Res->SetN( at*3 ); hasAtoms = true; }
        else if ( setup.Top()[at].Name() == BB_H_ )  { Res->SetH( at*3 ); hasAtoms = true; }
        else if ( setup.Top()[at].Name() == BB_CA_ ) { Res->SetCA( at*3 ); hasAtoms = true; }
      }
    }
    if (!hasAtoms) {
      if (nAtSelected > 0)
        mprintf("Warning: No atoms selected for res %s; skipping.\n",
                setup.Top().TruncResNameNum( Res->Num() ).c_str());
    } else {
      Res->SetSelected( true );
      ++nResSelected;
      // Determine previous residue
      if ( Res->N() != -1 ) {
        // Find C in previous res bonded to this N.
        prevresnum = FindResNum(setup.Top(), Res->N()/3, Res->Num(), BB_C_, "previous");
      }
      // Determine next residue
      if ( Res->C() != -1) {
        // Find N in next res bonded to this C.
        nextresnum = FindResNum(setup.Top(), Res->C()/3, Res->Num(), BB_N_, "next");
      }
#     ifdef DSSPDEBUG
      mprintf("\t %8i < %8i < %8i\n", prevresnum+1, Res->Num()+1, nextresnum+1);
#     endif
      // Set previous/next residue indices 
      if (prevresnum > -1) Res->SetPrevIdx( prevresnum );
      if (nextresnum > -1) Res->SetNextIdx( nextresnum );
      // Report if residue is missing atoms
      if (Res->IsMissingAtoms()) {
        mprintf("Warning: Res %s is missing atoms", setup.Top().TruncResNameNum( Res->Num() ).c_str());
        if (Res->C() == -1)  mprintf(" %s", *BB_C_);
        if (Res->O() == -1)  mprintf(" %s", *BB_N_);
        if (Res->N() == -1)  mprintf(" %s", *BB_O_);
        if (Res->H() == -1)  mprintf(" %s", *BB_H_);
        if (Res->CA() == -1) mprintf(" %s", *BB_CA_);
        mprintf("\n");
      }
      // Set up DataSet if necessary
      if (Res->Dset() == 0) {
        md.SetIdx( *ridx+1 );
        md.SetLegend( setup.Top().TruncResNameNum( *ridx ) );
        // Setup DataSet for this residue
        Res->SetDset( Init_.DSL().AddSet( dt, md ) );
        if (Res->Dset() == 0) {
          mprinterr("Error: Could not allocate DSSP data set for residue %i\n", *ridx+1);
          return Action::ERR;
        }
        if (outfile_ != 0) outfile_->AddDataSet( Res->Dset() );
      }
    } // END residue has selected atoms
  } // END loop over residues
  mprintf("\t%u of %i solute residues selected.\n", nResSelected, soluteRes.Size());

  // DEBUG - print each residue set up.
  if (debug_ > 0) {
    for (SSarrayType::const_iterator it = Residues_.begin(); it != Residues_.end(); ++it)
    {
      mprintf("    %8i", it->Num() + 1);
      PrintAtom(setup.Top(), BB_C_, it->C());
      PrintAtom(setup.Top(), BB_O_, it->O());
      PrintAtom(setup.Top(), BB_N_, it->N());
      PrintAtom(setup.Top(), BB_H_, it->H());
      PrintAtom(setup.Top(), BB_CA_, it->CA());
      mprintf(" Prev=%8i Next=%8i", it->PrevIdx(), it->NextIdx());
      mprintf("\n");
    }
  }

  return Action::OK;
}

/** Assign a bridge between the two specified residues. */
void Action_DSSP::AssignBridge(int idx1in, int idx2in, BridgeType btypeIn) {
  // By convention, always make idx1 the lower one
  int idx1, idx2;
  if (idx1in < idx2in) {
    idx1 = idx1in;
    idx2 = idx2in;
  } else {
    idx1 = idx2in;
    idx2 = idx1in;
  }

  SSres& Resi = Residues_[idx1];
  SSres& Resj = Residues_[idx2];
# ifdef DSSPDEBUG
  if (btypeIn == ANTIPARALLEL) {
    mprintf("\t\tAssignBridge %i to %i, Antiparallel\n", idx1+1, idx2+1);
  } else {
    mprintf("\t\tAssignBridge %i to %i, Parallel\n", idx1+1, idx2+1);
  }
# endif
  // Do not duplicate bridges
  if (Resi.IsBridgedWith(idx2)) {
#   ifdef DSSPDEBUG
    mprintf("\t\tAlready present.\n");
#   endif
    return;
  }

  Resi.SetBridge( idx2, btypeIn );
  Resj.SetBridge( idx1, btypeIn );
}

/*
void Action_DSSP::AssignBridge(int idx1in, int idx2in, BridgeType btypeIn, char& currentStrandChar) {
  // By convention, always make idx1 the lower one
  int idx1, idx2;
  if (idx1in < idx2in) {
    idx1 = idx1in;
    idx2 = idx2in;
  } else {
    idx1 = idx2in;
    idx2 = idx1in;
  }
  if (btypeIn == ANTIPARALLEL)
    mprintf("\t\tAssignBridge %i to %i, Antiparallel\n", idx1+1, idx2+1);
  else
    mprintf("\t\tAssignBridge %i to %i, Parallel\n", idx1+1, idx2+1);
  SSres& Resi = Residues_[idx1];
  SSres& Resj = Residues_[idx2];
  // Do not duplicate bridges
  if (Resi.IsBridgedWith(idx2)) {
    mprintf("\t\tAlready present.\n");
    return;
  }
  // Determine if we are already part of a ladder.
  char ladderChar = ' ';
  char resiLadderChar = Residues_[Resi.PrevIdx()].StrandChar();
  if (resiLadderChar == ' ')
    resiLadderChar = Residues_[Resi.NextIdx()].StrandChar();
  char resjLadderChar = Residues_[Resj.PrevIdx()].StrandChar();
  if (resjLadderChar == ' ')
    resjLadderChar = Residues_[Resj.NextIdx()].StrandChar();
  if (resiLadderChar == ' ' && resjLadderChar == ' ') {
    // If both are blank, new ladder.
    ladderChar = currentStrandChar;
    ++currentStrandChar;
  } else if (resiLadderChar != resjLadderChar) {
    // They do not match. New ladder.
    ladderChar = currentStrandChar;
    ++currentStrandChar;
  } else
    ladderChar = resiLadderChar;
  //else if (resiLadderChar != ' ')
  //  ladderChar = resiLadderChar;
  //else
  //  ladderChar = resjLadderChar;
  // Set the bridge; adjust character case if needed
  if (btypeIn == ANTIPARALLEL)
    ladderChar = toupper( ladderChar );
  mprintf("\t\tResi strand %c, resj strand %c, ladder char %c\n", resiLadderChar, resjLadderChar, ladderChar);
  Resi.SetBridge( idx2, ladderChar );
  Resj.SetBridge( idx1, ladderChar );
}
*/

/** Given that the current residue is in-between two bridging residues, check
  * if there is a bulge in the strand(s) it is bridging with. The largest
  * gap allowed in the other strand is 5 residues.
  */
void Action_DSSP::CheckBulge(int prevIdx, int currentIdx, int nextIdx) {
  if (prevIdx == -1 || nextIdx == -1) return;
  SSres const& prevRes = Residues_[prevIdx];
  SSres const& nextRes = Residues_[nextIdx];

  // Loop over indices bonded to the previous residue
  for (BridgeArray::const_iterator pr = prevRes.Bridges().begin(); pr != prevRes.Bridges().end(); ++pr)
  {
    // Loop over indices bonded to the next residue
    for (BridgeArray::const_iterator nr = nextRes.Bridges().begin(); nr != nextRes.Bridges().end(); ++nr)
    {
      int resGapSize = AbsResDelta( pr->Idx(), nr->Idx() );
#     ifdef DSSPDEBUG
      mprintf("DEBUG: Checking potential bulge for %i: %i type %i to %i type %i (%i)", currentIdx+1,
              pr->Idx()+1, (int)pr->Btype(), nr->Idx()+1, (int)nr->Btype(), resGapSize);
#     endif
      // Minimum allowed gap is 4 residues in between, so 5 residues total.
      // Types also need to match
      if (resGapSize < 6 && pr->Btype() == nr->Btype()) {
#       ifdef DSSPDEBUG
        mprintf(" Found!\n");
#       endif
        if (Residues_[prevIdx].SS() != ALPHA)
          Residues_[prevIdx].SetSS( EXTENDED );
        Residues_[currentIdx].SetSS( EXTENDED );
        if (Residues_[nextIdx].SS() != ALPHA)
          Residues_[nextIdx].SetSS( EXTENDED );
        // Set extended on other strand as well
        int sres0, sres1;
        if (pr->Idx() < nr->Idx()) {
          sres0 = pr->Idx();
          sres1 = nr->Idx();
        } else {
          sres0 = nr->Idx();
          sres1 = pr->Idx();
        }
        for (int sres = sres0; sres != sres1; sres++)
          if (Residues_[sres].SS() != ALPHA)
            Residues_[sres].SetSS( EXTENDED );
      }
#     ifdef DSSPDEBUG
      else
        mprintf("\n");
#     endif
    } // END loop over next residue bridging residues
  } // END loop over previous residue bridging residues
}

/** Determine CO-NH hbonds, loop over them to do SS assignment. */
int Action_DSSP::OverHbonds(int frameNum, ActionFrame& frm)
{
  t_total_.Start();
  t_calchb_.Start();
  // ----- Determine hydrogen bonding ------------
  int resi;
  int Nres = (int)Residues_.size();
#ifdef _OPENMP
  int mythread;
#pragma omp parallel private(resi, mythread)
{
  mythread = omp_get_thread_num();
  CO_NH_bondsArray_[mythread].clear();
#pragma omp for
#else /* _OPENMP */
  CO_NH_bonds_.clear();
#endif /* _OPENMP */
  for (resi = 0; resi < Nres; resi++)
  {
    Residues_[resi].Unassign();
    if (Residues_[resi].IsSelected())
    {
      SSres& ResCO = Residues_[resi];
      if (ResCO.HasCO())
      {
        const double* Cxyz = frm.Frm().CRD( ResCO.C() );
        const double* Oxyz = frm.Frm().CRD( ResCO.O() );
        for (int resj = 0; resj < Nres; resj++)
        {
          // resj is the N-H residue. resi is the C=O residue.
          // Ignore j-i = 0 (same residue) and (j-i = 1) peptide bond.
          int resDelta = resj - resi;
          if (resDelta < 0 || resDelta > 1) {
            SSres& ResNH = Residues_[resj];
            if (ResNH.IsSelected() && ResNH.HasNH())
            {
              const double* Nxyz = frm.Frm().CRD( ResNH.N() );
              const double* Hxyz = frm.Frm().CRD( ResNH.H() );

              double rON = 1.0/sqrt(DIST2_NoImage(Oxyz, Nxyz));
              double rCH = 1.0/sqrt(DIST2_NoImage(Cxyz, Hxyz));
              double rOH = 1.0/sqrt(DIST2_NoImage(Oxyz, Hxyz));
              double rCN = 1.0/sqrt(DIST2_NoImage(Cxyz, Nxyz));

              double E = DSSP_fac_ * (rON + rCH - rOH - rCN);
              if (E < DSSP_cut_) {
#               ifdef DSSPDEBUG
                mprintf("DBG: %i-CO --> %i-NH  E= %g\n", resi+1, resj+1, E);
#               endif
#               ifdef _OPENMP
                CO_NH_bondsArray_[mythread].insert( HbondPairType(resi, resj) );
#               else
                CO_NH_bonds_.insert( HbondPairType(resi, resj) );
#               endif
              }
//#             ifdef DSSPDEBUG
//              else if (resDelta < 6)
//                mprintf("DBG: No hbond %i-CO --> %i-NH  E= %g\n", resi+1, resj+1, E);
//#             endif
            } // END ResNH selected
          } // END residues spaced > 2 apart
        } // END inner loop over residues
      } // END has CO
    } // END ResCO selected
  } // END outer loop over residues
#ifdef _OPENMP
} // END pragma omp parallel
  for (unsigned int thread = 1; thread < CO_NH_bondsArray_.size(); thread++)
    for (HbondMapType::const_iterator hb = CO_NH_bondsArray_[thread].begin();
                                      hb != CO_NH_bondsArray_[thread].end(); ++hb)
      CO_NH_bondsArray_[0].insert( *hb );
  HbondMapType const& CO_NH_bonds_ = CO_NH_bondsArray_[0];
#endif
  t_calchb_.Stop();

  t_assign_.Start();
  // ----- Do basic assignment -------------------
  for (HbondMapType::const_iterator hb0 = CO_NH_bonds_.begin(); hb0 != CO_NH_bonds_.end(); ++hb0)
  {
    int riidx = hb0->first;
    int rjidx = hb0->second;
    SSres const& Resi = Residues_[riidx];
    SSres const& Resj = Residues_[rjidx];
#   ifdef DSSPDEBUG
    mprintf("Res %8i %c -- %8i %c", Resi.Num()+1, Resi.ResChar(),
                                    Resj.Num()+1, Resj.ResChar()); // DBG
#   endif
    // Spacing between residues i and j
    int resDelta = Resj.Num() - Resi.Num();
#   ifdef DSSPDEBUG
    mprintf("(%4i)\n", resDelta);
#   endif
    // Check for H bond from CO i to NH i+n
    if (resDelta == 3) {
      // 3-TURN
      Residues_[riidx  ].SetTurnBegin(T3);
      Residues_[riidx+1].SetTurn(T3);
      Residues_[riidx+2].SetTurn(T3);
      Residues_[riidx+3].SetTurnEnd(T3);
      Residues_[riidx+1].SetSS( TURN );
      Residues_[riidx+2].SetSS( TURN );
    } else if (resDelta == 4) {
      // 4-TURN
      Residues_[riidx  ].SetTurnBegin(T4);
      Residues_[riidx+1].SetTurn(T4);
      Residues_[riidx+2].SetTurn(T4);
      Residues_[riidx+3].SetTurn(T4);
      Residues_[riidx+4].SetTurnEnd(T4);
      Residues_[riidx+1].SetSS( TURN );
      Residues_[riidx+2].SetSS( TURN );
      Residues_[riidx+3].SetSS( TURN );
    } else if (resDelta == 5) {
      // 5-TURN
      Residues_[riidx  ].SetTurnBegin(T5);
      Residues_[riidx+1].SetTurn(T5);
      Residues_[riidx+2].SetTurn(T5);
      Residues_[riidx+3].SetTurn(T5);
      Residues_[riidx+4].SetTurn(T5);
      Residues_[riidx+5].SetTurnEnd(T5);
      Residues_[riidx+1].SetSS( TURN );
      Residues_[riidx+2].SetSS( TURN );
      Residues_[riidx+3].SetSS( TURN );
      Residues_[riidx+4].SetSS( TURN );
    }
    // Look for bridge. Start with the premise that this bond is part of one
    // of the 4 potential bridge patterns, then check if the compliment exists.
    HbondMapType::iterator hb;
    // Here we want absolute value of spacing.
    if (resDelta < 0) resDelta = -resDelta;
    if (resDelta > 2) {
      // Assume (i,j). Look for (j,i)
      hb = CO_NH_bonds_.find( HbondPairType(hb0->second, hb0->first) );
      if (hb != CO_NH_bonds_.end()) {
#       ifdef DSSPDEBUG
        mprintf("\t\t%i ANTI-PARALLELa with %i (%i to %i)\n", hb0->first+1, hb0->second+1, hb->first+1, hb->second+1);
#       endif
        AssignBridge(hb0->first, hb0->second, ANTIPARALLEL);
      }
    }
    resDelta = AbsResDelta(hb0->first+1, hb0->second-1);
    if (resDelta > 2) {
      // Assume (i-1, j+1). Look for (j-1, i+1)
      hb = CO_NH_bonds_.find( HbondPairType(hb0->second-2, hb0->first+2) );
      if (hb != CO_NH_bonds_.end()) {
#       ifdef DSSPDEBUG
        mprintf("\t\t%i ANTI-PARALLELb with %i (%i to %i)\n", hb0->first+2, hb0->second, hb->first+1, hb->second+1);
#       endif
        AssignBridge(hb0->first+1, hb0->second-1, ANTIPARALLEL);
      }
    }
    resDelta = AbsResDelta(hb0->first+1, hb0->second);
    if (resDelta > 2) {
      // Assume (i-1, j). Check for (j, i+1) PARALLEL
      hb = CO_NH_bonds_.find( HbondPairType(hb0->second, hb0->first+2) );
      if (hb != CO_NH_bonds_.end()) {
#       ifdef DSSPDEBUG
        mprintf("\t\t%i PARALLELa with %i (%i to %i)\n", hb0->first+2, hb0->second+1, hb->first+1, hb->second+1);
#       endif
        AssignBridge(hb0->first+1, hb0->second, PARALLEL);
      }
    }
    resDelta = AbsResDelta(hb0->second-1, hb0->first);
    if (resDelta > 2) {
      // Assume (j, i+1). Check for (i-1, j)
      hb = CO_NH_bonds_.find( HbondPairType(hb0->second-2, hb0->first) );
      if (hb != CO_NH_bonds_.end()) {
#       ifdef DSSPDEBUG
        mprintf("\t\t%i PARALLELb with %i (%i to %i)\n", hb0->second, hb0->first+1, hb->first+1, hb->second+1);
#       endif
        AssignBridge(hb0->second-1, hb0->first, PARALLEL);
      }
    }
  } // END loop over Hbonds

  // ----- Do SS assignment ----------------------
  // Priority is 'H', 'B', 'E', 'G', 'I', 'T', 'S' None
  //              8    7    6    5    4    3    2  1
  // First do alpha, extended, and bridge.
  for (resi = 0; resi < Nres; resi++)
  {
#   ifdef DSSPDEBUG
    mprintf("Residue %i\n", resi+1);
#   endif
    SSres& Resi = Residues_[resi];
    int prevIdx = Resi.PrevIdx();
    int nextIdx = Resi.NextIdx();
    int priority = Resi.SSpriority();
    if ( Resi.HasTurnStart(T4) && prevIdx != -1 && Residues_[prevIdx].HasTurnStart(T4) )
    {
      // Alpha helix.
#     ifdef DSSPDEBUG
      mprintf("ALPHA helix starting at %i\n", resi+1);
#     endif
      Residues_[resi  ].SetSS( ALPHA );
      Residues_[resi+1].SetSS( ALPHA );
      Residues_[resi+2].SetSS( ALPHA );
      Residues_[resi+3].SetSS( ALPHA );
    } else if (Resi.SS() != ALPHA) {
      if (priority < 6) {
        // Priority < 6 means not alpha or beta assigned yet.
        // Check for Beta structure
        bool prevHasBridge, prevHasExtended, nextHasBridge;
        if (prevIdx != -1) {
          prevHasBridge = Residues_[prevIdx].HasBridge();
          prevHasExtended = Residues_[prevIdx].SS() == EXTENDED;
        } else {
          prevHasBridge = false;
          prevHasExtended = false;
        }
        if (nextIdx != -1)
          nextHasBridge = Residues_[nextIdx].HasBridge();
        else
          nextHasBridge = false;
        if (Resi.HasBridge()) {
          // Regular Beta. Check if previous is assigned EXTENDED in case it 
          // was assigned via a Beta bulge.
          if ( prevHasBridge || nextHasBridge || prevHasExtended )
          {
#           ifdef DSSPDEBUG
            mprintf("Extended BETA bridge at %i\n", resi+1);
#           endif
            Resi.SetSS( EXTENDED );
          } else {
#           ifdef DSSPDEBUG
            mprintf("Isolated BETA bridge at %i.\n", resi+1);
#           endif
            Resi.SetSS( BRIDGE );
          }
        } else if (prevHasBridge && nextHasBridge) {
          // Potential Beta bulge. Check other strand.
          CheckBulge(prevIdx, resi, nextIdx);
        }
      } // END check for Beta structure
    } // END not alpha
  } // END loop over residues

  // Check for 3-10 helices. Do this separately so we dont assign regions
  // that are too small because other residues have already been assigned.
  for (resi = 1; resi < Nres-2; resi++) {
    if (Residues_[resi  ].SSpriority() < 6 &&
        Residues_[resi+1].SSpriority() < 6 &&
        Residues_[resi+2].SSpriority() < 6 &&
        Residues_[resi  ].HasTurnStart(T3) &&
        Residues_[resi-1].HasTurnStart(T3))
    {
      // 3-10 helix
#     ifdef DSSPDEBUG
      mprintf("3-10 helix starting at %i\n", resi+1);
#     endif
      Residues_[resi  ].SetSS( H3_10 );
      Residues_[resi+1].SetSS( H3_10 );
      Residues_[resi+2].SetSS( H3_10 );
    }
  }
  // Check for PI helices, similar to 3-10 helices.
  for (resi = 1; resi < Nres-4; resi++) {
    if (Residues_[resi  ].SSpriority() < 5 &&
        Residues_[resi+1].SSpriority() < 5 &&
        Residues_[resi+2].SSpriority() < 5 &&
        Residues_[resi+3].SSpriority() < 5 &&
        Residues_[resi+4].SSpriority() < 5 &&
        Residues_[resi  ].HasTurnStart(T5) &&
        Residues_[resi-1].HasTurnStart(T5))
    {
      // PI helix
#     ifdef DSSPDEBUG
      mprintf("PI helix starting at %i\n", resi+1);
#     endif
      Residues_[resi  ].SetSS( HPI );
      Residues_[resi+1].SetSS( HPI );
      Residues_[resi+2].SetSS( HPI );
      Residues_[resi+3].SetSS( HPI );
      Residues_[resi+4].SetSS( HPI );
    }
  }
  // Check for bends. Only do if no other assignment.
  for (resi = 0; resi < Nres; resi++) {
    if (Residues_[resi].IsSelected() && Residues_[resi].SS() == NONE)
    {
      int im2 = resi - 2;
      if (im2 > -1) {
        int ip2 = resi + 2;
        if (ip2 < Nres) {
          SSres& Resi = Residues_[resi];
          if (Residues_[im2].CA() != -1 && Resi.CA() != -1 && Residues_[ip2].CA() != -1) {
            const double* CAm2 = frm.Frm().CRD(Residues_[im2].CA());
            const double* CA0  = frm.Frm().CRD(Resi.CA());
            const double* CAp2 = frm.Frm().CRD(Residues_[ip2].CA());
            Vec3 CA1( CA0[0]-CAm2[0], CA0[1]-CAm2[1], CA0[2]-CAm2[2] );
            Vec3 CA2( CAp2[0]-CA0[0], CAp2[1]-CA0[1], CAp2[2]-CA0[2] );
            CA1.Normalize();
            CA2.Normalize();
            double bAngle = CA1.Angle(CA2);
#           ifdef DSSPDEBUG
            mprintf("DEBUG: Bend calc %i-%i-%i: %g deg.\n", resi-1, resi+1, resi+3, bAngle*Constants::RADDEG);
#           endif
            // 1.221730476 rad = 70 degrees
            if (bAngle > 1.221730476) {
              Resi.SetSS( BEND );
            }
          }
        }
      }
    } // END selected and no assignment 
  }
 
  // ----- Store data for each res. Get stats ----
  int totalSS[NSSTYPE_];
  std::fill( totalSS, totalSS + NSSTYPE_, 0 ); 
  int Nselected = 0;
  for (resi=0; resi < Nres; resi++) {
    SSres& Resi = Residues_[resi];
    if (Resi.IsSelected()) {
      if (betaDetail_ && (Resi.SS() == EXTENDED || Resi.SS() == BRIDGE))
      {
        BridgeType btype = Resi.DominantBridgeType();
        if (btype == ANTIPARALLEL)
          totalSS[BRIDGE]++;
        else if (btype == PARALLEL)
          totalSS[EXTENDED]++;
      } else
        totalSS[Residues_[resi].SS()]++;
      Residues_[resi].AccumulateData(frameNum, printString_, betaDetail_);
      Nselected++;
    }
  }
  for (int i = 0; i < NSSTYPE_; i++) {
    float fvar = (float)totalSS[i];
    fvar /= (float)Nselected;
    totalDS_[i]->Add(frameNum, &fvar);
  }


  t_assign_.Stop();
  t_total_.Stop();
  return 0;
}

Action::RetType Action_DSSP::DoAction(int frameNum, ActionFrame& frm)
{
  OverHbonds(frameNum, frm);
  Nframes_++;
  // DEBUG - Print basic assignment
  if (debug_ > 1) {
    for (SSarrayType::const_iterator it = Residues_.begin(); it != Residues_.end(); ++it)
      it->PrintSSchar();
  }
  return Action::OK;
}

#ifdef MPI
int Action_DSSP::SyncAction() {
  // Consolidate SScount data to master.
  for (SSarrayType::iterator res = Residues_.begin(); res != Residues_.end(); ++res)
    res->SyncToMaster( Init_.TrajComm() );

  // Calc total number of frames.
  int total_frames = 0;
  Init_.TrajComm().ReduceMaster( &total_frames, &Nframes_, 1, MPI_INT, MPI_SUM );
  if (Init_.TrajComm().Master())
    Nframes_ = total_frames;
  return 0;
}
#endif

// Action_DSSP::Print()
void Action_DSSP::Print() {
  if (dsetname_.empty()) return;
  t_total_.WriteTiming(1,"DSSP total");
  t_calchb_.WriteTiming(2, "Calc Hbonds", t_total_.Total());
  t_assign_.WriteTiming(2, "Assignment ", t_total_.Total());

  // Try not to print empty residues. Find the minimum and maximum residue
  // for which there is data. Output res nums start from 1.
  int min_res = -1;
  int max_res = -1;
  for (int resi = 0; resi != (int)Residues_.size(); resi++) {
    if (Residues_[resi].Dset() != 0) {
      if (min_res < 0) min_res = resi;
      if (resi > max_res) max_res = resi;
    }
  }
  if (min_res < 0 || max_res < min_res) {
    mprinterr("Error: No residues have SS data.\n");
    return;
  }
  // Calculate average of each SS type across all residues.
  if (dsspFile_ != 0) {
    std::vector<DataSet*> dsspData_(NSSTYPE_);
    Dimension Xdim( min_res + 1, 1, "Residue" );
    MetaData md(dsetname_, "avgss", MetaData::NOT_TS);
    // Set up a dataset for each SS type. TODO: NONE type?
    for (int ss = 1; ss < NSSTYPE_; ss++) {
      md.SetIdx(ss);
      md.SetLegend( SSname_[ss] );
      dsspData_[ss] = Init_.DSL().AddSet(DataSet::DOUBLE, md);
      dsspData_[ss]->SetDim(Dimension::X, Xdim);
      dsspFile_->AddDataSet( dsspData_[ss] ); 
    }
    
    // Calc the avg SS type for each residue that has data.
    int idx = 0;
    unsigned int norm = Nframes_;
    for (int resi = min_res; resi < max_res+1; resi++) {
      SSres& Resi = Residues_[resi];
      if (Resi.Dset() != 0) {
        //int Nframe = 0;
        //for (int ss = 0; ss < NSSTYPE_; ss++)
        //  Nframe += Resi.SScount((SStype)ss);
        //mprintf("DEBUG: Total frames for residue %i is %i\n", Resi.Num()+1, Nframe);
        for (int ss = 1; ss < NSSTYPE_; ss++) {
          double avg;
          if (betaDetail_ && (SStype)ss == EXTENDED)
            avg = (double)Resi.Bcount(PARALLEL);
          else if (betaDetail_ && (SStype)ss == BRIDGE)
            avg = (double)Resi.Bcount(ANTIPARALLEL);
          else
            avg = (double)Resi.SScount((SStype)ss);
          //mprintf("DEBUG:\t\tCount for type %i is %f\n", ss, avg);
          avg /= (double)norm;
          dsspData_[ss]->Add(idx, &avg);
        }
        ++idx;
      }
    }
  }
  // Print out SS assignment like PDB
  if (assignout_ != 0) {
      int total = 0;
      int startRes = -1;
      std::string resLine, ssLine;
      for (int resi = min_res; resi < max_res+1; resi++) {
        if (startRes == -1) startRes = resi;
        // Convert residue name.
        SSres& Resi = Residues_[resi];
        resLine += Resi.ResChar();
        // Figure out which SS element is dominant for res if selected
        if (Resi.Dset() != 0) {
          int dominantType = 0;
          int ssmax = 0;
          for (int ss = 0; ss < NSSTYPE_; ss++) {
            int sscount;
            if (betaDetail_ && (SStype)ss == EXTENDED)
              sscount = Resi.Bcount(PARALLEL);
            else if (betaDetail_ && (SStype)ss == BRIDGE)
              sscount = Resi.Bcount(ANTIPARALLEL);
            else
              sscount = Resi.SScount((SStype)ss);
            if ( sscount > ssmax ) {
              ssmax = sscount;
              dominantType = ss;
            }
          }
          ssLine += DSSP_char_[dominantType];
        } else
          ssLine += '-';
        total++;
        if ((total % 50) == 0 || resi == max_res) {
          assignout_->Printf("%-8i %s\n", startRes+1, resLine.c_str());
          assignout_->Printf("%8s %s\n\n", " ", ssLine.c_str());
          startRes = -1;
          resLine.clear();
          ssLine.clear();
        } else if ((total % 10) == 0) {
          resLine += ' '; 
          ssLine += ' ';
        }
      }
  }
}
