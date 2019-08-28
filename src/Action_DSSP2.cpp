#include <cmath> // sqrt
#include <algorithm> // std::fill, std::min, std::max
#include <set> // SET
#include "Action_DSSP2.h"
#include "CpptrajStdio.h"
#include "DataSet.h"
#include "DistRoutines.h"
#include "Timer.h"
#ifdef DSSPDEBUG
#include "Constants.h"
#endif

/** From the original Kabsch & Sander 1983 paper. Obtained via
  *   q1 = 0.42e
  *   q2 = 0.20e
  *   f = 332 (approximate conversion factor to get energies in kcal/mol)
  *   fac = q1*q2*f
  */
const double Action_DSSP2::DSSP_fac_ = 27.888;

const double Action_DSSP2::DSSP_cut_ = -0.5;

const char  Action_DSSP2::DSSP_char_[] = { ' ', 'E', 'B', 'G', 'H', 'I', 'T', 'S' };

const char* Action_DSSP2::SSchar_[]    = { "0", "E", "B", "G", "H", "I", "T", "S" };

const char* Action_DSSP2::SSname_[]={"None", "Extended", "Bridge", "3-10", "Alpha", "Pi", "Turn", "Bend"};

const std::string Action_DSSP2::SSzlabels_ = "zlabels None,Ext,Bridge,3-10,Alpha,Pi,Turn,Bend";

// ----- SSres -----------------------------------------------------------------
Action_DSSP2::SSres::SSres() :
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
  bridge1idx_(-1),
  b1type_(NO_BRIDGE),
  bridge2idx_(-1),
  b2type_(NO_BRIDGE),
  resChar_(' '),
  isSelected_(false)
{
  std::fill(SScount_, SScount_ + NSSTYPE_, 0);
  std::fill(Bcount_, Bcount_ + NBRIDGETYPE_, 0);
  std::fill(turnChar_, turnChar_ + NTURNTYPE_, ' ');
}

Action_DSSP2::SSres::SSres(SSres const& rhs) :
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
  bridge1idx_(rhs.bridge1idx_),
  b1type_(rhs.b1type_),
  bridge2idx_(rhs.bridge2idx_),
  b2type_(rhs.b2type_),
  resChar_(rhs.resChar_),
  isSelected_(rhs.isSelected_)
{
  std::copy(rhs.SScount_, rhs.SScount_ + NSSTYPE_, SScount_);
  std::copy(rhs.Bcount_, rhs.Bcount_ + NBRIDGETYPE_, Bcount_);
  std::copy(rhs.turnChar_, rhs.turnChar_ + NTURNTYPE_, turnChar_);
}
  
Action_DSSP2::SSres& Action_DSSP2::SSres::operator=(SSres const& rhs) {
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
  bridge1idx_ = rhs.bridge1idx_;
  b1type_ = rhs.b1type_;
  bridge2idx_ = rhs.bridge2idx_;
  b2type_ = rhs.b2type_;
  resChar_ = rhs.resChar_;
  isSelected_ = rhs.isSelected_;
  std::copy(rhs.SScount_, rhs.SScount_ + NSSTYPE_, SScount_);
  std::copy(rhs.Bcount_, rhs.Bcount_ + NBRIDGETYPE_, Bcount_);
  std::copy(rhs.turnChar_, rhs.turnChar_ + NTURNTYPE_, turnChar_);
  return *this;
}


void Action_DSSP2::SSres::Deselect() {
  isSelected_ = false;
  C_ = -1;
  H_ = -1;
  N_ = -1;
  O_ = -1;
  CA_ = -1;
}

void Action_DSSP2::SSres::Unassign() {
  sstype_ = NONE;
  bridge1idx_ = -1;
  b1type_ = NO_BRIDGE;
  bridge2idx_ = -1;
  b2type_ = NO_BRIDGE;
  std::fill(turnChar_, turnChar_ + NTURNTYPE_, ' ');
}

/** Accumulate SS data. */
void Action_DSSP2::SSres::AccumulateData(int frameNum, bool useString) {
  SScount_[sstype_]++;
  if (sstype_ == EXTENDED || sstype_ == BRIDGE) {
    if (b1type_ == ANTIPARALLEL || b2type_ == ANTIPARALLEL)
      Bcount_[ANTIPARALLEL]++;
    if (b1type_ == PARALLEL || b2type_ == PARALLEL)
      Bcount_[PARALLEL]++;
    // TODO bulge?
  } else
    Bcount_[NO_BRIDGE]++;
  if (useString)
    resDataSet_->Add(frameNum, SSchar_[sstype_]); 
  else
    resDataSet_->Add(frameNum, &sstype_);
}

/** Set turn beginning. */
void Action_DSSP2::SSres::SetTurnBegin(TurnType typeIn) {
  if (turnChar_[typeIn] == '<')
    turnChar_[typeIn] = 'X';
  else
    turnChar_[typeIn] = '>';
}

void Action_DSSP2::SSres::SetTurnEnd(TurnType typeIn) {
  if (turnChar_[typeIn] == '>')
    turnChar_[typeIn] = 'X';
  else
    turnChar_[typeIn] = '<';
}

void Action_DSSP2::SSres::SetTurn(TurnType typeIn) {
  // Do not overwrite an existing end character
  if (turnChar_[typeIn] == ' ') {
    if      (typeIn == T3) turnChar_[typeIn] = '3';
    else if (typeIn == T4) turnChar_[typeIn] = '4';
    else if (typeIn == T5) turnChar_[typeIn] = '5';
  }
}

int Action_DSSP2::SSres::ssPriority(SStype typeIn) {
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

int Action_DSSP2::SSres::SSpriority() const {
  return ssPriority(sstype_);
}

void Action_DSSP2::SSres::SetSS(SStype typeIn) {
  // TODO check if the priority check is necessary
  //if (ssPriority(typeIn) > ssPriority(sstype_))
    sstype_ = typeIn;
}

bool Action_DSSP2::SSres::HasTurnStart(TurnType typeIn) const {
  if (turnChar_[typeIn] == '>' ||
      turnChar_[typeIn] == 'X')
    return true;
  return false;
}

void Action_DSSP2::SSres::SetBridge(int idx, BridgeType btypeIn) {
  if (bridge1idx_ == -1) {
    bridge1idx_ = idx;
    b1type_ = btypeIn;
  } else if (bridge2idx_ == -1) {
    bridge2idx_ = idx;
    b2type_ = btypeIn;
  } else
    mprinterr("Error: Too many bridges for %i (to %i)\n", Num(), idx+1);
}

bool Action_DSSP2::SSres::HasBridge() const {
  if (bridge1idx_ != -1) return true;
  return false;
}

bool Action_DSSP2::SSres::IsBridgedWith(int idx2) const {
  if (bridge1idx_ == idx2) return true;
  if (bridge2idx_ == idx2) return true;
  return false;
}

/*char Action_DSSP2::SSres::StrandChar() const {
  // TODO ever b2?
  return ssChar_[B1];
}*/

void Action_DSSP2::SSres::PrintSSchar() const {
  static const char btypeChar[] = { ' ', 'p', 'A' };
  mprintf("\t%6i %c %c %c %c %c(%6i) %c(%6i) %6i %6i %6i %c\n", num_+1, resChar_,
          turnChar_[T3], turnChar_[T4], turnChar_[T5],
          btypeChar[b1type_], bridge1idx_+1, btypeChar[b2type_], bridge2idx_+1,
          Bcount_[NO_BRIDGE], Bcount_[PARALLEL], Bcount_[ANTIPARALLEL],
          DSSP_char_[sstype_]);
}

// ----- Action_DSSP2 ----------------------------------------------------------

Action_DSSP2::Action_DSSP2() :
  debug_(0),
  BB_N_("N"),
  BB_H_("H"),
  BB_C_("C"),
  BB_O_("O"),
  BB_CA_("CA"),
  outfile_(0),
  dsspFile_(0),
  assignout_(0)
{}

// Action_DSSP2::Help()
void Action_DSSP2::Help() const {
  mprintf("\t[<name>] [out <filename>] [<mask>] [sumout <filename>]\n"
          "\t[assignout <filename>] [totalout <filename>] [ptrajformat]\n"
          "\t[namen <N name>] [nameh <H name>] [nameca <CA name>]\n"
          "\t[namec <C name>] [nameo <O name>]\n"
          "  Calculate secondary structure content for residues in <mask>.\n"
          "  If sumout not specified, the filename specified by out is used with .sum suffix.\n");
}

// Action_DSSP2::Init()
Action::RetType Action_DSSP2::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  //Nframe_ = 0;
  // Get keywords
  outfile_ = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  std::string temp = actionArgs.GetStringKey("sumout");
  if (temp.empty() && outfile_ != 0) 
    temp = outfile_->DataFilename().Full() + ".sum";
  dsspFile_ = init.DFL().AddDataFile( temp );
  DataFile* totalout = init.DFL().AddDataFile( actionArgs.GetStringKey("totalout"), actionArgs );
  assignout_ = init.DFL().AddCpptrajFile(actionArgs.GetStringKey("assignout"), "SS assignment");
  printString_ = actionArgs.hasKey("ptrajformat");
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
  // Get masks
  if (Mask_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;

  // Set up the DSSP data set
  dsetname_ = actionArgs.GetStringNext();
  // Set default name if none specified
  if (dsetname_.empty())
    dsetname_ = init.DSL().GenerateDefaultName("DSSP");
  // Set up Z labels
  if (outfile_ != 0)
    outfile_->ProcessArgs(SSzlabels_);
  // Create data sets for total fraction SS vs time.
  for (int i = 0; i < NSSTYPE_; i++) {
    totalDS_[i] = init.DSL().AddSet(DataSet::FLOAT, MetaData(dsetname_, SSname_[i]));
    if (totalDS_[i] == 0) {
      mprinterr("Error: Could not create DSSP total frac v time data set.\n");
      return Action::ERR;
    }
    // For now dont add 'None' so colors match up.
    if (totalout != 0 && i > 0) totalout->AddDataSet( totalDS_[i] );
  }

  mprintf( "    SECSTRUCT: Calculating secondary structure using mask [%s]\n",Mask_.MaskString());
  if (outfile_ != 0) 
    mprintf("\tDumping results to %s\n", outfile_->DataFilename().full());
  if (dsspFile_ != 0)
    mprintf("\tSum results to %s\n", dsspFile_->DataFilename().full());
  if (printString_) { 
    mprintf("\tSS data for each residue will be stored as a string.\n");
    for (int i = 0; i < NSSTYPE_; i++)
      mprintf("\t\t%s = %s\n", SSchar_[i], SSname_[i]);
  } else {
    mprintf("\tSS data for each residue will be stored as integers.\n");
    for (int i = 0; i < NSSTYPE_; i++)
      mprintf("\t\t%i = %s\n", i, SSname_[i]);
  }
  if (assignout_ != 0)
    mprintf("\tOverall assigned SS will be written to %s\n", assignout_->Filename().full());
  mprintf("\tBackbone Atom Names: N=[%s]  H=[%s]  C=[%s]  O=[%s]  CA=[%s]\n",
          *BB_N_, *BB_H_, *BB_C_, *BB_O_, *BB_CA_ );
  mprintf("# Citation: Kabsch, W.; Sander, C.; \"Dictionary of Protein Secondary Structure:\n"
          "#           Pattern Recognition of Hydrogen-Bonded and Geometrical Features.\"\n"
          "#           Biopolymers (1983), V.22, pp.2577-2637.\n" );
  init.DSL().SetDataSetsPending(true);
  Init_ = init;
  return Action::OK;
}

static inline void PrintAtom(Topology const& top, NameType const& name, int idx) {
  if (idx > -1)
    mprintf(" '%s'=%-12s", *name, top.AtomMaskName(idx/3).c_str());
  else
    mprintf(" '%s'=%-12s", *name, "NONE");
}

// Action_DSSP2::Setup()
/** Set up secondary structure calculation for all residues selected by the
  * mask. A residue is selected if at least one atom in the residue is
  * selected by the mask. The coordinate indices (i.e. atom # * 3) for
  * the C, O, N, H, and CA atoms are set up if those atoms are present.
  * The residue is only initialized if it has not been previously selected
  * and set up by a prior call to setup.
  */
Action::RetType Action_DSSP2::Setup(ActionSetup& setup)
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
      // Determine the previous and next residues
      int prevresnum = -1;
      int nextresnum = -1;
      for (int at = thisRes.FirstAtom(); at != thisRes.LastAtom(); at++) {
        if ( setup.Top()[at].Element() != Atom::HYDROGEN ) {
          for (Atom::bond_iterator ib = setup.Top()[at].bondbegin();
                                   ib != setup.Top()[at].bondend(); ++ib)
          {
            if ( setup.Top()[*ib].ResNum() < *ridx ) {
              if (prevresnum != -1)
                mprintf("Warning: Multiple previous residues for res %i\n", *ridx+1);
              else
                prevresnum = setup.Top()[*ib].ResNum();
            } else if ( setup.Top()[*ib].ResNum() > *ridx ) {
              if (nextresnum != -1)
                mprintf("Warning: Multiple next residues for res %i\n", *ridx+1);
              else
                nextresnum = setup.Top()[*ib].ResNum();
            }
          }
        }
      }
      mprintf("\t %8i < %8i < %8i\n", prevresnum+1, *ridx+1, nextresnum+1);
      // Here we assume that residues are sequential!
      if (prevresnum > -1) Res->SetPrevIdx( Res-Residues_.begin()-1 );
      if (nextresnum > -1) Res->SetNextIdx( Res-Residues_.begin()+1 );
    }
    // Determine if this residue is selected
    if (Mask_.AtomsInCharMask(thisRes.FirstAtom(), thisRes.LastAtom())) {
      Res->SetSelected( true );
      ++nResSelected;
      // Determine atom indices
      for (int at = thisRes.FirstAtom(); at != thisRes.LastAtom(); at++)
      {
        if      ( setup.Top()[at].Name() == BB_C_ )  Res->SetC( at*3 );
        else if ( setup.Top()[at].Name() == BB_O_ )  Res->SetO( at*3 );
        else if ( setup.Top()[at].Name() == BB_N_ )  Res->SetN( at*3 );
        else if ( setup.Top()[at].Name() == BB_H_ )  Res->SetH( at*3 );
        else if ( setup.Top()[at].Name() == BB_CA_ ) Res->SetCA( at*3 );
      }
      // Check if residue is missing atoms
      if (Res->IsMissingAtoms()) {
        mprintf("Warning: Res %s is missing atoms", setup.Top().TruncResNameNum( *ridx ).c_str());
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
    } // END residue is selected
  }
  mprintf("\t%u of %zu solute residues selected.\n", nResSelected, soluteRes.Size());

  // DEBUG - print each residue set up.
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
  

  return Action::OK;
}

void Action_DSSP2::AssignBridge(int idx1in, int idx2in, BridgeType btypeIn) {
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

  if (btypeIn == ANTIPARALLEL) {
    mprintf("\t\tAssignBridge %i to %i, Antiparallel\n", idx1+1, idx2+1);
    //bchar = 'A';
  } else {
    mprintf("\t\tAssignBridge %i to %i, Parallel\n", idx1+1, idx2+1);
    //bchar = 'p';
  }

  // Do not duplicate bridges
  if (Resi.IsBridgedWith(idx2)) {
    mprintf("\t\tAlready present.\n");
    return;
  }

  Resi.SetBridge( idx2, btypeIn );
  Resj.SetBridge( idx1, btypeIn );
}

/*
void Action_DSSP2::AssignBridge(int idx1in, int idx2in, BridgeType btypeIn, char& currentStrandChar) {
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

// TODO use Num()?
static inline int AbsResDelta(int idx1, int idx2) {
  int resDelta =  idx1 - idx2;
  if (resDelta < 0) resDelta = -resDelta;
  return resDelta;
}

static inline void SetMin(int& resGapSize, int& sres0, int& sres1, int Nres, int Pres)
{
  int itmp = AbsResDelta(Nres, Pres);
  if (itmp < resGapSize) {
    resGapSize = itmp;
    sres0 = std::min(Pres, Nres);
    sres1 = std::max(Pres, Nres);
  }
}


int Action_DSSP2::OverHbonds(int frameNum, ActionFrame& frm)
{
  Timer t_overhbonds;
  t_overhbonds.Start();
  // ----- Determine hydrogen bonding ------------
  typedef std::pair<int,int> HbondPairType;
  typedef std::set<HbondPairType> HbondMapType;
  /// Map resi (CO) to resj (NH) potential bridge hbonds
  HbondMapType CO_NH_bonds;

  Timer t_calchb;
  t_calchb.Start();
  int resi;
  int Nres = (int)Residues_.size();
#ifdef _OPENMP
#pragma omp parallel private(resi)
{
#pragma omp for
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
          // Need to consider adjacent residues since delta (2 res) turns possible 
          if (resi != resj) {
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
                CO_NH_bonds.insert( HbondPairType(resi, resj) );
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
#endif
  t_calchb.Stop();

  Timer t_assign;
  t_assign.Start();
  // ----- Do basic assignment -------------------
  for (HbondMapType::const_iterator hb0 = CO_NH_bonds.begin(); hb0 != CO_NH_bonds.end(); ++hb0)
  {
    int riidx = hb0->first;
    int rjidx = hb0->second;
    SSres const& Resi = Residues_[riidx];
    SSres const& Resj = Residues_[rjidx];
    mprintf("Res %8i %c -- %8i %c", Resi.Num()+1, Resi.ResChar(),
                                    Resj.Num()+1, Resj.ResChar()); // DBG
    // Spacing between residues i and j
    int resDelta = Resj.Num() - Resi.Num();
    mprintf("(%4i)\n", resDelta);
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
      hb = CO_NH_bonds.find( HbondPairType(hb0->second, hb0->first) );
      if (hb != CO_NH_bonds.end()) {
        mprintf("\t\t%i ANTI-PARALLELa with %i (%i)\n", hb0->first+1, hb0->second+1, hb->first+1);
        AssignBridge(hb0->first, hb0->second, ANTIPARALLEL);
      }
    }
    resDelta = AbsResDelta(hb0->first+1, hb0->second-1);
    if (resDelta > 2) {
      // Assume (i-1, j+1). Look for (j-1, i+1)
      hb = CO_NH_bonds.find( HbondPairType(hb0->second-2, hb0->first+2) );
      if (hb != CO_NH_bonds.end()) {
        mprintf("\t\t%i ANTI-PARALLELb with %i (%i)\n", hb0->first+2, hb0->second, hb->first+1);
        AssignBridge(hb0->first+1, hb0->second-1, ANTIPARALLEL);
      }
    }
    resDelta = AbsResDelta(hb0->first+1, hb0->second);
    if (resDelta > 2) {
      // Assume (i-1, j). Check for (j, i+1) PARALLEL
      hb = CO_NH_bonds.find( HbondPairType(hb0->second, hb0->first+2) );
      if (hb != CO_NH_bonds.end()) {
        mprintf("\t\t%i PARALLELa with %i (%i)\n", hb0->first+2, hb0->second+1, hb->first+1);
        AssignBridge(hb0->first+1, hb0->second, PARALLEL);
      }
    }
    resDelta = AbsResDelta(hb0->second-1, hb0->first);
    if (resDelta > 2) {
      // Assume (j, i+1). Check for (i-1, j)
      hb = CO_NH_bonds.find( HbondPairType(hb0->second-2, hb0->first) );
      if (hb != CO_NH_bonds.end()) {
        mprintf("\t\t%i PARALLELb with %i (%i)\n", hb0->second, hb0->first+1, hb->first+1);
        AssignBridge(hb0->second-1, hb0->first, PARALLEL);
      }
    }
  } // END loop over Hbonds

  // ----- Do SS assignment ----------------------
  // Priority is 'H', 'B', 'E', 'G', 'I', 'T', 'S' None
  //              8    7    6    5    4    3    2  1
  for (resi = 0; resi < Nres; resi++)
  {
    mprintf("Residue %i\n", resi+1);
    SSres& Resi = Residues_[resi];
    int prevRes = resi - 1;
    int nextRes = resi + 1;
    int priority = Resi.SSpriority();
    if ( Resi.HasTurnStart(T4) && prevRes > -1 && Residues_[prevRes].HasTurnStart(T4) )
    {
      // Alpha helix.
      mprintf("ALPHA helix starting at %i\n", resi+1);
      Residues_[resi  ].SetSS( ALPHA );
      Residues_[resi+1].SetSS( ALPHA );
      Residues_[resi+2].SetSS( ALPHA );
      Residues_[resi+3].SetSS( ALPHA );
    } else if (Resi.SS() != ALPHA) {
      if (priority < 6) {
        // Priority < 6 means not alpha or beta assigned yet.
        // Check for Beta structure
        bool prevHasBridge = (prevRes > -1   && Residues_[prevRes].HasBridge());
        bool nextHasBridge = (nextRes < Nres && Residues_[nextRes].HasBridge());
        if (Resi.HasBridge()) {
          // Regular Beta. Check if previous is assigned EXTENDED in case it 
          // was assigned via a Beta bulge.
          if ( prevHasBridge || nextHasBridge || Residues_[prevRes].SS() == EXTENDED )
          {
            mprintf("Extended BETA bridge at %i\n", resi+1);
            Resi.SetSS( EXTENDED );
          } else {
            mprintf("Isolated BETA bridge at %i.\n", resi+1);
            Resi.SetSS( BRIDGE );
          }
        } else if (prevHasBridge && nextHasBridge) {
          // Potential Beta bulge. Check other strand.
          int presb1 = Residues_[prevRes].Bridge1Idx();
          int presb2 = Residues_[prevRes].Bridge2Idx();
          int nresb1 = Residues_[nextRes].Bridge1Idx();
          int nresb2 = Residues_[nextRes].Bridge2Idx();
          mprintf("Potential bulge? Prev res bridge res: %i %i  Next res bridge res: %i %i\n", presb1, presb2, nresb1, nresb2);
          // The largest allowed gap in the other strand is 5 residues.
          // Since we know that next res and previous res both have at 
          // least 1 bridge, Next B1 - Prev B1 can be the gap to beat.
          int resGapSize = AbsResDelta(nresb1, presb1);
          int sres0 = std::min(presb1, nresb1);
          int sres1 = std::max(presb1, nresb1);
          if (presb2 != -1)
            SetMin( resGapSize, sres0, sres1, nresb1, presb2 );
          if (nresb2 != -1) {
            SetMin( resGapSize, sres0, sres1, nresb2, presb1 ); 
            if (presb2 != -1)
              SetMin( resGapSize, sres0, sres1, nresb2, presb2 );
          }
          mprintf("Min res gap size on other strand = %i (%i to %i)\n", resGapSize, sres0, sres1);
          // Minimum allowed gap is 4 residues in between, so 5 residues total.
          if (resGapSize < 6) {
            mprintf("Beta bulge.\n");
            if (Residues_[prevRes].SS() != ALPHA)
              Residues_[prevRes].SetSS( EXTENDED );
            Resi.SetSS( EXTENDED );
            if (Residues_[nextRes].SS() != ALPHA)
              Residues_[nextRes].SetSS( EXTENDED );
            // Set extended on other strand as well
            for (int sres = sres0; sres != sres1; sres++)
              if (Residues_[sres].SS() != ALPHA)
                Residues_[sres].SetSS( EXTENDED );
          }
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
      mprintf("3-10 helix starting at %i\n", resi+1);
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
      mprintf("PI helix starting at %i\n", resi+1);
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
    if (Residues_[resi].IsSelected()) {
      Residues_[resi].AccumulateData(frameNum, printString_);
      Nselected++;
    }
  }
  for (int i = 0; i < NSSTYPE_; i++) {
    float fvar = (float)totalSS[i];
    fvar /= (float)Nselected;
    totalDS_[i]->Add(frameNum, &fvar);
  }


  t_assign.Stop();

  t_overhbonds.Stop();
  t_calchb.WriteTiming(1, "Calc Hbonds", t_overhbonds.Total());
  t_assign.WriteTiming(1, "Assignment ", t_overhbonds.Total());
  t_overhbonds.WriteTiming(0,"Over Hbonds");
  return 0;
}

Action::RetType Action_DSSP2::DoAction(int frameNum, ActionFrame& frm)
{
  OverHbonds(frameNum, frm);
  // DEBUG - Print basic assignment
  for (SSarrayType::const_iterator it = Residues_.begin(); it != Residues_.end(); ++it)
    it->PrintSSchar();
  return Action::OK;
}

// Action_DSSP::Print()
void Action_DSSP2::Print() {
  if (dsetname_.empty()) return;
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
    for (int resi = min_res; resi < max_res+1; resi++) {
      if (Residues_[resi].Dset() != 0) {
        int Nframe = 0;
        for (int ss = 0; ss < NSSTYPE_; ss++)
          Nframe += Residues_[resi].SScount((SStype)ss);
        for (int ss = 1; ss < NSSTYPE_; ss++) {
          double avg = (double)Residues_[resi].SScount((SStype)ss);
          avg /= (double)Nframe;
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
        resLine += Residues_[resi].ResChar();
        // Figure out which SS element is dominant for res if selected
        if (Residues_[resi].Dset() != 0) {
          int dominantType = 0;
          int ssmax = 0;
          for (int ss = 0; ss < NSSTYPE_; ss++) {
            if ( Residues_[resi].SScount((SStype)ss) > ssmax ) {
              ssmax = Residues_[resi].SScount((SStype)ss);
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
