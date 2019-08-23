#include <cmath> // sqrt
#include <algorithm> // std::fill
#include "Action_DSSP2.h"
#include "CpptrajStdio.h"
#include "DataSet.h"
#include "DistRoutines.h"

// ----- SSres -----------------------------------------------------------------
Action_DSSP2::SSres::SSres() :
  resDataSet_(0),
  chirality_(0),
//  pattern_(NOHBOND),
  sstype_(NONE),
  idx_(-1),
  C_(-1),
  O_(-1),
  N_(-1),
  H_(-1),
  CA_(-1),
  bridge1idx_(-1),
  bridge2idx_(-1),
  resChar_(' '),
  TURN3_(' '),
  TURN4_(' '),
  TURN5_(' '),
  BRIDGE1_(' '),
  BRIDGE2_(' '),
  SHEET_(' '),
  isSelected_(false)
{
  std::fill(SScount_, SScount_ + NSSTYPE_, 0); 
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
  CO_HN_Hbonds_.clear();
//  pattern_ = NOHBOND;
  sstype_ = NONE;
  bridge1idx_ = -1;
  bridge2idx_ = -1;
  TURN3_ = ' ';
  TURN4_ = ' ';
  TURN5_ = ' ';
  BRIDGE1_ = ' ';
  BRIDGE2_ = ' ';
  SHEET_ = ' ';
}

// ----- Action_DSSP2 ----------------------------------------------------------

/** From the original Kabsch & Sander 1983 paper. Obtained via
  *   q1 = 0.42e
  *   q2 = 0.20e
  *   f = 332 (approximate conversion factor to get energies in kcal/mol)
  *   fac = q1*q2*f
  */
const double Action_DSSP2::DSSP_fac_ = 27.888;

const double Action_DSSP2::DSSP_cut_ = -0.5;

const char  Action_DSSP2::DSSP_char_[] = { ' ', 'E', 'B', 'G', 'H', 'I', 'T', 'S' };

const char* Action_DSSP2::SSchar_[]    = { "0", "b", "B", "G", "H", "I", "T", "S" };

const char* Action_DSSP2::SSname_[]={"None", "Extended", "Bridge", "3-10", "Alpha", "Pi", "Turn", "Bend"};

const std::string Action_DSSP2::SSzlabels_ = "zlabels None,Para,Anti,3-10,Alpha,Pi,Turn,Bend";

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
    if (Res->Idx() != -1) {
      // Residue has previously been set up. Check that indices match.
      if (Res->Idx() != *ridx) {
        mprinterr("Error: Solute residue index %i does not match previously setup\n"
                  "Error: index %i\n", *ridx, Res->Idx());
        return Action::ERR;
      }
    } else {
      // Set up Residue. TODO also molecule index?
      Res->SetIdx( *ridx );
      Res->SetResChar( setup.Top().Res(*ridx).SingleCharName() );
    }
    Residue const& thisRes = setup.Top().Res( *ridx );
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
    mprintf("    %8i", it->Idx() + 1);
    PrintAtom(setup.Top(), BB_C_, it->C());
    PrintAtom(setup.Top(), BB_O_, it->O());
    PrintAtom(setup.Top(), BB_N_, it->N());
    PrintAtom(setup.Top(), BB_H_, it->H());
    PrintAtom(setup.Top(), BB_CA_, it->CA());
    mprintf("\n");
  }
  

  return Action::OK;
}

// Action_DSSP2::DoAction()
Action::RetType Action_DSSP2::DoAction(int frameNum, ActionFrame& frm)
{
  int resi;
  int Nres = (int)Residues_.size();
  // The first pass is used to determine hydrogen bonding
#ifdef _OPENMP
#pragma omp parallel private(resi)
{
#pragma omp for
#endif /* _OPENMP */
  for (resi = 0; resi < Nres; resi++)
  {
    if (Residues_[resi].IsSelected())
    {
      SSres& ResCO = Residues_[resi];
      ResCO.Unassign();
      if (ResCO.HasCO())
      {
        const double* Cxyz = frm.Frm().CRD( ResCO.C() );
        const double* Oxyz = frm.Frm().CRD( ResCO.O() );
        for (int resj = 0; resj < Nres; resj++)
        {
          // Only consider residues spaced more than 2 apart.
          int resDelta = resi - resj;
          if (resDelta < 0) resDelta = -resDelta;
          if (resDelta > 2) {
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
                ResCO.AddHbond( resj );
              }
            } // END ResNH selected
          } // END residues spaced > 2 apart
        } // END inner loop over residues
      } // END has CO
    } // END ResCO selected
  } // END outer loop over residues
#ifdef _OPENMP
} // END pragma omp parallel
#endif
  // Do basic assignment
  for (SSarrayType::const_iterator Resi = Residues_.begin(); Resi != Residues_.end(); ++Resi)
  {
    mprintf("Res %8i %c", Resi->Idx()+1, Resi->ResChar()); // DBG
    for (SSres::const_iterator rjidx = Resi->begin(); rjidx != Resi->end(); ++rjidx)
    {
      SSres& Resj = Residues_[*rjidx];
      mprintf(" %8i", Resj.Idx()+1); // DBG
      int resDelta = Resj.Idx() - Resi->Idx();
      if (resDelta < 0) resDelta = -resDelta;
      mprintf("(%4i)", resDelta);
    }
    mprintf("\n"); // DBG
  }
  return Action::OK;
}
