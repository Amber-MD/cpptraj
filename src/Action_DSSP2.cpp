#include <algorithm> // std::fill
#include "Action_DSSP2.h"
#include "CpptrajStdio.h"

// ----- SSres -----------------------------------------------------------------
Action_DSSP2::SSres::SSres() :
  resDataSet_(0),
  chirality_(0),
  pattern_(NOHBOND),
  sstype_(NONE),
  idx_(-1),
  C_(-1),
  O_(-1),
  N_(-1),
  H_(-1),
  CA_(-1),
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
    }
    Residue const& thisRes = setup.Top().Res( *ridx );
    // Determine if this residue is selected
    if (Mask_.AtomsInCharMask(thisRes.FirstAtom(), thisRes.LastAtom())) {
      Res->SetSelected( true );
      Res->SetC(  setup.Top().FindAtomInResidue(*ridx, BB_C_)  );
      Res->SetO(  setup.Top().FindAtomInResidue(*ridx, BB_O_)  );
      Res->SetN(  setup.Top().FindAtomInResidue(*ridx, BB_N_)  );
      Res->SetH(  setup.Top().FindAtomInResidue(*ridx, BB_H_)  );
      Res->SetCA( setup.Top().FindAtomInResidue(*ridx, BB_CA_) );
    }
  }

  return Action::OK;
}

// Action_DSSP2::DoAction()
Action::RetType Action_DSSP2::DoAction(int frameNum, ActionFrame& frm)
{

}
