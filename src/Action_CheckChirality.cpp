#include "Action_CheckChirality.h"
#include "CpptrajStdio.h"
#include "TorsionRoutines.h"
#include "DataSet_Mesh.h"

void Action_CheckChirality::Help() const {
  mprintf("\t[<name>] [<mask1>] [out <filename>]\n"
          "  Check the chirality of AA residues in <mask1>.\n");
}

// Action_CheckChirality::Init()
Action::RetType Action_CheckChirality::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get keywords
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  // Get Masks
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );
  // Set up DataSets
  setname_ = actionArgs.GetStringNext();
  if (setname_.empty())
    setname_ = init.DSL().GenerateDefaultName("CHIRAL");
  MetaData md(setname_, "L", MetaData::NOT_TS);
  data_L_ = init.DSL().AddSet( DataSet::XYMESH, md);
  md.SetAspect("D");
  data_D_ = init.DSL().AddSet( DataSet::XYMESH, md); 
  if (data_L_ == 0 || data_D_ == 0) return Action::ERR;
  data_L_->SetupFormat().SetFormatWidthPrecision( 8, 0 );
  data_D_->SetupFormat().SetFormatWidthPrecision( 8, 0 );
  if (outfile != 0) {
    outfile->AddDataSet( data_L_ );
    outfile->AddDataSet( data_D_ );
  }
# ifdef MPI
  data_L_->SetNeedsSync( false );
  data_D_->SetNeedsSync( false );
# endif
  mprintf("    CHECKCHIRALITY: Check chirality for AA residues in mask '%s'\n",
          Mask1_.MaskString());
  if (outfile != 0)
    mprintf("\tOutput to file %s\n", outfile->DataFilename().full());
  if (!setname_.empty())
    mprintf("\tData set name: %s\n", setname_.c_str());
  Init_ = init;

  return Action::OK;
}

// Action_CheckChirality::Setup()
/** Set angle up for this parmtop. Get masks etc. */
Action::RetType Action_CheckChirality::Setup(ActionSetup& setup) {
  if (setup.Top().SetupCharMask(Mask1_)) return Action::ERR;
  if (Mask1_.None()) {
    mprinterr("Warning: Mask '%s' selects no atoms.\n", Mask1_.MaskString());
    return Action::SKIP;
  }
  // Reset any existing ResidueInfos to inactive
  for (Rarray::iterator ri = resInfo_.begin(); ri != resInfo_.end(); ++ri)
    ri->isActive_ = 0;
  int Nactive = 0;
  // Loop over all non-solvent residues
  int resnum = 0;
  std::vector<std::string> NotFound;
  for (Topology::res_iterator res = setup.Top().ResStart();
                              res != setup.Top().ResEnd(); ++res, ++resnum)
  {
    int firstAtom = res->FirstAtom();
    int molnum = setup.Top()[firstAtom].MolNum();
    // Skip solvent
    if (!setup.Top().Mol(molnum).IsSolvent())
    {
      if (Mask1_.AtomsInCharMask( firstAtom, res->LastAtom() ))
      { 
        int n_atom = setup.Top().FindAtomInResidue(resnum, "N");
        int ca_atom = setup.Top().FindAtomInResidue(resnum, "CA");
        int c_atom = setup.Top().FindAtomInResidue(resnum, "C");
        int cb_atom = setup.Top().FindAtomInResidue(resnum, "CB");
        if (n_atom == -1 || ca_atom == -1 || c_atom == -1 || cb_atom == -1)
        {
          NotFound.push_back( setup.Top().TruncResNameNum(resnum) );
          continue;
        }
        Nactive++;
        // See if a data set is already present
        Rarray::iterator rinfo = resInfo_.end();
        for (Rarray::iterator ri = resInfo_.begin(); ri != resInfo_.end(); ++ri)
          if (resnum == ri->num_) {
            rinfo = ri;
            break;
          }
        if (rinfo != resInfo_.end()) {
          // Update coord indices in ResidueInfo
          rinfo->isActive_ = 1;
          rinfo->n_ = n_atom * 3;
          rinfo->ca_ = ca_atom * 3;
          rinfo->c_ = c_atom * 3;
          rinfo->cb_ = cb_atom * 3;
        } else {
          // New ResidueInfo
          ResidueInfo RI;
          RI.num_ = resnum;
          RI.isActive_ = 1;
          RI.n_ = n_atom * 3;
          RI.ca_ = ca_atom * 3;
          RI.c_ = c_atom * 3;
          RI.cb_ = cb_atom * 3;
          RI.N_L_ = 0;
          RI.N_D_ = 0;
          resInfo_.push_back( RI );
        }
      } // END atoms in mask
    } // END not solvent
  } // END loop over residues

  if (Nactive == 0) {
    mprintf("Warning: No valid residues selected from '%s'\n", setup.Top().c_str());
    return Action::SKIP;
  }
  mprintf("\tChecking chirality for %i residues\n", Nactive);
  if (!NotFound.empty()) {
    mprintf("\tSome atoms not found for %zu residues (this is expected for e.g. GLY)\n\t",
            NotFound.size());
    for (std::vector<std::string>::const_iterator rn = NotFound.begin();
                                                  rn != NotFound.end(); ++rn)
      mprintf(" %s", rn->c_str());
    mprintf("\n");
  }
  return Action::OK;
}

// Action_CheckChirality::DoAction()
Action::RetType Action_CheckChirality::DoAction(int frameNum, ActionFrame& frm) {
  for (Rarray::iterator ri = resInfo_.begin(); ri != resInfo_.end(); ++ri)
  {
    double torsion = Torsion( frm.Frm().CRD(ri->n_),
                              frm.Frm().CRD(ri->ca_),
                              frm.Frm().CRD(ri->c_),
                              frm.Frm().CRD(ri->cb_) );
    if (torsion < 0.0)
      ri->N_L_++;
    else
      ri->N_D_++;
  }

  return Action::OK;
}

#ifdef MPI
int Action_CheckChirality::SyncAction() {
  int total_L, total_D;
  for (Rarray::iterator ri = resInfo_.begin(); ri != resInfo_.end(); ++ri) {
    Init_.TrajComm().Reduce( &total_L, &(ri->N_L_), 1, MPI_INT, MPI_SUM );
    Init_.TrajComm().Reduce( &total_D, &(ri->N_D_), 1, MPI_INT, MPI_SUM );
    if (Init_.TrajComm().Master()) {
      ri->N_L_ = total_L;
      ri->N_D_ = total_D;
    }
  }
  return 0;
}
#endif

void Action_CheckChirality::Print() {
  data_L_->ModifyDim(Dimension::X).SetLabel("Res");
  data_D_->ModifyDim(Dimension::X).SetLabel("Res");
  DataSet_Mesh& dsetL = static_cast<DataSet_Mesh&>( *data_L_ );
  DataSet_Mesh& dsetD = static_cast<DataSet_Mesh&>( *data_D_ );
  for (Rarray::const_iterator ri = resInfo_.begin(); ri != resInfo_.end(); ++ri) {
    dsetL.AddXY( ri->num_+1, ri->N_L_ );
    dsetD.AddXY( ri->num_+1, ri->N_D_ );
  }
}
