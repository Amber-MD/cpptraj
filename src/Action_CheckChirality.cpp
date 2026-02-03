#include "Action_CheckChirality.h"
#include "CpptrajStdio.h"
#include "TorsionRoutines.h"
#include "DataSet_Mesh.h"
//#incl ude "Constants.h"
#include "Structure/Chirality.h"

void Action_CheckChirality::Help() const {
  //mprintf("\t[<name>] [<mask1>] [out <filename>] [byatom]\n"
  mprintf("\t[<name>] [<mask1>] [out <filename>]\n"
          "  Check the chirality of amino acid residues in <mask1>.\n");
//          "  If 'byatom' is specified, check the chirality of any atoms in\n"
//          "  <mask1> that are chiral centers.\n");
}

/** CONSTRUCTOR */
Action_CheckChirality::Action_CheckChirality() :
  data_L_(0),
  data_D_(0),
  outfile_(0),
  byatom_(false),
  currentTop_(0),
  debug_(0)
{}

// Action_CheckChirality::Init()
Action::RetType Action_CheckChirality::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  // Get keywords
  outfile_ = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  //byatom_ = actionArgs.hasKey("byatom");i //FIXME - not ready for primetime, needs more testing
  byatom_ = false;
  // Get Masks
  if (Mask1_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;
  // Set up DataSets
  setname_ = actionArgs.GetStringNext();
  if (setname_.empty())
    setname_ = init.DSL().GenerateDefaultName("CHIRAL");
  if (!byatom_) {
    MetaData md(setname_, "L", MetaData::NOT_TS);
    data_L_ = init.DSL().AddSet( DataSet::XYMESH, md);
    md.SetAspect("D");
    data_D_ = init.DSL().AddSet( DataSet::XYMESH, md); 
    if (data_L_ == 0 || data_D_ == 0) return Action::ERR;
    data_L_->SetupFormat().SetFormatWidthPrecision( 8, 0 );
    data_D_->SetupFormat().SetFormatWidthPrecision( 8, 0 );
    if (outfile_ != 0) {
      outfile_->AddDataSet( data_L_ );
      outfile_->AddDataSet( data_D_ );
    }
#   ifdef MPI
    data_L_->SetNeedsSync( false );
    data_D_->SetNeedsSync( false );
#   endif
  }
  if (!byatom_) {
    mprintf("    CHECKCHIRALITY: Check chirality for AA residues in mask '%s'\n",
            Mask1_.MaskString());
  } else {
    mprintf("    CHECKCHIRALITY: Check chirality for atoms in mask '%s'\n",
            Mask1_.MaskString());
  }
  if (outfile_ != 0)
    mprintf("\tOutput to file %s\n", outfile_->DataFilename().full());
  if (!setname_.empty())
    mprintf("\tData set name: %s\n", setname_.c_str());
  Init_ = init;

  return Action::OK;
}

// Action_CheckChirality::Setup()
/** Set angle up for this parmtop. Get masks etc. */
Action::RetType Action_CheckChirality::Setup(ActionSetup& setup) {
//  top_ = setup.TopAddress(); // DEBUG
  if (setup.Top().SetupCharMask(Mask1_)) return Action::ERR;
  if (Mask1_.None()) {
    mprinterr("Warning: Mask '%s' selects no atoms.\n", Mask1_.MaskString());
    return Action::SKIP;
  }
  if (byatom_) {
    currentTop_ = setup.TopAddress();
    if (atomChiralitySets_.empty()) {
      // First time
      atomChiralitySets_.assign( setup.Top().Natom(), 0 );
    } else if ((unsigned int)setup.Top().Natom() > atomChiralitySets_.size()) {
      atomChiralitySets_.resize( setup.Top().Natom(), 0 );
    }
    // For each selected atom, create a set
    for (int at = 0; at != setup.Top().Natom(); ++at) {
      if (Mask1_.AtomInCharMask(at)) {
        std::string legend = setup.Top().AtomMaskName(at);
        if ( atomChiralitySets_[at] == 0 ) {
          if (setup.Top()[at].Nbonds() < 3) {
            mprintf("Warning: Atom %s has < 3 bonds, not checking for chirality.\n", legend.c_str());
          } else {
            // Initial set up
            //atomChiralitySets_[at] = Init_.DSL().AddSet( DataSet::FLOAT, MetaData(setname_, at+1) );
            atomChiralitySets_[at] = Init_.DSL().AddSet( DataSet::STRING, MetaData(setname_, at+1) );
            if (atomChiralitySets_[at] == 0) {
              mprinterr("Error: Could not allocate set for atom '%s'\n", legend.c_str());
              return Action::ERR;
            }
            atomChiralitySets_[at]->SetLegend( legend );
            if (outfile_ != 0) outfile_->AddDataSet( atomChiralitySets_[at] );
          }
        } else {
          // Check an already set up dataset
          if (legend != atomChiralitySets_[at]->Meta().Legend()) {
            mprintf("Warning: After topology change, atom %i has changed from %s to %s\n",
                    at+1, atomChiralitySets_[at]->legend(), legend.c_str());
          }
        }
      } // END atom is selected
    } // END loop over atoms
    mprintf("\tChecking chirality for %i atoms.\n", Mask1_.Nselected());
  } else {
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
  } // END residue setup
  return Action::OK;
}

// Action_CheckChirality::DoAction()
Action::RetType Action_CheckChirality::DoAction(int frameNum, ActionFrame& frm) {
  if (byatom_) {
    for (std::vector<DataSet*>::const_iterator it = atomChiralitySets_.begin();
                                               it != atomChiralitySets_.end(); ++it)
    {
      if (*it != 0) {
        DataSet* ds = *it;
        int iat = ds->Meta().Idx() - 1;
        Cpptraj::Structure::ChiralType chirality =
          Cpptraj::Structure::DetermineChirality(iat, *currentTop_, frm.Frm(), debug_);
        if (chirality == Cpptraj::Structure::CHIRALITY_ERR)
          return Action::ERR;
        //double tval_degrees = AtomC.TorsionVal() * Constants::RADDEG;
        //float ftval = (float)tval_degrees;
        //ds->Add( frameNum, &ftval );
        const char* chiralstr = Cpptraj::Structure::chiralStr( chirality );
        ds->Add( frameNum, chiralstr );
      }
    }
  } else {
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
//    DEBUG
//    int at = ri->ca_ / 3;
//    Cpptraj::Structure::ChiralType chirality = Cpptraj::Structure::IS_UNKNOWN_CHIRALITY;
//    if (top_->Atoms()[at].Nbonds() > 2)
//      chirality = Cpptraj::Structure::DetermineChirality(at, *top_, frm.Frm(), 0);
//    mprintf("Chirality around %s is %s\n", top_->AtomMaskName(at).c_str(), Cpptraj::Structure::chiralStr(chirality));
//    DEBUG
    }
  }

  return Action::OK;
}

#ifdef MPI
int Action_CheckChirality::SyncAction() {
  if (byatom_) return 0;
  int total_L, total_D;
  for (Rarray::iterator ri = resInfo_.begin(); ri != resInfo_.end(); ++ri) {
    Init_.TrajComm().ReduceMaster( &total_L, &(ri->N_L_), 1, MPI_INT, MPI_SUM );
    Init_.TrajComm().ReduceMaster( &total_D, &(ri->N_D_), 1, MPI_INT, MPI_SUM );
    if (Init_.TrajComm().Master()) {
      ri->N_L_ = total_L;
      ri->N_D_ = total_D;
    }
  }
  return 0;
}
#endif

void Action_CheckChirality::Print() {
  if (byatom_) return;
  data_L_->ModifyDim(Dimension::X).SetLabel("Res");
  data_D_->ModifyDim(Dimension::X).SetLabel("Res");
  DataSet_Mesh& dsetL = static_cast<DataSet_Mesh&>( *data_L_ );
  DataSet_Mesh& dsetD = static_cast<DataSet_Mesh&>( *data_D_ );
  for (Rarray::const_iterator ri = resInfo_.begin(); ri != resInfo_.end(); ++ri) {
    dsetL.AddXY( ri->num_+1, ri->N_L_ );
    dsetD.AddXY( ri->num_+1, ri->N_D_ );
  }
}
