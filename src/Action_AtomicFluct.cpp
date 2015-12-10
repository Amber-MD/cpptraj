#include <cmath> // sqrt
#include "Action_AtomicFluct.h"
#include "CpptrajStdio.h"
#include "Constants.h" // PI
#include "DataSet_Mesh.h"
#include "PDBfile.h"

// CONSTRUCTOR
Action_AtomicFluct::Action_AtomicFluct() :
  sets_(0),
  bfactor_(false),
  calc_adp_(false),
  adpoutfile_(0),
  fluctParm_(0),
  outtype_(BYATOM),
  dataout_(0)
{}

void Action_AtomicFluct::Help() const {
  mprintf("\t[out <filename>] [<mask>] [byres | byatom | bymask] [bfactor]\n"
          "\t[calcadp [adpout <file>]]\n"
          "\t%s\n"
          "  Calculate atomic fluctuations of atoms in <mask>\n", ActionFrameCounter::HelpText);
}

// Action_AtomicFluct::Init()
Action::RetType Action_AtomicFluct::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get frame # keywords
  if (InitFrameCounter(actionArgs)) return Action::ERR;
  // Get other keywords
  bfactor_ = actionArgs.hasKey("bfactor");
  calc_adp_ = actionArgs.hasKey("calcadp");
  adpoutfile_ = init.DFL().AddCpptrajFile(actionArgs.GetStringKey("adpout"), "PDB w/ADP",
                                    DataFileList::PDB);;
  if (adpoutfile_!=0) calc_adp_ = true; // adpout implies calcadp
  if (calc_adp_ && !bfactor_) bfactor_ = true;
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs ); 
  if (actionArgs.hasKey("byres"))
    outtype_ = BYRES;
  else if (actionArgs.hasKey("bymask"))
    outtype_ = BYMASK;
  else if (actionArgs.hasKey("byatom") || actionArgs.hasKey("byatm"))
    outtype_ = BYATOM;
  // Get Mask
  Mask_.SetMaskString( actionArgs.GetMaskNext()  );
  // Get DataSet name
  std::string setname = actionArgs.GetStringNext();
  // Add output dataset
  MetaData md( setname );
  md.SetTimeSeries( MetaData::NOT_TS );
  if (bfactor_)
    md.SetLegend("B-factors");
  else
    md.SetLegend("AtomicFlx");
  dataout_ = init.DSL().AddSet( DataSet::XYMESH, md, "Fluct" );
  if (dataout_ == 0) {
    mprinterr("Error: AtomicFluct: Could not allocate dataset for output.\n");
    return Action::ERR; 
  }
  if (outfile != 0) 
    outfile->AddDataSet( dataout_ );

  mprintf("    ATOMICFLUCT: calculating");
  if (bfactor_)
    mprintf(" B factors");
  else
    mprintf(" atomic positional fluctuations");
  if (outfile != 0)
    mprintf(", output to file %s", outfile->DataFilename().full());
  mprintf("\n                 Atom mask: [%s]\n",Mask_.MaskString());
  FrameCounterInfo();
  if (calc_adp_) {
    mprintf("\tCalculating anisotropic displacement parameters.\n");
    if (adpoutfile_!=0) mprintf("\tWriting PDB with ADP to '%s'\n", adpoutfile_->Filename().full());
  }
  if (!setname.empty())
    mprintf("\tData will be saved to set named %s\n", setname.c_str());

  return Action::OK;
}

// Action_AtomicFluct::Setup()
Action::RetType Action_AtomicFluct::Setup(ActionSetup& setup) {

  if (SumCoords_.Natom()==0) {
    // Set up frames if not already set
    SumCoords_.SetupFrame(setup.Top().Natom());
    SumCoords2_.SetupFrame(setup.Top().Natom());
    SumCoords_.ZeroCoords();
    SumCoords2_.ZeroCoords();
    if (calc_adp_) {
      Cross_.SetupFrame(setup.Top().Natom());
      Cross_.ZeroCoords();
    }
    // This is the parm that will be used for this calc
    fluctParm_ = setup.TopAddress();
    // Set up atom mask
    if (setup.Top().SetupCharMask( Mask_ )) {
      mprinterr("Error: Could not set up mask [%s]\n",Mask_.MaskString());
      return Action::ERR;
    }
    Mask_.MaskInfo();
    if (Mask_.None()) {
      mprintf("Warning: No atoms selected [%s]\n",Mask_.MaskString());
      return Action::SKIP;
    }
  } else if (setup.Top().Natom() != SumCoords_.Natom()) {
    // Check that current #atoms matches
    mprinterr("Error: AtomicFluct not yet supported for mulitple topologies with different\n");
    mprinterr("       #s of atoms (set up for %i, this topology has %i\n",
              SumCoords_.Natom(), setup.Top().Natom());
    return Action::ERR;
  } 
  // NOTE: Print warning here when setting up multiple topologies?
  return Action::OK;
}

// Action_AtomicFluct::DoAction()
Action::RetType Action_AtomicFluct::DoAction(int frameNum, ActionFrame& frm) {
  if ( CheckFrameCounter( frameNum ) ) return Action::OK;
  SumCoords_ += frm.Frm();
  SumCoords2_ += ( frm.Frm() * frm.Frm() );
  if (calc_adp_) {
    for (int i = 0; i < SumCoords_.size(); i+=3) {
      Cross_[i  ] += frm.Frm()[i  ] * frm.Frm()[i+1]; // U12
      Cross_[i+1] += frm.Frm()[i  ] * frm.Frm()[i+2]; // U13
      Cross_[i+2] += frm.Frm()[i+1] * frm.Frm()[i+2]; // U23
    }
  }
  ++sets_;
  return Action::OK;
}

// Action_AtomicFluct::Print() 
void Action_AtomicFluct::Print() {
  mprintf("    ATOMICFLUCT: Calculating fluctuations for %i sets.\n",sets_);

  double Nsets = (double)sets_;
  // SumCoords will hold the average: <R>
  SumCoords_.Divide(Nsets);
  // SumCoords2 will hold the variance: <R^2> - <R>^2
  SumCoords2_.Divide(Nsets);
  SumCoords2_ = SumCoords2_ - (SumCoords_ * SumCoords_);
  // Cross terms: XY, XZ, YZ
  if (calc_adp_)
    Cross_.Divide(Nsets);

  // Hold fluctuation results - initialize to 0
  std::vector<double> Results( SumCoords2_.Natom(), 0 );
  std::vector<double>::iterator result = Results.begin();

  if (bfactor_) {
    // Set up b factor normalization
    // B-factors are (8/3)*PI*PI * <r>**2 hence we do not sqrt the fluctuations
    double bfac = (8.0/3.0)*Constants::PI*Constants::PI;
    for (int i = 0; i < SumCoords2_.size(); i+=3) {
      double fluct = SumCoords2_[i] + SumCoords2_[i+1] + SumCoords2_[i+2];
      if (fluct > 0) 
        *result = bfac * fluct;
      ++result;
      if (calc_adp_) {
        int atom = (i/3);
        if (Mask_.AtomInCharMask(atom)) {
          int resnum = (*fluctParm_)[atom].ResNum();
          int u11 = (int)(SumCoords2_[i  ] * 10000);
          int u22 = (int)(SumCoords2_[i+1] * 10000);
          int u33 = (int)(SumCoords2_[i+2] * 10000);
          // Calculate covariance: <XY> - <X><Y> etc.
          int u12 = (int)((Cross_[i  ] - SumCoords_[i  ] * SumCoords_[i+1]) * 10000);
          int u13 = (int)((Cross_[i+1] - SumCoords_[i  ] * SumCoords_[i+2]) * 10000);
          int u23 = (int)((Cross_[i+2] - SumCoords_[i+1] * SumCoords_[i+2]) * 10000);
          PDBfile& adpout = static_cast<PDBfile&>( *adpoutfile_ );
          adpout.WriteANISOU(
            atom+1, (*fluctParm_)[atom].c_str(), fluctParm_->Res(resnum).c_str(),
            fluctParm_->Res(resnum).ChainID(), fluctParm_->Res(resnum).OriginalResNum(),
            u11, u22, u33, u12, u13, u23, (*fluctParm_)[atom].ElementName(), 0 );
        }
      }
    }
  } else {
    // Atomic fluctuations
    for (int i = 0; i < SumCoords2_.size(); i+=3) {
      double fluct = SumCoords2_[i] + SumCoords2_[i+1] + SumCoords2_[i+2];
      if (fluct > 0)
        *result = sqrt(fluct);
      ++result;
    }
  }

  DataSet_Mesh& dset = static_cast<DataSet_Mesh&>( *dataout_ );
  if (outtype_ == BYATOM) {
    // By atom output
    dset.ModifyDim(Dimension::X).SetLabel("Atom");
    for (int atom = 0; atom < (int)Results.size(); atom++ ) {
      if (Mask_.AtomInCharMask(atom))
        dset.AddXY( atom+1, Results[atom] );
    }
  } else if (outtype_ == BYRES) { 
    // By residue output
    dset.ModifyDim(Dimension::X).SetLabel("Res");
    for (Topology::res_iterator residue = fluctParm_->ResStart();
                                residue != fluctParm_->ResEnd(); ++residue) {
      double xi = 0.0;
      double fluct = 0.0;
      for (int atom = (*residue).FirstAtom(); atom < (*residue).LastAtom(); atom++) {
        if ( Mask_.AtomInCharMask(atom) ) {
          double mass = (*fluctParm_)[atom].Mass(); 
          xi += mass;
          fluct += Results[atom] * mass;
        }
      }
      if (xi > Constants::SMALL) 
        dset.AddXY( residue - fluctParm_->ResStart() + 1, fluct / xi );
    }
  } else if (outtype_ == BYMASK) {
    // By mask output
    dset.ModifyDim(Dimension::X).SetLabel( Mask_.MaskExpression() );
    double xi = 0.0;
    double fluct = 0.0;
    for (int atom = 0; atom < (int)Results.size(); atom++) {
      if (Mask_.AtomInCharMask(atom)) {
        double mass = (*fluctParm_)[atom].Mass();
        xi += mass;
        fluct += Results[atom] * mass;
      }
    }
    if (xi > Constants::SMALL) 
      dset.AddXY( 1, fluct / xi );
  }
}
