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
  usePdbRes_(false),
  adpoutfile_(0),
  fluctParm_(0),
  outtype_(BYATOM),
  dataout_(0)
{}

void Action_AtomicFluct::Help() const {
  mprintf("\t[out <filename>] [<mask>] [byres [pdbres] | byatom | bymask] [bfactor]\n"
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
  if (actionArgs.hasKey("byres")) {
    outtype_ = BYRES;
    usePdbRes_ = actionArgs.hasKey("pdbres");
  } else if (actionArgs.hasKey("bymask"))
    outtype_ = BYMASK;
  else if (actionArgs.hasKey("byatom") || actionArgs.hasKey("byatm"))
    outtype_ = BYATOM;
  // Get Mask
  if (Mask_.SetMaskString( actionArgs.GetMaskNext()  )) return Action::ERR;
  // Get DataSet name
  std::string setname = actionArgs.GetStringNext();
  // Add output dataset
  MetaData md( setname, "", MetaData::NOT_TS );
  if (setname.empty()) {
    // Only overwrite legend if no name specified.
    if (bfactor_)
      md.SetLegend("B-factors");
    else
      md.SetLegend("AtomicFlx");
  }
  dataout_ = init.DSL().AddSet( DataSet::XYMESH, md, "Fluct" );
  if (dataout_ == 0) {
    mprinterr("Error: AtomicFluct: Could not allocate dataset for output.\n");
    return Action::ERR; 
  }
# ifdef MPI
  dataout_->SetNeedsSync( false ); // Not a time series
  trajComm_ = init.TrajComm();
# endif
  if (outfile != 0) 
    outfile->AddDataSet( dataout_ );

  mprintf("    ATOMICFLUCT: calculating");
  if (bfactor_)
    mprintf(" B factors");
  else
    mprintf(" atomic positional fluctuations");
  switch (outtype_) {
    case BYATOM: mprintf(" for atoms.\n"); break;
    case BYRES : mprintf(" over residues.\n"); break;
    case BYMASK: mprintf(" over entire atom mask.\n"); break;
  }
  if (usePdbRes_) mprintf("\tUsing PDB residue numbers if present in topology.\n");
  if (outfile != 0)
    mprintf("\tOutput to file %s\n", outfile->DataFilename().full());
  mprintf("\tAtom mask: [%s]\n",Mask_.MaskString());
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
  // Set up atom mask
  if (setup.Top().SetupIntegerMask( Mask_ )) return Action::ERR;
  Mask_.MaskInfo();
  if (Mask_.None()) {
    mprintf("Warning: No atoms selected for mask [%s]\n", Mask_.MaskString());
    return Action::SKIP;
  }
  if (SumCoords_.Natom()==0) {
    // First time setup
    SumCoords_.SetupFrame(Mask_.Nselected());
    SumCoords2_.SetupFrame(Mask_.Nselected());
    SumCoords_.ZeroCoords();
    SumCoords2_.ZeroCoords();
    if (calc_adp_) {
      Cross_.SetupFrame(Mask_.Nselected());
      Cross_.ZeroCoords();
    }
    // This is the parm that will be used for this calc
    fluctParm_ = setup.TopAddress();
  } else if (Mask_.Nselected() != SumCoords_.Natom()) {
    // Check that current #atoms matches
    mprinterr("Error: AtomicFluct not yet supported for mulitple topologies with different\n");
    mprinterr("       #s of atoms (set up for %i, this topology has %i\n",
              SumCoords_.Natom(), Mask_.Nselected());
    return Action::ERR;
  } else {
    if (fluctParm_ != setup.TopAddress())
      mprintf("Warning: Topology is changing. Will base output only using topology '%s'.\n",
              fluctParm_->c_str());
  }
  return Action::OK;
}

// Action_AtomicFluct::DoAction()
Action::RetType Action_AtomicFluct::DoAction(int frameNum, ActionFrame& frm) {
  if ( CheckFrameCounter( frm.TrajoutNum() ) ) return Action::OK;
  int sidx = 0;
  for (AtomMask::const_iterator atm = Mask_.begin(); atm != Mask_.end(); ++atm, sidx += 3)
  {
    int fidx = *atm * 3;
    SumCoords_[sidx  ]  +=  frm.Frm()[fidx  ];
    SumCoords2_[sidx  ] += (frm.Frm()[fidx  ]*frm.Frm()[fidx  ]);
    SumCoords_[sidx+1]  +=  frm.Frm()[fidx+1];
    SumCoords2_[sidx+1] += (frm.Frm()[fidx+1]*frm.Frm()[fidx+1]);
    SumCoords_[sidx+2]  +=  frm.Frm()[fidx+2];
    SumCoords2_[sidx+2] += (frm.Frm()[fidx+2]*frm.Frm()[fidx+2]);
  }
  if (calc_adp_) {
    sidx = 0;
    for (AtomMask::const_iterator atm = Mask_.begin(); atm != Mask_.end(); ++atm, sidx += 3)
    {
      int fidx = *atm * 3;
      Cross_[sidx  ] += frm.Frm()[fidx  ] * frm.Frm()[fidx+1]; // U12
      Cross_[sidx+1] += frm.Frm()[fidx  ] * frm.Frm()[fidx+2]; // U13
      Cross_[sidx+2] += frm.Frm()[fidx+1] * frm.Frm()[fidx+2]; // U23
    }
  }
  ++sets_;
  return Action::OK;
}

#ifdef MPI
int Action_AtomicFluct::SyncAction() {
  int total_frames = 0;
  trajComm_.ReduceMaster( &total_frames, &sets_, 1, MPI_INT, MPI_SUM );
  if (trajComm_.Master())
    sets_ = total_frames;
  SumCoords_.SumToMaster(trajComm_);
  SumCoords2_.SumToMaster(trajComm_);
  Cross_.SumToMaster(trajComm_);
  return 0;
}
#endif

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
        int idx = (i/3);
        int atom = Mask_[idx];
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
    for (int idx = 0; idx < (int)Results.size(); idx++)
      dset.AddXY( Mask_[idx]+1, Results[idx] );
  } else if (outtype_ == BYRES) {
    // By residue output
    dset.ModifyDim(Dimension::X).SetLabel("Res");
    int lastidx = (int)Results.size() - 1;
    double fluct = 0.0;
    double xi = 0.0;
    for (int idx = 0; idx < (int)Results.size(); idx++) {
      int atom = Mask_[idx];
      double mass = (*fluctParm_)[atom].Mass();
      fluct += Results[idx] * mass;
      xi += mass;
      int currentres = (*fluctParm_)[atom].ResNum();
      int nextres;
      if (idx != lastidx) {
        int nextatom = Mask_[idx+1];
        nextres = (*fluctParm_)[nextatom].ResNum();
      } else
        nextres = -1;
      if (nextres != currentres) {
        int resnum;
        if (usePdbRes_)
          resnum = fluctParm_->Res(currentres).OriginalResNum();
        else
          resnum = currentres + 1;
        dset.AddXY( resnum, fluct / xi );
        xi = 0.0;
        fluct = 0.0;
      }
    }
  } else if (outtype_ == BYMASK) {
    // By mask output
    dset.ModifyDim(Dimension::X).SetLabel( Mask_.MaskExpression() );
    double xi = 0.0;
    double fluct = 0.0;
    for (int idx = 0; idx < (int)Results.size(); idx++) {
      int atom = Mask_[idx];
      double mass = (*fluctParm_)[atom].Mass();
      xi += mass;
      fluct += Results[idx] * mass;
    }
    if (xi > Constants::SMALL) 
      dset.AddXY( 1, fluct / xi );
  }
}
