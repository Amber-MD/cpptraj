#include <cmath> // sqrt
#include "Action_AtomicFluct.h"
#include "CpptrajStdio.h"
#include "Constants.h" // PI
#include "DataSet_double.h"

// CONSTRUCTOR
Action_AtomicFluct::Action_AtomicFluct() :
  sets_(0),
  bfactor_(false),
  fluctParm_(0),
  outtype_(BYATOM),
  dataout_(0),
  outfile_(0)
{}

void Action_AtomicFluct::Help() {
  mprintf("\t[out <filename>] [<mask>] [byres | byatom | bymask] [bfactor]\n");
  mprintf("\t%s\n", ActionFrameCounter::HelpText);
  mprintf("\tCalculate atomic fluctuations of atoms in <mask>\n");
}

// Action_AtomicFluct::init()
Action::RetType Action_AtomicFluct::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get frame # keywords
  if (InitFrameCounter(actionArgs)) return Action::ERR;
  // Get other keywords
  bfactor_ = actionArgs.hasKey("bfactor");
  outfile_ = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs ); 
  if (actionArgs.hasKey("byres"))
    outtype_ = BYRES;
  else if (actionArgs.hasKey("bymask"))
    outtype_ = BYMASK;
  else if (actionArgs.hasKey("byatom") || actionArgs.hasKey("byatm"))
    outtype_ = BYATOM;
  // Get Mask
  Mask.SetMaskString( actionArgs.GetMaskNext()  );
  // Get DataSet name
  std::string setname = actionArgs.GetStringNext();
  // Add output dataset
  dataout_ = DSL->AddSet( DataSet::DOUBLE, setname, "Fluct" );
  if (dataout_ == 0) {
    mprinterr("Error: AtomicFluct: Could not allocate dataset for output.\n");
    return Action::ERR; 
  }
  if (outfile_ != 0) 
    outfile_->AddSet( dataout_ );
  mprintf("    ATOMICFLUCT: calculating");
  if (bfactor_)
    mprintf(" B factors");
  else
    mprintf(" atomic positional fluctuations");
  if (outfile_ != 0)
    mprintf(", output to file %s",outfile_->DataFilename().base());
  mprintf("\n                 Atom mask: [%s]\n",Mask.MaskString());
  FrameCounterInfo();
  if (!setname.empty())
    mprintf("\tData will be saved to set named %s\n", setname.c_str());

  return Action::OK;
}

// Action_AtomicFluct::setup()
Action::RetType Action_AtomicFluct::Setup(Topology* currentParm, Topology** parmAddress) {

  if (SumCoords_.Natom()==0) {
    // Set up frames if not already set
    SumCoords_.SetupFrame(currentParm->Natom());
    SumCoords2_.SetupFrame(currentParm->Natom());
    SumCoords_.ZeroCoords();
    SumCoords2_.ZeroCoords();
    // This is the parm that will be used for this calc
    fluctParm_ = currentParm;
    // Set up atom mask
    if (currentParm->SetupCharMask( Mask )) {
      mprinterr("Error: AtomicFluct: Could not set up mask [%s]\n",Mask.MaskString());
      return Action::ERR;
    }
    if (Mask.None()) {
      mprinterr("Error: AtomicFluct: No atoms selected [%s]\n",Mask.MaskString());
      return Action::ERR;
    }
  } else if (currentParm->Natom() != SumCoords_.Natom()) {
    // Check that current #atoms matches
    mprinterr("Error: AtomicFluct not yet supported for mulitple topologies with different\n");
    mprinterr("       #s of atoms (set up for %i, this topology has %i\n",
              SumCoords_.Natom(), currentParm->Natom());
    return Action::ERR;
  } 
  // NOTE: Print warning here when setting up multiple topologies?
  return Action::OK;
}

// Action_AtomicFluct::action()
Action::RetType Action_AtomicFluct::DoAction(int frameNum, Frame* currentFrame, 
                                             Frame** frameAddress) 
{
  if ( CheckFrameCounter( frameNum ) ) return Action::OK;
  SumCoords_ += *currentFrame;
  SumCoords2_ += ( (*currentFrame) * (*currentFrame) ) ;
  ++sets_;
  return Action::OK;
}

// Action_AtomicFluct::print() 
void Action_AtomicFluct::Print() {
  mprintf("    ATOMICFLUCT: Calculating fluctuations for %i sets.\n",sets_);

  double Nsets = (double)sets_;
  SumCoords_.Divide(Nsets);
  SumCoords2_.Divide(Nsets);
  //SumCoords2_ = SumCoords2_ - (SumCoords_ * SumCoords);
  SumCoords_ *= SumCoords_;
  SumCoords2_ -= SumCoords_;

  // DEBUG
  //mprintf("DEBUG: Converting to double: Original first coord:\n");
  //SumCoords2_.printAtomCoord(0);
  //mprintf("       Converted first coord: %lf %lf %lf\n",XYZ[0],XYZ[1],XYZ[2]);

  // Hold fluctuation results - initialize to 0
  std::vector<double> Results( SumCoords2_.Natom(), 0 );
  std::vector<double>::iterator result = Results.begin();

  if (bfactor_) {
    // Set up b factor normalization
    // B-factors are (8/3)*PI*PI * <r>**2 hence we do not sqrt the fluctuations
    outfile_->Dim(Dimension::Y).SetLabel("B-factors");
    double bfac = (8.0/3.0)*PI*PI;
    for (int i = 0; i < SumCoords2_.size(); i+=3) {
      double fluct = SumCoords2_[i] + SumCoords2_[i+1] + SumCoords2_[i+2];
      if (fluct > 0) 
        *result = bfac * fluct;
      ++result;
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

  int minElt = -1;
  DataSet_double& dset = static_cast<DataSet_double&>( *dataout_ );
  if (outtype_ == BYATOM) {
    // By atom output
    outfile_->Dim(Dimension::X).SetLabel("Atom");
    for (int atom = 0; atom < (int)Results.size(); atom++ ) {
      if (Mask.AtomInCharMask(atom)) {
        if (minElt == -1) minElt = atom;
        dset.AddElement( Results[atom] );
      }
    }
  } else if (outtype_ == BYRES) { 
    // By residue output
    outfile_->Dim(Dimension::X).SetLabel("Res");
    for (Topology::res_iterator residue = fluctParm_->ResStart();
                                residue != fluctParm_->ResEnd(); ++residue) {
      double xi = 0.0;
      double fluct = 0.0;
      for (int atom = (*residue).FirstAtom(); atom < (*residue).LastAtom(); atom++) {
        if ( Mask.AtomInCharMask(atom) ) {
          double mass = (*fluctParm_)[atom].Mass(); 
          xi += mass;
          fluct += Results[atom] * mass;
        }
      }
      if (xi > SMALL) { 
        dset.AddElement( fluct / xi );
        if (minElt == -1) minElt = (int)(residue - fluctParm_->ResStart());
      }
    }
  } else if (outtype_ == BYMASK) {
    // By mask output
    minElt = 0;
    outfile_->Dim(Dimension::X).SetLabel( Mask.MaskExpression() );
    double xi = 0.0;
    double fluct = 0.0;
    for (int atom = 0; atom < (int)Results.size(); atom++) {
      if (Mask.AtomInCharMask(atom)) {
        double mass = (*fluctParm_)[atom].Mass();
        xi += mass;
        fluct += Results[atom] * mass;
      }
    }
    if (xi > SMALL) 
      dset.AddElement( fluct / xi );
  }
  if (minElt > -1)
    outfile_->Dim(Dimension::X).SetMin( minElt+1 );
}

