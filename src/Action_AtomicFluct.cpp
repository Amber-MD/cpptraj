#include <cmath> // sqrt
#include "Action_AtomicFluct.h"
#include "CpptrajStdio.h"
#include "Constants.h" // PI
#include "DataSet_double.h"

// CONSTRUCTOR
Action_AtomicFluct::Action_AtomicFluct() :
  sets_(0),
  bfactor_(false),
  calc_adp_(false),
  fluctParm_(0),
  outtype_(BYATOM),
  dataout_(0),
  outfile_(0)
{}

void Action_AtomicFluct::Help() {
  mprintf("\t[out <filename>] [<mask>] [byres | byatom | bymask] [bfactor]\n"
          "\t[calcadp [adpout <file>]]\n"
          "\t%s\n\tCalculate atomic fluctuations of atoms in <mask>\n",
          ActionFrameCounter::HelpText);
}

// Action_AtomicFluct::Init()
Action::RetType Action_AtomicFluct::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get frame # keywords
  if (InitFrameCounter(actionArgs)) return Action::ERR;
  // Get other keywords
  bfactor_ = actionArgs.hasKey("bfactor");
  calc_adp_ = actionArgs.hasKey("calcadp");
  if (calc_adp_) {
     adpoutname_ = actionArgs.GetStringKey("adpout");
     if (!bfactor_) bfactor_ = true;
  }
  outfile_ = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs ); 
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
  mprintf("\n                 Atom mask: [%s]\n",Mask_.MaskString());
  FrameCounterInfo();
  if (calc_adp_)
    mprintf("\tCalculating anisotropic displacement parameters.\n");
  if (!setname.empty())
    mprintf("\tData will be saved to set named %s\n", setname.c_str());

  return Action::OK;
}

// Action_AtomicFluct::Setup()
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
    if (currentParm->SetupCharMask( Mask_ )) {
      mprinterr("Error: Could not set up mask [%s]\n",Mask_.MaskString());
      return Action::ERR;
    }
    Mask_.MaskInfo();
    if (Mask_.None()) {
      mprinterr("Error: AtomicFluct: No atoms selected [%s]\n",Mask_.MaskString());
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

// Action_AtomicFluct::DoAction()
Action::RetType Action_AtomicFluct::DoAction(int frameNum, Frame* currentFrame, 
                                             Frame** frameAddress) 
{
  if ( CheckFrameCounter( frameNum ) ) return Action::OK;
  SumCoords_ += *currentFrame;
  SumCoords2_ += ( (*currentFrame) * (*currentFrame) ) ;
  ++sets_;
  return Action::OK;
}

// Action_AtomicFluct::Print() 
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
    CpptrajFile adpout;
    if (calc_adp_) adpout.OpenWrite( adpoutname_ );
    // Set up b factor normalization
    // B-factors are (8/3)*PI*PI * <r>**2 hence we do not sqrt the fluctuations
    // TODO: Set Y axis label in DataFile
    //outfile_->Dim(Dimension::Y).SetLabel("B-factors");
    double bfac = (8.0/3.0)*PI*PI;
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
          int u12 = (int)((SumCoords2_[i  ] + SumCoords2_[i+1]) * 10000);
          int u13 = (int)((SumCoords2_[i  ] + SumCoords2_[i+2]) * 10000);
          int u23 = (int)((SumCoords2_[i+1] + SumCoords2_[i+2]) * 10000);
          adpout.Printf("ANISOU%5i %4s%4s %c%4i%c %7i%7i%7i%7i%7i%7i      %2s%2i\n",
                        atom+1, (*fluctParm_)[atom].c_str(), fluctParm_->Res(resnum).c_str(),
                        ' ', resnum+1, ' ', u11, u22, u33, u12, u13, u23,
                        (*fluctParm_)[atom].ElementName(), 0);
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

  int minElt = -1;
  DataSet_double& dset = static_cast<DataSet_double&>( *dataout_ );
  if (outtype_ == BYATOM) {
    // By atom output
    dataout_->Dim(Dimension::X).SetLabel("Atom");
    for (int atom = 0; atom < (int)Results.size(); atom++ ) {
      if (Mask_.AtomInCharMask(atom)) {
        if (minElt == -1) minElt = atom;
        dset.AddElement( Results[atom] );
      }
    }
  } else if (outtype_ == BYRES) { 
    // By residue output
    dataout_->Dim(Dimension::X).SetLabel("Res");
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
      if (xi > SMALL) { 
        dset.AddElement( fluct / xi );
        if (minElt == -1) minElt = (int)(residue - fluctParm_->ResStart());
      }
    }
  } else if (outtype_ == BYMASK) {
    // By mask output
    minElt = 0;
    dataout_->Dim(Dimension::X).SetLabel( Mask_.MaskExpression() );
    double xi = 0.0;
    double fluct = 0.0;
    for (int atom = 0; atom < (int)Results.size(); atom++) {
      if (Mask_.AtomInCharMask(atom)) {
        double mass = (*fluctParm_)[atom].Mass();
        xi += mass;
        fluct += Results[atom] * mass;
      }
    }
    if (xi > SMALL) 
      dset.AddElement( fluct / xi );
  }
  if (minElt > -1)
    dataout_->Dim(Dimension::X).SetMin( minElt+1 );
}
