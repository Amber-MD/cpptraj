#include <cmath>
#include "Action_Volume.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"

// CONSTRUCTOR
Action_Volume::Action_Volume() : vol_(0) {} 

// void Action_Volume::Help()
void Action_Volume::Help() const {
  mprintf("\t[<name>] [out <filename>]\n"
          "  Calculate unit cell volume in Ang^3.\n");
}

// Action_Volume::Init()
Action::RetType Action_Volume::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  image_.InitImaging( true );
  // Get keywords
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  // Dataset
  vol_ = init.DSL().AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"Vol");
  if (vol_==0) return Action::ERR;
  // Add dataset to data file list
  if (outfile != 0) outfile->AddDataSet( vol_ );

  mprintf("    VOLUME: Calculating unit cell volume in Ang^3.\n");
  if (outfile != 0)
    mprintf("\tOutput to '%s'.\n", outfile->DataFilename().full());

  return Action::OK;
}

// Action_Volume::Setup()
Action::RetType Action_Volume::Setup(ActionSetup& setup) {
  image_.SetupImaging( setup.CoordInfo().TrajBox().Type() );
  if (!image_.ImagingEnabled()) {
    mprintf("Warning: No unit cell information, volume cannot be calculated for '%s'\n",
            setup.Top().c_str());
    return Action::SKIP;
  }

  return Action::OK;  
}

// Action_Volume::DoAction()
Action::RetType Action_Volume::DoAction(int frameNum, ActionFrame& frm) {
  Matrix_3x3 ucell, recip;
  double volume = 0.0;
  if (image_.ImageType() == ORTHO)
    volume = frm.Frm().BoxCrd().BoxX() *
             frm.Frm().BoxCrd().BoxY() *
             frm.Frm().BoxCrd().BoxZ();
  else if (image_.ImageType() == NONORTHO)
    volume = frm.Frm().BoxCrd().ToRecip( ucell, recip );
  vol_->Add(frameNum, &volume);

  return Action::OK;
}

// Action_Volume::Print()
void Action_Volume::Print() {
  if (vol_ != 0 && vol_->Size() > 0) {
    double stdev;
    double avg = ((DataSet_1D*)vol_)->Avg(stdev); 
    mprintf("    VOLUME: Avg= %g  Stdev= %g (%zu elements), Ang^3\n",
            avg, stdev, vol_->Size());
  }
}
