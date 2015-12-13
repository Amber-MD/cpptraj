#include <cmath>
#include "Action_Volume.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Volume::Action_Volume() :
  vol_(0), sum_(0.0), sum2_(0.0), nframes_(0)
{ } 

void Action_Volume::Help() const {
  mprintf("\t[<name>] [out <filename>]\n  Calculate unit cell volume.\n");
}

// Action_Volume::Init()
Action::RetType Action_Volume::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  image_.InitImaging( true );
  // Get keywords
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  sum_ = 0.0;
  sum2_ = 0.0;
  nframes_ = 0;
  // Dataset
  vol_ = init.DSL().AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"Vol");
  if (vol_==0) return Action::ERR;
  // Add dataset to data file list
  if (outfile != 0) outfile->AddDataSet( vol_ );

  mprintf("    VOLUME:");
  if (outfile != 0)
    mprintf(" Output to '%s'.", outfile->DataFilename().full());
  mprintf("\n");

  return Action::OK;
}

// Action_Volume::Setup()
/** Set angle up for this parmtop. Get masks etc.
  */
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
  sum_ += volume;
  sum2_ += (volume * volume);
  nframes_++;

  return Action::OK;
}

void Action_Volume::Print() {
  double avg = 0.0, stdev = 0.0;
  if (nframes_ > 0) {
    avg = sum_ / (double)nframes_;
    stdev = (sum2_ / (double)nframes_) - (avg * avg);
    if (stdev > 0.0)
      stdev = sqrt(stdev);
    else
      stdev = 0.0;
  }
  mprintf("    VOLUME: Avg= %g  Stdev= %g (%i frames)\n", avg, stdev, nframes_);
} 
