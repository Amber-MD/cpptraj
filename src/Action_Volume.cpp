#include <cmath>
#include "Action_Volume.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Volume::Action_Volume() :
  vol_(0), sum_(0.0), sum2_(0.0), nframes_(0)
{ } 

void Action_Volume::Help() {
  mprintf("\t[<name>] [out <filename>]\n  Calculate unit cell volume.\n");
}

// Action_Volume::Init()
Action::RetType Action_Volume::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  image_.InitImaging( true );
  // Get keywords
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  sum_ = 0.0;
  sum2_ = 0.0;
  nframes_ = 0;
  // Dataset
  vol_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"Vol");
  if (vol_==0) return Action::ERR;
  // Add dataset to data file list
  if (outfile != 0) outfile->AddSet( vol_ );

  mprintf("    VOLUME:");
  if (outfile != 0)
    mprintf(" Output to '%s'.", outfile->DataFilename().full());
  mprintf("\n");

  return Action::OK;
}

// Action_Volume::Setup()
/** Set angle up for this parmtop. Get masks etc.
  */
Action::RetType Action_Volume::Setup(Topology* currentParm, Topology** parmAddress) {
  image_.SetupImaging( currentParm->BoxType() );
  if (!image_.ImagingEnabled()) {
    mprintf("Warning: No unit cell information, volume cannot be calculated for '%s'\n",
            currentParm->c_str());
    return Action::ERR;
  }

  return Action::OK;  
}

// Action_Volume::DoAction()
Action::RetType Action_Volume::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  Matrix_3x3 ucell, recip;
  double volume = 0.0;
  if (image_.ImageType() == ORTHO)
    volume = currentFrame->BoxCrd().BoxX() *
             currentFrame->BoxCrd().BoxY() *
             currentFrame->BoxCrd().BoxZ();
  else if (image_.ImageType() == NONORTHO)
    volume = currentFrame->BoxCrd().ToRecip( ucell, recip );
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
