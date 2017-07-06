#include <cmath>
#include "Action_Volume.h"
#include "CpptrajStdio.h"
#include "Constants.h" // NA

// CONSTRUCTOR
Action_Volume::Action_Volume() :
  vol_(0), sum_(0.0), sum2_(0.0), nframes_(0), calcDensity_(false)
{ } 

// void Action_Volume::Help()
void Action_Volume::Help() const {
  mprintf("\t[<name>] [out <filename>] [density]\n"
          "  Calculate unit cell volume in Ang^3. If 'density' is specified\n"
          "  the system density in g/cm^3 instead.\n");
}

/** Convert units from amu/Ang^3 to g/cm^3 */
const double Action_Volume::AMU_ANG_TO_G_CM3 = Constants::NA * 1E-24;

// Action_Volume::Init()
Action::RetType Action_Volume::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  image_.InitImaging( true );
  // Get keywords
  calcDensity_ = actionArgs.hasKey("density");
  if (calcDensity_)
    mask_.SetMaskString( "*" );
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  sum_ = 0.0;
  sum2_ = 0.0;
  nframes_ = 0;
  // Dataset
  if (!calcDensity_)
    vol_ = init.DSL().AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"Vol");
  else
    vol_ = init.DSL().AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"Density");
  if (vol_==0) return Action::ERR;
  // Add dataset to data file list
  if (outfile != 0) outfile->AddDataSet( vol_ );

  mprintf("    VOLUME:");
  if (calcDensity_)
    //mprintf(" Calculating density (in g/cm^3) of atoms in mask '%s'.\n", mask_.MaskString());
    mprintf(" Calculating system density (in g/cm^3).\n");
  else
    mprintf(" Calculating unit cell volume in Ang^3.\n");
  if (outfile != 0)
    mprintf("\tOutput to '%s'.", outfile->DataFilename().full());
  mprintf("\n");

  return Action::OK;
}

// Action_Volume::Setup()
/** Set angle up for this parmtop. Get masks etc.
  */
Action::RetType Action_Volume::Setup(ActionSetup& setup) {
  if (calcDensity_) {
    if (setup.Top().SetupIntegerMask( mask_ )) return Action::ERR;
    //mask_.MaskInfo();
    if (mask_.None()) {
      mprintf("Warning: No atoms selected.\n");
      return Action::SKIP;
    }
  }
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
  if (calcDensity_) {
    double totalMass = 0.0;
    for (AtomMask::const_iterator at = mask_.begin(); at != mask_.end(); ++at)
      totalMass += frm.Frm().Mass(*at);
    volume = totalMass / (volume * AMU_ANG_TO_G_CM3);
  }
  vol_->Add(frameNum, &volume);
  sum_ += volume;
  sum2_ += (volume * volume);
  nframes_++;

  return Action::OK;
}

// Action_Volume::Print()
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
  const char* dstr = "";
  const char* units = "Ang^3";
  if (calcDensity_) {
    dstr = "Density ";
    units = "g/cm^3";
  }
  mprintf("    VOLUME: %sAvg= %g  Stdev= %g (%i frames), %s\n",
          dstr, avg, stdev, nframes_, units);
}
