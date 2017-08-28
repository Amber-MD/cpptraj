#include "Action_InfraredSpectrum.h"
#include "CpptrajStdio.h"
#include "Constants.h"

Action_InfraredSpectrum::Action_InfraredSpectrum() :
  currentTop_(0),
  maxLag_(-1),
  previousNselected_(-1),
  useFFT_(true)
{}

// Action_InfraredSpectrum::Help()
void Action_InfraredSpectrum::Help() const {

}

// Action_InfraredSpectrum::Init()
Action::RetType Action_InfraredSpectrum::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  DataFile* outfile =  init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  useFFT_ = !actionArgs.hasKey("direct");
  maxLag_ = actionArgs.getKeyInt("maxlag", -1);
  mask_.SetMaskString( actionArgs.GetMaskNext() );
  previousNselected_ = -1;
  // DataSet
  Vel_ = (DataSet_Vector*)
         init.DSL().AddSet(DataSet::VECTOR,
                           MetaData(actionArgs.GetStringNext(), "raw"), "IR");
  if (Vel_ == 0) return Action::ERR;
  if (outfile != 0) outfile->AddDataSet( Vel_ );

  mprintf("    INFRARED SPECTRUM:\n");
  mprintf("\tFor atoms in mask '%s'\n", mask_.MaskString());
  if (maxLag_ < 1)
    mprintf("\tMaximum lag will be half total # of frames");
  else
    mprintf("\tMaximum lag is %i frames", maxLag_);
  if (useFFT_)
    mprintf("\tUsing FFT to calculate autocorrelation function.\n");
  else
    mprintf("\tUsing direct method to calculate autocorrelation function.\n");
  if (outfile != 0)
    mprintf("\tOutput to '%s'\n", outfile->DataFilename().full());
  mprintf("\tRaw velocity*charge data in '%s'\n", Vel_->Meta().PrintName().c_str());
  return Action::OK;
}

// Action_InfraredSpectrum::Setup()
Action::RetType Action_InfraredSpectrum::Setup(ActionSetup& setup)
{
  if (!setup.CoordInfo().HasVel()) {
    mprinterr("Error: No velocity info present in frames.\n");
    return Action::ERR;
  }
  if (setup.Top().SetupIntegerMask( mask_ )) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprintf("Warning: No atoms selected by mask.\n");
    return Action::SKIP;
  }
  if (previousNselected_ != -1 && mask_.Nselected() != previousNselected_)
    mprintf("Warning: Selected # atoms has changed; was %i, now is %i\n",
            previousNselected_, mask_.Nselected());
  previousNselected_ = mask_.Nselected();
  // TODO: Cache charges?
  currentTop_ = setup.TopAddress();

  return Action::OK;
}

// Action_InfraredSpectrum::DoAction()
Action::RetType Action_InfraredSpectrum::DoAction(int frameNum, ActionFrame& frm)
{
  Vec3 sum(0.0);
  for (AtomMask::const_iterator atm = mask_.begin(); atm != mask_.end(); ++atm)
    sum += Vec3(frm.Frm().VelXYZ(*atm)) * Constants::AMBERTIME_TO_PS * (*currentTop_)[*atm].Charge();
  Vel_->AddVxyz( sum );
  return Action::OK;
}

//void Action_InfraredSpectrum
