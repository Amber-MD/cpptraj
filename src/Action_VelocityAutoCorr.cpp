#include "Action_VelocityAutoCorr.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_VelocityAutoCorr::Action_VelocityAutoCorr() :
  velAC_(0), useVelInfo_(false), frameIdx_(0) {}

void Action_VelocityAutoCorr::Help() {
  mprintf("\t[<set name>] [<mask>] [usevelocity] [out <filename>]\n"
          "  Calculate velocity auto-correlation function for atoms in <mask>\n");
}

// Action_VelocityAutoCorr::Init()
Action::RetType Action_VelocityAutoCorr::Init(ArgList& actionArgs, TopologyList* PFL,
                                              FrameList* FL, DataSetList* DSL, 
                                              DataFileList* DFL, int debugIn)
{
  useVelInfo_ = actionArgs.hasKey("usevelocity");
  mask_.SetMaskString( actionArgs.GetMaskNext() );
  DataFile* outfile =  DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  frameIdx_ = 0;
  // Set up output data set
  velAC_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "VAC");
  if (velAC_ == 0) return Action::ERR;
  if (outfile != 0) outfile->AddSet( velAC_ );

  mprintf("    VELOCITYAUTOCORR:\n"
          "\tCalculate velocity auto-correlation function for atoms in mask '%s'\n",
          mask_.MaskString());
  if (useVelInfo_)
    mprintf("\tUsing velocity information present in frames.\n");
  else
    mprintf("\tCalculating velocities between consecutive frames.\n");
  if (outfile != 0)
    mprintf("\tOutput data set '%s' to '%s'\n", velAC_->Legend().c_str(), 
            outfile->DataFilename().full());
  return Action::OK;
}

// Action_VelocityAutoCorr::Setup()
/** For this to be valid the same # of atoms should be selected each time. */
Action::RetType Action_VelocityAutoCorr::Setup(Topology* currentParm,
                                               Topology** parmAddress)
{
  if (currentParm->SetupIntegerMask( mask_ )) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprintf("Warning: No atoms selected by mask.\n");
    return Action::ERR;
  }
  // If using velocity info, check that it is present.
  if (useVelInfo_ && !currentParm->HasVelInfo()) {
    mprinterr("Error: 'usevelocity' specified but no velocity info assocated with %s\n",
              currentParm->c_str());
    return Action::ERR;
  }
  // If velocities at time 0 have been set, check if the current number of 
  // selected atoms has remained the same.
  if (!Velocity0_.empty()) {
    if ((int)Velocity0_.size() != mask_.Nselected()) {
      mprinterr("Error: # of selected atoms %i has changed (previously %zu)\n",
                mask_.Nselected(), Velocity0_.size());
      return Action::ERR;
    }
  } else
    Velocity0_.resize( mask_.Nselected(), Vec3(0.0) );
  return Action::OK;
}

// Action_VelocityAutoCorr::DoAction()
Action::RetType Action_VelocityAutoCorr::DoAction(int frameNum, 
                                                  Frame* currentFrame,
                                                  Frame** frameAddress)
{
  double sum = 0.0;
  if (!useVelInfo_) {
    if (frameIdx_ == 1) {
      // This is the first frame which we can calculate pseudo-velocity.
      Varray::iterator v0 = Velocity0_.begin();
      for (AtomMask::const_iterator atom = mask_.begin();
                                    atom != mask_.end(); 
                                  ++atom, ++v0)
      {
        *v0 = Vec3( currentFrame->XYZ(*atom) ) - Vec3( previousFrame_.XYZ(*atom) );
        v0->Normalize();
        sum += (*v0) * (*v0);
      }
      sum /= (double)mask_.Nselected();
      velAC_->Add( frameNum - 1, &sum );
    } else if ( frameIdx_ > 1 ) {
      Varray::const_iterator v0 = Velocity0_.begin();
      for (AtomMask::const_iterator atom = mask_.begin();
                                    atom != mask_.end(); 
                                  ++atom, ++v0)
      {
        Vec3 vt = Vec3( currentFrame->XYZ(*atom) ) - Vec3( previousFrame_.XYZ(*atom) );
        vt.Normalize();
        sum += *v0 * vt;
        //sum += *v0 * (Vec3( currentFrame->XYZ(*atom) ) - Vec3( previousFrame_.XYZ(*atom) ));
      }
      sum /= (double)mask_.Nselected();
      velAC_->Add( frameNum - 1, &sum );
    }
    previousFrame_ = *currentFrame;
  } else {
    if (frameIdx_ == 0) {
      Varray::iterator v0 = Velocity0_.begin();
      for (AtomMask::const_iterator atom = mask_.begin();
                                    atom != mask_.end(); 
                                  ++atom, ++v0)
      {
        *v0 = Vec3( currentFrame->Vel( *atom ) );
        sum += (*v0) * (*v0);
      }
    } else { // frameIdx_ > 0
      Varray::const_iterator v0 = Velocity0_.begin();
      for (AtomMask::const_iterator atom = mask_.begin();
                                    atom != mask_.end();
                                  ++atom, ++v0)
        sum += *v0 * Vec3( currentFrame->Vel( *atom ) );
    }
    sum /= (double)mask_.Nselected();
    velAC_->Add( frameNum, &sum );
  }
  frameIdx_++;
  return Action::OK;
}
