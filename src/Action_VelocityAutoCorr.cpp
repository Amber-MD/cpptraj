#include "Action_VelocityAutoCorr.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"
#include "DataSet_Mesh.h"
#include "DataSet_double.h"
#include "Corr.h"
#ifdef _OPENMP
#  include "omp.h"
#endif

// CONSTRUCTOR
Action_VelocityAutoCorr::Action_VelocityAutoCorr() :
  useVelInfo_(false), useFFT_(false), VAC_(0), tstep_(0.0), maxLag_(0) {}

void Action_VelocityAutoCorr::Help() {
  mprintf("\t[<set name>] [<mask>] [usevelocity] [out <filename>]\n"
          "\t[maxlag <time>] [tstep <timestep>] [usefft]\n"
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
  maxLag_ = actionArgs.getKeyInt("maxlag", -1);
  tstep_ = actionArgs.getKeyDouble("tstep", 1.0);
  useFFT_ = actionArgs.hasKey("usefft");
  // Set up output data set
  VAC_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "VAC");
  if (VAC_ == 0) return Action::ERR;
  VAC_->Dim(0).SetStep( tstep_ );
  if (outfile != 0) outfile->AddSet( VAC_ ); 

  mprintf("    VELOCITYAUTOCORR:\n"
          "\tCalculate velocity auto-correlation function for atoms in mask '%s'\n",
          mask_.MaskString());
  if (useVelInfo_)
    mprintf("\tUsing velocity information present in frames.\n");
  else
    mprintf("\tCalculating velocities between consecutive frames.\n");
  if (outfile != 0)
    mprintf("\tOutput data set '%s' to '%s'\n", VAC_->Legend().c_str(), 
            outfile->DataFilename().full());
  if (maxLag_ < 1)
    mprintf("\tMaximum lag will be half total # of frames");
  else
    mprintf("\tMaximum lag is %i frames", maxLag_);
  mprintf(", time step is %f ps\n", tstep_);
  if (useFFT_)
    mprintf("\tUsing FFT to calculate autocorrelation function.\n");
  else
    mprintf("\tUsing direct method to calculate autocorrelation function.\n");
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
  // If we have already started recording velocities, check that the current 
  // number of selected atoms has remained the same.
  if (!Vel_.empty()) {
    if ((int)Vel_.size() != mask_.Nselected()) {
      mprinterr("Error: # of selected atoms %i has changed (previously %zu)\n",
                mask_.Nselected(), Vel_.size());
      return Action::ERR;
    }
  } else
    Vel_.resize( mask_.Nselected() );
  return Action::OK;
}

// Action_VelocityAutoCorr::DoAction()
Action::RetType Action_VelocityAutoCorr::DoAction(int frameNum, 
                                                  Frame* currentFrame,
                                                  Frame** frameAddress)
{
  if (!useVelInfo_) {
    // Calculate pseudo-velocities between frames.
    if (!previousFrame_.empty()) {
      // This is the first frame which we can calculate pseudo-velocity.
      VelArray::iterator vel = Vel_.begin();
      for (AtomMask::const_iterator atom = mask_.begin();
                                    atom != mask_.end(); 
                                  ++atom, ++vel)
        vel->AddVxyz( Vec3(currentFrame->XYZ(*atom)) - Vec3(previousFrame_.XYZ(*atom)) );
    }
    previousFrame_ = *currentFrame;
  } else {
    // Use velocity information in the frame.
    VelArray::iterator vel = Vel_.begin();
    for (AtomMask::const_iterator atom = mask_.begin();
                                  atom != mask_.end(); 
                                ++atom, ++vel)
      vel->AddVxyz( Vec3(currentFrame->Vel( *atom )) );
  }
  return Action::OK;
}

// Action_VelocityAutoCorr::Print()
void Action_VelocityAutoCorr::Print() {
  if (Vel_.empty()) return;
  mprintf("    VELOCITYAUTOCORR:\n");
  mprintf("\t%zu vectors have been saved, total length of each = %zu\n",
          Vel_.size(), Vel_[0].Size());
  unsigned int maxlag;
  if (maxLag_ <= 0) {
    maxlag = Vel_[0].Size() / 2;
    mprintf("\tSetting maximum lag to 1/2 total time (%u)\n", maxlag);
  } else if (maxLag_ > (int)Vel_[0].Size()) {
    maxlag = Vel_[0].Size();
    mprintf("\tSpecified maximum lag > total length, setting to %u\n", maxlag);
  } else
    maxlag = (unsigned int)maxLag_;
  // Allocate space for output correlation function values.
  DataSet_double& Ct = static_cast<DataSet_double&>( *VAC_ );
  Ct.Resize( maxlag );
  if (!useFFT_) {
    // DIRECT METHOD 
    ParallelProgress progress( maxlag );
    unsigned int t, dtmax, dt;
#   ifdef _OPENMP
#   pragma omp parallel private(t, dtmax, dt) firstprivate(progress)
    {
      progress.SetThread(omp_get_thread_num());
#     pragma omp for schedule(dynamic)
#   endif
      for (t = 0; t < maxlag; ++t)
      {
        progress.Update( t );
        dtmax = Vel_[0].Size() - t;
        for (dt = 0; dt < dtmax; ++dt)
        {
          for (VelArray::const_iterator vel = Vel_.begin(); vel != Vel_.end(); ++vel)
            Ct[t] += (*vel)[dt] * (*vel)[dt + t];
        }
        Ct[t] /= (double)(dtmax * Vel_.size());
        //mprintf("\t%i %f\n", t, sum);
      }
#   ifdef _OPENMP
    } // END pragma omp parallel
#   endif
    progress.Finish();
    double norm = 1.0 / Ct[0];
    for (unsigned int t = 0; t < maxlag; ++t)
      Ct[t] *= norm;
  } else {
    // FFT METHOD
    CorrF_FFT pubfft( Vel_[0].Size() );
    ComplexArray data1 = pubfft.Array();
    ProgressBar progress( Vel_.size() );
    unsigned int nvel = 0;
    for (VelArray::iterator vel = Vel_.begin(); vel != Vel_.end(); ++vel, ++nvel)
    {
      progress.Update( nvel );
      vel->CalcSphericalHarmonics( 1 );
      for (int midx = -1; midx < 2; ++midx)
      {
        data1.Assign( vel->SphericalHarmonics( midx ) );
        data1.PadWithZero( Vel_[0].Size() );
        pubfft.AutoCorr( data1 );
        for (unsigned int i = 0; i < maxlag; i++)
          Ct[i] += data1[i*2];
      }
    }
    double norm = DataSet_Vector::SphericalHarmonicsNorm( 1 ) / (double)Vel_.size();
    for (unsigned int i = 0; i < maxlag; ++i)
      Ct[i] *= (norm / (Vel_[0].Size() - i));
  }
  // Integration to get diffusion coefficient.
  mprintf("\tIntegrating data set %s, step is %f\n", VAC_->Legend().c_str(), VAC_->Dim(0).Step());
  DataSet_Mesh mesh;
  mesh.SetMeshXY( static_cast<DataSet_1D const&>(*VAC_) );
  double total = mesh.Integrate_Trapezoid();
  mprintf("\tIntegral= %g  Integral/3= %g\n", total, total / 3.0);
}
