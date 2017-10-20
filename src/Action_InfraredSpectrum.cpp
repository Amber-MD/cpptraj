#include "Action_InfraredSpectrum.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"
#include "Constants.h"
#include "DataSet_double.h"
#include "Corr.h"
#ifdef _OPENMP
#  include <omp.h>
#endif

Action_InfraredSpectrum::Action_InfraredSpectrum() : Action(HIDDEN),
  Vel_(0),
  VAC_(0),
  currentTop_(0),
  tstep_(0.0),
  maxLag_(-1),
  previousNselected_(-1),
  useFFT_(true)
{}

// Action_InfraredSpectrum::Help()
void Action_InfraredSpectrum::Help() const {
  mprintf("\t[<name>] [<mask>] [out <file>] [direct] [maxlag <lag>] [tstep <step>]\n"
          "\t[rawout <raw vector>]\n");
}

// Action_InfraredSpectrum::Init()
Action::RetType Action_InfraredSpectrum::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  DataFile* rawfile = init.DFL().AddDataFile( actionArgs.GetStringKey("rawout"), actionArgs );
  useFFT_ = !actionArgs.hasKey("direct");
  maxLag_ = actionArgs.getKeyInt("maxlag", -1);
  tstep_ = actionArgs.getKeyDouble("tstep", 1.0);
  mask_.SetMaskString( actionArgs.GetMaskNext() );
  previousNselected_ = -1;
  // DataSet
  Vel_ = (DataSet_Vector*)
         init.DSL().AddSet(DataSet::VECTOR,
                           MetaData(actionArgs.GetStringNext(), "raw"), "IR");
  if (Vel_ == 0) return Action::ERR;
  if (rawfile != 0) rawfile->AddDataSet( Vel_ );
  VAC_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(Vel_->Meta().Name(), "ac"));
  if (VAC_ == 0) return Action::ERR;
  if (outfile != 0) outfile->AddDataSet( VAC_ );
  Vel_->SetDim(Dimension::X, Dimension(0.0, tstep_, "Time (ps)"));
  VAC_->SetDim(Dimension::X, Dimension(0.0, tstep_, "Time (ps)"));

  mprintf("    INFRARED SPECTRUM:\n");
  mprintf("\tFor atoms in mask '%s'\n", mask_.MaskString());
  if (maxLag_ < 1)
    mprintf("\tMaximum lag will be half total # of frames");
  else
    mprintf("\tMaximum lag is %i frames", maxLag_);
  mprintf(", time step between frames is %f ps\n", tstep_);
  if (useFFT_)
    mprintf("\tUsing FFT to calculate autocorrelation function.\n");
  else
    mprintf("\tUsing direct method to calculate autocorrelation function.\n");
  mprintf("\tAutocorrelation function data in '%s'\n", VAC_->Meta().PrintName().c_str());
  if (outfile != 0)
    mprintf("\tAutocorrelation function output to '%s'\n", outfile->DataFilename().full());
  mprintf("\tRaw velocity*charge vector data in '%s'\n", Vel_->Meta().PrintName().c_str());
  if (rawfile != 0)
    mprintf("\tRaw velocity*charge vector output to '%s'\n", rawfile->DataFilename().full());
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

// Action_InfraredSpectrum::Print()
void Action_InfraredSpectrum::Print() {
  if (Vel_ == 0 || Vel_->Size() < 1) return;
  mprintf("    INFRARED SPECTRUM:\n");
  int maxlag;
  if (maxLag_ <= 0) {
    maxlag = (int)Vel_->Size() / 2;
    mprintf("\tSetting maximum lag to 1/2 total frames (%i)\n", maxlag);
  } else if (maxLag_ > (int)Vel_->Size()) {
    maxlag = (int)Vel_->Size();
    mprintf("\tSpecified maximum lag > total length, setting to %i\n", maxlag);
  } else
    maxlag = maxLag_;
  // Allocate space for output correlation function values.
  DataSet_double& Ct = static_cast<DataSet_double&>( *VAC_ );
  Ct.Resize( maxlag );
  if (!useFFT_) {
    // DIRECT METHOD 
    ParallelProgress progress( maxlag );
    int t;
    unsigned int dtmax, dt;
#   ifdef _OPENMP
#   pragma omp parallel private(t, dtmax, dt) firstprivate(progress)
    {
      progress.SetThread(omp_get_thread_num());
#     pragma omp for schedule(dynamic)
#   endif
      for (t = 0; t < maxlag; ++t)
      {
        progress.Update( t );
        dtmax = Vel_->Size() - t;
        for (dt = 0; dt < dtmax; ++dt)
          Ct[t] += (*Vel_)[dt] * (*Vel_)[dt + t];
        Ct[t] /= (double)dtmax;
        //mprintf("\tCt[%i]= %f\n", t, Ct[t]); // DEBUG
      }
#   ifdef _OPENMP
    } // END pragma omp parallel
#   endif
    progress.Finish();
  } else {
    // FFT METHOD
    // Since FFT is cyclic, unroll vectors into a 1D array; in the resulting
    // transformed array after FFT, every 3rd value will be the correlation
    // via dot products that we want (once it is normalized).
    unsigned int total_length = Vel_->Size() * 3;
    CorrF_FFT pubfft;
    pubfft.CorrSetup( total_length );
    ComplexArray data1 = pubfft.Array();
    //mprintf("Complex Array Size is %i (%i actual)\n", data1.size(), data1.size()*2);
    // Place vector from each frame into 1D array
    unsigned int nd = 0; // Will be used to index complex data
    for (DataSet_Vector::const_iterator vec = Vel_->begin(); vec != Vel_->end(); ++vec, nd+=6)
    {
      data1[nd  ] = (*vec)[0]; data1[nd+1] = 0.0;
      data1[nd+2] = (*vec)[1]; data1[nd+3] = 0.0;
      data1[nd+4] = (*vec)[2]; data1[nd+5] = 0.0;
    }
    data1.PadWithZero( total_length );
    pubfft.AutoCorr( data1 );
    // Normalization
    nd = 0;
    for (int t = 0; t < maxlag; t++, nd += 3)
      Ct[t] = data1[nd*2] * ( 3.0 / (double)((total_length - nd)) );
  }
}
