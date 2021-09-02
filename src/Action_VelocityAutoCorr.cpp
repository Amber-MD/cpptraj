#include "Action_VelocityAutoCorr.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"
#include "DataSet_double.h"
#include "Constants.h"
#include "Corr.h"
#ifdef _OPENMP
#  include <omp.h>
#endif

// CONSTRUCTOR
Action_VelocityAutoCorr::Action_VelocityAutoCorr() :
  diffout_(0),
  VAC_(0),
  diffConst_(0),
  tstep_(0.0),
  maxLag_(0),
  useVelInfo_(false),
  useFFT_(false),
  normalize_(false)
{}

void Action_VelocityAutoCorr::Help() const {
  mprintf("\t[<set name>] [<mask>] [usecoords] [out <filename>] [diffout <file>]\n"
          "\t[maxlag <frames>] [tstep <timestep>] [direct] [norm]\n"
          "  Calculate velocity auto-correlation function for atoms in <mask>. If\n"
          "  'diffout' specified write calculated diffusion constants to <file>,\n"
          "  otherwise they will be written to STDOUT.\n");
}

// Action_VelocityAutoCorr::Init()
Action::RetType Action_VelocityAutoCorr::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  if (actionArgs.hasKey("usevelocity")) {
    mprinterr("Error: The 'usevelocity' keyword is deprecated. Velocity information\n"
              "Error:   is now used by default if present. To force cpptraj to use\n"
              "Error:   coordinates to estimate velocities (not recommended) use the\n"
              "Error:   'usecoords' keyword.\n");
    return Action::ERR;
  }
  useVelInfo_ = !actionArgs.hasKey("usecoords");
  if (mask_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;
  DataFile* outfile =  init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  diffout_ = init.DFL().AddCpptrajFile( actionArgs.GetStringKey("diffout"),
                                        "VAC diffusion constants", DataFileList::TEXT, true );
  maxLag_ = actionArgs.getKeyInt("maxlag", -1);
  tstep_ = actionArgs.getKeyDouble("tstep", 1.0);
  useFFT_ = !actionArgs.hasKey("direct");
  normalize_ = actionArgs.hasKey("norm");
  // Set up output data set
  VAC_ = init.DSL().AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "VAC");
  if (VAC_ == 0) return Action::ERR;
  // TODO: This should just be a scalar
  diffConst_ = init.DSL().AddSet(DataSet::DOUBLE,
                                 MetaData(VAC_->Meta().Name(), "D", MetaData::NOT_TS));
  if (diffConst_ == 0) return Action::ERR;
  if (outfile != 0) outfile->AddDataSet( VAC_ );
# ifdef MPI
  trajComm_ = init.TrajComm(); 
  if (trajComm_.Size() > 1 && !useVelInfo_)
    mprintf("\nWarning: When calculating velocities between consecutive frames,\n"
            "\nWarning:   'velocityautocorr' in parallel will not work correctly if\n"
            "\nWarning:   coordinates have been modified by previous actions (e.g. 'rms').\n\n");
  diffConst_->SetNeedsSync( false );
# endif
  mprintf("    VELOCITYAUTOCORR:\n"
          "\tCalculate velocity auto-correlation function for atoms in mask '%s'\n",
          mask_.MaskString());
  if (useVelInfo_)
    mprintf("\tUsing velocity information present in frames.\n");
  else
    mprintf("\tCalculating velocities between consecutive frames from coordinates.\n");
  if (outfile != 0)
    mprintf("\tOutput velocity autocorrelation function '%s' to '%s'\n", VAC_->legend(), 
            outfile->DataFilename().full());
  mprintf("\tWriting diffusion constants to '%s'\n", diffout_->Filename().full());
  if (maxLag_ < 1)
    mprintf("\tMaximum lag will be half total # of frames");
  else
    mprintf("\tMaximum lag is %i frames", maxLag_);
  mprintf(", time step between frames is %f ps\n", tstep_);
  if (useFFT_)
    mprintf("\tUsing FFT to calculate autocorrelation function.\n");
  else
    mprintf("\tUsing direct method to calculate autocorrelation function.\n");
  if (normalize_)
    mprintf("\tNormalizing autocorrelation function to 1.0\n");
  return Action::OK;
}

// Action_VelocityAutoCorr::Setup()
/** For this to be valid the same # of atoms should be selected each time. */
Action::RetType Action_VelocityAutoCorr::Setup(ActionSetup& setup) {
  if (setup.Top().SetupIntegerMask( mask_ )) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprintf("Warning: No atoms selected by mask.\n");
    return Action::SKIP;
  }
  // If using velocity info, check that it is present.
  if (useVelInfo_ && !setup.CoordInfo().HasVel()) {
    mprinterr("Error: No velocity info present in frames.\n");
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
Action::RetType Action_VelocityAutoCorr::DoAction(int frameNum, ActionFrame& frm) {
  if (!useVelInfo_) {
    // Calculate pseudo-velocities between frames.
    if (!previousFrame_.empty()) {
      // This is the first frame which we can calculate pseudo-velocity.
      VelArray::iterator vel = Vel_.begin();
      for (AtomMask::const_iterator atom = mask_.begin();
                                    atom != mask_.end(); 
                                  ++atom, ++vel)
        vel->AddVxyz( (Vec3(frm.Frm().XYZ(*atom)) - Vec3(previousFrame_.XYZ(*atom))) / tstep_ );
    }
    previousFrame_ = frm.Frm();
  } else {
    // Use velocity information in the frame.
    // Assume V is in Amber units.
    VelArray::iterator vel = Vel_.begin();
    for (AtomMask::const_iterator atom = mask_.begin();
                                  atom != mask_.end(); 
                                ++atom, ++vel)
      vel->AddVxyz( Vec3(frm.Frm().VelXYZ( *atom )) * Constants::AMBERTIME_TO_PS );
  }
  return Action::OK;
}

#ifdef MPI
// Action_VelocityAutoCorr::ParallelPreviousFramesRequired()
int Action_VelocityAutoCorr::ParallelPreviousFramesRequired() const {
  if (useVelInfo_)
    return 0;
  else
    return 1;
}

// Action_VelocityAutoCorr::ParallelPreloadFrames()
int Action_VelocityAutoCorr::ParallelPreloadFrames(FArray const& preload_frames) {
  if (!useVelInfo_) {
    unsigned int idx = preload_frames.size() - 1;
    previousFrame_ = preload_frames[idx];
  }
  return 0;
}

// Action_VelocityAutoCorr::SyncAction()
int Action_VelocityAutoCorr::SyncAction() {
  if (Vel_.empty()) return 0;
  // Get total number of frames. Assume same # vectors in each thread.
  int nframes = (int)Vel_[0].Size();
  std::vector<int> rank_frames( trajComm_.Size() );
  trajComm_.GatherMaster( &nframes, 1, MPI_INT, &rank_frames[0] );
  int total_frames = 0;
  if (trajComm_.Master())
    for (int rank = 0; rank < trajComm_.Size(); rank++)
      total_frames += rank_frames[rank];
  // Sync each vector set.
  for (VelArray::iterator vel = Vel_.begin(); vel != Vel_.end(); ++vel)
    vel->Sync( total_frames, rank_frames, trajComm_ );
  return 0;
}
#endif

// Action_VelocityAutoCorr::Print()
void Action_VelocityAutoCorr::Print() {
  if (Vel_.empty()) return;
  mprintf("    VELOCITYAUTOCORR:\n");
  mprintf("\t%zu vectors have been saved, total length of each = %zu frames.\n",
          Vel_.size(), Vel_[0].Size());
  int maxlag;
  if (maxLag_ <= 0) {
    maxlag = (int)Vel_[0].Size() / 2;
    mprintf("\tSetting maximum lag to 1/2 total frames (%i)\n", maxlag);
  } else if (maxLag_ > (int)Vel_[0].Size()) {
    maxlag = (int)Vel_[0].Size();
    mprintf("\tSpecified maximum lag > total length, setting to %i\n", maxlag);
  } else
    maxlag = maxLag_;
  // DEBUG
  //for (VelArray::iterator vel = Vel_.begin(); vel != Vel_.end(); ++vel) {
  //  mprintf("Vector %u:\n", vel - Vel_.begin());
  //  for (DataSet_Vector::iterator vec = vel->begin(); vec != vel->end(); ++vec)
  //    mprintf("\t%u {%f %f %f}\n", vec - vel->begin(), (*vec)[0], (*vec)[1], (*vec)[2]);
  //}
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
        dtmax = Vel_[0].Size() - t;
        for (dt = 0; dt < dtmax; ++dt)
        {
          for (VelArray::const_iterator vel = Vel_.begin(); vel != Vel_.end(); ++vel)
            Ct[t] += (*vel)[dt] * (*vel)[dt + t];
        }
        Ct[t] /= (double)(dtmax * Vel_.size());
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
    unsigned int total_length = Vel_[0].Size() * 3;
    CorrF_FFT pubfft;
    pubfft.CorrSetup( total_length );
    ComplexArray data1 = pubfft.Array();
    //mprintf("Complex Array Size is %i (%i actual)\n", data1.size(), data1.size()*2);
    ProgressBar progress( Vel_.size() );
    unsigned int nvel = 0;
    for (VelArray::iterator vel = Vel_.begin(); vel != Vel_.end(); ++vel, ++nvel)
    {
      //mprintf("Vector %u\n", vel - Vel_.begin()); // DEBUG
      progress.Update( nvel );
      // Place vector from each frame into 1D array
      unsigned int nd = 0; // Will be used to index complex data
      for (DataSet_Vector::const_iterator vec = vel->begin(); vec != vel->end(); ++vec, nd+=6)
      {
        //mprintf("\tFrame %u assigned to complex array %u\n", vec - vel->begin(), nd); // DEBUG
        data1[nd  ] = (*vec)[0]; data1[nd+1] = 0.0;
        data1[nd+2] = (*vec)[1]; data1[nd+3] = 0.0;
        data1[nd+4] = (*vec)[2]; data1[nd+5] = 0.0;
        //mprintf("\t  Complex[%u]= %f  Complex[%u]= %f  Complex[%u]= %f\n", // DEBUG
        //        nd, data1[nd], nd+2, data1[nd+2], nd+4, data1[nd+4]);
      }
      data1.PadWithZero( total_length );
      //mprintf("Before FFT:\n"); // DEBUG
      //for (unsigned int cd = 0; cd < (unsigned int)data1.size()*2; cd += 2)
      //  mprintf("\t%u: %f + %fi\n", cd/2, data1[cd], data1[cd+1]);
      pubfft.AutoCorr( data1 );
      //mprintf("Results of FFT:\n"); // DEBUG
      //for (unsigned int cd = 0; cd < (unsigned int)data1.size()*2; cd += 2)
      //  mprintf("\t%u: %f + %fi\n", cd/2, data1[cd], data1[cd+1]);
      // Increment nd by 3 here so it can be used for normalization.
      nd = 0;
      for (int t = 0; t < maxlag; t++, nd += 3) {
        //mprintf("\tdata1[%u] = %f", nd*2, data1[nd*2]); // DEBUG
        //mprintf("  norm= %u", total_length - nd); // DEBUG
        //Ct[t] += data1[nd*2] * 3.0 / (double)(total_length - nd);
        Ct[t] += data1[nd*2];
        //mprintf("  Ct[%i]= %f\n", t, Ct[t]); // DEBUG
      }
    }
    // Normalization
    for (int t = 0, nd = 0; t < maxlag; t++, nd += 3)
      Ct[t] *= ( 3.0 / (double)((total_length - nd) * Vel_.size()) );
      //Ct[t] /= (double)Vel_.size();
  }
  // Integration to get diffusion coefficient.
  VAC_->SetDim(Dimension::X, Dimension(0.0, tstep_, "Time (ps)"));
  mprintf("\tIntegrating data set %s, step is %f\n", VAC_->legend(), VAC_->Dim(0).Step());
  double total = Ct.Integrate( DataSet_1D::TRAPEZOID );
  const double ANG2_PS_TO_CM2_S = 10.0; // Convert Ang^2/ps to 1E-5 cm^2/s
  const char* tab = "\t";
  if (!diffout_->IsStream()) {
    mprintf("\tDiffusion constants output to '%s'\n", diffout_->Filename().full());
    diffout_->Printf("# Diffusion constants from VAC for atoms in '%s'\n", mask_.MaskString());
    tab = "";
  } 
  diffout_->Printf("%s3D= %g Ang.^2/ps, %g x10^-5 cm^2/s\n", tab, total,
                   total * ANG2_PS_TO_CM2_S);
  double diffusionConst = total * ANG2_PS_TO_CM2_S / 3.0;
  diffout_->Printf("%s D= %g Ang.^2/ps, %g x10^-5 cm^2/s\n", tab, total/3.0,
                   diffusionConst);
  diffConst_->Add(0, &diffusionConst);
  diffout_->Printf("%s6D= %g Ang.^2/ps, %g x10^-5 cm^2/s\n", tab, total*2.0,
                   total * ANG2_PS_TO_CM2_S * 2.0);
  if (normalize_) {
    // Normalize VAC fn to 1.0
    mprintf("\tNormalizing VAC function to 1.0, C[0]= %g\n", Ct[0]);
    double norm = 1.0 / Ct[0];
    for (int t = 0; t < maxlag; ++t)
      Ct[t] *= norm;
  }
}
