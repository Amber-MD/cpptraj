#include "Analysis_CalcDiffusion.h"
#include "CpptrajStdio.h"
#include "DataSet_double.h"
#include "ProgressBar.h"
#include "OnlineVarT.h"
#include <algorithm> // std::min, std::copy
#include <cmath> // sqrt
#ifdef _OPENMP
# include <omp.h>
#endif
#ifdef MPI
#include "DataSet_Coords_CRD.h"
#endif

/** CONSTRUCTOR */
Analysis_CalcDiffusion::Analysis_CalcDiffusion() :
  TgtTraj_(0),
  debug_(0),
  maxlag_(0),
  time_(0),
  avg_x_(0),
  avg_y_(0),
  avg_z_(0),
  avg_r_(0),
  avg_a_(0)
{}

// Analysis_CalcDiffusion::Help()
void Analysis_CalcDiffusion::Help() const {
  mprintf("\t[crdset <coords set>] [maxlag <maxlag>] [<mask>] [time <dt>]\n"
          "\t[<name>] [out <file>] [diffout <file>]\n");
}

// Analysis_CalcDiffusion::Setup()
Analysis::RetType Analysis_CalcDiffusion::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
# ifdef MPI
  trajComm_ = setup.TrajComm();
# endif
  debug_ = debugIn;
  // Attempt to get coords dataset from datasetlist
  std::string setname = analyzeArgs.GetStringKey("crdset");
  TgtTraj_ = (DataSet_Coords*)setup.DSL().FindCoordsSet( setname );
  if (TgtTraj_ == 0) {
    mprinterr("Error: Could not locate COORDS set corresponding to '%s'\n",
              setname.c_str());
    Help();
    return Analysis::ERR;
  }
# ifdef MPI
  if (trajComm_.Size() > 1) {
    if (TgtTraj_->Type() != DataSet::COORDS && 
        TgtTraj_->Type() != DataSet::TRAJ)
    {
      mprinterr("Error: Parallel calcdiffusion only works with COORDS and TRAJ sets currently.\n");
      return Analysis::ERR;
    }
  }
# endif
  maxlag_ = analyzeArgs.getKeyInt("maxlag", -1);
  time_ = analyzeArgs.getKeyDouble("time", 1.0);
  DataFile* outfile = setup.DFL().AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs );
  results_.AddDiffOut(setup.DFL(), analyzeArgs.GetStringKey("diffout"));
  // Mask
  if (mask1_.SetMaskString( analyzeArgs.GetMaskNext() )) {
    mprinterr("Error: Could not set mask string.\n");
    return Analysis::ERR;
  }
  // Add DataSets
  std::string dsname_ = analyzeArgs.GetStringNext();
  if (dsname_.empty())
    dsname_ = setup.DSL().GenerateDefaultName("Diff");
  avg_x_ = setup.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "X"));
  avg_y_ = setup.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "Y"));
  avg_z_ = setup.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "Z"));
  avg_r_ = setup.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "R"));
  avg_a_ = setup.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "A"));
  if (avg_x_ == 0 || avg_y_ == 0 || avg_z_ == 0 || avg_r_ == 0 || avg_a_ == 0) {
    mprinterr("Error: Could not allocate one or more average diffusion sets.\n");
    return Analysis::ERR;
  }
  if (outfile != 0) {
    outfile->AddDataSet( avg_r_ );
    outfile->AddDataSet( avg_x_ );
    outfile->AddDataSet( avg_y_ );
    outfile->AddDataSet( avg_z_ );
    outfile->AddDataSet( avg_a_ );
  }
  // Set X dim
  Dimension Xdim_ = Dimension(0.0, time_, "Time");
  avg_x_->SetDim(Dimension::X, Xdim_);
  avg_y_->SetDim(Dimension::X, Xdim_);
  avg_z_->SetDim(Dimension::X, Xdim_);
  avg_r_->SetDim(Dimension::X, Xdim_);
  avg_a_->SetDim(Dimension::X, Xdim_);
  // Set up diffusion sets
  if (results_.CreateDiffusionSets(setup.DSL(), dsname_))
    return Analysis::ERR;

  mprintf("    CALCDIFFUSION: Calculating diffusion from COORDS set '%s'\n", TgtTraj_->legend());
  if (maxlag_ > 0)
    mprintf("\tMaximum lag is %i frames.\n", maxlag_);
  mprintf("\tUsing atoms selected by mask '%s'\n", mask1_.MaskString());
  mprintf("\tData set name: %s\n", dsname_.c_str());
  if (outfile != 0)
    mprintf("\tOutput to '%s'\n", outfile->DataFilename().full());
  results_.Info();
# ifdef MPI
  mprintf("\tDividing frames among %i processes.\n", trajComm_.Size());
# endif
  //rprintf("DEBUG: Setup finished successfully.\n");

  return Analysis::OK;
}

# ifdef MPI
/** Flatten a Stats<double> array. */
static inline void flattenArray(std::vector<double>& out, std::vector< Stats<double> > const& in)
{
  out.clear();
  out.reserve( in.size()*3 );
  for (std::vector< Stats<double> >::const_iterator it = in.begin(); it != in.end(); ++it) {
    out.push_back( it->nData() );
    out.push_back( it->mean() );
    out.push_back( it->M2() );
  }
}
  
/** Combine a Stats<double> array on master. */
static inline void sumToMaster(std::vector<double>& dbuf, std::vector< Stats<double> >& AX,
                               Parallel::Comm const& trajComm, int maxlag)
{
  dbuf.resize(maxlag*3);
  if (trajComm.Master()) {
    for (int rank = 1; rank < trajComm.Size(); rank++) {
      // Receive flattened array from rank
      trajComm.SendMaster( &dbuf[0], dbuf.size(), rank, MPI_DOUBLE);
      // Combine elements
      for (int idx = 0; idx < maxlag; idx++) {
        int i3 = idx * 3;
        AX[idx].Combine( Stats<double>( dbuf[i3], dbuf[i3+1], dbuf[i3+2] ) );
      }
    }
  } else {
    // Flatten array
    flattenArray( dbuf, AX );
    // Send to master
    trajComm.SendMaster( &dbuf[0], dbuf.size(), trajComm.Rank(), MPI_DOUBLE );
  }
}
# endif

// Analysis_CalcDiffusion::Analyze()
Analysis::RetType Analysis_CalcDiffusion::Analyze() {
# ifdef MPI
  // Need to make sure each process has access to the frames it needs.
  // FIXME do a better job determining what frames are actually needed on each process.
  if (TgtTraj_->Type() == DataSet::COORDS) {
    DataSet_Coords_CRD& crd = static_cast<DataSet_Coords_CRD&>( *TgtTraj_ );
    if (crd.Bcast( trajComm_ )) {
      rprinterr("Error: COORDS broadcast failed.\n");
      return Analysis::ERR;
    }
  }
  //rprintf("DEBUG: COORDS set has %zu frames.\n", TgtTraj_->Size());
# endif
  if (TgtTraj_->Size() < 1) {
    mprinterr("Error: COORDS set '%s' is empty.\n", TgtTraj_->legend());
    return Analysis::ERR;
  }
  if (TgtTraj_->Size() == 1) {
    mprinterr("Error: COORDS set '%s' has only 1 frame.\n");
    return Analysis::ERR;
  }
  if (TgtTraj_->Top().SetupIntegerMask( mask1_ )) {
    mprinterr("Error: Could not set up mask '%s'\n", mask1_.MaskString());
    return Analysis::ERR;
  }
  mask1_.MaskInfo();
  if (mask1_.None()) {
    mprinterr("Error: Nothing selected by mask '%s'\n", mask1_.MaskString());
    return Analysis::ERR;
  }

  int maxframes = (int)TgtTraj_->Size();
  if (maxlag_ < 1) {
    mprintf("\t'maxlag' not specified, using half the number of input frames.\n");
    maxlag_ = (maxframes / 2);
  } else {
    int halfmax = (maxframes / 2);
    if (maxlag_ > halfmax)
      mprintf("Warning: Lag %i is more than half the number of input frames (N/2=%i)\n", maxlag_, halfmax);
  }
  mprintf("\tMax lag is %i frames.\n", maxlag_);
  int stopframe = maxframes - maxlag_;

  mprintf("\tCalculating diffusion from set '%s' for atoms in mask '%s' from t=0 to %g ps.\n",
          TgtTraj_->legend(), mask1_.MaskString(), (double)(maxlag_-1) * time_);
  mprintf("\tUsing frames 1 to %i as time origins.\n", stopframe+1);

  if (stopframe < 1) {
    mprinterr("Error: Stop frame is less than 1.\n");
    return Analysis::ERR;
  }

  // Allocate temp. sets
  typedef std::vector< Stats<double> > Darray;
  typedef std::vector< Darray > Tarray;
  int nthreads = 1;
# ifdef _OPENMP
# pragma omp parallel
  {
# pragma omp master
  nthreads = omp_get_num_threads();
  }
# endif
  if (nthreads > 1) mprintf("\tParallelizing calculation with %i OpenMP threads.\n", nthreads);
  Tarray thread_X_( nthreads );
  Tarray thread_Y_( nthreads );
  Tarray thread_Z_( nthreads );
  Tarray thread_A_( nthreads );
  Tarray thread_R_( nthreads );
  for (int t = 0; t < nthreads; t++) {
    thread_X_[t].resize( maxlag_ );
    thread_Y_[t].resize( maxlag_ );
    thread_Z_[t].resize( maxlag_ );
    thread_A_[t].resize( maxlag_ );
    thread_R_[t].resize( maxlag_ );
  }

  int idx0;
  Frame frm0 = TgtTraj_->AllocateFrame();
  Frame frm1 = frm0;

  int my_start, my_stop;
#ifdef MPI
  int my_frames = trajComm_.DivideAmongProcesses( my_start, my_stop, stopframe + 1 );
  std::vector<int> all_frames( trajComm_.Size() );
  std::vector<int> all_start( trajComm_.Size() );
  std::vector<int> all_stop( trajComm_.Size() );
  trajComm_.GatherMaster(&my_frames, 1, MPI_INT, &all_frames[0]);
  trajComm_.GatherMaster(&my_start, 1, MPI_INT, &all_start[0]);
  trajComm_.GatherMaster(&my_stop, 1, MPI_INT, &all_stop[0]);
  if (trajComm_.Master()) {
    for (int rank = 0; rank < trajComm_.Size(); rank++)
      mprintf("\tRank %i processing %i frames (%i to %i).\n", rank, all_frames[rank], all_start[rank]+1, all_stop[rank]);
  }
  if (my_frames < 1)
    rprintf("Warning: Rank is processing less than 1 frame. Should use fewer processes.\n");
#else /* MPI */
  my_start = 0; 
  my_stop = stopframe + 1;
#endif /* MPI */

  ParallelProgress progress( my_stop );

  int mythread = 0;
# ifdef _OPENMP
# pragma omp parallel private(idx0, mythread) firstprivate(frm0, frm1, progress)
  {
  mythread = omp_get_thread_num();
  progress.SetThread( mythread );
# pragma omp for schedule(dynamic)
# endif /* _OPENMP */
  // LOOP OVER FRAMES
  for (idx0 = my_start; idx0 < my_stop; idx0++)
  {
    progress.Update( idx0 );
//    mprintf("DEBUG: (t=%g) %i to", (double)idx0*time_, idx0);
    TgtTraj_->GetFrame(idx0, frm0);
    int endidx = std::min(idx0 + maxlag_, maxframes);
    int tidx = 0;
    for (int idx1 = idx0; idx1 < endidx; idx1++, tidx++)
    {
//      mprintf(" %i", idx1);
      // TODO for idx1==idx0 this is the same frame
      TgtTraj_->GetFrame(idx1, frm1);
      // Loop over atoms
      for (AtomMask::const_iterator at = mask1_.begin(); at != mask1_.end(); ++at)
      {
        const double* xyz0 = frm0.XYZ( *at );
        const double* xyz1 = frm1.XYZ( *at );
        double delx = xyz1[0] - xyz0[0];
        double dely = xyz1[1] - xyz0[1];
        double delz = xyz1[2] - xyz0[2];
        // Calc distances for this atom
        double distx = delx * delx;
        double disty = dely * dely;
        double distz = delz * delz;
        double dist2 = distx + disty + distz;
//        mprintf("DEBUG: At=%i  frm %i to %i  t=%g  d2=%g\n", *at+1, idx0+1, idx1+1, (double)tidx*time_, dist2);
        // Accumulate distances
        thread_X_[mythread][tidx].accumulate( distx );
        thread_Y_[mythread][tidx].accumulate( disty );
        thread_Z_[mythread][tidx].accumulate( distz );
        thread_R_[mythread][tidx].accumulate( dist2 );
        thread_A_[mythread][tidx].accumulate( sqrt(dist2) );
      } // END loop over atoms
    } // END inner loop
//    mprintf("\n");
  } // END outer loop
# ifdef _OPENMP
  } // END omp parallel
  // Sum into thread 0 array
  for (int t = 1; t < nthreads; t++) {
    for (idx0 = 0; idx0 < maxlag_; idx0++) {
      thread_X_[0][idx0].Combine( thread_X_[t][idx0] );
      thread_Y_[0][idx0].Combine( thread_Y_[t][idx0] );
      thread_Z_[0][idx0].Combine( thread_Z_[t][idx0] );
      thread_R_[0][idx0].Combine( thread_R_[t][idx0] );
      thread_A_[0][idx0].Combine( thread_A_[t][idx0] );
    }
  }
# endif
  progress.Finish();

# ifdef MPI
  // Sum arrays back down to the master
  std::vector<double> dbuf( maxlag_*3, 0 );
  sumToMaster( dbuf, thread_X_[0], trajComm_, maxlag_ );
  sumToMaster( dbuf, thread_Y_[0], trajComm_, maxlag_ );
  sumToMaster( dbuf, thread_Z_[0], trajComm_, maxlag_ );
  sumToMaster( dbuf, thread_R_[0], trajComm_, maxlag_ );
  sumToMaster( dbuf, thread_A_[0], trajComm_, maxlag_ );
# endif /* MPI */
  // Store averages
  DataSet_double& AX = static_cast<DataSet_double&>( *avg_x_ );
  DataSet_double& AY = static_cast<DataSet_double&>( *avg_y_ );
  DataSet_double& AZ = static_cast<DataSet_double&>( *avg_z_ );
  DataSet_double& AA = static_cast<DataSet_double&>( *avg_a_ );
  DataSet_double& AR = static_cast<DataSet_double&>( *avg_r_ );
  AX.Resize( maxlag_ );
  AY.Resize( maxlag_ );
  AZ.Resize( maxlag_ );
  AA.Resize( maxlag_ );
  AR.Resize( maxlag_ );

  unsigned int maxcount = thread_X_[0][0].nData();
  unsigned int mincount = maxcount;
  for (idx0 = 0; idx0 < maxlag_; idx0++) {
    maxcount = std::max( maxcount, (unsigned int)thread_X_[0][idx0].nData() );
    mincount = std::min( mincount, (unsigned int)thread_X_[0][idx0].nData() );
    if (debug_ > 0)
      mprintf("DEBUG: Average at t=%g is from %g data points.\n", (double)idx0*time_, thread_X_[0][idx0].nData());
    AX[idx0] = thread_X_[0][idx0].mean();
    AY[idx0] = thread_Y_[0][idx0].mean();
    AZ[idx0] = thread_Z_[0][idx0].mean();
    AR[idx0] = thread_R_[0][idx0].mean();
    AA[idx0] = thread_A_[0][idx0].mean();
  }
  if (maxcount == mincount)
    mprintf("\t%u data points contributed to each average.\n", maxcount);
  else
    mprintf("\tMax # points averaged = %u, min # points averaged = %u\n", maxcount, mincount);

  // Calculate diffusion constants
  std::string const& name = avg_r_->Meta().Name();
  unsigned int set = 0;
  results_.CalcDiffusionConst( set, avg_r_, 3, name + "_AvgDr" );
  results_.CalcDiffusionConst( set, avg_x_, 1, name + "_AvgDx" );
  results_.CalcDiffusionConst( set, avg_y_, 1, name + "_AvgDy" );
  results_.CalcDiffusionConst( set, avg_z_, 1, name + "_AvgDz" );
  //rprintf("DEBUG: Analysis finished successfully.\n");

  return Analysis::OK;
}
