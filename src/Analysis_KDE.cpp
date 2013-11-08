#include <cmath> // pow
#include "Analysis_KDE.h"
#include "CpptrajStdio.h"
#include "DataSet_double.h"
#include "Constants.h" // TWOPI
#ifdef _OPENMP
#  include "omp.h"
#endif

// CONSTRUCTOR
Analysis_KDE::Analysis_KDE() : 
  data_(0), q_data_(0), bandwidth_(0.0), output_(0), kldiv_(0),
  amddata_(0), calcFreeE_(false),
  Kernel_(&Analysis_KDE::GaussianKernel) {}

void Analysis_KDE::Help() {
  mprintf("\t<dataset> [bandwidth <bw>] [out <file>] [name <dsname>]\n"
          "\t[min <min>] [max <max] [step <step>] [bins <bins>]\n"
          "\t[kldiv <dsname2> [klout <outfile>]] [amd <amdboost_data>]\n"
          "\tHistogram 1D data set using a kernel density estimator.\n");
}

// Analysis_KDE::Setup()
Analysis::RetType Analysis_KDE::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  Dimension Xdim;
  if (analyzeArgs.Contains("min"))
    Xdim.SetMin( analyzeArgs.getKeyDouble("min", 0.0) );
  if (analyzeArgs.Contains("max"))
    Xdim.SetMax( analyzeArgs.getKeyDouble("max", 0.0) );
  Xdim.SetStep( analyzeArgs.getKeyDouble("step", -1.0) );
  Xdim.SetBins( analyzeArgs.getKeyInt("bins", -1) );
  if (Xdim.Step() < 0.0 && Xdim.Bins() < 0) {
    mprinterr("Error: Must set either bins or step.\n");
    return Analysis::ERR;
  }
  std::string setname = analyzeArgs.GetStringKey("name");
  bandwidth_ = analyzeArgs.getKeyDouble("bandwidth", -1.0);
  DataFile* outfile = DFLin->AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs );
  DataFile* klOutfile = 0;
  // Get second data set for KL divergence calc.
  std::string q_dsname = analyzeArgs.GetStringKey("kldiv");
  if (!q_dsname.empty()) {
    q_data_ = datasetlist->GetDataSet( q_dsname );
    if (q_data_ == 0) {
      mprinterr("Error: Data set %s not found.\n", q_dsname.c_str());
      return Analysis::ERR;
    }
    if (q_data_->Ndim() != 1) {
      mprinterr("Error: Only 1D data sets supported.\n");
      return Analysis::ERR;
    }
    klOutfile = DFLin->AddDataFile( analyzeArgs.GetStringKey("klout"), analyzeArgs );
  } else {
    q_data_ = 0;
    kldiv_ = 0;
  }
  // Get AMD boost data set
  std::string amdname = analyzeArgs.GetStringKey("amd");
  if (!amdname.empty()) {
    amddata_ = datasetlist->GetDataSet( amdname );
    if (amddata_ == 0) {
      mprinterr("Error: AMD data set %s not found.\n", amdname.c_str());
      return Analysis::ERR;
    }
    if (amddata_->Ndim() != 1) {
      mprinterr("Error: AMD data set must be 1D.\n");
      return Analysis::ERR;
    }
  } else
    amddata_ = 0;
  calcFreeE_ = analyzeArgs.hasKey("free");

  // Get data set
  data_ = datasetlist->GetDataSet( analyzeArgs.GetStringNext() );
  if (data_ == 0) {
    mprinterr("Error: No data set or invalid data set name specified\n");
    return Analysis::ERR;
  }
  if (data_->Ndim() != 1) {
    mprinterr("Error: Only 1D data sets supported.\n");
    return Analysis::ERR;
  }
  
  // Output data set
  output_ = datasetlist->AddSet(DataSet::DOUBLE, setname, "kde");
  output_->SetDim(Dimension::X, Xdim);
  if (outfile != 0) outfile->AddSet( output_ );
  // Output for KL divergence calc.
  if ( q_data_ != 0 ) {
    kldiv_ = datasetlist->AddSetAspect(DataSet::DOUBLE, output_->Name(), "kld");
    if (klOutfile != 0) klOutfile->AddSet( kldiv_ );
  }

  mprintf("    KDE: Using gaussian KDE to histogram set \"%s\"\n", data_->Legend().c_str());
  if (amddata_!=0)
    mprintf("\tPopulating bins using AMD boost from data set %s\n",
            amddata_->Legend().c_str());
  if (q_data_ != 0) {
    mprintf("\tCalculating Kullback-Leibler divergence with set \"%s\"\n", 
            q_data_->Legend().c_str());
  }
  if (bandwidth_ < 0.0)
    mprintf("\tBandwidth will be estimated.\n");
  else
    mprintf("\tBandwidth= %f\n", bandwidth_);
  return Analysis::OK;
}

const double Analysis_KDE::ONE_OVER_ROOT_TWOPI = 1.0 / sqrt( TWOPI );

double Analysis_KDE::GaussianKernel(double u) const {
  return ( ONE_OVER_ROOT_TWOPI * exp( -0.5 * u * u ) );
}

Analysis::RetType Analysis_KDE::Analyze() {
  // Declare variables here in case of OPENMP
  unsigned int frame, bin, validPoint;
  double val, increment;
  bool Pzero, Qzero;
# ifdef _OPENMP
  int numthreads, mythread, nt;
  double** out_thread;
# endif

  Dimension& Xdim = output_->Dim(0);
  DataSet_1D const& In = static_cast<DataSet_1D const&>( *data_ );
  size_t inSize = In.Size();
  // Set output set dimensions from input set if necessary.
  if (!Xdim.MinIsSet())
    Xdim.SetMin( In.Min() );
  if (!Xdim.MaxIsSet())
    Xdim.SetMax( In.Max() );
  if (Xdim.CalcBinsOrStep()) return Analysis::ERR;
  Xdim.PrintDim();
  // TODO: When calculating KL divergence check q_data_ dimension.

  // Allocate output set
  DataSet_double& Out = static_cast<DataSet_double&>( *output_ );
  Out.Resize( Xdim.Bins() );
  size_t outSize = Out.Size();

  // Estimate bandwidth from normal distribution approximation if necessary.
  if (bandwidth_ < 0.0) {
    double stdev;
    In.Avg( stdev );
    double N_to_1_over_5 = pow( (double)inSize, (-1.0/5.0) );
    bandwidth_ = 1.06 * stdev * N_to_1_over_5;
    mprintf("\tDetermined bandwidth from normal distribution approximation: %f\n", bandwidth_);
  }

  // Set up increments
  std::vector<double> Increments(inSize, 1.0);
  if (amddata_ != 0) {
    DataSet_1D& AMD = static_cast<DataSet_1D&>( *amddata_ );
    if (AMD.Size() != inSize) {
      if (AMD.Size() < inSize) {
        mprinterr("Error: Size of AMD data set %zu < input data set %zu\n",
                  AMD.Size(), inSize);
        return Analysis::ERR;
      } else {
        mprintf("Warning: Size of AMD data set %zu > input data set %zu\n",
                AMD.Size(), inSize);
      }
    }
    for (unsigned int i = 0; i < inSize; i++)
      Increments[i] = exp( AMD.Dval(i) );
  }
 
  double total = 0.0;
  if (q_data_ == 0) {
    // Calculate KDE, loop over input data
#ifdef _OPENMP
#   pragma omp parallel private(frame, bin, val, increment, mythread) reduction(+:total)
    {
      mythread = omp_get_thread_num();
      // Prevent race conditions by giving each thread its own histogram.
#     pragma omp master
      {
        numthreads = omp_get_num_threads();
        mprintf("\tParallelizing calculation with %i threads\n", numthreads);
        out_thread = new double*[ numthreads ];
        for (nt = 0; nt < numthreads; nt++) {
          out_thread[nt] = new double[ outSize ]; 
          std::fill(out_thread[nt], out_thread[nt] + outSize, 0.0);
        }
      }
#     pragma omp barrier
#     pragma omp for
# endif
      for (frame = 0; frame < inSize; frame++) {
        val = In.Dval(frame);
        increment = Increments[frame];
        total += increment;
        // Apply kernel across histogram
        for (bin = 0; bin < outSize; bin++)
#         ifdef _OPENMP
          out_thread[mythread][bin] +=
#         else
          Out[bin] += 
#         endif
            (increment * (this->*Kernel_)( (Xdim.Coord(bin) - val) / bandwidth_ ));
      }
#ifdef _OPENMP
    } // END parallel block
    // Combine results from each thread histogram into Out
    for (int i = 0; i < numthreads; i++) {
      for (unsigned int j = 0; j < outSize; j++)
        Out[j] += out_thread[i][j];
      delete[] out_thread[i];
    }
    delete[] out_thread;
#endif
  } else {
    // Calculate Kullback-Leibler divergence vs time
    DataSet_1D const& Qdata = static_cast<DataSet_1D const&>( *q_data_ );
    if (inSize != Qdata.Size()) {
      mprintf("Warning: Size of %s (%zu) != size of %s (%zu)\n",
                In.Legend().c_str(), In.Size(), Qdata.Legend().c_str(), Qdata.Size());
      inSize = std::min( inSize, Qdata.Size() );
      mprintf("Warning:  Only using %zu data points.\n", inSize);
    }
    DataSet_double& klOut = static_cast<DataSet_double&>( *kldiv_ );
    std::vector<double> Qhist( Xdim.Bins(), 0.0 ); // Raw Q histogram.
    klOut.Resize( inSize ); // Hold KL div vs time
    double val_p, val_q, KL, norm, xcrd, Pnorm, Qnorm;
    // Loop over input P and Q data
    unsigned int nInvalid = 0;
    for (frame = 0; frame < inSize; frame++) {
      increment = Increments[frame];
      total += increment;
      // Apply kernel across P and Q, calculate KL divergence as we go. 
      val_p = In.Dval(frame);
      val_q = Qdata.Dval(frame);
      KL = 0.0;
      validPoint = 0; // 0 in this context means true
      // NOTE: This normalization is subject to precision loss.
      //       For a calculation between 2 sets of ~850000 frames the max
      //       precision loss is on the order of 0.001, and the avg. deviation
      //       is ~2E-05.
      norm = Xdim.Step() / (total * bandwidth_);
#     ifdef _OPENMP
#     pragma omp parallel private(bin, xcrd, Pnorm, Qnorm, Pzero, Qzero) reduction(+:KL, validPoint)
{
#     pragma omp for
#     endif
      for (bin = 0; bin < outSize; bin++) {
        xcrd = Xdim.Coord(bin);
        Out[bin]   += (increment * (this->*Kernel_)( (xcrd - val_p) / bandwidth_ ));
        Qhist[bin] += (increment * (this->*Kernel_)( (xcrd - val_q) / bandwidth_ ));
        // KL only defined when Q and P are non-zero, or both zero.
        if (validPoint == 0) {
          // Normalize for this frame
          Pnorm = Out[bin] * norm;
          Qnorm = Qhist[bin] * norm;
          Pzero = (Pnorm == 0.0);
          Qzero = (Qnorm == 0.0);
          if (!Pzero && !Qzero)
            KL += ( log( Pnorm / Qnorm ) * Pnorm );
          else if ( Pzero != Qzero )
            validPoint++;
        }
      }
#ifdef _OPENMP
}
#endif
      if (validPoint == 0) {
        klOut[frame] = KL;
      } else {
        //mprintf("Warning:\tKullback-Leibler divergence is undefined for frame %u\n", i+1);
        nInvalid++;
      }
    }
    if (nInvalid > 0)
      mprintf("Warning:\tKullback-Leibler divergence was undefined for %u frames.\n", nInvalid);
  }

  // Normalize
  for (unsigned int j = 0; j < Out.Size(); j++)
    Out[j] /= (total * bandwidth_);

  // Calc free E
  if (calcFreeE_) {
    double minFreeE = 0.0;
    for (unsigned int j = 0; j < Out.Size(); j++) {
      Out[j] = -log( Out[j] );
      if (j == 0)
        minFreeE = Out[j];
      else if (Out[j] < minFreeE)
        minFreeE = Out[j];
    }
    for (unsigned int j = 0; j < Out.Size(); j++)
      Out[j] -= minFreeE;
  }

  return Analysis::OK;
}
