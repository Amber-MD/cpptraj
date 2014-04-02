#include <cmath> // pow
#include "Analysis_KDE.h"
#include "CpptrajStdio.h"
#include "DataSet_double.h"
#include "Constants.h" // TWOPI, GASK_KCAL
#include <limits> // Minimum double val for checking zero
#ifdef _OPENMP
#  include "omp.h"
#endif

// CONSTRUCTOR
Analysis_KDE::Analysis_KDE() : 
  data_(0), q_data_(0), bandwidth_(0.0), output_(0), kldiv_(0),
  amddata_(0), calcFreeE_(false), Temp_(0.0),
  Kernel_(&Analysis_KDE::GaussianKernel) {}

void Analysis_KDE::Help() {
  mprintf("\t<dataset> [bandwidth <bw>] [out <file>] [name <dsname>]\n"
          "\t[min <min>] [max <max] [step <step>] [bins <bins>] [free]\n"
          "\t[kldiv <dsname2> [klout <outfile>]] [amd <amdboost_data>]\n"
          "  Histogram 1D data set using a kernel density estimator.\n");
}

Analysis::RetType Analysis_KDE::Setup(DataSet_1D* dsIn, std::string const& histname,
                                       std::string const& outfilenameIn,
                                       bool minArgSetIn, double minIn,
                                       bool maxArgSetIn, double maxIn,
                                       double stepIn, int binsIn, double tempIn,
                                       DataSetList& datasetlist, DataFileList& DFLin)
{
  if (dsIn == 0) return Analysis::ERR;
  data_ = dsIn;
  q_data_ = 0;
  kldiv_ = 0;
  amddata_ = 0;
  bandwidth_ = -1.0;
  Dimension Xdim;
  if (minArgSetIn)
    Xdim.SetMin( minIn );
  if (maxArgSetIn)
    Xdim.SetMax( maxIn );
  Xdim.SetStep( stepIn );
  Xdim.SetBins( binsIn );
  if (Xdim.Step() < 0.0 && Xdim.Bins() < 0) {
    mprinterr("Error: Must set either bins or step.\n");
    return Analysis::ERR;
  }
  Temp_ = tempIn;
  if (Temp_ != -1.0)
    calcFreeE_ = true;
  else
    calcFreeE_ = false;
  std::string setname = histname;
  if (histname.empty()) {
    if (calcFreeE_)
      setname="FreeE_";
    else
      setname="KDE_";
  }
  setname += dsIn->Legend();
  DataFile* outfile = DFLin.AddDataFile( outfilenameIn );
  output_ = datasetlist.AddSetIdxAspect(DataSet::DOUBLE, setname, dsIn->Idx(), dsIn->Aspect());
  if (output_ == 0) return Analysis::ERR;
  output_->SetDim(Dimension::X, Xdim);
  if (outfile != 0) outfile->AddSet( output_ );
  return Analysis::OK;
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
  Temp_ = analyzeArgs.getKeyDouble("free",-1.0);
  if (Temp_!=-1.0) 
    calcFreeE_ = true;
  else
    calcFreeE_ = false;
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
  if (output_ == 0) return Analysis::ERR;
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
  if (calcFreeE_)
    mprintf("\tFree energy in kcal/mol will be calculated from bin populations at %f K.\n",Temp_);
  return Analysis::OK;
}

const double Analysis_KDE::ONE_OVER_ROOT_TWOPI = 1.0 / sqrt( Constants::TWOPI );

double Analysis_KDE::GaussianKernel(double u) const {
  return ( ONE_OVER_ROOT_TWOPI * exp( -0.5 * u * u ) );
}

// Analysis_KDE::Analyze()
Analysis::RetType Analysis_KDE::Analyze() {
  Dimension& Xdim = output_->Dim(0);
  DataSet_1D const& Pdata = static_cast<DataSet_1D const&>( *data_ );
  size_t inSize = Pdata.Size();
  // Set output set dimensions from input set if necessary.
  if (!Xdim.MinIsSet())
    Xdim.SetMin( Pdata.Min() );
  if (!Xdim.MaxIsSet())
    Xdim.SetMax( Pdata.Max() );
  if (Xdim.CalcBinsOrStep()) return Analysis::ERR;
  Xdim.PrintDim();
  // TODO: When calculating KL divergence check q_data_ dimension.

  // Allocate output set
  DataSet_double& P_hist = static_cast<DataSet_double&>( *output_ );
  P_hist.Resize( Xdim.Bins() );
  size_t outSize = P_hist.Size();

  // Estimate bandwidth from normal distribution approximation if necessary.
  if (bandwidth_ < 0.0) {
    double stdev;
    Pdata.Avg( stdev );
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

# ifdef _OPENMP
  // Prevent race conditions by giving each thread its own histogram
  int numthreads, mythread;
  double **P_thread, **Q_thread;
# pragma omp parallel
  {
#   pragma omp master
    {
      numthreads = omp_get_num_threads();
    }
  }
  mprintf("\tParallelizing calculation with %i threads\n", numthreads);
  P_thread = new double*[ numthreads ];
  for (int nt = 0; nt < numthreads; nt++) {
    P_thread[nt] = new double[ outSize ];
    std::fill(P_thread[nt], P_thread[nt] + outSize, 0.0);
  }
  if (q_data_ != 0) {
    Q_thread = new double*[ numthreads ];
    for (int nt = 0; nt < numthreads; nt++) {
     Q_thread[nt] = new double[ outSize ];
      std::fill(Q_thread[nt], Q_thread[nt] + outSize, 0.0);
    }
  }
# endif
 
  unsigned int frame, bin;
  double increment;
  double total = 0.0;
  if (q_data_ == 0) {
    double val;
    // Calculate KDE, loop over input data
#   ifdef _OPENMP
#   pragma omp parallel private(frame, bin, val, increment, mythread) reduction(+:total)
    {
      mythread = omp_get_thread_num();
#     pragma omp for
#   endif
      for (frame = 0; frame < inSize; frame++) {
        val = Pdata.Dval(frame);
        increment = Increments[frame];
        total += increment;
        // Apply kernel across histogram
        for (bin = 0; bin < outSize; bin++)
#         ifdef _OPENMP
          P_thread[mythread][bin] +=
#         else
          P_hist[bin] += 
#         endif
            (increment * (this->*Kernel_)( (Xdim.Coord(bin) - val) / bandwidth_ ));
      }
#   ifdef _OPENMP
    } // END parallel block
#   endif
  } else {
    // Calculate Kullback-Leibler divergence vs time
    DataSet_1D const& Qdata = static_cast<DataSet_1D const&>( *q_data_ );
    if (inSize != Qdata.Size()) {
      mprintf("Warning: Size of %s (%zu) != size of %s (%zu)\n",
                Pdata.Legend().c_str(), Pdata.Size(), Qdata.Legend().c_str(), Qdata.Size());
      inSize = std::min( inSize, Qdata.Size() );
      mprintf("Warning:  Only using %zu data points.\n", inSize);
    }
    DataSet_double& klOut = static_cast<DataSet_double&>( *kldiv_ );
    std::vector<double> Q_hist( Xdim.Bins(), 0.0 ); // Raw Q histogram.
    klOut.Resize( inSize ); // Hold KL div vs time
    double val_p, val_q, KL, xcrd, Pnorm, Qnorm, normP, normQ, Pbin, Qbin;
    bool Pzero, Qzero, validPoint;
    // Loop over input P and Q data
    unsigned int nInvalid = 0;
#   ifdef _OPENMP
#   pragma omp parallel private(frame, bin, validPoint, increment, val_p, val_q, KL, xcrd, Pnorm, Qnorm, Pzero, Qzero, normP, normQ, Pbin, Qbin) reduction(+:total,nInvalid)
    {
      mythread = omp_get_thread_num();
#     pragma omp for
#   endif
      for (frame = 0; frame < inSize; frame++) {
        increment = Increments[frame];
        total += increment;
        // Apply kernel across P and Q, calculate KL divergence as we go. 
        val_p = Pdata.Dval(frame);
        val_q = Qdata.Dval(frame);
        normP = 0.0;
        normQ = 0.0;
        validPoint = true;
        for (bin = 0; bin < outSize; bin++) {
          xcrd = Xdim.Coord(bin);
          Pbin = (increment * (this->*Kernel_)( (xcrd - val_p) / bandwidth_ ));
          Qbin = (increment * (this->*Kernel_)( (xcrd - val_q) / bandwidth_ ));
          normP += Pbin;
          normQ += Qbin;
#         ifdef _OPENMP
          P_thread[mythread][bin] += Pbin;
          Q_thread[mythread][bin] += Qbin;
#         else
          P_hist[bin] += Pbin;
          Q_hist[bin] += Qbin;
#         endif
        }
        normP = 1.0 / normP;
        normQ = 1.0 / normQ;
        KL = 0.0;
        for (bin = 0; bin < outSize; bin++) {
          // KL only defined when Q and P are non-zero, or both zero.
          // Normalize for this frame
#         ifdef _OPENMP
          Pnorm = P_thread[mythread][bin] * normP;
          Qnorm = Q_thread[mythread][bin] * normQ;
#         else
          Pnorm = P_hist[bin] * normP;
          Qnorm = Q_hist[bin] * normQ;
#         endif
          Pzero = (Pnorm <= std::numeric_limits<double>::min());
          Qzero = (Qnorm <= std::numeric_limits<double>::min());
          if (!Pzero && !Qzero)
            KL += ( log( Pnorm / Qnorm ) * Pnorm );
          else if ( Pzero != Qzero ) {
            validPoint = false;
            break;
          }
        }
        if (validPoint) {
          klOut[frame] = KL;
        } else {
          //mprintf("Warning:\tKullback-Leibler divergence is undefined for frame %u\n", frame+1);
          nInvalid++;
        }
      } // END KL divergence calc loop over frames
#   ifdef _OPENMP
    }
#   endif 
    if (nInvalid > 0)
      mprintf("Warning:\tKullback-Leibler divergence was undefined for %u frames.\n", nInvalid);
  }
# ifdef _OPENMP
  // Combine results from each thread histogram into P_hist
  for (int i = 0; i < numthreads; i++) {
    for (unsigned int j = 0; j < outSize; j++)
      P_hist[j] += P_thread[i][j];
    delete[] P_thread[i];
  }
  delete[] P_thread;
  // Clean up Q thread memory
  if (q_data_ != 0) {
    for (int i = 0; i < numthreads; i++)
      delete[] Q_thread[i];
    delete[] Q_thread;
  }
# endif

  // Normalize
  for (unsigned int j = 0; j < P_hist.Size(); j++)
    P_hist[j] /= (total * bandwidth_);

  // Calc free E
  if (calcFreeE_) {
    double KT = (-Constants::GASK_KCAL * Temp_);
    double minFreeE = 0.0;
    for (unsigned int j = 0; j < P_hist.Size(); j++) {
      P_hist[j] = log( P_hist[j] ) * KT;
      if (j == 0)
        minFreeE = P_hist[j];
      else if (P_hist[j] < minFreeE)
        minFreeE = P_hist[j];
    }
    for (unsigned int j = 0; j < P_hist.Size(); j++)
      P_hist[j] -= minFreeE;
  }

  return Analysis::OK;
}
