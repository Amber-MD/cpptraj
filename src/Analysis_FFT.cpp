#include <cstring> // memset
#include <cmath> // sqrt
#include "Analysis_FFT.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // doubleToString
#include "PubFFT.h"

// CONSTRUCTOR
Analysis_FFT::Analysis_FFT() :
  outfile_(0),
  maxsize_(0),
  dt_(0.0)
{}

void Analysis_FFT::Help() {
  mprintf("\t<dset0> [<dset1> ...] [out <outfile>] [name <outsetname>] [dt <samp_int>]\n");
  mprintf("\tPerform fast-Fourier transformation of data set(s)\n");
}

// Analysis_FFT::Setup()
Analysis::RetType Analysis_FFT::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  std::string setname = analyzeArgs.GetStringKey("name");
  outfile_ = DFLin->AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  dt_ = analyzeArgs.getKeyDouble("dt",1.0);
  // Select datasets from remaining args
  ArgList dsetArgs = analyzeArgs.RemainingArgs();
  for (ArgList::const_iterator dsa = dsetArgs.begin(); dsa != dsetArgs.end(); ++dsa)
    input_dsets_ += datasetlist->GetMultipleSets( *dsa );
  if (input_dsets_.empty()) {
    mprinterr("Error: FFT: No data sets selected.\n");
    return Analysis::ERR;
  }
  // If setname is empty generate a default name
  if (setname.empty())
    setname = datasetlist->GenerateDefaultName( "FFT" );
  // Setup output datasets. Also ensure all DataSets have the same # of points. 
  int idx = 0;
  maxsize_ = 0;
  if ( input_dsets_.size() == 1 )
    idx = -1; // Only one input set, no need to refer to it by index
  for ( DataSetList::const_iterator DS = input_dsets_.begin(); 
                                    DS != input_dsets_.end(); ++DS) 
  {
    // Check for empty set
    if ( (*DS)->Empty() ) {
      mprintf("Warning: FFT: Set %s is empty, skipping.\n", (*DS)->Legend().c_str() );
      continue;
    }
    if ( maxsize_ == 0 ) 
      maxsize_ = (*DS)->Size();
    else if ( (*DS)->Size() != maxsize_ ) {
      mprintf("Warning: FFT: Set %s does not have same size (%i) as initial set (%i). Skipping.\n",
              (*DS)->Legend().c_str(), (*DS)->Size(), maxsize_ );
      continue;
    }
    DataSet* dsout = datasetlist->AddSetIdx( DataSet::DOUBLE, setname, idx++ );
    if (dsout==0) return Analysis::ERR;
    dsout->SetLegend( (*DS)->Legend() );
    output_dsets_.push_back( dsout );
    if (outfile_ != 0) outfile_->AddSet( dsout );
  }

  mprintf("    FFT: Calculating FFT for %i data sets (of size %i):\n", 
          input_dsets_.size(), maxsize_ );
  mprintf("\tTime step: %f\n", dt_);
  if ( !setname.empty() )
    mprintf("\tSet name: %s\n", setname.c_str() );
  if ( outfile_ != 0 )
    mprintf("\tOutfile name: %s\n", outfile_->DataFilename().base());

  return Analysis::OK;
}

/** Calculate FFT for input DataSets. FFT magnitude is reported. Magnitude
  * is normalized by N / 2. Only data up to the Nyquist frequency is used. 
  */
//TODO: Deal with vectors
Analysis::RetType Analysis_FFT::Analyze() {
  //PubFFT pubfft( maxsize_ );
  PubFFT pubfft;
  pubfft.SetupFFTforN( maxsize_ );
  mprintf("DEBUG: FFT size is %i\n",pubfft.size());
  int ndata = pubfft.size() * 2; // space for (real + img) per datapoint
  double *data1 = new double[ ndata ];

  double sr = 1.0 / dt_;              // 1 / sampling interval, sampling rate (freq)
  double fnyquist = sr / 2.0;         // Nyquist frequency
  double total_time = dt_ * maxsize_; // Total time (fundamental period)
  double f0 = 1.0 / total_time;       // Fundamental frequency (first harmonic)
  if (outfile_ != 0) 
    outfile_->ProcessArgs("xlabel Freq. xmin 0 xstep " + doubleToString( f0 ));
  double norm = (double)(maxsize_ / 2);

  std::vector<DataSet*>::iterator dsout = output_dsets_.begin();
  for (DataSetList::const_iterator DS = input_dsets_.begin(); 
                                   DS != input_dsets_.end(); ++DS)
  {
    mprintf("\t\tCalculating FFT for set %s\n", (*DS)->Legend().c_str());
    // Reset data1 so it is padded with zeros
    memset( data1, 0, ndata*sizeof(double) );
    // Place data from DS in real spots in data1
    int datasize =  (*DS)->Size();
    mprintf("\t\t\tDT=%f ps, SR= %f ps^-1, FC= %f ps^-1, total time=%f ps, f0=%f ps^-1\n",
            dt_, sr, fnyquist, total_time, f0);
    for (int i = 0; i < datasize; ++i)
      data1[i*2] = (*DS)->Dval(i);
    // DEBUG
    //for (int i = 0; i < pubfft.size()*2; i+=2)
    //  mprintf("\t\t\t%i FFTinR=%f  FFTinI=%f\n",i/2,data1[i],data1[i+1]);
    // Perform FFT
    pubfft.Forward( data1 );
    // Place real data from FFT in output Data up to the Nyquist frequency
    int i2 = 0;
    for (int i1 = 0; i1 < datasize; ++i1) {
     double freq = i1 * f0;
     if (freq > fnyquist) break;
     double magnitude = sqrt(data1[i2]*data1[i2] + data1[i2+1]*data1[i2+1]);
     magnitude /= norm;
     //mprintf("\t\t\tReal=%f  Img=%f  Mag=%f\n",data1[i2],data1[i2+1],magnitude);
     (*dsout)->Add( i1, &magnitude );
     i2 += 2;
    }
    ++dsout;
  }
  delete[] data1;
  return Analysis::OK;
}
