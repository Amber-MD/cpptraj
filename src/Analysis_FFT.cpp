#include <cstring> // memset
#include <cmath> // sqrt
#include "Analysis_FFT.h"
#include "CpptrajStdio.h"
#include "PubFFT.h"

// CONSTRUCTOR
Analysis_FFT::Analysis_FFT() :
  maxsize_(0),
  dt_(0.0),
  f0_(0.0)
{}

// Analysis_FFT::Setup()
/** Usage: fft <sets arg> [out <outfile>] [name <outsetname>] [d] */
int Analysis_FFT::Setup(DataSetList* datasetlist) {
  std::string setname_ = analyzeArgs_.GetStringKey("name");
  outfilename_ = analyzeArgs_.GetStringKey("out");
  dt_ = analyzeArgs_.getKeyDouble("dt",1.0);
  // Select datasets
  input_dsets_ = datasetlist->GetMultipleSets( analyzeArgs_.GetStringNext() );
  if (input_dsets_.empty()) {
    mprinterr("Error: FFT: No data sets selected.\n");
    return 1;
  }
  // If setname is empty generate a default name
  if (setname_.empty())
    setname_ = datasetlist->GenerateDefaultName( "FFT" );
  // Setup output datasets - Also determine max input DataSet size
  int idx = 0;
  maxsize_ = 0;
  if ( input_dsets_.size() == 1 )
    idx = -1; // Only one input set, no need to refer to it by index
  for ( DataSetList::const_iterator DS = input_dsets_.begin(); 
                                    DS != input_dsets_.end(); ++DS) 
  {
    DataSet* dsout = datasetlist->AddSetIdx( DataSet::DOUBLE, setname_, idx++ );
    if (dsout==NULL) return 1;
    if ( (*DS)->Size() > maxsize_ ) maxsize_ = (*DS)->Size();
    output_dsets_.push_back( dsout );
  }

  mprintf("    FFT: Calculating FFT for %i data sets (max size %i):\n", 
          input_dsets_.size(), maxsize_ );
  mprintf("\tTime step: %f\n", dt_);
  if ( !setname_.empty() )
    mprintf("\tSet name: %s\n", setname_.c_str() );
  if ( !outfilename_.empty() )
    mprintf("\tOutfile name: %s\n", outfilename_.c_str());

  return 0;
}

int Analysis_FFT::Analyze() {
  PubFFT pubfft( maxsize_ );
  //PubFFT pubfft;
  //pubfft.SetupFFTforN( maxsize_ );
  int ndata = pubfft.size() * 2; // space for (real + img) per datapoint
  double *data1 = new double[ ndata ];

  double sr = 1.0 / dt_; // 1 / sampling interval, sampling rate (freq)

  std::vector<DataSet*>::iterator dsout = output_dsets_.begin();
  for (DataSetList::const_iterator DS = input_dsets_.begin(); 
                                   DS != input_dsets_.end(); ++DS)
  {
    mprintf("\t\tCalculating FFT for set %s\n", (*DS)->Legend().c_str());
    // Reset data1 so it is padded with zeros
    memset( data1, 0, ndata*sizeof(double) );
    // Place data from DS in real spots in data1
    int datasize =  (*DS)->Size();
    double total_time = dt_ * datasize;
    f0_ = 1.0 / total_time;
    mprintf("\t\t\tDT=%f ps, SR= %f ps^-1, total time=%f ps, f0=%f ps^-1\n",dt_,sr,
            total_time, f0_);
    for (int i = 0; i < datasize; ++i)
      data1[i*2] = (*DS)->Dval(i);
    // DEBUG
    //for (int i = 0; i < pubfft.size()*2; i+=2)
    //  mprintf("\t\t\t%i FFTinR=%f  FFTinI=%f\n",i/2,data1[i],data1[i+1]);
    // Perform FFT
    pubfft.Forward( data1 );
    // Place real data from FFT in output Data
    int i2 = 0;
    for (int i1 = 0; i1 < datasize; ++i1) {
     double magnitude = sqrt(data1[i2]*data1[i2] + data1[i2+1]*data1[i2+1]);
     //mprintf("\t\t\tReal=%f  Img=%f  Mag=%f\n",data1[i2],data1[i2+1],magnitude);
     (*dsout)->Add( i1, &magnitude );
     i2 += 2;
    }
    ++dsout;
  }
  delete[] data1;
  return 0;
}

void Analysis_FFT::Print( DataFileList* datafilelist ) {
  if (!outfilename_.empty()) {
    for (std::vector<DataSet*>::iterator dsout = output_dsets_.begin();
                                         dsout != output_dsets_.end(); ++dsout)
      datafilelist->Add( outfilename_.c_str(), *dsout );
    DataFile* DF = datafilelist->GetDataFile( outfilename_ );
    if (DF != NULL) 
      DF->ProcessArgs("xlabel Freq. xmin 0 xstep " + doubleToString( f0_ ));
  }
}
