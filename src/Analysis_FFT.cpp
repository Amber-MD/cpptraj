#include <cstring> // memset
#include "Analysis_FFT.h"
#include "CpptrajStdio.h"
#include "PubFFT.h"

// CONSTRUCTOR
Analysis_FFT::Analysis_FFT() :
  maxsize_(0),
  powerspectrum_(true)
{}

// Analysis_FFT::Setup()
int Analysis_FFT::Setup(DataSetList* datasetlist) {
  std::string setname_ = analyzeArgs_.GetStringKey("name");
  outfilename_ = analyzeArgs_.GetStringKey("out");
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
    DataSet* dsout = datasetlist->AddSetIdx( DataSet::DOUBLE, setname_, idx );
    if (dsout==NULL) return 1;
    dsout->SetLegend( (*DS)->Legend() );
    if ( (*DS)->Size() > maxsize_ ) maxsize_ = (*DS)->Size();
    output_dsets_.push_back( dsout );
    if (!powerspectrum_) {
      // If not calculating power spectrum create space to hold imaginary components
      DataSet* imout = datasetlist->AddSetIdxAspect( DataSet::DOUBLE, setname_, idx, "img" );
      img_dsets_.push_back( imout );
    }
    ++idx;
  }

  if (powerspectrum_) 
    mprintf("    FFT: Calculating FFT power spectrum");
  else 
    mprintf("    FFT: Calculating FFT");
  mprintf(" for %i data sets (max size %i):\n", input_dsets_.size(), maxsize_ );
  if ( !setname_.empty() )
    mprintf("\tSet name: %s\n", setname_.c_str() );
  if ( !outfilename_.empty() )
    mprintf("\tOutfile name: %s\n", outfilename_.c_str());

  return 0;
}

int Analysis_FFT::Analyze() {
  PubFFT pubfft( maxsize_ );
  int ndata = pubfft.size() * 2; // space for (real + img) per datapoint
  double *data1 = new double[ ndata ];

  std::vector<DataSet*>::iterator dsout = output_dsets_.begin();
  std::vector<DataSet*>::iterator imout = img_dsets_.begin();
  for (DataSetList::const_iterator DS = input_dsets_.begin(); 
                                   DS != input_dsets_.end(); ++DS)
  {
    mprintf("\t\tCalculating FFT for set %s\n", (*DS)->Legend().c_str());
    // Reset data1 so it is padded with zeros
    memset( data1, 0, ndata*sizeof(double) );
    // Place data from DS in real spots in data1
    int datasize =  (*DS)->Size();
    for (int i = 0; i < datasize; ++i)
      data1[i*2] = (*DS)->Dval(i);
    // Perform FFT
    pubfft.Forward( data1 );
    // Place real data from FFT in output Data
    int i2 = 0;
    if (powerspectrum_) {
      for (int i1 = 0; i1 < datasize; ++i1) {
        data1[i2  ] /= (double)datasize;
        data1[i2+1] /= (double)datasize;
        double dval = data1[i2] * data1[i2] + data1[i2+1] * data1[i2+1];
        (*dsout)->Add( i1, &dval );
        i2 += 2;
      }
    } else {
      for (int i1 = 0; i1 < datasize; ++i1) {
        data1[i2  ] /= (double)datasize;
        (*dsout)->Add( i1, data1 + (i2++) );
        data1[i2  ] /= (double)datasize;
        (*imout)->Add( i1, data1 + (i2++) );
      }
    }
    ++dsout;
    ++imout;
  }
  delete[] data1;
  return 0;
}

void Analysis_FFT::Print( DataFileList* datafilelist ) {
  if (!outfilename_.empty()) {
    std::vector<DataSet*>::iterator imout = img_dsets_.begin();
    for (std::vector<DataSet*>::iterator dsout = output_dsets_.begin();
                                         dsout != output_dsets_.end(); ++dsout)
    {
      datafilelist->Add( outfilename_.c_str(), *dsout );
      if (!powerspectrum_)
        datafilelist->Add( outfilename_.c_str(), *(imout++) );
    }
    //DataFile* DF = datafilelist->GetDataFile( outfilename_.c_str());
    //if (DF != NULL) 
    //  DF->ProcessArgs("xlabel DataSets");
  }
}
