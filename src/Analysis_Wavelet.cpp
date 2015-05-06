#include <cmath> // pow
#include "Analysis_Wavelet.h"
#include "CpptrajStdio.h"
#include "Matrix.h"
#include "DistRoutines.h"
#include "Constants.h"

// CONSTRUCTOR
Analysis_Wavelet::Analysis_Wavelet() :
  coords_(0),
  S0_(0.0),
  ds_(0.0),
  correction_(0.0),
  chival_(0.0),
  wavelet_type_(W_NONE),
  nb_(0)
{}

// Wavelet functions: takes array of precomputed prefactors and scaling factor.
ComplexArray Analysis_Wavelet::F_Morlet(std::vector<int> const& K, double S) const {
  unsigned int N = K.size();
  ComplexArray out( N );
  double scale = sqrt(S);

  for (unsigned int i = 0; i != N; i++) {
    double t = (double)K[i] / S;
    unsigned int idx = ((i+N/2+1)%N) * 2;
    out[idx  ]=((double) pow(Constants::PI,-0.25)*exp(-1*t*t/2)*cos(2*Constants::PI*t)) / scale;
    out[idx+1]=((double) pow(Constants::PI,-0.25)*exp(-1*t*t/2)*sin(2*Constants::PI*t)) / scale;
  }
  return out;
}

ComplexArray Analysis_Wavelet::F_Paul(std::vector<int> const& K, double S) const {
  static const double q_paul = 8*sqrt(2/(35*Constants::PI));
  unsigned int N = K.size();
  ComplexArray out( N );

  for (unsigned int i = 0; i != N; i++) {
    double t = (double)K[i] / S;
    unsigned int idx = ((i+N/2+1)%N) * 2;
    double scale = pow((1+t*t),5);
    out[idx  ]=((double) q_paul * (1-10*pow(t,2)+5*pow(t,4))  ) / scale;
    out[idx+1]=((double) q_paul * (5*t-10*pow(t,3)+5*pow(t,5))) / scale;
  }
  return out;
}

const Analysis_Wavelet::WaveletToken Analysis_Wavelet::Tokens_[] = {
  {"morlet", "Morlet"}, // W_MORLET 
  {"paul",   "Paul"}    // W_PAUL
};

/// Provide keywords
void Analysis_Wavelet::Help() {
  mprintf("\t[crdset <set name>] nb <n scaling vals> [s0 <s0>] [ds <ds>]\n"
          "\t[correction <correction>] [chival <chival>] [type <wavelet>]\n");
}

// Analysis_Wavelet::Setup
Analysis::RetType Analysis_Wavelet::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  // Attempt to get COORDS DataSet from DataSetList. If none specified the
  // default COORDS set will be used.
  std::string setname = analyzeArgs.GetStringKey("crdset");
  coords_ = (DataSet_Coords*)datasetlist->FindCoordsSet( setname );
  if (coords_ == 0) {
    mprinterr("Error: Could not locate COORDS set corresponding to %s\n", setname.c_str());
    return Analysis::ERR;
  }
  // Get keywords
  // TODO: Check defaults
  nb_ = analyzeArgs.getKeyInt("nb", 0); // FIXME: Should be more descriptive? nscale?
  if (nb_ < 1) {
    mprinterr("Error: Scaling number must be > 0\n");
    return Analysis::ERR;
  }
  S0_ = analyzeArgs.getKeyDouble("s0", 0.2);
  ds_ = analyzeArgs.getKeyDouble("ds", 1.0/3.0);
  correction_ = analyzeArgs.getKeyDouble("correction", 1.01);
  chival_ = analyzeArgs.getKeyDouble("chival", 0.2231);
  // Wavelet type: default to Morlet
  std::string wavelet_name = analyzeArgs.GetStringKey("type");
  if (wavelet_name.empty())
    wavelet_type_ = W_MORLET;
  else {
    wavelet_type_ = W_NONE;
    for (int itoken = 0; itoken != (int)W_NONE; itoken++)
      if (wavelet_name.compare(Tokens_[itoken].key_) == 0) {
        wavelet_type_ = (WaveletType)itoken;
        break;
      } 
    if (wavelet_type_ == W_NONE) {
      mprinterr("Error: Unrecognized wavelet type: %s\n", wavelet_name.c_str());
      return Analysis::ERR;
    }
  }
  // Atom mask
  mask_.SetMaskString( analyzeArgs.GetMaskNext() );

  mprintf("    WAVELET: Using COORDS set '%s', wavelet type %s\n",
          coords_->legend(), Tokens_[wavelet_type_].description_);
  mprintf("\tCalculating for atoms in mask '%s'\n", mask_.MaskString());
  mprintf("\tScaling wavelet %i times starting from %g with delta of %g\n", 
          nb_, S0_, ds_);
  mprintf("\tCorrection: %g\n", correction_);
  mprintf("\tChiVal:     %g\n", chival_);


  return Analysis::OK;
}

// Analysis_Wavelet::Analyze()
Analysis::RetType Analysis_Wavelet::Analyze() {
  // Step 1 - Create a matrix that is #atoms rows by #frames - 1 cols,
  //          where matrix(frame, atom) is the distance that atom has
  //          travelled from the previous frame.
  // TODO: Implement this in Action_Matrix()?
  mprintf("    WAVELET:\n");
  // First set up atom mask.
  if (coords_->Top().SetupIntegerMask( mask_ )) return Analysis::ERR;
  mask_.MaskInfo();
  int natoms = mask_.Nselected();
  int nframes = (int)coords_->Size();
  if (natoms < 1 || nframes < 2) {
    mprinterr("Error: Not enough frames (%i) or atoms (%i) in '%s'\n",
              nframes, natoms, coords_->legend());
    return Analysis::ERR;
  }
  mprintf("\t%i frames, %i atoms, distance matrix will require %.2f MB\n",
          nframes, natoms, (double)(nframes * natoms * sizeof(double)) / (1024 * 1024));
  Matrix<double> d_matrix;
  d_matrix.resize(nframes, natoms);
  // Get initial frame.
  Frame currentFrame, lastFrame;
  currentFrame.SetupFrameFromMask( mask_, coords_->Top().Atoms() );
  lastFrame = currentFrame;
  coords_->GetFrame( 0, lastFrame, mask_ );
  // Iterate over frames
  for (int frm = 1; frm != nframes; frm++) {
    coords_->GetFrame( frm, currentFrame, mask_ );
    int idx = frm; // Position in distance matrix; start at column 'frame'
    for (int at = 0; at != natoms; at++, idx += nframes)
      // Distance of atom at in current frame from last frame.
      d_matrix[idx] = DIST2_NoImage( currentFrame.XYZ(at), lastFrame.XYZ(at) );
    lastFrame = currentFrame;
  }
  // DEBUG: Write matrix to file.
  CpptrajFile dmatrixOut;
  dmatrixOut.OpenWrite("dmatrix.dat");
  Matrix<double>::iterator mval = d_matrix.begin();
  for (int row = 0; row != natoms; row++) {
    for (int col = 0; col != nframes; col++)
      dmatrixOut.Printf("%g ", *(mval++));
    dmatrixOut.Printf("\n");
  }
  dmatrixOut.CloseFile();

  // Calculate scale factors
  typedef std::vector<double> Darray;
  Darray scaleVector;
  scaleVector.reserve( nb_ );
  for (int iscale = 0; iscale != nb_; iscale++)
    scaleVector.push_back( S0_ * pow(2.0, iscale * ds_) );
  mprintf("DEBUG: Scaling factors:");
  for (Darray::const_iterator sval = scaleVector.begin(); sval != scaleVector.end(); ++sval)
    mprintf(" %g", *sval);
  mprintf("\n");

  // Precompute some factors for calculating scaled wavelets.
  std::vector<int> arrayK( nframes );
  arrayK[0] = -1 * (nframes/2);
  for (int i = 1; i != nframes; i++)
    arrayK[i] = arrayK[i-1] + 1;
  mprintf("DEBUG: K:");
  for (std::vector<int>::const_iterator kval = arrayK.begin(); kval != arrayK.end(); ++kval)
    mprintf(" %i", *kval);
  mprintf("\n");

  // Step 2 - Get FFT of wavelet for each scale
  mprintf("\tMemory required for scaled wavelet array: %.2f MB\n",
          (double)(2 * nframes * nb_ * sizeof(double)) / (1024 * 1024));
  typedef std::vector<ComplexArray> WaveletArray;
  WaveletArray FFT_of_Scaled_Wavelets;
  FFT_of_Scaled_Wavelets.reserve( nb_ );
  for (Darray::const_iterator sval = scaleVector.begin(); sval != scaleVector.end(); ++sval)
  {
    ComplexArray scaledWavelet;
    switch (wavelet_type_) {
      case W_MORLET: scaledWavelet = F_Morlet(arrayK, *sval); break;
      case W_PAUL  : scaledWavelet = F_Paul(arrayK, *sval); break;
      case W_NONE  : return Analysis::ERR; // Sanity check
    }
  }
  //  wavelet_fft_scale.push_back( ScaledWaveletFFT(nframes, *sval) );
  
      

  return Analysis::OK;
}
