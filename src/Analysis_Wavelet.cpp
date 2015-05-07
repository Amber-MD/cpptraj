#include <cmath> // pow
#include "Analysis_Wavelet.h"
#include "CpptrajStdio.h"
#include "Matrix.h"
#include "DistRoutines.h"
#include "Constants.h"
#include "PubFFT.h"

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

static inline void PrintComplex(const char* title, ComplexArray const& C) {
  if (title != 0) mprintf("DEBUG: %s:", title);
  for (ComplexArray::iterator cval = C.begin(); cval != C.end(); cval += 2)
    mprintf(" (%g,%g)", *cval, *(cval+1));
  mprintf("\n");
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
    //const double* curr = currentFrame.XYZ(0);
    //const double* last = lastFrame.XYZ(0);
    //mprintf("DEBUG: last frame: %g %g %g  this frame: %g %g %g\n",
    //        last[0], last[1], last[2],
    //        curr[0], curr[1], curr[2]);
    int idx = frm; // Position in distance matrix; start at column 'frame'
    for (int at = 0; at != natoms; at++, idx += nframes) {
      // Distance of atom at in current frame from last frame.
      d_matrix[idx] = sqrt(DIST2_NoImage( currentFrame.XYZ(at), lastFrame.XYZ(at) ));
      //mprintf("DEBUG: Col(frm) %i Row(at) %i Idx %i distance: %g\n",
      //        frm, at, idx, d_matrix[idx]);
    }
    //lastFrame = currentFrame; // TODO: Re-enable?
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

  // Precompute some factors for calculating scaled wavelets.
  double one_over_sqrt_N = 1.0 / sqrt( nframes );
  std::vector<int> arrayK( nframes );
  arrayK[0] = -1 * (nframes/2);
  for (int i = 1; i != nframes; i++)
    arrayK[i] = arrayK[i-1] + 1;
  mprintf("DEBUG: K:");
  for (std::vector<int>::const_iterator kval = arrayK.begin(); kval != arrayK.end(); ++kval)
    mprintf(" %i", *kval);
  mprintf("\n");

  // Step 2 - Get FFT of wavelet for each scale.
  PubFFT pubfft;
  pubfft.SetupFFTforN( nframes );
  mprintf("\tMemory required for scaled wavelet array: %.2f MB\n",
          (double)(2 * nframes * nb_ * sizeof(double)) / (1024 * 1024));
  typedef std::vector<ComplexArray> WaveletArray;
  WaveletArray FFT_of_Scaled_Wavelets;
  FFT_of_Scaled_Wavelets.reserve( nb_ );
  typedef std::vector<double> Darray;
  Darray scaleVector;
  scaleVector.reserve( nb_ );
  Darray MIN( nb_ * 2 );
  for (int iscale = 0; iscale != nb_; iscale++)
  {
    // Calculate and store scale factor.
    scaleVector.push_back( S0_ * pow(2.0, iscale * ds_) );
    // Populate MIN array
    MIN[iscale    ] = (0.00647*pow((correction_*scaleVector.back()),1.41344)+19.7527)*chival_;
    MIN[iscale+nb_] = correction_*scaleVector.back();
    // Calculate scalved wavelet
    ComplexArray scaledWavelet;
    switch (wavelet_type_) {
      case W_MORLET: scaledWavelet = F_Morlet(arrayK, scaleVector.back()); break;
      case W_PAUL  : scaledWavelet = F_Paul(arrayK, scaleVector.back()); break;
      case W_NONE  : return Analysis::ERR; // Sanity check
    }
    PrintComplex("wavelet_before_fft", scaledWavelet);
    // Perform FFT
    pubfft.Forward( scaledWavelet );
    // Normalize
    scaledWavelet.Normalize( one_over_sqrt_N );
    PrintComplex("wavelet_after_fft", scaledWavelet);
    FFT_of_Scaled_Wavelets.push_back( scaledWavelet );
  }
  mprintf("DEBUG: Scaling factors:");
  for (Darray::const_iterator sval = scaleVector.begin(); sval != scaleVector.end(); ++sval)
    mprintf(" %g", *sval);
  mprintf("\n");
  mprintf("DEBUG: MIN:");
  for (int i = 0; i != nb_; i++)
    mprintf(" %g", MIN[i]);
  mprintf("\n");

  // Step 3 - For each atom, calculate the convolution of scaled wavelets
  //          with rows (atom distance vs frame) via dot product of the 
  //          frequency domains, i.e. Fourier-transformed, followed by an
  //          inverse FT.
  Matrix<double> MAX;
  MAX.resize( nframes, natoms );
  Darray output( nframes ); // Scratch space
  for (int at = 0; at != natoms; at++) {
    ComplexArray AtomSignal( nframes ); // Initializes to zero
    // Calculate the distance variance for this atom and populate the array.
    int midx = at * nframes; // Index into d_matrix
    int cidx = 0;            // Index into AtomSignal
    double d_avg = 0.0;
    double d_var = 0.0;
    for (int frm = 0; frm != nframes; frm++, cidx += 2, midx++) {
      d_avg += d_matrix[midx];
      d_var += (d_matrix[midx] * d_matrix[midx]);
      AtomSignal[cidx] = d_matrix[midx];
    }
    d_var = (d_var - ((d_avg * d_avg) / (double)nframes)) / ((double)(nframes - 1));
    mprintf("VARIANCE: %g\n", d_var);
    double var_norm = 1.0 / d_var;
    // Calculate FT of atom signal
    pubfft.Forward( AtomSignal );
    PrintComplex("AtomSignal", AtomSignal);
    // Normalize
    AtomSignal.Normalize( one_over_sqrt_N );
    // Calculate dot product of atom signal with each scaled FT wavelet
    for (int iscale = 0; iscale != nb_; iscale++) {
      ComplexArray dot = AtomSignal.TimesComplexConj( FFT_of_Scaled_Wavelets[iscale] );
      // Inverse FT of dot product
      pubfft.Back( dot );
      PrintComplex("InverseFT_Dot", dot);
      // Chi-squared testing
      midx = at * nframes;
      cidx = 0;
      for (int frm = 0; frm != nframes; frm++, cidx += 2, midx++) {
        output[frm] = (dot[cidx]*dot[cidx] + dot[cidx+1]*dot[cidx+1]) * var_norm;
        if (output[frm] < MIN[iscale])
          output[frm] = 0.0;
        if (output[frm] > MAX[midx]) {
          MAX[midx] = output[frm];
          //Indices[midx] = iscale
        }
      }
      mprintf("DEBUG: AbsoluteValue:");
      for (Darray::const_iterator dval = output.begin(); dval != output.end(); ++dval)
        mprintf(" %g", *dval);
      mprintf("\n");
    } // END loop over scales
  } // END loop over atoms
      
  return Analysis::OK;
}
