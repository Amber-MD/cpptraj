#include <cmath> // pow
#include <algorithm> // sort, unique
#include "Analysis_Wavelet.h"
#include "CpptrajStdio.h"
#include "Matrix.h"
#include "DistRoutines.h"
#include "Constants.h"
#include "PubFFT.h"
#include "StringRoutines.h"
#include "ProgressBar.h"
#ifdef _OPENMP
# include <omp.h>
#endif

// CONSTRUCTOR
Analysis_Wavelet::Analysis_Wavelet() :
  coords_(0),
  S0_(0.0),
  ds_(0.0),
  correction_(0.0),
  chival_(0.0),
  wavelet_type_(W_NONE),
  nb_(0),
  clustermap_(0),
  epsilon_(0.0),
  epsilon2_(0.0),
  Avg_(0.0),
  minPoints_(0),
  nClusters_(0)
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
void Analysis_Wavelet::Help() const {
  mprintf("\t[crdset <set name>] nb <n scaling vals> [s0 <s0>] [ds <ds>]\n"
          "\t[correction <correction>] [chival <chival>] [type <wavelet>]\n"
          "\t[out <filename>] [name <setname>]\n"
          "    <wavelet>: morlet, paul\n");
}

// Analysis_Wavelet::Setup
Analysis::RetType Analysis_Wavelet::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  // Attempt to get COORDS DataSet from DataSetList. If none specified the
  // default COORDS set will be used.
  std::string setname = analyzeArgs.GetStringKey("crdset");
  coords_ = (DataSet_Coords*)setup.DSL().FindCoordsSet( setname );
  if (coords_ == 0) {
    mprinterr("Error: Could not locate COORDS set corresponding to %s\n", setname.c_str());
    return Analysis::ERR;
  }
  // Get keywords
  DataFile* outfile = setup.DFL().AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs );
  setname = analyzeArgs.GetStringKey("name");
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
  // Wavelet map clustering
  doClustering_ = false;
  DataFile* clustermapout = 0;
  DataFile* clusterout = 0;
  if (analyzeArgs.hasKey("cluster")) {
    doClustering_ = true;
    minPoints_ = analyzeArgs.getKeyInt("minpoints", 4);
    epsilon_ = analyzeArgs.getKeyDouble("epsilon", 10.0);
    epsilon2_ = epsilon_ * epsilon_;
    clustermapout = setup.DFL().AddDataFile( analyzeArgs.GetStringKey("clustermapout"),
                                             analyzeArgs );
    clusterout = setup.DFL().AddDataFile( analyzeArgs.GetStringKey("clusterout"),
                                          analyzeArgs );
#   ifdef _OPENMP
    int numthreads = 0;
#   pragma omp parallel
    {
    mythread_ = omp_get_thread_num();
    if (mythread_ == 0)
      numthreads = omp_get_num_threads();
    }
    mprintf("\tParallelizing calc with %i threads.\n", numthreads);
    thread_neighbors_.resize( numthreads );
#   endif
  }
  // Atom mask
  mask_.SetMaskString( analyzeArgs.GetMaskNext() );
  // Set up output data set
  output_ = setup.DSL().AddSet( DataSet::MATRIX_FLT, setname, "WAVELET" );
  if (output_ == 0) return Analysis::ERR;
  if (outfile != 0) outfile->AddDataSet( output_ );
  if (doClustering_) {
    std::string const& dname = output_->Meta().Name();
    // Set up wavelet map clustering output set
    clustermap_ = setup.DSL().AddSet( DataSet::MATRIX_FLT, MetaData(dname, "clustermap") );
    if (clustermap_ == 0) return Analysis::ERR;
    if (clustermapout != 0) clustermapout->AddDataSet( clustermap_ );
    Dimension Cdim( 0, 1, "Cluster" );
    c_points_ = setup.DSL().AddSet(DataSet::INTEGER, MetaData(dname, "points"));
    c_minatm_ = setup.DSL().AddSet(DataSet::INTEGER, MetaData(dname, "minatm"));
    c_maxatm_ = setup.DSL().AddSet(DataSet::INTEGER, MetaData(dname, "maxatm"));
    c_minfrm_ = setup.DSL().AddSet(DataSet::INTEGER, MetaData(dname, "minfrm"));
    c_maxfrm_ = setup.DSL().AddSet(DataSet::INTEGER, MetaData(dname, "maxfrm"));
    c_avgval_ = setup.DSL().AddSet(DataSet::FLOAT,   MetaData(dname, "avgval"));
    if (c_points_==0 || c_minatm_==0 || c_maxatm_==0 ||
        c_minfrm_==0 || c_maxfrm_==0 || c_avgval_==0)
      return Analysis::ERR;
    c_points_->SetDim(0, Cdim);
    c_minatm_->SetDim(0, Cdim);
    c_maxatm_->SetDim(0, Cdim);
    c_minfrm_->SetDim(0, Cdim);
    c_maxfrm_->SetDim(0, Cdim);
    c_avgval_->SetDim(0, Cdim);
    if (clusterout != 0) {
      clusterout->AddDataSet( c_points_ );
      clusterout->AddDataSet( c_minatm_ );
      clusterout->AddDataSet( c_maxatm_ );
      clusterout->AddDataSet( c_minfrm_ );
      clusterout->AddDataSet( c_maxfrm_ );
      clusterout->AddDataSet( c_avgval_ );
    }
  }

  mprintf("    WAVELET: Using COORDS set '%s', wavelet type %s\n",
          coords_->legend(), Tokens_[wavelet_type_].description_);
  mprintf("\tCalculating for atoms in mask '%s'\n", mask_.MaskString());
  mprintf("\tScaling wavelet %i times starting from %g with delta of %g\n", 
          nb_, S0_, ds_);
  mprintf("\tCorrection: %g\n", correction_);
  mprintf("\tChiVal:     %g\n", chival_);
  if (outfile != 0) mprintf("\tOutput to '%s'\n", outfile->DataFilename().full());
  if (doClustering_) {
    mprintf("\tPerforming regional clustering on resulting wavelet map.\n");
    mprintf("\t  Wavelet map cluster set: '%s'\n", clustermap_->legend());
    if (clustermapout != 0)
      mprintf("\t  Wavelet map cluster output to '%s'\n", clustermapout->DataFilename().full());
    mprintf("\t  minpoints= %i, epsilon= %f\n", minPoints_, epsilon_);
  }
  return Analysis::OK;
}

#ifdef DEBUG_WAVELET
static inline void PrintComplex(const char* title, ComplexArray const& C) {
  if (title != 0) mprintf("DEBUG: %s:", title);
  for (ComplexArray::iterator cval = C.begin(); cval != C.end(); cval += 2)
    mprintf(" (%g,%g)", *cval, *(cval+1));
  mprintf("\n");
}
#endif

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
  Matrix<double> d_matrix;
  mprintf("\t%i frames, %i atoms, distance matrix will require %s\n", nframes, natoms,
          ByteString(d_matrix.sizeInBytes(nframes, natoms), BYTE_DECIMAL).c_str());
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
      // Distance of atom in currentFrame from its position in lastFrame.
      d_matrix[idx] = sqrt(DIST2_NoImage( currentFrame.XYZ(at), lastFrame.XYZ(at) ));
    //lastFrame = currentFrame; // TODO: Re-enable?
  }
# ifdef DEBUG_WAVELET
  // DEBUG: Write matrix to file.
  CpptrajFile dmatrixOut; // DEBUG
  dmatrixOut.OpenWrite("dmatrix.dat");
  Matrix<double>::iterator mval = d_matrix.begin();
  for (int row = 0; row != natoms; row++) {
    for (int col = 0; col != nframes; col++)
      dmatrixOut.Printf("%g ", *(mval++));
    dmatrixOut.Printf("\n");
  }
  dmatrixOut.CloseFile();
# endif

  // Precompute some factors for calculating scaled wavelets.
  const double one_over_sqrt_N = 1.0 / sqrt(static_cast<double>( nframes ));
  std::vector<int> arrayK( nframes );
  arrayK[0] = -1 * (nframes/2);
  for (int i = 1; i != nframes; i++)
    arrayK[i] = arrayK[i-1] + 1;
# ifdef DEBUG_WAVELET
  mprintf("DEBUG: K:");
  for (std::vector<int>::const_iterator kval = arrayK.begin(); kval != arrayK.end(); ++kval)
    mprintf(" %i", *kval);
  mprintf("\n");
# endif

  // Step 2 - Get FFT of wavelet for each scale.
  PubFFT pubfft;
  pubfft.SetupFFTforN( nframes );
  mprintf("\tMemory required for scaled wavelet array: %s\n",
          ByteString(2 * nframes * nb_ * sizeof(double), BYTE_DECIMAL).c_str());
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
#   ifdef DEBUG_WAVELET
    PrintComplex("wavelet_before_fft", scaledWavelet);
#   endif
    // Perform FFT
    pubfft.Forward( scaledWavelet );
    // Normalize
    scaledWavelet.Normalize( one_over_sqrt_N );
#   ifdef DEBUG_WAVELET
    PrintComplex("wavelet_after_fft", scaledWavelet);
#   endif
    FFT_of_Scaled_Wavelets.push_back( scaledWavelet );
  }
# ifdef DEBUG_WAVELET
  mprintf("DEBUG: Scaling factors:");
  for (Darray::const_iterator sval = scaleVector.begin(); sval != scaleVector.end(); ++sval)
    mprintf(" %g", *sval);
  mprintf("\n");
  mprintf("DEBUG: MIN:");
  for (int i = 0; i != nb_; i++)
    mprintf(" %g", MIN[i]);
  mprintf("\n");
# endif

  // Step 3 - For each atom, calculate the convolution of scaled wavelets
  //          with rows (atom distance vs frame) via dot product of the 
  //          frequency domains, i.e. Fourier-transformed, followed by an
  //          inverse FT.
  DataSet_MatrixFlt& OUT = static_cast<DataSet_MatrixFlt&>( *output_ );
  mprintf("\tMemory required for output matrix: %s\n",
          ByteString(Matrix<float>::sizeInBytes(nframes, natoms), BYTE_DECIMAL).c_str());
  OUT.Allocate2D( nframes, natoms ); // Should initialize to zero
  OUT.SetDim(Dimension::X, Dimension(1, 1, "Frame"));
  OUT.SetDim(Dimension::Y, Dimension(1, 1, "Atom"));
  Matrix<double> MAX;
  mprintf("\tMemory required for Max array: %s\n",
          ByteString(MAX.sizeInBytes(nframes, natoms), BYTE_DECIMAL).c_str());
  MAX.resize( nframes, natoms );
  Darray magnitude( nframes ); // Scratch space for calculating magnitude across rows
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
#   ifdef DEBUG_WAVELET
    mprintf("VARIANCE: %g\n", d_var);
#   endif
    double var_norm = 1.0 / d_var;
    // Calculate FT of atom signal
    pubfft.Forward( AtomSignal );
#   ifdef DEBUG_WAVELET
    PrintComplex("AtomSignal", AtomSignal);
#   endif
    // Normalize
    AtomSignal.Normalize( one_over_sqrt_N );
    // Calculate dot product of atom signal with each scaled FT wavelet
    for (int iscale = 0; iscale != nb_; iscale++) {
      ComplexArray dot = AtomSignal.TimesComplexConj( FFT_of_Scaled_Wavelets[iscale] );
      // Inverse FT of dot product
      pubfft.Back( dot );
#     ifdef DEBUG_WAVELET
      PrintComplex("InverseFT_Dot", dot);
#     endif
      // Chi-squared testing
      midx = at * nframes;
      cidx = 0;
      for (int frm = 0; frm != nframes; frm++, cidx += 2, midx++) {
        magnitude[frm] = (dot[cidx]*dot[cidx] + dot[cidx+1]*dot[cidx+1]) * var_norm;
        if (magnitude[frm] < MIN[iscale])
          magnitude[frm] = 0.0;
        if (magnitude[frm] > MAX[midx]) {
          MAX[midx] = magnitude[frm];
          //Indices[midx] = iscale
          OUT[midx] = (float)(correction_ * scaleVector[iscale]);
        }
      }
#     ifdef DEBUG_WAVELET
      mprintf("DEBUG: AbsoluteValue:");
      for (Darray::const_iterator dval = magnitude.begin(); dval != magnitude.end(); ++dval)
        mprintf(" %g", *dval);
      mprintf("\n");
#     endif
    } // END loop over scales
  } // END loop over atoms
# ifdef DEBUG_WAVELET 
  // DEBUG: Print MAX
  CpptrajFile maxmatrixOut; // DEBUG
  maxmatrixOut.OpenWrite("maxmatrix.dat");
  for (int col = 0; col != nframes; col++) {
    for (int row = 0; row != natoms; row++)
      maxmatrixOut.Printf("%g ", MAX.element(col, row));
    maxmatrixOut.Printf("\n");
  }
  maxmatrixOut.CloseFile();
# endif

  // Step 4 - Wavelet map clustering. Try to automatically identify regions
  //          where important motions are occurring.
  if (doClustering_) {
    if (ClusterMap( OUT )) return Analysis::ERR;
  }
      
  return Analysis::OK;
}

// -----------------------------------------------------------------------------
static inline void IdxToColRow(int idx, int ncols, int& col, int& row) {
  col = idx % ncols;
  row = idx / ncols;
}

// Analysis_Wavelet::ClusterMap()
int Analysis_Wavelet::ClusterMap(DataSet_MatrixFlt const& matrix) {
  // Set up output cluster map
  DataSet_MatrixFlt& outmap = static_cast<DataSet_MatrixFlt&>( *clustermap_ );
  outmap.Allocate2D( matrix.Ncols(), matrix.Nrows() );
  std::fill( outmap.begin(), outmap.end(), -1.0 );

  // Go through the set, calculate the average. Also determine the max.
  double maxVal = matrix.GetElement(0);
  unsigned int maxIdx = 0;
  Avg_ = 0.0;
  for (unsigned int i = 1; i != matrix.Size(); i++) {
    double val = matrix.GetElement(i);
    if (val > maxVal) {
      maxVal = val;
      maxIdx = i;
    }
    Avg_ += val;
  }
  Avg_ /= (double)matrix.Size();
  int maxCol, maxRow;
  IdxToColRow( maxIdx, matrix.Ncols(), maxCol, maxRow );
  mprintf("\t%zu elements, max= %f at index %u (frame %i, atom %i), Avg= %f\n",
          matrix.Size(), maxVal, maxIdx, maxCol+1, maxRow+1, Avg_);
  mprintf("\tPoints below %f will be treated as noise.\n", Avg_);

  // Use DBSCAN-style algorithm to cluster points. Any point less than the
  // average will be considered noise.
  std::vector<bool> Visited( matrix.Size(), false );
  const char UNASSIGNED = 'U';
  const char NOISE = 'N';
  const char INCLUSTER = 'C';
  std::vector<char> Status( matrix.Size(), UNASSIGNED );
  // Loop over all points in matrix
  Iarray NeighborPts;    // Will hold all neighbors of the current point
  Iarray Npts2;          // Will hold neighbors of a neighbor
  Iarray cluster_frames; // Hold indices of current cluster
  ProgressBar progress(matrix.Size());
  int iterations = 0;
# ifdef TIMER
  t_overall_.Start();
# endif
  for (int point = 0; point != (int)matrix.Size(); point++)
  {
    if (!Visited[point])
    {
      progress.Update(iterations++);
      // Mark this point as visited.
      Visited[point] = true;
      double val = matrix.GetElement(point);
      // If this point is less than the cutoff just mark it as noise.
      if (val < Avg_)
        Status[point] = NOISE;
      else
      {
        // Determine how many other points are near this point.
#       ifdef TIMER
        t_query1_.Start();
#       endif
        RegionQuery( NeighborPts, val, point, matrix );
#       ifdef TIMER
        t_query1_.Stop();
#       endif
#       ifdef DEBUG_CLUSTERMAP
        mprintf("\tPoint %u\n", point);
        mprintf("\t\t%u neighbors:\n", NeighborPts.size());
#       endif
        // If # of neighbors less than cutoff, noise; otherwise cluster.
        if ((int)NeighborPts.size() < minPoints_)
        {
#         ifdef DEBUG_CLUSTERMAP
          mprintf(" NOISE\n");
#         endif
          Status[point] = NOISE;
        }
        else
        {
          // Expand cluster
          cluster_frames.clear();
          cluster_frames.push_back( point );
          // NOTE: Use index instead of iterator since NeighborPts may be
          //       modified inside this loop.
          unsigned int endidx = NeighborPts.size();
          for (unsigned int idx = 0; idx < endidx; ++idx)
          {
            int neighbor_pt = NeighborPts[idx];
            double neighbor_val = matrix.GetElement(neighbor_pt);
            if (!Visited[neighbor_pt])
            {
              progress.Update(iterations++);
#             ifdef DEBUG_CLUSTERMAP
              mprintf(" %i", neighbor_pt + 1);
#             endif
              // Mark this neighbor as visited
              Visited[neighbor_pt] = true;
              // Determine how many other points are near this neighbor
#             ifdef TIMER
              t_query2_.Start();
#             endif
              RegionQuery( Npts2, neighbor_val, neighbor_pt, matrix );
#             ifdef TIMER
              t_query2_.Stop();
#             endif
              if ((int)Npts2.size() >= minPoints_)
              {
                // Add other points to current neighbor list
                NeighborPts.insert( NeighborPts.end(), Npts2.begin(), Npts2.end() );
                endidx = NeighborPts.size();
              }
            }
            // If neighbor is not yet part of a cluster, add it to this one.
            if (Status[neighbor_pt] != INCLUSTER)
            {
              cluster_frames.push_back( neighbor_pt );
              Status[neighbor_pt] = INCLUSTER;
            }
          }
          // Remove duplicate frames
          // TODO: Take care of this in Renumber?
          std::sort(cluster_frames.begin(), cluster_frames.end());
          Iarray::iterator it = std::unique(cluster_frames.begin(),
                                            cluster_frames.end());
          cluster_frames.resize( std::distance(cluster_frames.begin(),it) );
          // Add cluster to the list
#         ifdef DEBUG_CLUSTERMAP
          mprintf("\n");
#         endif
          AddCluster( cluster_frames, matrix );
        }
      } // END value > cutoff
    } // END if not visited
  } // END loop over matrix points
# ifdef TIMER
  t_overall_.Stop();
  t_query1_.WriteTiming(2, "Region Query 1:", t_overall_.Total());
  t_query2_.WriteTiming(2, "Region Query 2:", t_overall_.Total());
  t_overall_.WriteTiming(1, "Overall:");
# endif
  mprintf("\t%zu clusters:\n", clusters_.size());

  // Sort by number of points
  std::sort(clusters_.begin(), clusters_.end());
  // Renumber clusters
  int cnum = 0;
  for (Carray::iterator CL = clusters_.begin(); CL != clusters_.end(); ++CL)
  {
    CL->SetCnum( cnum );
    for (Iarray::const_iterator pt = CL->Points().begin(); pt != CL->Points().end(); ++pt)
      outmap[*pt] = cnum;
    // Save cluster data to sets
    int ival = (int)CL->Points().size();
    c_points_->Add( cnum, &ival );
    ival = CL->MinRow()+1;
    c_minatm_->Add( cnum, &ival );
    ival = CL->MaxRow()+1;
    c_maxatm_->Add( cnum, &ival );
    ival = CL->MinCol()+1;
    c_minfrm_->Add( cnum, &ival );
    ival = CL->MaxCol()+1;
    c_maxfrm_->Add( cnum, &ival );
    float fval = (float)CL->Avg();
    c_avgval_->Add( cnum, &fval );
    cnum++;
  }
  //  mprintf("\t %i: %zu points, atoms %i-%i, frames %i-%i, avg= %f\n",
  //          CL->Cnum(), CL->Points().size(),
  //          CL->MinRow()+1, CL->MaxRow()+1,
  //          CL->MinCol()+1, CL->MaxCol()+1, CL->Avg());

  return 0;
}

// Analysis_Wavelet::RegionQuery()
void Analysis_Wavelet::RegionQuery(Iarray& NeighborPts, double val, int point,
                                   DataSet_2D const& matrix)
{
  NeighborPts.clear();
# ifdef _OPENMP
  for (std::vector<Iarray>::iterator it = thread_neighbors_.begin();
                                     it != thread_neighbors_.end(); ++it)
    it->clear();
# endif
  int point_col, point_row;
  IdxToColRow( point, matrix.Ncols(), point_col, point_row );
  int msize = (int)matrix.Size();
  int otherpoint;

# ifdef _OPENMP
# pragma omp parallel private(otherpoint)
  {
# pragma omp for
# endif
  for (otherpoint = 0; otherpoint < msize; otherpoint++)
  {
    if (point != otherpoint) {
      // Distance calculation.
      int other_col, other_row;
      IdxToColRow( otherpoint, matrix.Ncols(), other_col, other_row );
      double other_val = matrix.GetElement( otherpoint );
      double dv = val - other_val;
      double dr = (double)(point_row - other_row);
      double dc = (double)(point_col - other_col);
      double dist2 = dv*dv + dr*dr + dc*dc;

      if ( dist2 < epsilon2_ )
#       ifdef _OPENMP
        thread_neighbors_[mythread_].push_back( otherpoint );
#       else
        NeighborPts.push_back( otherpoint );
#       endif
    }
  }
# ifdef _OPENMP
  } // END pragma omp parallel
  for (std::vector<Iarray>::const_iterator it = thread_neighbors_.begin();
                                           it != thread_neighbors_.end(); ++it)
    for (Iarray::const_iterator pt = it->begin(); pt != it->end(); ++pt)
      NeighborPts.push_back( *pt );
# endif
}

// Analysis_Wavelet::AddCluster()
void Analysis_Wavelet::AddCluster(Iarray const& points, DataSet_2D const& matrix)
{
# ifdef DEBUG_CLUSTERMAP
  mprintf("Cluster %i (%zu):", nClusters_, points.size());
# endif
  int ncols = (int)matrix.Ncols();
  int min_col, min_row;
  IdxToColRow(points.front(), ncols, min_col, min_row);
  int max_col = min_col;
  int max_row = min_row;
  int row, col;
  double cavg = 0.0;
  for (Iarray::const_iterator pt = points.begin(); pt != points.end(); ++pt) {
#   ifdef DEBUG_CLUSTERMAP
    mprintf(" %i", *pt);
#   endif
    // Determine min/max row/col and average of all points
    IdxToColRow( *pt, ncols, col, row );
    min_col = std::min(col, min_col);
    max_col = std::max(col, max_col);
    min_row = std::min(row, min_row);
    max_row = std::max(row, max_row);
    cavg += matrix.GetElement( *pt );
  }
  cavg /= (double)points.size();
  clusters_.push_back( Cluster(points, cavg, nClusters_, min_col, max_col, min_row, max_row) );
# ifdef DEBUG_CLUSTERMAP
  mprintf("\n");
# endif
  nClusters_++;
}
