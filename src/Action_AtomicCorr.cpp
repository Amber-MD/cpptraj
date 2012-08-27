#include <cmath> // sqrt
#include "Action_AtomicCorr.h"
#include "CpptrajStdio.h"
#include "TriangleMatrix.h"
#include "vectormath.h"
#ifdef _OPENMP
#  include "omp.h"
#endif

// CONSTRUCTOR
Action_AtomicCorr::Action_AtomicCorr() {}

/** Usage: atomiccorr [<mask>] out <filename> */
int Action_AtomicCorr::init() {
  outname_ = actionArgs.GetStringKey("out");
  if (outname_.empty()) {
    mprinterr("Error: atomiccorr: No output filename specified [out <filename>]\n");
    return 1;
  }
  mask_.SetMaskString( actionArgs.GetMaskNext() );

  mprintf("    ATOMICCORR: Correlation of atomic motions will be calculated for\n");
  mprintf("\tatoms in mask [%s], output to file %s\n", mask_.MaskString(), outname_.c_str());

  return 0;
}

int Action_AtomicCorr::setup() {
  if (currentParm->SetupIntegerMask( mask_ )) return 1;
  mask_.MaskInfo();
  if (mask_.None()) return 1;
  // Setup output array
  if ( mask_.Nselected() > (int)atom_vectors_.size() )
    atom_vectors_.resize( mask_.Nselected() );
  return 0;
}

int Action_AtomicCorr::action() {
  // On first pass through refframe will be empty and first frame will become ref.
  if (!refframe_.empty()) {
    // For each atom in mask, calc delta position.
    ACvector::iterator atom_vector = atom_vectors_.begin();
    for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom) 
    {
      const double* tgtxyz = currentFrame->XYZ( *atom );
      const double* refxyz = refframe_.XYZ( *atom );
      (*atom_vector).push_back( (float)(tgtxyz[0] - refxyz[0]) );
      (*atom_vector).push_back( (float)(tgtxyz[1] - refxyz[1]) );
      (*atom_vector).push_back( (float)(tgtxyz[2] - refxyz[2]) );
      ++atom_vector;
    }
  }
  // Store this frame as new reference frame
  refframe_ = *currentFrame;
  return 0;
}

void Action_AtomicCorr::print() {
  unsigned int idx, idx3;
  const float *v1, *v2;
  double V1[3], V2[3], corr_coeff;
  if (atom_vectors_.empty()) {
    mprinterr("Error: atomiccorr: No atomic vectors calcd.\n");
    return;
  }
  TriangleMatrix* tmatrix = (TriangleMatrix*)DSL->AddSet( DataSet::TRIMATRIX,
                                                          "", "ACorr" );
  if (tmatrix == NULL) {
    mprinterr("Error: atomiccorr: Could not allocate output dataset.\n");
    return;
  }
  tmatrix->Setup( atom_vectors_.size() );
  // Calculate correlation coefficient of each atomic vector to each other
  ACvector::iterator av_end = atom_vectors_.end();
  ACvector::iterator av_end1 = av_end - 1;
  for (ACvector::iterator vec1 = atom_vectors_.begin(); vec1 != av_end1; ++vec1)
  {
    for (ACvector::iterator vec2 = vec1 + 1; vec2 != av_end; ++vec2)
    {
      // Make sure same # of frames in each and not empty
      if ( (*vec1).empty() || (*vec2).empty() ) {
        mprintf("Warning: autocorr: A vector is empty: Vec%zu=%zu, Vec%zu=%zu\n",
                vec1 - atom_vectors_.begin(), (*vec1).size(),
                vec2 - atom_vectors_.begin(), (*vec2).size());
        tmatrix->AddElement((float)0.0);
      } else if ( (*vec1).size() != (*vec2).size() ) {
        mprintf("Warning: atomiccorr: Vec %zu and Vec %zu do not have same # of frames.\n",
                vec1 - atom_vectors_.begin(), vec2 - atom_vectors_.begin());
        tmatrix->AddElement((float)0.0);
      } else {
        corr_coeff = 0.0;
        unsigned int vec1size = (*vec1).size() / 3;
#ifdef _OPENMP
#pragma omp parallel private(idx, idx3, v1, v2, V1, V2) reduction(+: corr_coeff)
{
#pragma omp for
#endif
        for (idx = 0; idx < vec1size; ++idx) {
          idx3 = idx * 3;
          v1 = &((*vec1)[idx3]);
          v2 = &((*vec2)[idx3]);
          V1[0] = (double)v1[0];
          V1[1] = (double)v1[1];
          V1[2] = (double)v1[2];
          V2[0] = (double)v2[0];
          V2[1] = (double)v2[1];
          V2[2] = (double)v2[2];
          normalize( V1 );
          normalize( V2 );
          corr_coeff += dot_product( V1, V2 );
        }
#ifdef _OPENMP
} // END pragma omp parallel
#endif
        corr_coeff /= (double)vec1size;
        tmatrix->AddElement( corr_coeff );
      }
    } // END inner loop
  } // END outer loop

  // Add dataset to output file
  DataFile* outfile = DFL->AddSetToFile( outname_, (DataSet*)tmatrix );
  outfile->ProcessArgs("xlabel Atom");
}

