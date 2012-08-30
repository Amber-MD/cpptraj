#include <cmath> // sqrt
#include <map>
#include "Action_AtomicCorr.h"
#include "CpptrajStdio.h"
#include "TriangleMatrix.h"
#include "vectormath.h"
#include "ProgressBar.h"
#ifdef _OPENMP
#  include "omp.h"
#endif

// CONSTRUCTOR
Action_AtomicCorr::Action_AtomicCorr() :
  acorr_mode_(ATOM),
  cut_(0.0),
  min_(0)
{}

const char Action_AtomicCorr::ModeString[2][8] = {"atom", "residue"};

/** Usage: atomiccorr [<mask>] out <filename> */
int Action_AtomicCorr::init() {
  outname_ = actionArgs.GetStringKey("out");
  if (outname_.empty()) {
    mprinterr("Error: atomiccorr: No output filename specified [out <filename>]\n");
    return 1;
  }
  cut_ = actionArgs.getKeyDouble("cut", 0.0);
  if (cut_ < 0.0 || cut_ > 1.0) {
    mprinterr("Error: atomiccorr: cut value must be between 0 and 1.\n");
    return 1;
  }
  min_ = actionArgs.getKeyInt("min",0);
  if (actionArgs.hasKey("byatom"))
    acorr_mode_ = ATOM;
  else if (actionArgs.hasKey("byres"))
    acorr_mode_ = RES;
  mask_.SetMaskString( actionArgs.GetMaskNext() );

  mprintf("    ATOMICCORR: Correlation of %s motions will be calculated for\n",
          ModeString[acorr_mode_]);
  mprintf("\tatoms in mask [%s], output to file %s\n", mask_.MaskString(), outname_.c_str());
  if (cut_ != 0)
    mprintf("\tOnly correlations greater than %.2f or less than -%.2f will be printed.\n",
            cut_,cut_);
  if (min_!=0)
    mprintf("\tOnly correlations for %ss > %i apart will be calculated.\n",
            ModeString[acorr_mode_],min_);

  return 0;
}

int Action_AtomicCorr::setup() {
  if (currentParm->SetupIntegerMask( mask_ )) return 1;
  mask_.MaskInfo();
  if (mask_.None()) return 1;
  if (acorr_mode_ == ATOM) {
    // Setup output array; labels and index
    atom_vectors_.clear();
    for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
      atom_vectors_.push_back( AtomVector(integerToString( *atom ), *atom) );
  } else {
    std::map<int,AtomMask> rmaskmap;
    // Find which residues selected atoms belong to.
    for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom) 
    {
      int current_res = (*currentParm)[*atom].ResNum();
      std::map<int,AtomMask>::iterator rmask = rmaskmap.find( current_res );
      if ( rmask == rmaskmap.end() ) {
        // Residue not yet in map.
        AtomMask newmask;
        newmask.AddAtom( *atom );
        rmaskmap.insert( std::pair<int,AtomMask>( current_res, newmask ) );
      } else {
        // Residue is already in map. Add this atom.
        (*rmask).second.AddAtom( *atom );
      }
    }
    // Place selected residues in mask vector and setup output array; labels and index.
    resmasks_.clear();
    atom_vectors_.clear();
    for (std::map<int,AtomMask>::iterator rmask = rmaskmap.begin();
                                          rmask != rmaskmap.end(); ++rmask)
    {
      if (debug > 0)
        mprintf("DBG:\tRes mask for %i has %i atoms\n",(*rmask).first,(*rmask).second.Nselected());
      resmasks_.push_back( (*rmask).second );
      atom_vectors_.push_back( AtomVector( currentParm->TruncResNameNum( (*rmask).first ),
                                           (*rmask).first ) );
    }
    mprintf("\tSelected %zu residues.\n", resmasks_.size());
  }
  return 0;
}

int Action_AtomicCorr::action() {
  double RXYZ[3], CXYZ[3];
  // On first pass through refframe will be empty and first frame will become ref.
  if (!refframe_.empty()) {
    ACvector::iterator atom_vector = atom_vectors_.begin();
    if (acorr_mode_ == ATOM) {
      // For each atom in mask, calc delta position.
      for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom) 
      {
        const double* tgtxyz = currentFrame->XYZ( *atom );
        const double* refxyz = refframe_.XYZ( *atom );
        (*atom_vector).push_back( (float)(tgtxyz[0] - refxyz[0]) );
        (*atom_vector).push_back( (float)(tgtxyz[1] - refxyz[1]) );
        (*atom_vector).push_back( (float)(tgtxyz[2] - refxyz[2]) );
        ++atom_vector;
      }
    } else {
      for (std::vector<AtomMask>::iterator rmask = resmasks_.begin();
                                           rmask != resmasks_.end(); ++rmask)
      {
        currentFrame->GeometricCenter( CXYZ, *rmask );
        refframe_.GeometricCenter( RXYZ, *rmask );
        (*atom_vector).push_back( (float)(CXYZ[0] - RXYZ[0]) );
        (*atom_vector).push_back( (float)(CXYZ[1] - RXYZ[1]) );
        (*atom_vector).push_back( (float)(CXYZ[2] - RXYZ[2]) );
        ++atom_vector;
      } 
    }
  }
  // Store this frame as new reference frame
  refframe_ = *currentFrame;
  return 0;
}

void Action_AtomicCorr::print() {
  int idx, idx3, vec1size;
  double V1[3], V2[3], corr_coeff;
  mprintf("    ATOMICCORR: Calculating correlations between %s vectors:\n",
          ModeString[acorr_mode_]);
  if (atom_vectors_.empty()) {
    mprinterr("Error: atomiccorr: No vectors calcd.\n");
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
  ProgressBar progress( tmatrix->Size() );
  int iprogress = 0;
  for (ACvector::iterator vec1 = atom_vectors_.begin(); vec1 != av_end1; ++vec1)
  {
    for (ACvector::iterator vec2 = vec1 + 1; vec2 != av_end; ++vec2)
    {
      corr_coeff = 0.0;
      progress.Update( iprogress++ );
      // If vectors are too close, skip. vec2 always > vec1
      if ( (*vec2) - (*vec1) > min_ ) {
        // Make sure same # of frames in each and not empty
        if ( (*vec1).empty() || (*vec2).empty() ) {
          mprintf("Warning: autocorr: A vector is empty: Vec%zu=%zu, Vec%zu=%zu\n",
                  vec1 - atom_vectors_.begin(), (*vec1).size(),
                  vec2 - atom_vectors_.begin(), (*vec2).size());
        } else if ( (*vec1).size() != (*vec2).size() ) {
          mprintf("Warning: atomiccorr: Vec %zu and Vec %zu do not have same # of frames.\n",
                  vec1 - atom_vectors_.begin(), vec2 - atom_vectors_.begin());
        } else {
          vec1size = (int)((*vec1).size() / 3);
#ifdef _OPENMP
#pragma omp parallel private(idx, idx3, V1, V2) reduction(+: corr_coeff)
{
#pragma omp for
#endif
          for (idx = 0; idx < vec1size; ++idx) {
            idx3 = idx * 3;
            (*vec1).XYZ(V1, idx3);
            (*vec2).XYZ(V2, idx3);
            normalize( V1 );
            normalize( V2 );
            corr_coeff += dot_product( V1, V2 );
          }
#ifdef _OPENMP
} // END pragma omp parallel
#endif
          corr_coeff /= (double)vec1size;
          if (fabs(corr_coeff) <= cut_)
            corr_coeff = 0.0;
        }
      }
      tmatrix->AddElement( corr_coeff );
    } // END inner loop
  } // END outer loop

  // Add dataset to output file
  DataFile* outfile = DFL->AddSetToFile( outname_, (DataSet*)tmatrix );
  if (acorr_mode_ == ATOM)
    outfile->ProcessArgs("xlabel Atom");
  else
    outfile->ProcessArgs("xlabel Residue");
  std::string ylabels = "ylabels ";
  for (ACvector::iterator atom = atom_vectors_.begin();
                          atom != atom_vectors_.end(); ++atom)
    ylabels += ( (*atom).Label() + "," );
  outfile->ProcessArgs( ylabels );
}

