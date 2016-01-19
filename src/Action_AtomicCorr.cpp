#include <cmath> // sqrt
#include <map>
#include "Action_AtomicCorr.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString
#include "DataSet_MatrixFlt.h"
#include "ProgressBar.h"
#ifdef _OPENMP
#  include "omp.h"
#endif

// CONSTRUCTOR
Action_AtomicCorr::Action_AtomicCorr() :
  acorr_mode_(ATOM),
  cut_(0.0),
  min_(0),
  debug_(0),
  dset_(0),
  outfile_(0)
{}

void Action_AtomicCorr::Help() const {
  mprintf("\t[<mask>] [out <filename>] [cut <cutoff>] [min <min spacing>]\n"
          "\t[byatom | byres]\n"
          "  Calculate average correlations between the motion of atoms in <mask>.\n");
}

const char* Action_AtomicCorr::ModeString[] = {"atom", "residue"};

// Action_AtomicCorr::Init()
Action::RetType Action_AtomicCorr::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  outfile_ = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  cut_ = actionArgs.getKeyDouble("cut", 0.0);
  if (cut_ < 0.0 || cut_ > 1.0) {
    mprinterr("Error: cut value must be between 0 and 1.\n");
    return Action::ERR;
  }
  min_ = actionArgs.getKeyInt("min",0);
  if (actionArgs.hasKey("byatom"))
    acorr_mode_ = ATOM;
  else if (actionArgs.hasKey("byres"))
    acorr_mode_ = RES;
  mask_.SetMaskString( actionArgs.GetMaskNext() );
  
  // Set up DataSet
  dset_ = init.DSL().AddSet( DataSet::MATRIX_FLT, actionArgs.GetStringNext(), "ACorr" );
  if (dset_ == 0) {
    mprinterr("Error: Could not allocate output data set.\n");
    return Action::ERR;
  }
  // Add DataSet to output file
  if (outfile_ != 0) outfile_->AddDataSet( dset_ );

  mprintf("    ATOMICCORR: Correlation of %s motions will be calculated for\n",
          ModeString[acorr_mode_]);
  mprintf("\tatoms in mask [%s]", mask_.MaskString());
  if (outfile_ != 0) mprintf(", output to file %s", outfile_->DataFilename().full());
  mprintf("\n\tData saved in set '%s'\n", dset_->legend());
  if (cut_ != 0)
    mprintf("\tOnly correlations greater than %.2f or less than -%.2f will be printed.\n",
            cut_,cut_);
  if (min_!=0)
    mprintf("\tOnly correlations for %ss > %i apart will be calculated.\n",
            ModeString[acorr_mode_],min_);
# ifdef MPI
  trajComm_ = init.TrajComm();
  if (trajComm_.Size() > 1)
    mprintf("\nWarning: 'atomiccorr' in parallel will not work correctly if coordinates have\n"
              "Warning:   been modified by previous actions (e.g. 'rms').\n\n");
# endif
  return Action::OK;
}

Action::RetType Action_AtomicCorr::Setup(ActionSetup& setup) {
  if (setup.Top().SetupIntegerMask( mask_ )) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) return Action::SKIP;
  if (acorr_mode_ == ATOM) {
    // Setup output array; labels and index
    atom_vectors_.clear();
    for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
      atom_vectors_.push_back( AtomVector(integerToString( *atom + 1 ), *atom) );
  } else {
    std::map<int,AtomMask> rmaskmap;
    // Find which residues selected atoms belong to.
    for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom) 
    {
      int current_res = setup.Top()[*atom].ResNum();
      std::map<int,AtomMask>::iterator rmask = rmaskmap.find( current_res );
      if ( rmask == rmaskmap.end() ) {
        // Residue not yet in map.
        AtomMask newmask;
        newmask.AddAtom( *atom );
        rmaskmap.insert( std::pair<int,AtomMask>( current_res, newmask ) );
      } else {
        // Residue is already in map. Add this atom.
        rmask->second.AddAtom( *atom );
      }
    }
    // Place selected residues in mask vector and setup output array; labels and index.
    resmasks_.clear();
    atom_vectors_.clear();
    for (std::map<int,AtomMask>::const_iterator rmask = rmaskmap.begin();
                                                rmask != rmaskmap.end(); ++rmask)
    {
      if (debug_ > 0)
        mprintf("DBG:\tRes mask for %i has %i atoms\n", rmask->first, rmask->second.Nselected());
      resmasks_.push_back( rmask->second );
      atom_vectors_.push_back( AtomVector( setup.Top().TruncResNameNum( rmask->first ),
                                           rmask->first ) );
    }
    mprintf("\tSelected %zu residues.\n", resmasks_.size());
  }
  return Action::OK;
}

#ifdef MPI
int Action_AtomicCorr::ParallelPreloadFrames(FArray const& preload_frames) {
  unsigned int idx = preload_frames.size() - 1;
  previousFrame_ = preload_frames[idx];
  return 0;
}

int Action_AtomicCorr::SyncAction() {
  if (trajComm_.Size() < 2 || atom_vectors_.empty()) return 0;
  // atom_vectors_ should be same size on all ranks, and length of all
  // should be the same. Need to know length on each vector.
  int vlength = (int)atom_vectors_[0].size();
  if (trajComm_.Master()) {
    // MASTER
    int* rank_sizes = new int[ trajComm_.Size() ];
    trajComm_.GatherMaster( &vlength, 1, MPI_INT, rank_sizes );
    // Compute total number of elements 
    int total_elements = rank_sizes[0];
    for (int rank = 1; rank < trajComm_.Size(); rank++)
      total_elements += rank_sizes[rank];
    mprintf("DEBUG: Total atomiccorr elements: %i\n", total_elements);
    int idx = 0; // For tag
    for (ACvector::iterator av = atom_vectors_.begin();
                            av != atom_vectors_.end(); ++av, ++idx)
    {
      // Resize atomic vector to hold all incoming data from ranks.
      av->resize( total_elements );
      float* ptr = av->Fptr() + rank_sizes[0];
      // Receive data from all ranks for this vector
      for (int rank = 1; rank < trajComm_.Size(); rank++) {
        trajComm_.Recv( ptr, rank_sizes[rank], MPI_FLOAT, rank, 1400 + idx );
        ptr += rank_sizes[rank];
      } 
    }
    delete[] rank_sizes;
  } else {
    // RANK
    trajComm_.GatherMaster( &vlength, 1, MPI_INT, 0 );
    int idx = 0; // For tag
    for (ACvector::iterator av = atom_vectors_.begin(); 
                            av != atom_vectors_.end(); ++av, ++idx) 
      trajComm_.Send( av->Fptr(), av->size(), MPI_FLOAT, 0, 1400 + idx );
  }
  return 0;
}
#endif

Action::RetType Action_AtomicCorr::DoAction(int frameNum, ActionFrame& frm) {
  // On first pass through refframe will be empty and first frame will become ref.
  if (!previousFrame_.empty()) {
    ACvector::iterator atom_vector = atom_vectors_.begin();
    if (acorr_mode_ == ATOM) {
      // For each atom in mask, calc delta position.
      for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom) 
      {
        const double* tgtxyz = frm.Frm().XYZ( *atom );
        const double* refxyz = previousFrame_.XYZ( *atom );
        atom_vector->push_back( (float)(tgtxyz[0] - refxyz[0]) );
        atom_vector->push_back( (float)(tgtxyz[1] - refxyz[1]) );
        atom_vector->push_back( (float)(tgtxyz[2] - refxyz[2]) );
        ++atom_vector;
      }
    } else {
      for (std::vector<AtomMask>::const_iterator rmask = resmasks_.begin();
                                                 rmask != resmasks_.end(); ++rmask)
      {
        Vec3 CXYZ = frm.Frm().VGeometricCenter( *rmask );
        Vec3 RXYZ = previousFrame_.VGeometricCenter( *rmask );
        atom_vector->push_back( (float)(CXYZ[0] - RXYZ[0]) );
        atom_vector->push_back( (float)(CXYZ[1] - RXYZ[1]) );
        atom_vector->push_back( (float)(CXYZ[2] - RXYZ[2]) );
        ++atom_vector;
      } 
    }
  }
  // Store this frame as new reference frame
  previousFrame_ = frm.Frm();
  return Action::OK;
}

void Action_AtomicCorr::Print() {
  int idx, idx3, vec1size;
  Vec3 V1, V2;
  mprintf("    ATOMICCORR: Calculating correlations between %s vectors:\n",
          ModeString[acorr_mode_]);
  if (atom_vectors_.empty()) {
    mprinterr("Error: No vectors calculated.\n");
    return;
  }
  DataSet_MatrixFlt* tmatrix = static_cast<DataSet_MatrixFlt*>( dset_ );
  tmatrix->AllocateTriangle( atom_vectors_.size() );
  // Calculate correlation coefficient of each atomic vector to each other
  ACvector::iterator av_end = atom_vectors_.end();
  ACvector::iterator av_end1 = av_end - 1;
  ProgressBar progress( tmatrix->Size() );
  int iprogress = 0;
  for (ACvector::const_iterator vec1 = atom_vectors_.begin(); vec1 != av_end1; ++vec1)
  {
    for (ACvector::const_iterator vec2 = vec1 + 1; vec2 != av_end; ++vec2)
    {
      double corr_coeff = 0.0;
      progress.Update( iprogress++ );
      // If vectors are too close, skip. vec2 always > vec1
      if ( (*vec2) - (*vec1) > min_ ) {
        // Make sure same # of frames in each and not empty
        if ( vec1->empty() || vec2->empty() ) {
          mprintf("Warning: A vector is empty: Vec%zu=%zu, Vec%zu=%zu\n",
                  vec1 - atom_vectors_.begin(), vec1->size(),
                  vec2 - atom_vectors_.begin(), vec2->size());
        } else if ( vec1->size() != vec2->size() ) {
          mprintf("Warning: Vec %zu and Vec %zu do not have same # of frames.\n",
                  vec1 - atom_vectors_.begin(), vec2 - atom_vectors_.begin());
        } else {
          vec1size = (int)(vec1->size() / 3);
#ifdef _OPENMP
#pragma omp parallel private(idx, idx3, V1, V2) reduction(+: corr_coeff)
{
#pragma omp for
#endif
          for (idx = 0; idx < vec1size; ++idx) {
            idx3 = idx * 3;
            V1 = vec1->VXYZ(idx3);
            V2 = vec2->VXYZ(idx3);
            V1.Normalize();
            V2.Normalize();
            corr_coeff += (V1 * V2);
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

  if (acorr_mode_ == ATOM) {
    Dimension dim(1.0, 1.0, "Atom");
    dset_->SetDim(Dimension::X, dim);
    dset_->SetDim(Dimension::Y, dim);
  } else {
    Dimension dim(1.0, 1.0, "Residue");
    dset_->SetDim(Dimension::X, dim);
    dset_->SetDim(Dimension::Y, dim);
  }
  std::string labels;
  for (ACvector::const_iterator atom = atom_vectors_.begin();
                                atom != atom_vectors_.end(); ++atom)
    labels += ( atom->Label() + "," );
  if (outfile_ != 0)
    outfile_->ProcessArgs( "xlabels " + labels + " ylabels " + labels );
}
