#include <cmath> // sqrt
#include "Action_Matrix.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"

// CONSTRUCTOR
Action_Matrix::Action_Matrix() :
  Mat_(0),
  outfile_(0),
  outtype_(BYATOM),
  type_(DataSet_Matrix::NO_OP),
  snap_(0),
  order_(2),
  useMask2_(false),
  useMass_(false),
  CurrentParm_(0)
{}

void Action_Matrix::Help() {
  mprintf("\t[out <filename>] %s\n", ActionFrameCounter::HelpText);
  mprintf("\t[name <name>] [ byatom | byres [mass] | bymask [mass] ]\n");
  mprintf("\t[ ired [order <#>] ]\n");
  mprintf("\t[ {distcovar | idea} <mask1> ]\n");
  mprintf("\t[ {dist | correl | covar | mwcovar} <mask1> [<mask2>]\n");
  mprintf("\tCalculate a matrix of the specified type from input coordinates.\n");
}

// Action_Matrix::init()
Action::RetType Action_Matrix::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get Keywords
  filename_ = actionArgs.GetStringKey("out");
  // Get start/stop/offset
  if (InitFrameCounter(actionArgs)) return Action::ERR;
  // Determine matrix type
  type_ = DataSet_Matrix::TypeFromArg(actionArgs);
  // Output type
  if (actionArgs.hasKey("byres"))
    outtype_ = BYRESIDUE;
  else if (actionArgs.hasKey("bymask"))
    outtype_ = BYMASK;
  else if (actionArgs.hasKey("byatom"))
    outtype_ = BYATOM;
  else
    outtype_ = BYATOM;
  // Check if output type is valid for matrix type
  if ( outtype_ != BYATOM && (type_ == DataSet_Matrix::COVAR || 
                              type_ == DataSet_Matrix::MWCOVAR || 
                              type_ == DataSet_Matrix::IRED ) )
  {
    mprinterr("Error: matrix: for COVAR, MWCOVAR, or IRED matrix only byatom output possible\n");
    return Action::ERR;
  }
  // Get matrix name
  std::string name = actionArgs.GetStringKey("name");
  // UseMass
  useMass_ = actionArgs.hasKey("mass");

  if (type_ != DataSet_Matrix::IRED) {
    // Get masks if not IRED
    mask1_.SetMaskString( actionArgs.GetMaskNext() );
    std::string maskexpr = actionArgs.GetMaskNext();
    if (!maskexpr.empty()) useMask2_ = true;
    if ( useMask2_ && (type_ == DataSet_Matrix::IDEA || type_ == DataSet_Matrix::DISTCOVAR) )
    {
      mprinterr("Error: Mask 2 [%s] specified but not used for %s\n",
                maskexpr.c_str(), DataSet_Matrix::MatrixTypeString[type_]);
      useMask2_ = false;
      return Action::ERR;
    }
    if (useMask2_)
      mask2_.SetMaskString( maskexpr );
  } else {
    // Setup IRED vectors and determine Legendre order
    order_ = actionArgs.getKeyInt("order",1);
    if (order_ <= 0) {
      mprinterr("Error: matrix: order parameter <= 0, ignoring command\n");
      return Action::ERR;
    }
    for ( DataSetList::const_iterator DS = DSL->begin(); DS != DSL->end(); ++DS) {
      if ( (*DS)->Type() == DataSet::VECTOR ) {
        DataSet_Vector* Vtmp = (DataSet_Vector*)(*DS);
        if (Vtmp->IsIred())
          IredVectors_.push_back( Vtmp );
      }
    }
    if (IredVectors_.empty()) {
      mprinterr("Error: matrix: no vectors defined for IRED\n");
      return Action::ERR;
    }
  }
 
  // Set up matrix DataSet and type
  Mat_ = (DataSet_Matrix*)DSL->AddSet(DataSet::MATRIX, name, "Mat");
  if (Mat_ == 0) return Action::ERR;
  Mat_->SetType( type_ );
  // Add set to output file if doing BYATOM output
  if (outtype_ == BYATOM)
    outfile_ = DFL->AddSetToFile(filename_, Mat_);

  mprintf("    MATRIX: Calculating %s, output is", DataSet_Matrix::MatrixTypeString[type_]);
  switch (outtype_) {
    case BYATOM:    mprintf(" by atom"); break;
    case BYRESIDUE: mprintf(" by residue"); break;
    case BYMASK:    mprintf(" by mask"); break;
  }
  if (outtype_ != BYATOM) {
    if (useMass_)
      mprintf(" using mass weighting");
    else
      mprintf(" using no mass weighting");
  }
  mprintf("\n");
  if (type_ == DataSet_Matrix::IRED)
    mprintf("            %u IRED vecs, Order of Legendre polynomials: %i\n",
            IredVectors_.size(), order_);
  if (!filename_.empty()) {
    mprintf("            Printing to file %s\n",filename_.c_str());
    if (outtype_ != BYATOM) {
      mprintf("Warning: Output type is not 'byatom'. File will not be stored on internal\n");
      mprintf("Warning: DataFile stack, only basic formatting is possible.\n");
    }
  }
  if (!name.empty())
    mprintf("            Storing matrix on internal stack with name: %s\n", 
            Mat_->Legend().c_str());
  FrameCounterInfo();
  if (type_ != DataSet_Matrix::IRED) {
    mprintf("            Mask1: %s\n",mask1_.MaskString());
    if (useMask2_)
      mprintf("            Mask2: %s\n",mask2_.MaskString());
  }

  mprintf("\n");

  return Action::OK;
}

// Action_Matrix::FillMassArray()
void Action_Matrix::FillMassArray(Topology* currentParm, std::vector<double>& mass, AtomMask& mask)
{
  mass.clear();
  if (type_ == DataSet_Matrix::MWCOVAR) {
    mass.reserve( mask.Nselected() );
    for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom) 
      mass.push_back( (*currentParm)[ *atom ].Mass() );
  }
}

// Action_Matrix::setup()
Action::RetType Action_Matrix::Setup(Topology* currentParm, Topology** parmAddress) {
  size_t mask1tot = 0; // Will be # of columns
  size_t mask2tot = 0; // Will be # of rows if not symmetric matrix

  // Set up masks. Store mass info for masks if MWCOVAR
  if (type_ != DataSet_Matrix::IRED) {
    if (currentParm->SetupIntegerMask(mask1_)) return Action::ERR;
    mask1_.MaskInfo();
    if (mask1_.None()) {
      mprinterr("Error: No atoms selected for mask1.\n");
      return Action::ERR;
    }
    FillMassArray(currentParm, mass1_, mask1_); // MWCOVAR only
    if (useMask2_) {
      if (currentParm->SetupIntegerMask(mask2_)) return Action::ERR;
      mask2_.MaskInfo(); 
      if (mask2_.None()) {
        mprinterr("Error: No atoms selected for mask2.\n");
        return Action::ERR;
      }
      FillMassArray(currentParm, mass2_, mask2_); // MWCOVAR only
    }
    mask1tot = (size_t)mask1_.Nselected();
    mask2tot = (size_t)mask2_.Nselected();
  } else {
    // IRED - matrix # cols = # of IRED vectors
    mask1tot = IredVectors_.size();
  }
  if (mask1tot < mask2tot) {
    mprinterr("Error: DataSet_Matrix: # of atoms in mask1 < # of atoms in mask2\n");
    return Action::ERR;
  }

  // Allocate vector memory
  size_t vectsize = 0;
  switch( type_ ) {
    case DataSet_Matrix::DIST     : break;
    case DataSet_Matrix::DISTCOVAR: vectsize = mask1tot * (mask1tot - 1) / 2; break;
    case DataSet_Matrix::CORREL   :
    case DataSet_Matrix::COVAR    :
    case DataSet_Matrix::MWCOVAR  : vectsize = (mask1tot + mask2tot) * 3; break;
    default                       : vectsize = mask1tot + mask2tot;
  }
  if (Mat_->AllocateVectors( vectsize )) return Action::ERR;

  // Allocate matrix memory.
  size_t matrixSize = 0;
  size_t nrows = 0;
  size_t ncols = 0;
  if (mask2tot == 0) {
    // "Upper right half" matrix, including main diagonal.
    switch( type_ ) {
      case DataSet_Matrix::DISTCOVAR:
        ncols = vectsize;
        matrixSize = mask1tot * (mask1tot - 1) * (ncols + 1) / 4;
        break;
      case DataSet_Matrix::COVAR:
      case DataSet_Matrix::MWCOVAR:
        ncols = mask1tot * 3;
        matrixSize = 3 * ncols * (mask1tot + 1) / 2;
        break;
      default:
        ncols = mask1tot;
        matrixSize = ncols * (mask1tot + 1) / 2;
    }
  } else {
    // Full matrix - no DISTCOVAR, IDEA, or IRED possible
    switch( type_ ) {
      case DataSet_Matrix::DISTCOVAR:
      case DataSet_Matrix::IDEA:
      case DataSet_Matrix::IRED:
        mprinterr("Error: Second mask (full matrix) not supported for DISTCOVAR,\n");
        mprinterr("Error: IDEA, or IRED matrix.\n");
        return Action::ERR;
        break;
      case DataSet_Matrix::COVAR:
      case DataSet_Matrix::MWCOVAR:
        matrixSize = 9 * mask1tot * mask2tot;
        ncols = mask1tot * 3;
        nrows = mask2tot * 3;
        break;
      default:
        matrixSize = mask1tot * mask2tot;
        ncols = mask1tot;
        nrows = mask2tot;
    }
  }
  if (Mat_->AllocateMatrix(ncols, nrows, mask1tot, matrixSize)) return Action::ERR;

  // Mass info needed for MWCOVAR analysis, store in matrix dataset.
  if (type_ == DataSet_Matrix::MWCOVAR) 
    Mat_->StoreMass( mass1_ );
  CurrentParm_ = currentParm;
  return Action::OK;
}

// -----------------------------------------------------------------------------
// LegendrePoly()
/** Calculate Legendre Polynomial. Used only by IRED. */
static double LegendrePoly(int order, double val) {
  if (order == 0)
    return 1.0;
  else if (order == 1)
    return val;

  double pNminus1 = 1.0;
  double pN = val;
  double twox = 2.0 * val;
  double f2 = val;
  double d = 1.0;

  for(int i=2; i<=order; i++){
    double f1 = d++;
    f2 += twox;
    double pNplus1 = (f2 * pN - f1 * pNminus1) / d;
    pNminus1 = pN;
    pN = pNplus1;
  }
  return pN;
}

/** Calc isotropic reorientational eigenmode dynamics.
  * See JACS 2002, 124, 4522, eq. A14 
  * CAVEAT: omegaK-omegaL is not "just" the intra molecular angle there.
  */
void Action_Matrix::CalcIredMatrix() {
  DataSet_Matrix::iterator v2idx1 = Mat_->v2begin();
  // Store length of IRED vectors in vect2
  for (std::vector<DataSet_Vector*>::iterator Vtmp = IredVectors_.begin();
                                              Vtmp != IredVectors_.end(); ++Vtmp)
    *(v2idx1++) = sqrt( (*Vtmp)->Dot( *(*Vtmp) ) );

  // Loop over all pairs of IRED vectors.
  DataSet_Matrix::iterator mat = Mat_->begin();
  DataSet_Matrix::iterator v1idx = Mat_->v1begin();
  v2idx1 = Mat_->v2begin();
  for (std::vector<DataSet_Vector*>::iterator Vtmp = IredVectors_.begin();
                                              Vtmp != IredVectors_.end(); ++Vtmp)
  {
    double len1 = *v2idx1;
    DataSet_Matrix::iterator v2idx2 = v2idx1;
    for (std::vector<DataSet_Vector*>::iterator Vtmp2 = Vtmp;
                                                Vtmp2 != IredVectors_.end(); ++Vtmp2)
    {
      double len2 = *(v2idx2++);
      double legendre = LegendrePoly(order_, (*Vtmp)->Dot( *(*Vtmp2) ) / (len1 * len2) );
      *(mat++) += legendre;
      if (Vtmp == Vtmp2)
        *(v1idx++) += legendre;
    }
    ++v2idx1;
  }
}

/** Calc Distance Matrix */
void Action_Matrix::CalcDistanceMatrix(Frame* currentFrame) {
  DataSet_Matrix::iterator mat = Mat_->begin();
  if (!useMask2_) {
    // Upper Triangle
    for (AtomMask::const_iterator atom2 = mask1_.begin(); atom2 != mask1_.end(); ++atom2)
      for (AtomMask::const_iterator atom1 = atom2; atom1 != mask1_.end(); ++atom1)
        *(mat++) += sqrt(DIST2_NoImage(currentFrame->XYZ(*atom2), currentFrame->XYZ(*atom1)));
  } else {
    // Full matrix
    for (AtomMask::const_iterator atom2 = mask2_.begin(); atom2 != mask2_.end(); ++atom2)
      for (AtomMask::const_iterator atom1 = mask1_.begin(); atom1 != mask1_.end(); ++atom1)
        *(mat++) += sqrt(DIST2_NoImage(currentFrame->XYZ(*atom2), currentFrame->XYZ(*atom1)));
  }
}

// Action_Matrix::StoreVec()
void Action_Matrix::StoreVec(DataSet_Matrix::iterator& v1, DataSet_Matrix::iterator& v2,
                             const double* XYZ) 
{
  *(v1++) += XYZ[0];
  *(v2++) += (XYZ[0] * XYZ[0]);
  *(v1++) += XYZ[1];
  *(v2++) += (XYZ[1] * XYZ[1]);
  *(v1++) += XYZ[2];
  *(v2++) += (XYZ[2] * XYZ[2]);
}

/** Calc Covariance Matrix */
void Action_Matrix::CalcCovarianceMatrix(Frame* currentFrame) {
  DataSet_Matrix::iterator mat = Mat_->begin();
  DataSet_Matrix::iterator v1idx1 = Mat_->v1begin();
  DataSet_Matrix::iterator v2idx1 = Mat_->v2begin();
  if (useMask2_) {
    // Full Matrix
    // Position for mask2 halfway through vect/vect2
    DataSet_Matrix::iterator v1idx2 = Mat_->v1begin() + (mask1_.Nselected() * 3); 
    DataSet_Matrix::iterator v2idx2 = Mat_->v2begin() + (mask1_.Nselected() * 3); 
    bool storeVecj = true; // Only store vecj|vecj^2 first time through inner loop
    // OUTER LOOP
    for (AtomMask::const_iterator atom2 = mask2_.begin(); atom2 != mask2_.end(); ++atom2)
    {
      const double* XYZi = currentFrame->XYZ( *atom2 );
      // Store veci and veci^2
      StoreVec(v1idx2, v2idx2, XYZi);
      // Loop over X, Y, and Z of veci
      for (int iidx = 0; iidx < 3; ++iidx) {
        double Vi = XYZi[iidx];
        // INNER LOOP
        for (AtomMask::const_iterator atom1 = mask1_.begin(); atom1 != mask1_.end(); ++atom1)
        {
          const double* XYZj = currentFrame->XYZ( *atom1 );
          // Store vecj and vecj^2, first time through only
          if (storeVecj) 
            StoreVec(v1idx1, v2idx1, XYZj);
          *(mat++) += Vi * XYZj[0];
          *(mat++) += Vi * XYZj[1];
          *(mat++) += Vi * XYZj[2];
        }
        storeVecj = false;
      }
    }
  } else {
    // Half Matrix
    // OUTER LOOP
    for (AtomMask::const_iterator atom2 = mask1_.begin(); atom2 != mask1_.end(); ++atom2)
    {
      const double* XYZi = currentFrame->XYZ( *atom2 );
      // Store veci and veci^2
      StoreVec(v1idx1, v2idx1, XYZi);
      // Loop over X, Y, and Z of veci
      for (int iidx = 0; iidx < 3; ++iidx) {
        double Vi = XYZi[iidx];
        // INNER LOOP
        for (AtomMask::const_iterator atom1 = atom2; atom1 != mask1_.end(); ++atom1)
        {
          const double* XYZj = currentFrame->XYZ( *atom1 );
          if ( atom1 == atom2 ) { // TODO: This saves 3 doubles / frame. Worth it?
            for (int jidx = iidx; jidx < 3; ++jidx)
              *(mat++) += Vi * XYZj[jidx]; // Vi * j{0,1,2}, Vi * j{1,2}, Vi * j{2}
          } else {
            *(mat++) += Vi * XYZj[0];
            *(mat++) += Vi * XYZj[1];
            *(mat++) += Vi * XYZj[2];
          }
        }
      }
    }
  }
}

/** Calc Isotropically distributed ensemble matrix.
  * See Proteins 2002, 46, 177; eq. 7 
  */
void Action_Matrix::CalcIdeaMatrix(Frame* currentFrame) {
  DataSet_Matrix::iterator mat = Mat_->begin();
  DataSet_Matrix::iterator v1idx1 = Mat_->v1begin();
  DataSet_Matrix::iterator v2idx1 = Mat_->v2begin();
  // Get COM
  Vec3 COM = currentFrame->VCenterOfMass( mask1_ );
  // Get ri, rj, and calc ri*rj
  // Matrix IDEA only uses 1 mask.
  for (AtomMask::const_iterator atomi = mask1_.begin(); atomi != mask1_.end(); ++atomi)
  {
    Vec3 ri = currentFrame->XYZ(*atomi);
    ri -= COM;
    for (AtomMask::const_iterator atomj = atomi; atomj != mask1_.end(); ++atomj)
    {
      Vec3 rj = currentFrame->XYZ(*atomj);
      rj -= COM;
      double val = ri * rj;
      *(mat++) += val;
      if (atomj == atomi) {
        *(v1idx1++) += val;
        *(v2idx1++) += (val * val);
      }
    }
  }
}

/** Calc correlation matrix. */
void Action_Matrix::CalcCorrelationMatrix(Frame* currentFrame) {
  DataSet_Matrix::iterator mat = Mat_->begin();
  DataSet_Matrix::iterator v1idx1 = Mat_->v1begin();
  DataSet_Matrix::iterator v2idx1 = Mat_->v2begin();
  if (useMask2_) {
    // Full Matrix
    // Position for mask2 halfway through vect/vect2
    DataSet_Matrix::iterator v1idx2 = Mat_->v1begin() + (mask1_.Nselected() * 3);
    DataSet_Matrix::iterator v2idx2 = Mat_->v2begin() + (mask1_.Nselected() * 3);
    bool storeVecj = true; // Only store vecj|vecj^2 first time through inner loop
    // OUTER LOOP
    for (AtomMask::const_iterator atom2 = mask2_.begin(); atom2 != mask2_.end(); ++atom2)
    {
      const double* XYZi = currentFrame->XYZ( *atom2 );
      // Store veci and veci^2
      StoreVec(v1idx2, v2idx2, XYZi);
      // INNER LOOP
      for (AtomMask::const_iterator atom1 = mask1_.begin(); atom1 != mask1_.end(); ++atom1)
      {
        const double* XYZj = currentFrame->XYZ( *atom1 );
        // Store vecj and vecj^2, first time through only
        if (storeVecj)
          StoreVec(v1idx1, v2idx1, XYZj);
        *(mat++) += (XYZi[0]*XYZj[0] + XYZi[1]*XYZj[1] + XYZi[2]*XYZj[2]);
      }
      storeVecj = false;
    }
  } else {
    // Half Matrix
    // OUTER LOOP
    for (AtomMask::const_iterator atom2 = mask1_.begin(); atom2 != mask1_.end(); ++atom2)
    {
      const double* XYZi = currentFrame->XYZ( *atom2 );
      // Store veci and veci^2
      StoreVec(v1idx1, v2idx1, XYZi);
      for (AtomMask::const_iterator atom1 = atom2; atom1 != mask1_.end(); ++atom1)
      {
        const double* XYZj = currentFrame->XYZ( *atom1 );
        *(mat++) += (XYZi[0]*XYZj[0] + XYZi[1]*XYZj[1] + XYZi[2]*XYZj[2]);
      }
    }
  }
}

/** Calculate distance covariance matrix. */
void Action_Matrix::CalcDistanceCovarianceMatrix(Frame* currentFrame) {
  // Calculate all distance pairs for mask 1
  DataSet_Matrix::iterator pair_j = Mat_->v2begin();
  AtomMask::const_iterator mask1end = mask1_.end() - 1;
  for (AtomMask::const_iterator atom1 = mask1_.begin(); atom1 != mask1end; ++atom1)
    for (AtomMask::const_iterator atom2 = atom1 + 1; atom2 != mask1_.end(); ++atom2)
      *(pair_j++) = sqrt(DIST2_NoImage(currentFrame->XYZ(*atom1), currentFrame->XYZ(*atom2)));
  // Create matrix from all distance pairs
  DataSet_Matrix::iterator mat = Mat_->begin();
  DataSet_Matrix::iterator v1idx = Mat_->v1begin();
  for (DataSet_Matrix::iterator pair_i = Mat_->v2begin();
                                pair_i != Mat_->v2end(); ++pair_i)
  {
    for (pair_j = pair_i; pair_j != Mat_->v2end(); ++pair_j) {
      *(mat++) += ( (*pair_i) * (*pair_j) );
      if (pair_i == pair_j)
        *(v1idx++) += (*pair_i);
    }
  }
}

// Action_Matrix::action()
Action::RetType Action_Matrix::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  // Check if this frame should be processed
  if ( CheckFrameCounter( frameNum ) ) return Action::OK;
  // Increment number of snapshots
  ++snap_; 

  switch (type_) {
    case DataSet_Matrix::DIST     : CalcDistanceMatrix(currentFrame); break;
    case DataSet_Matrix::COVAR    :
    case DataSet_Matrix::MWCOVAR  : CalcCovarianceMatrix(currentFrame); break;
    case DataSet_Matrix::CORREL   : CalcCorrelationMatrix(currentFrame); break;
    case DataSet_Matrix::DISTCOVAR: CalcDistanceCovarianceMatrix(currentFrame); break;
    case DataSet_Matrix::IDEA     : CalcIdeaMatrix(currentFrame); break;
    case DataSet_Matrix::IRED     : CalcIredMatrix(); break;
    default: return Action::ERR;
  }

  return Action::OK;
}

// -----------------------------------------------------------------------------
// Action_Matrix::FinishCovariance()
void Action_Matrix::FinishCovariance() {
  double Mass = 1.0;
  double mass2 = 1.0;
  DataSet_Matrix::iterator mat = Mat_->begin();
  // Calc <riri> - <ri><ri>
  Mat_->Vect2MinusVect(); // TODO: Is this ever used??
  // Calc <rirj> - <ri><rj>
  if (useMask2_) {
    // Full Matrix
    std::vector<double>::iterator m2 = mass2_.begin();
    // Position for mask2 halfway through vect/vect2
    DataSet_Matrix::iterator v1idx2begin = Mat_->v1begin() + Mat_->Ncols();
    for (DataSet_Matrix::iterator v1idx2 = v1idx2begin; 
                                  v1idx2 != Mat_->v1end(); v1idx2 += 3)
    {
      if (type_ == DataSet_Matrix::MWCOVAR)
        mass2 = *(m2++);
      for (int iidx = 0; iidx < 3; ++iidx) {
        std::vector<double>::iterator m1 = mass1_.begin();
        double Vi = *(v1idx2 + iidx);
        for (DataSet_Matrix::iterator v1idx1 = Mat_->v1begin();
                                      v1idx1 != v1idx2begin; v1idx1 += 3)
        {
          if (type_ == DataSet_Matrix::MWCOVAR)
            Mass = sqrt( mass2 * *(m1++) );
          *mat = (*mat - (Vi * *(v1idx1  ))) * Mass;
          ++mat;
          *mat = (*mat - (Vi * *(v1idx1+1))) * Mass;
          ++mat;
          *mat = (*mat - (Vi * *(v1idx1+2))) * Mass;
          ++mat;
          //*(mat++) -= Vi * *(v1idx1  );
          //*(mat++) -= Vi * *(v1idx1+1);
          //*(mat++) -= Vi * *(v1idx1+2);
        }
      }
    }
  } else {
    // Half Matrix
    std::vector<double>::iterator m2 = mass1_.begin();
    for (DataSet_Matrix::iterator v1idx2 = Mat_->v1begin();
                                  v1idx2 != Mat_->v1end(); v1idx2 += 3)
    {
      if (type_ == DataSet_Matrix::MWCOVAR)
        mass2 = *m2;
      for (int iidx = 0; iidx < 3; ++iidx) {
        std::vector<double>::iterator m1 = m2;
        double Vi = *(v1idx2 + iidx);
        for (DataSet_Matrix::iterator v1idx1 = v1idx2;
                                      v1idx1 != Mat_->v1end(); v1idx1 += 3)
        {
          if (type_ == DataSet_Matrix::MWCOVAR)
            Mass = sqrt( mass2 * *(m1++) );
          if ( v1idx1 == v1idx2 ) {
            for (int jidx = iidx; jidx < 3; ++jidx) {
              *mat = (*mat - (Vi * *(v1idx1 + jidx))) * Mass;
              ++mat;
              //*(mat++) -= Vi * *(v1idx1 + jidx); // Vi * j{0,1,2}, Vi * j{1,2}, Vi * j{2}
            }
          } else {
            *mat = (*mat - (Vi * *(v1idx1  ))) * Mass;
            ++mat;
            *mat = (*mat - (Vi * *(v1idx1+1))) * Mass;
            ++mat;
            *mat = (*mat - (Vi * *(v1idx1+2))) * Mass;
            ++mat;
          }
        }
      }
      ++m2;
    }
  }
}

// DotProdAndNorm()
static inline void DotProdAndNorm(DataSet_Matrix::iterator& mat,
                                  DataSet_Matrix::iterator& vecti, 
                                  DataSet_Matrix::iterator& vectj,
                                  DataSet_Matrix::iterator& vect2i,
                                  DataSet_Matrix::iterator& vect2j)
{
  *(mat) -= ( *(vectj  ) * *(vecti  ) +
              *(vectj+1) * *(vecti+1) +
              *(vectj+2) * *(vecti+2)   );
  // Normalize
  *(mat++) /= sqrt( (*(vect2j) + *(vect2j+1) + *(vect2j+2)) *
                    (*(vect2i) + *(vect2i+1) + *(vect2i+2))   );
}

// Action_Matrix::FinishCorrelation()
void Action_Matrix::FinishCorrelation() {
  DataSet_Matrix::iterator mat = Mat_->begin();
  // Calc <ri * ri> - <ri> * <ri>
  Mat_->Vect2MinusVect(); // TODO: Is this only for Correlation? 
  // Calc <ri * rj> - <ri> * <rj>
  if (useMask2_) {
    // Full Matrix
    // Position for mask2 halfway through vect/vect2
    // Vect has 3 entries per atom, but matrix only has 1 (as opposed to 
    // COVAR/MWCOVAR which has 3 for vect and matrix).
    DataSet_Matrix::iterator v1idx2begin = Mat_->v1begin() + Mat_->Ncols() * 3;
    DataSet_Matrix::iterator v2idx2      = Mat_->v2begin() + Mat_->Ncols() * 3;
    for (DataSet_Matrix::iterator v1idx2 = v1idx2begin;
                                  v1idx2 != Mat_->v1end(); v1idx2 += 3)
    {
      DataSet_Matrix::iterator v2idx1 = Mat_->v2begin();
      for (DataSet_Matrix::iterator v1idx1 = Mat_->v1begin();
                                    v1idx1 != v1idx2begin; v1idx1 += 3)
      {
        DotProdAndNorm( mat, v1idx2, v1idx1, v2idx2, v2idx1 );
        v2idx1 += 3;
      }
      v2idx2 += 3;
    }
  } else {
    // Half Matrix
    DataSet_Matrix::iterator v2idx2 = Mat_->v2begin();
    for (DataSet_Matrix::iterator v1idx2 = Mat_->v1begin();
                                  v1idx2 != Mat_->v1end(); v1idx2 += 3)
    {
      DataSet_Matrix::iterator v2idx1 = v2idx2;
      for (DataSet_Matrix::iterator v1idx1 = v1idx2;
                                    v1idx1 != Mat_->v1end(); v1idx1 += 3)
      {
        DotProdAndNorm( mat, v1idx2, v1idx1, v2idx2, v2idx1 );
        v2idx1 += 3;
      }
      v2idx2 += 3;
    }
  }
}

// Action_Matrix::FinishDistanceCovariance()
void Action_Matrix::FinishDistanceCovariance() {
  DataSet_Matrix::iterator mat = Mat_->begin();
  for (DataSet_Matrix::iterator pair_i = Mat_->v1begin();
                                pair_i != Mat_->v1end(); ++pair_i)
    for (DataSet_Matrix::iterator pair_j = pair_i; 
                                  pair_j != Mat_->v1end(); ++pair_j) 
      *(mat++) -= ((*pair_i) * (*pair_j));
}

// Action_Matrix::print()
void Action_Matrix::Print() {
  // ---------- Calculate average over number of sets ------
  Mat_->DivideBy((double)snap_);

  switch (type_) {
    case DataSet_Matrix::COVAR    :
    case DataSet_Matrix::MWCOVAR  : FinishCovariance(); break;
    case DataSet_Matrix::CORREL   : FinishCorrelation(); break;
    case DataSet_Matrix::DISTCOVAR: FinishDistanceCovariance(); break;
    case DataSet_Matrix::IDEA     : Mat_->DivideBy(3.0); break;
    default: break; 
  }

  // If byres/bymask output is desired, write the current matrix now since 
  // it is not stored in DataFileList.
  // TODO: Convert byres to new matrix so it can be output formatted.
  // TODO: Check that currentParm is still valid?
  if (!filename_.empty() && outtype_ == BYRESIDUE) {
    // ---------- Print out BYRESIDUE
    CpptrajFile outfile;
    outfile.OpenWrite(filename_);
    // Convert masks to char masks in order to check whether an atom
    // is selected.
    mask1_.ConvertToCharMask();
    if (useMask2_)
      mask2_.ConvertToCharMask();
    else
      mask2_ = mask1_;
    // Loop over residue pairs
    int crow = 0;
    double mass = 1.0;
    for (Topology::res_iterator resi = CurrentParm_->ResStart(); 
                                resi != CurrentParm_->ResEnd(); ++resi) { // Row
      bool printnewline = false;
      int crowold = crow;
      int ccol = 0;
      for (Topology::res_iterator resj = CurrentParm_->ResStart();
                                  resj != CurrentParm_->ResEnd(); ++resj) { // Column
        bool printval = false;
        double val = 0;
        double valnorm = 0;
        crow = crowold;
        int ccolold = ccol;
        for (int atomi = (*resi).FirstAtom(); atomi < (*resi).LastAtom(); ++atomi)
        {
          if ( mask2_.AtomInCharMask(atomi) ) {
            ccol = ccolold;
            for (int atomj = (*resj).FirstAtom(); atomj < (*resj).LastAtom(); ++atomj)
            {
              if ( mask1_.AtomInCharMask(atomj) ) {
                if (useMass_)
                  mass = (*CurrentParm_)[atomi].Mass() * (*CurrentParm_)[atomj].Mass();
                valnorm += mass;
                printval = printnewline = true;
                //mprintf("Res %i-%i row=%i col=%i\n",resi+1,resj+1,crow,ccol);
                val += (Mat_->GetElement( crow, ccol ) * mass);
                ++ccol;
              }
            }
            ++crow;
          }
        }
        if (printval) outfile.Printf("%6.2f ",val / valnorm);
      }
      if (printnewline) outfile.Printf("\n");
    }
    outfile.CloseFile();
  } else if (!filename_.empty() && outtype_ == BYMASK) {
    // ---------- Print out BYMASK
    CpptrajFile outfile;
    outfile.OpenWrite(filename_);
    // If only 1 mask, internal average over mask1, otherwise
    //   i==0: mask1/mask1 
    //   i==1: mask1/mask2 
    //   i==2: mask2/mask2
    int iend;
    if (!useMask2_)
      iend = 1;
    else
      iend = 3;
    AtomMask::const_iterator maskAbegin = mask1_.begin();
    AtomMask::const_iterator maskAend = mask1_.end();
    AtomMask::const_iterator maskBbegin = mask1_.begin();
    AtomMask::const_iterator maskBend = mask1_.end();
    double mass = 1.0;
    for (int i = 0; i < iend; ++i) {
      if (i > 0) {
        maskAbegin = maskBbegin;
        maskAend = maskBend;
        maskBbegin = mask2_.begin();
        maskBend = mask2_.end();
      }
      double val = 0;
      double valnorm = 0;
      int crow = 0;
      for (AtomMask::const_iterator atomj = maskBbegin; atomj != maskBend; ++atomj)
      {
        int ccol = 0;
        for (AtomMask::const_iterator atomi = maskAbegin; atomi != maskAend; ++atomi)
        {
          if (useMass_)
            mass = (*CurrentParm_)[*atomj].Mass() * (*CurrentParm_)[*atomi].Mass();
          valnorm += mass;
          val += (Mat_->GetElement( crow, ccol ) * mass);
          ++ccol;
        }
        ++crow;
      }
      outfile.Printf("%6.2f ", val / valnorm);
    }
    outfile.Printf("\n");
    outfile.CloseFile();
  }

  // Process output file args
  if (outfile_ != 0) {
    outfile_->ProcessArgs("xlabel Atom");
    outfile_->ProcessArgs("square2d noxcol noheader");
  }
  return;
}
