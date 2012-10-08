#include <cmath> // sqrt
#include "Action_Matrix.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"

// CONSTRUCTOR
Action_Matrix::Action_Matrix() :
  Mat_(0),
  outfile_(0),
  outtype_(BYATOM),
  type_(DataSet_Matrix::MATRIX_NULL),
  start_(0),
  stop_(-1),
  offset_(1),
  order_(2),
  useMask2_(false)
{}

const char Action_Matrix::MatrixModeString[8][27] = {
  "UNDEFINED",
  "distance matrix",
  "covar matrix",
  "mass weighted covar matrix",
  "correlation matrix",
  "distance covar matrix",
  "idea matrix",
  "ired matrix"
};

// Action_Matrix::init()
int Action_Matrix::init() {
  // Get Keywords
  filename_ = actionArgs.GetStringKey("out");
  start_ = actionArgs.getKeyInt("start", 1);
  stop_ = actionArgs.getKeyInt("stop", -1);
  if (stop_ == -1)
    stop_ = actionArgs.getKeyInt("end", -1);
  offset_ = actionArgs.getKeyInt("offset", 1);
  // Actual start frame arg should be from 0
  --start_;
  // Determine matrix type
  type_ = DataSet_Matrix::MATRIX_DIST;
  if (actionArgs.hasKey("distcovar"))
    type_ = DataSet_Matrix::MATRIX_DISTCOVAR;
  else if (actionArgs.hasKey("mwcovar"))
    type_ = DataSet_Matrix::MATRIX_MWCOVAR;
  else if (actionArgs.hasKey("dist"))
    type_ = DataSet_Matrix::MATRIX_DIST;
  else if (actionArgs.hasKey("covar"))
    type_ = DataSet_Matrix::MATRIX_COVAR;
  else if (actionArgs.hasKey("correl"))
    type_ = DataSet_Matrix::MATRIX_CORREL;
  else if (actionArgs.hasKey("idea"))
    type_ = DataSet_Matrix::MATRIX_IDEA;
  else if (actionArgs.hasKey("ired"))
    type_ = DataSet_Matrix::MATRIX_IRED;
  // Output type
  if (actionArgs.hasKey("byres"))
    outtype_ = BYRESIDUE;
  else if (actionArgs.hasKey("bymask"))
    outtype_ = BYMASK;
  else if (actionArgs.hasKey("byatom"))
    outtype_ = BYATOM;
  else
    outtype_ = BYATOM;
  // Get matrix name
  std::string name = actionArgs.GetStringKey("name");
  // UseMass
  useMass_ = actionArgs.hasKey("mass");

  if (type_ != DataSet_Matrix::MATRIX_IRED) {
    // Get masks if not IRED
    mask1_.SetMaskString( actionArgs.GetMaskNext() );
    std::string maskexpr = actionArgs.GetMaskNext();
    if (!maskexpr.empty()) useMask2_ = true;
    if ( type_ == DataSet_Matrix::MATRIX_IDEA ||
         type_ == DataSet_Matrix::MATRIX_DISTCOVAR )
    {
      mprintf("Warning: Mask 2 specified but not used for %s\n",MatrixModeString[type_]);
      useMask2_ = false;
    }
    if (useMask2_)
      mask2_.SetMaskString( maskexpr );
  } else {
    // Setup IRED vectors and determine Legendre order
    order_ = actionArgs.getKeyInt("order",1);
    if (order_ <= 0) {
      mprinterr("Error: matrix: order parameter <= 0, ignoring command\n");
      return 1;
    }
    DataSet_Vector* Vtmp;
    DSL->VectorBegin();
    while ( (Vtmp = (DataSet_Vector*)DSL->NextVector()) != 0 ) {
      if (Vtmp->IsIred()) 
        IredVectors_.push_back( Vtmp );
    }
    if (IredVectors_.empty()) {
      mprinterr("Error: matrix: no vectors defined for IRED\n");
      return 1;
    }
  }
 
  // Set up matrix DataSet and type
  Mat_ = (DataSet_Matrix*)DSL->AddSet(DataSet::MATRIX, name, "Mat");
  if (Mat_ == 0) return 1;
  Mat_->SetType( type_ );
  // Add set to output file if not doing ptraj-compatible output
  //if (!ptrajoutput_)
    outfile_ = DFL->AddSetToFile(filename_, Mat_);

  mprintf("    MATRIX: Calculating %s", MatrixModeString[ type_ ]);
  switch (outtype_) {
    case BYATOM:    mprintf(" by atom"); break;
    case BYRESIDUE: mprintf(" by residue"); break;
    case BYMASK:    mprintf(" by mask"); break;
  }
  if (useMass_)
    mprintf(", using mass weighting\n");
  else
    mprintf(", using no mass weighting\n");

  if (type_==DataSet_Matrix::MATRIX_IRED)
    mprintf("            Order of Legendre polynomials: %i\n",order_);
  if (!filename_.empty())
    mprintf("            Printing to file %s\n",filename_.c_str());
  if (!name.empty())
    mprintf("            Storing matrix on internal stack with name: %s\n", 
            Mat_->Legend().c_str());
  if (start_!=0 || stop_!=-1 || offset_!=1) {
    mprintf("            start: %i  stop:",start_);
    if (stop_==-1)
      mprintf(" Final frame");
    else
      mprintf(" %i",stop_);
    mprintf("  offset: %i\n",offset_);
  }
  if (type_!=DataSet_Matrix::MATRIX_IRED) {
    mprintf("            Mask1: %s\n",mask1_.MaskString());
    if (useMask2_)
      mprintf("            Mask2: %s\n",mask2_.MaskString());
  }

  mprintf("\n");

  return 0;
}

// Action_Matrix::FillMassArray()
void Action_Matrix::FillMassArray(std::vector<double>& mass, AtomMask& mask) {
  mass.clear();
  if (type_ == DataSet_Matrix::MATRIX_MWCOVAR) {
    mass.reserve( mask.Nselected() );
    for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom) 
      mass.push_back( (*currentParm)[ *atom ].Mass() );
  }
}

// Action_Matrix::setup()
int Action_Matrix::setup() {
  // Set up masks. Store mass info for masks if MWCOVAR
  if (type_!=DataSet_Matrix::MATRIX_IRED) {
    if (currentParm->SetupIntegerMask(mask1_)) return 1;
    mask1_.MaskInfo();
    if (mask1_.None()) {
      mprinterr("Error: No atoms selected for mask1.\n");
      return 1;
    }
    FillMassArray(mass1_, mask1_); // MWCOVAR only
    if (useMask2_) {
      if (currentParm->SetupIntegerMask(mask2_)) return 1;
      mask2_.MaskInfo(); 
      if (mask2_.None()) {
        mprinterr("Error: No atoms selected for mask2.\n");
        return 1;
      }
      FillMassArray(mass2_, mask2_); // MWCOVAR only
    }
  }

  // Allocate matrix memory.
  if (Mat_->MatrixAlloc(mask1_, mask2_, currentParm->Atoms())) {
    mprinterr("Error: matrix: Could not allocate memory.\n");
    return 1;
  }

  return 0;
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

/** Calc Distance Matrix */
void Action_Matrix::CalcDistanceMatrix() {
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
void Action_Matrix::CalcCovarianceMatrix() {
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
void Action_Matrix::CalcIdeaMatrix() {
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
void Action_Matrix::CalcCorrelationMatrix() {
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

// Action_Matrix::action()
int Action_Matrix::action() {
  // If the current frame is less than start exit
  if (frameNum < start_) return 0;
  // If the current frame is greater than stop exit
  if (stop_!=-1 && frameNum >= stop_) return 0;
  // Increment number of snapshots, update next target frame
  Mat_->IncrementSnap(); 
  start_ += offset_;

  switch (type_) {
    case DataSet_Matrix::MATRIX_DIST: CalcDistanceMatrix(); break;
    case DataSet_Matrix::MATRIX_COVAR:
    case DataSet_Matrix::MATRIX_MWCOVAR: CalcCovarianceMatrix(); break;
    case DataSet_Matrix::MATRIX_CORREL: CalcCorrelationMatrix(); break;
    case DataSet_Matrix::MATRIX_IDEA: CalcIdeaMatrix(); break;

    default: return 1;
  }

  return 0;
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
      if (type_ == DataSet_Matrix::MATRIX_MWCOVAR)
        mass2 = *(m2++);
      for (int iidx = 0; iidx < 3; ++iidx) {
        std::vector<double>::iterator m1 = mass1_.begin();
        double Vi = *(v1idx2 + iidx);
        for (DataSet_Matrix::iterator v1idx1 = Mat_->v1begin();
                                      v1idx1 != v1idx2begin; v1idx1 += 3)
        {
          if (type_ == DataSet_Matrix::MATRIX_MWCOVAR)
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
      if (type_ == DataSet_Matrix::MATRIX_MWCOVAR)
        mass2 = *m2;
      for (int iidx = 0; iidx < 3; ++iidx) {
        std::vector<double>::iterator m1 = m2;
        double Vi = *(v1idx2 + iidx);
        for (DataSet_Matrix::iterator v1idx1 = v1idx2;
                                      v1idx1 != Mat_->v1end(); v1idx1 += 3)
        {
          if (type_ == DataSet_Matrix::MATRIX_MWCOVAR)
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

// Action_Matrix::print()
void Action_Matrix::print() {
  // ---------- Calculate average over number of sets ------
  Mat_->DivideBy((double)Mat_->Nsnap());

  switch (type_) {
    case DataSet_Matrix::MATRIX_COVAR:
    case DataSet_Matrix::MATRIX_MWCOVAR: FinishCovariance(); break;
    case DataSet_Matrix::MATRIX_CORREL: FinishCorrelation(); break;
    case DataSet_Matrix::MATRIX_IDEA: Mat_->DivideBy(3.0); break;

    default: break; 
  }

  // Process output file args
  if (outfile_ != 0) {
    outfile_->ProcessArgs("xlabel Atom");
    outfile_->ProcessArgs("square2d noxcol noheader");
  }
  return;
}
