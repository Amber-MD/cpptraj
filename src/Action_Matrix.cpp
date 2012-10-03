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
    case BYATOM:    mprintf("by atom"); break;
    case BYRESIDUE: mprintf("by residue"); break;
    case BYMASK:    mprintf("by mask"); break;
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

// Action_Matrix::setup()
int Action_Matrix::setup() {
  // Set up masks
  if (type_!=DataSet_Matrix::MATRIX_IRED) {
    if (currentParm->SetupIntegerMask(mask1_)) return 1;
    mask1_.MaskInfo();
    if (mask1_.None()) {
      mprinterr("Error: No atoms selected for mask1.\n");
      return 1;
    }
    if (useMask2_) {
      if (currentParm->SetupIntegerMask(mask2_)) return 1;
      mask2_.MaskInfo(); 
      if (mask2_.None()) {
        mprinterr("Error: No atoms selected for mask2.\n");
        return 1;
      }
    }
  }

  // Allocate matrix memory.
  if (Mat_->MatrixAlloc(mask1_, mask2_, currentParm->Atoms())) {
    mprinterr("Error: matrix: Could not allocate memory.\n");
    return 1;
  }

  return 0;
}

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

// Action_Matrix::action()
int Action_Matrix::action() {
  // If the current frame is less than start exit
  if (frameNum < start_) return 0;
  // If the current frame is greater than stop exit
  if (stop_!=-1 && frameNum >= stop_) return 0;
  // Increment number of snapshots, update next target frame
  Mat_->IncrementSnap(); 
  start_ += offset_;

  DataSet_Matrix::iterator mat = Mat_->begin();
  switch (type_) {
    // ---------------------------------------------
    // ** Calc Distance Matrix **
    case DataSet_Matrix::MATRIX_DIST:
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
      break;
    default: return 1;
  }

  return 0;
}

void Action_Matrix::print() {
  // ---------- Calculate average over number of sets ------
  Mat_->AverageOverSnapshots();

  if (outfile_ != 0) {
    outfile_->ProcessArgs("xlabel Atom");
  }
  return;
}
