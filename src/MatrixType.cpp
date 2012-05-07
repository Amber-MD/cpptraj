#include "MatrixType.h"
#include "VectorType.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
MatrixType::MatrixType() : 
  type_(MATRIX_NULL),
  vect_(0),
  vect2_(0),
  mat_(0),
  vectsize_(0),
  matsize_(0),
  mask1tot_(0),
  mask2tot_(0),
  snap_(0),
  start_(0),
  stop_(0),
  offset_(0),
  order_(0),
  outtype_(BYATOM)
{}

// DESTRUCTOR
MatrixType::~MatrixType() {
  if (vect_!=0) delete[] vect_;
  if (vect2_!=0) delete[] vect2_;
  if (mat_!=0) delete[] mat_;
}

const char MatrixType::MatrixModeString[8][27] = { 
  "UNDEFINED",
  "distance matrix",
  "covar matrix",
  "mass weighted covar matrix",
  "correlation matrix",
  "distance covar matrix",
  "idea matrix",
  "ired matrix"
};

/** 
  *  matrix dist|covar|mwcovar|distcovar|correl|idea|ired
  *                                             [name <name>] [order <order>]
  *                                             [<mask1>] [<mask2>] [out <filename>] 
  *                                             [start <start>] [stop <stop>] [offset <offset>]
  *                                             [byatom|byres|bymask] [mass]
  * 
  *  - If MATRIX_IRED, mask1 and mask2 are ignored and the number of matrix elements
  *    is determined by the number of vector definitions given PRIOR to the
  *    matrix command. Here, only the "upper right half" of the matrix is allocated.
  *  - Otherwise:
  *    - Upon input, ||mask1|| >= ||mask2||; this is checked below.
  *    - If only mask1 (or none) is given, only the "upper right half" of the matrix
  *      is allocated, including the main diagonal.
  *      Non-squared elements ii are contained in "vect", squared are in "vect2".
  *      This is done to be consistent if mask1 and mask2 is given and mask1 != mask2.
  *      In the case of MATRIX_DISTCOVAR and MATRIX_IRED, only (mask1tot * (mask1tot - 1)/2) 
  *      resp. mask1tot elements of "vect" are used; "vect2" acts as a temporary array to store 
  *      distances resp. vector lengths for each snapshot.
  *    - If both mask1 and mask2 are given, the full matrix is allocated, assuming that atoms
  *      in both masks do not necessarily correspond. (To generate full, symmetric matrices, 
  *      call the function with mask1 == mask2 upon input.)
  *    - The matrix will be stored internally with the name "name" on the matrixStack for later
  *      processing (w/ the "analyze matrix" command) ONLY if mask1 (or none) is given.
  * 
  *  - For "covar, mwcovar, distcovar, idea, ired", only "byatom" output may be chosen.
  *  - Since "distcovar, idea, ired" is mainly intended for subsequent analysis with 
  *    "analyze matrix", only input of mask1 (or none) is possible.
  */
int MatrixType::init(  ) {
  filename_ = actionArgs.GetStringKey("out");

  start_ = actionArgs.getKeyInt("start", 1);
  stop_ = actionArgs.getKeyInt("stop", -1);
  if (stop_ == -1)
    stop_ = actionArgs.getKeyInt("end", -1);
  offset_ = actionArgs.getKeyInt("offset", 1);
  // Actual start arg should be from 0
  --start_;

  order_ = actionArgs.getKeyInt("order",1);
  if (order_ <= 0) {
    mprinterr("Error: matrix: order parameter <= 0, ignoring command\n");
    return 1;
  }

  if (actionArgs.hasKey("byres"))
    outtype_ = BYRESIDUE;
  else if (actionArgs.hasKey("bymask"))
    outtype_ = BYMASK;
  else if (actionArgs.hasKey("byatom"))
    outtype_ = BYATOM;
  else
    outtype_ = BYATOM;

  useMass_ = actionArgs.hasKey("mass");

  // Determine matrix type
  if (actionArgs.hasKey("distcovar"))
    type_ = MATRIX_DISTCOVAR;
  else if (actionArgs.hasKey("mwcovar"))
    type_ = MATRIX_MWCOVAR;
  else if (actionArgs.hasKey("dist"))
    type_ = MATRIX_DIST;
  else if (actionArgs.hasKey("covar"))
    type_ = MATRIX_COVAR;
  else if (actionArgs.hasKey("correl"))
    type_ = MATRIX_CORREL;
  else if (actionArgs.hasKey("idea"))
    type_ = MATRIX_IDEA;
  else if (actionArgs.hasKey("ired"))
    type_ = MATRIX_IRED;
  else
    type_ = MATRIX_DIST;

  // Sanity check
  if (type_==MATRIX_NULL) return 1;

  // Get Masks
  mask1tot_ = 0;
  char* maskexpr = actionArgs.getNextMask();
  // IRED Setup
  if (type_ == MATRIX_IRED) {
    if (maskexpr!=NULL) {
      mprinterr("Error: matrix: mask input does not work with ired\n");
      return 1;
    }
    // Count the number of previously defined IRED vectors.
    VectorType *Vtmp;
    DSL->VectorBegin();
    while ( (Vtmp = (VectorType*)DSL->NextVector()) != 0 ) {
      if (Vtmp->Mode() == VectorType::VECTOR_IRED)
        ++mask1tot_;
    }
    if (mask1tot_==0) {
      mprinterr("Error: matrix: no vector defined for IRED\n");
      return 1;
    }
  } else {
    mask1_.SetMaskString( maskexpr );
    maskexpr = actionArgs.getNextMask();
    if (maskexpr!=NULL)
      mask2_.SetMaskString( maskexpr );
  } 

  // Get matrix name
  // NOTE: Unlike ptraj where this was done after the 'mass' keyword check,
  //       here it is done after all other keywords are processed. 
  name_ = actionArgs.GetStringNext();

  // Check arguments
  if ( !name_.empty() && mask2_.MaskStringSet() ) {
    mprintf("Error: matrix: matrix only stored if no mask2\n");
    return 1;
  }

  if ( (type_ == MATRIX_COVAR ||
        type_ == MATRIX_MWCOVAR ||
        type_ == MATRIX_IRED)      && outtype_ != BYATOM)
  {
    mprinterr("Error: matrix: for COVAR, MWCOVAR, or IRED matrix only byatom output possible\n");
    return 1;
  }

  if ( (type_ == MATRIX_DISTCOVAR || type_ == MATRIX_IDEA) &&
       (mask2_.MaskStringSet() || outtype_ != BYATOM) )
  {
    mprinterr(
      "Error: matrix: DISTCOVAR or IDEA matrix only generated if no mask2 and byatom output\n");
    return 1;
  }

  // If the name is not empty, add Matrix to datasetlist
  // TODO: Check for name conflicts
  if (!name_.empty()) {
    DSL->AddDataSet( (DataSet*)this );
    // Since this now exists in the DataSetList and ActionList,
    // set the noDelete flag.
    SetNoDelete();
  }

  Info();

  return 0;
}

void MatrixType::Info() {
  mprintf("    MATRIX: Calculating %s ", MatrixModeString[type_]);
  switch (outtype_) {
    case BYATOM:    mprintf("by atom"); break;
    case BYRESIDUE: mprintf("by residue"); break;
    case BYMASK:    mprintf("by mask"); break;
  }
  if (!filename_.empty())
    mprintf(", dumping to file %s",filename_.c_str());
  if (useMass_)
    mprintf(", using mass weighting\n");
  else
    mprintf(", using no mass weighting\n");

  if (type_==MATRIX_IRED)
    mprintf("            Order of Legendre polynomials: %i\n",order_);
  if (!name_.empty())
    mprintf("            Storing matrix on internal stack with name: %s\n", name_.c_str());
  if (start_!=0 || stop_!=-1 || offset_!=1) {
    mprintf("            start: %i  stop:",start_);
    if (stop_==-1)
      mprintf(" Final frame");
    else
      mprintf(" %i",stop_);
    mprintf("  offset: %i\n",offset_);
  }
  if (type_!=MATRIX_IRED) {
    mprintf("            Mask1: %s\n",mask1_.MaskString());
    if (mask2_.MaskStringSet())
    mprintf("            Mask2: %s\n",mask2_.MaskString());
  }
}

int MatrixType::setup() {
  // For now, only allow 1 prmtop
  
  // Set up masks
  if (type_ != MATRIX_IRED) {
    if (currentParm->SetupIntegerMask(mask1_)) return 1;
    mprintf("\tMask1[%s]= %i atoms.\n",mask1_.MaskString(), mask1_.Nselected());
    if (mask1_.None()) {
      mprinterr("Error: No atoms selected for mask1.\n");
      return 1;
    }
    if (mask1tot_ > 0 && mask1tot_ != mask1_.Nselected()) {
      mprinterr("Error: # of atoms in mask 1 has changed. This is currently not\n");
      mprinterr("       supported.\n");
      return 1;
    }
    mask1tot_ = mask1_.Nselected();

    if (mask2_.MaskStringSet()) {
      if (currentParm->SetupIntegerMask(mask2_)) return 1;
      mprintf("\tMask2[%s]= %i atoms.\n",mask2_.MaskString(), mask2_.Nselected());
      if (mask2_.None()) {
        mprinterr("Error: No atoms selected for mask2.\n");
        return 1;
      }
      if (mask2tot_ > 0 && mask2tot_ != mask2_.Nselected()) {
        mprinterr("Error: # of atoms in mask 2 has changed. This is currently not\n");
        mprinterr("       supported.\n");
        return 1;
      }
    }
    mask2tot_ = mask2_.Nselected();

    if (mask1tot_ < mask2tot_) {
      mprinterr("Error: matrix: # of atoms in mask1 < # of atoms in mask2\n");
      return 1;
    }
  }

  // Allocate vector memory
  if (vect_==0 && vect2_==0) {
    if (type_ != MATRIX_DIST) {
      // No vector necessary for distance matrix
      if (type_ == MATRIX_DISTCOVAR)
        vectsize_ = mask1tot_ * (mask1tot_ - 1) / 2;
      else
        vectsize_ = mask1tot_ + mask2tot_;
      vect_ = new double[ vectsize_ * 3];
      vect2_ = new double[ vectsize_ * 3];
      for (int i = 0; i < vectsize_*3; ++i) {
        vect_[i] = 0;
        vect2_[i] = 0;
      }
    }
  }

  // Allocate matrix memory
  if (mat_==0) {
    if (mask2tot_ == 0) {
      // "upper right half" matrix, including main diagonal
      if (type_ == MATRIX_DISTCOVAR)
        matsize_ = mask1tot_ * (mask1tot_ - 1) * (mask1tot_ * (mask1tot_ - 1) / 2 + 1) / 4;
      else if (type_ == MATRIX_COVAR || 
               type_ == MATRIX_MWCOVAR)
        matsize_ = 9 * mask1tot_ * (mask1tot_ + 1) / 2;
      else // MATRIX_DIST || MATRIX_CORREL || MATRIX_IDEA || MATRIX_IRED
        matsize_ = mask1tot_ * (mask1tot_ + 1) / 2;
    } else {
      // full matrix -> no MATRIX_DISTCOVAR, MATRIX_IDEA, or MATRIX_IRED possible 
      if (type_ == MATRIX_COVAR || 
          type_ == MATRIX_MWCOVAR)
        matsize_ = 9 * mask1tot_ * mask2tot_;
      else // MATRIX_DIST || MATRIX_CORREL 
        matsize_ = mask1tot_ * mask2tot_;
    }
    mat_ = new double[ matsize_ ];
    for (int i = 0; i < matsize_; ++i)
      mat_[i] = 0;
  }
  
  return 0;
}
