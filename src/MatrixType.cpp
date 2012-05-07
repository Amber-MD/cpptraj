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
  std::string filename_ = actionArgs.GetStringKey("out");

  start_ = actionArgs.getKeyInt("start", 1);
  stop_ = actionArgs.getKeyInt("stop", -1);
  if (stop_ == -1)
    stop_ = actionArgs.getKeyInt("end", -1);
  offset_ = actionArgs.getKeyInt("offset", 1);

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
    type_ = MATRIX_MWCOVAR;
  else if (actionArgs.hasKey("covar"))
    type_ = MATRIX_MWCOVAR;
  else if (actionArgs.hasKey("correl"))
    type_ = MATRIX_MWCOVAR;
  else if (actionArgs.hasKey("idea"))
    type_ = MATRIX_MWCOVAR;
  else if (actionArgs.hasKey("ired"))
    type_ = MATRIX_MWCOVAR;
  else
    type_ = MATRIX_DIST;

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
  } 

  // Get matrix name
  // NOTE: Unlike ptraj where this was done after the 'mass' keyword check
  //       this is done at the very end here.

  return 0;
}
