#include <cmath> // sqrt
#include "MatrixType.h"
#include "VectorType.h"
#include "CpptrajStdio.h"
#include "vectormath.h" // vector_sub, dot_product

// CONSTRUCTOR
MatrixType::MatrixType() : 
  type_(MATRIX_NULL),
  vect_(0),
  vect2_(0),
  mat_(0),
  vectsize_(0),
  matsize_(0),
  mask2expr_(NULL),
  mask1tot_(0),
  mask2tot_(0),
  Nelt_(0),
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

// MatrixType::init()
/** matrix dist|covar|mwcovar|distcovar|correl|idea|ired
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
int MatrixType::init() {
  filename_ = actionArgs.GetStringKey("out");

  start_ = actionArgs.getKeyInt("start", 1);
  stop_ = actionArgs.getKeyInt("stop", -1);
  if (stop_ == -1)
    stop_ = actionArgs.getKeyInt("end", -1);
  offset_ = actionArgs.getKeyInt("offset", 1);
  // Actual start frame arg should be from 0
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
    mask2expr_ = actionArgs.getNextMask();
  } 

  // Get matrix name
  // NOTE: Unlike ptraj where this was done after the 'mass' keyword check,
  //       here it is done after all other keywords are processed. 
  name_ = actionArgs.GetStringNext();

  // Check arguments
  if ( !name_.empty() && mask2expr_!=NULL ) {
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
       (mask2expr_!=NULL || outtype_ != BYATOM) )
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

// MatrixType::Info()
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
    if (mask2expr_!=NULL)
      mprintf("            Mask2: %s\n",mask2expr_);
  }
}

// MatrixType::setup()
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

    if (mask2expr_!=NULL) {
      mask2_.SetMaskString(mask2expr_);
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
    } else {
      mask2_ = mask1_;
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
    if (mask2expr_ == NULL) {
      // "upper right half" matrix, including main diagonal.
      // Nelt is set for use in indexing half matrix with HalfMatrixIndex.
      if (type_ == MATRIX_DISTCOVAR) {
        matsize_ = mask1tot_ * (mask1tot_ - 1) * (mask1tot_ * (mask1tot_ - 1) / 2 + 1) / 4;
        Nelt_ = mask1tot_ * (mask1tot_ - 1) / 2;
      } else if (type_ == MATRIX_COVAR || 
                 type_ == MATRIX_MWCOVAR) {
        matsize_ = 9 * mask1tot_ * (mask1tot_ + 1) / 2;
        Nelt_ = mask1tot_ * 3;
      } else { // MATRIX_DIST || MATRIX_CORREL || MATRIX_IDEA || MATRIX_IRED
        matsize_ = mask1tot_ * (mask1tot_ + 1) / 2;
        Nelt_ = mask1tot_; 
      }
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

// MatrixType::action()
int MatrixType::action() {
  // If the current frame is less than start exit
  if (frameNum < start_) return 0;
  // If the current frame is greater than stop exit
  if (stop_!=-1 && frameNum >= stop_) return 0;
  // Increment number of snapshots, update next target frame
  ++snap_;
  start_ += offset_;

  // ---------- Calc Distance Matrix -------------
  if (type_ == MATRIX_DIST) {
    int idx = 0;
    if (mask2expr_==NULL) {
      // Upper triangle
      for (AtomMask::const_iterator atom2 = mask1_.begin();
                                    atom2 != mask1_.end(); ++atom2)
      { 
        for (AtomMask::const_iterator atom1 = atom2; 
                                      atom1 != mask1_.end(); ++atom1)
        {
          mat_[idx] += currentFrame->DIST(*atom2, *atom1);
          ++idx;
        }
      }
    } else {
      // Full matrix
      for (AtomMask::const_iterator atom2 = mask2_.begin();
                                    atom2 != mask2_.end(); ++atom2)
      {
        for (AtomMask::const_iterator atom1 = mask1_.begin();
                                      atom1 != mask1_.end(); ++atom1)
        {
          mat_[idx] += currentFrame->DIST(*atom2, *atom1);
          ++idx;
        }
      }
    }

  // ---------- Calc covariance or correlation ---
  } else if ( type_ == MATRIX_COVAR ||
              type_ == MATRIX_MWCOVAR ||
              type_ == MATRIX_CORREL )
  {
    int midx = 0;
    int vidx1 = 0;  // Index into vector arrays
    int vidx2 = 0;
    double XYZi[3]; // For retrieving atom coordinates
    double XYZj[3];
    int lend;
    bool iscovariance;
    if ( type_ == MATRIX_CORREL ) {
      lend = 1;
      iscovariance = false;
    } else {
      lend = 3;
      iscovariance = true;
    }
    int kidx = 1;
    // BEGIN OUTER LOOP 
    for (AtomMask::const_iterator atom2 = mask2_.begin();
                                  atom2 != mask2_.end(); ++atom2)
    {
      currentFrame->GetAtomXYZ(XYZi, *atom2);
      if (mask2expr_!=NULL) {
        for (int ixyz = 0; ixyz < 3; ++ixyz) {
          vect_[vidx1] += XYZi[ixyz];
          vect2_[vidx1++] += (XYZi[ixyz]*XYZi[ixyz]);
        }
      }
      for (int l = 0; l < lend; ++l) {
        // BEGIN INNER LOOP
        for (AtomMask::const_iterator atom1 = mask1_.begin();
                                      atom1 != mask1_.end(); ++atom1)
        {
          currentFrame->GetAtomXYZ(XYZj, *atom1);
          if (kidx==1) {
            for (int ixyz = 0; ixyz < 3; ++ixyz) {
              vect_[vidx2] += XYZj[ixyz];
              vect2_[vidx2++] += (XYZj[ixyz]*XYZj[ixyz]);
            }
          }
          if ( (mask2expr_==NULL && *atom1 >= *atom2) || mask2expr_!=NULL ) {
            if ( iscovariance ) {
              if (mask2expr_==NULL && *atom2 == *atom1) {
                // Half matrix or diagonal
                if (l==0) {
                  mat_[midx++] += XYZi[0]*XYZj[0];
                  mat_[midx++] += XYZi[0]*XYZj[1];
                  mat_[midx++] += XYZi[0]*XYZj[2];
                } else if (l==1) {
                  mat_[midx++] += XYZi[1]*XYZj[1];
                  mat_[midx++] += XYZi[1]*XYZj[2];
                } else  // l==2
                  mat_[midx++] += XYZi[2]*XYZj[2];
              } else if ((mask2expr_==NULL && *atom2 < *atom1) || mask2expr_!=NULL) {
                // Half matrix and atom2 < atom1, or Full matrix
                if (l==0) {
                  mat_[midx++] += XYZi[0]*XYZj[0];
                  mat_[midx++] += XYZi[0]*XYZj[1];
                  mat_[midx++] += XYZi[0]*XYZj[2];
                } else if (l==1) {
                  mat_[midx++] += XYZi[1]*XYZj[0];
                  mat_[midx++] += XYZi[1]*XYZj[1];
                  mat_[midx++] += XYZi[1]*XYZj[2];
                } else { // l==2
                  mat_[midx++] += XYZi[2]*XYZj[0];
                  mat_[midx++] += XYZi[2]*XYZj[1];
                  mat_[midx++] += XYZi[2]*XYZj[2];
                }
              }
            } else {
              mat_[midx++] += XYZi[0]*XYZj[0] + XYZi[1]*XYZj[1] + XYZi[2]*XYZj[2];
            }
          }
        } // END INNER LOOP
        kidx = 0;
      } // END LOOP OVER l
    } // END OUTER LOOP

  // ---------- Calc Isotropically distributed ---
  //            ensemble matrix.
  // See Proteins 2002, 46, 177; eq. 7
  } else if (type_ == MATRIX_IDEA) {
    // Get COM
    double COM[3];
    currentFrame->CenterOfMass( &mask1_, COM );
    // Get ri, rj, and calc ri*rj
    int midx = 0;
    int vidx = 0;
    double ri[3];
    double rj[3];
    // Matrix IDEA only uses 1 mask.
    for (AtomMask::const_iterator atomi = mask1_.begin();
                                  atomi != mask1_.end(); ++atomi)
    {
      currentFrame->GetAtomXYZ(ri, *atomi);
      vector_sub(ri, ri, COM);
      for (AtomMask::const_iterator atomj = atomi; 
                                    atomj != mask1_.end(); ++atomj)
      {
        currentFrame->GetAtomXYZ(rj, *atomj);
        vector_sub(rj, rj, COM);
        double val = dot_product(ri, rj);
        mat_[midx++] += val; 
        if (*atomj == *atomi) {
          vect_[vidx] += val;
          vect2_[vidx] += (val * val);
          ++vidx;
        }
      }
    }
  } 

  return 0;
}

// MatrixType::HalfMatrixIndex()
int MatrixType::HalfMatrixIndex(int row, int col) {
  int i, j;
  if (row<=col) {
    i = row;
    j = col;
  } else {
    i = col;
    j = row;
  }
  //mprintf("CDBG:\ti=%i j=%i N=%i idx=%i\n",i,j,Nelt_,
  //        (i * Nelt_ - (i * (i-1) / 2) + (j - i)));
  return (i * Nelt_ - (i * (i-1) / 2) + (j - i));
}

// MatrixType::print()
void MatrixType::print() {
  int err = 0;
  CpptrajFile outfile;
  // ---------- Calculate average over number of sets
  double dsnap = (double)snap_;
  if (vect_!=0) {
    for (int i = 0; i < vectsize_*3; ++i) {
      vect_[i] /= dsnap;
      vect2_[i] /= dsnap;
    }
  }
  for (int i = 0; i < matsize_; ++i)
    mat_[i] /= dsnap;

  // ---------- Isotropically distributed ensmble matrix
  if (type_ == MATRIX_IDEA) {
    for (int i = 0; i < vectsize_*3; ++i) {
      vect_[i] /= 3;
      vect2_[i] /= 3;
    }
    for (int i = 0; i < matsize_; ++i)
      mat_[i] /= 3.0;

  // ---------- Calc covariance or correlation matrix
  } else if (type_==MATRIX_COVAR || type_==MATRIX_MWCOVAR || type_==MATRIX_CORREL) {
    // Calc <riri> - <ri><ri>
    for (int i = 0; i < vectsize_*3; ++i)
      vect2_[i] -= (vect_[i]*vect_[i]);
    // Calc <rirj> - <ri><rj>
    int idx = 0;
    int vidx; 
    if (mask2expr_==NULL)
      vidx = 0;
    else
      vidx = mask1tot_;
    if (type_==MATRIX_COVAR || type_==MATRIX_MWCOVAR) {
      for (int i = 0; i < mask2tot_; ++i) {
        for (int l = 0; l < 3; ++l) {
          for (int j = 0; j < mask1tot_; ++j) {
            if ( (mask2expr_==NULL && j>=i) || mask2expr_!=NULL ) {
              if (mask2expr_==NULL && i==j) { 
                if (l==0) {
                  mat_[idx++] -= vect_[(vidx + i)*3  ] * vect_[j*3  ];
                  mat_[idx++] -= vect_[(vidx + i)*3  ] * vect_[j*3+1];
                  mat_[idx++] -= vect_[(vidx + i)*3  ] * vect_[j*3+2];
                } else if (l==1) {
                  mat_[idx++] -= vect_[(vidx + i)*3+1] * vect_[j*3+1];
                  mat_[idx++] -= vect_[(vidx + i)*3+1] * vect_[j*3+2];
                } else // l==2
                  mat_[idx++] -= vect_[(vidx + i)*3+2] * vect_[j*3+2];
              } else if ((mask2expr_==NULL && i<j) || mask2expr_!=NULL) {
                if (l==0) {
                  mat_[idx++] -= vect_[(vidx + i)*3  ] * vect_[j*3  ];
                  mat_[idx++] -= vect_[(vidx + i)*3  ] * vect_[j*3+1];
                  mat_[idx++] -= vect_[(vidx + i)*3  ] * vect_[j*3+2];
                } else if (l==1) {
                  mat_[idx++] -= vect_[(vidx + i)*3+1] * vect_[j*3  ];
                  mat_[idx++] -= vect_[(vidx + i)*3+1] * vect_[j*3+1];
                  mat_[idx++] -= vect_[(vidx + i)*3+1] * vect_[j*3+2];
                } else { // l==2
                  mat_[idx++] -= vect_[(vidx + i)*3+2] * vect_[j*3  ];
                  mat_[idx++] -= vect_[(vidx + i)*3+2] * vect_[j*3+1];
                  mat_[idx++] -= vect_[(vidx + i)*3+2] * vect_[j*3+2];
                } 
              }
            }
          }
        }
      }
      // Add in mass information for MWCOVAR
      if (type_ == MATRIX_MWCOVAR) {
        idx = 0;
        int row3 = 0;
        double mass;
        for (AtomMask::const_iterator atomi = mask2_.begin();
                                      atomi != mask2_.end(); ++atomi)
        {
          for (int k = 0; k < 3; ++k) {
            int col3 = 0;
            for (AtomMask::const_iterator atomj = mask1_.begin();
                                          atomj != mask1_.end(); ++atomj)
            {
              if (*atomi == *atomj) {
                vect2_[col3  ] *= mass;
                vect2_[col3+1] *= mass;
                vect2_[col3+2] *= mass;
              }
              if (mask2expr_==NULL && *atomj >= *atomi) { // half-matrix
                mass = sqrt( (*currentParm)[*atomi].Mass() * (*currentParm)[*atomj].Mass() );
                if (row3+k <= col3) {
                  idx = HalfMatrixIndex( row3+k, col3  );
                  mat_[idx] *= mass;
                }
                if (row3+k <= col3+1) {
                  idx = HalfMatrixIndex( row3+k, col3+1);
                  mat_[idx] *= mass;
                }
                if (row3+k <= col3+2) {
                  idx = HalfMatrixIndex( row3+k, col3+2);
                  mat_[idx] *= mass;
                }
              } else if (mask2expr_!=NULL) { // full matrix
                mass = sqrt( (*currentParm)[*atomi].Mass() * (*currentParm)[*atomj].Mass() );
                mat_[idx++] *= mass;
                mat_[idx++] *= mass;
                mat_[idx++] *= mass;
              }
              col3 += 3;
            } // END LOOP OVER MASK1
          }
          row3 += 3;
        } // END LOOP OVER MASK2
      }
    } else { // MATRIX_CORREL
      for (int i = 0; i < mask2tot_; ++i) {
        for (int j = 0; j < mask1tot_; ++j) {
          if ( (mask2expr_==NULL && j>=i) || mask2expr_!=NULL ) {
            mat_[idx] -= (vect_[j*3  ] * vect_[(vidx + i)*3  ] +
                          vect_[j*3+1] * vect_[(vidx + i)*3+1] +
                          vect_[j*3+2] * vect_[(vidx + i)*3+2]);
            // Normalize
            mat_[idx] /= sqrt((vect2_[j*3] + vect2_[j*3+1] + vect2_[j*3+2]) *
                              (vect2_[(vidx + i)*3  ] + 
                               vect2_[(vidx + i)*3+1] +
                               vect2_[(vidx + i)*3+2]));
            ++idx;
          }
        }
      }
    } 
  }

  // Open file
  if (filename_.empty())
    err = outfile.SetupWrite(NULL, debug);
  else
    err = outfile.SetupWrite(filename_.c_str(), debug);
  if (err!=0) return;
  if (outfile.OpenFile()) return;

  // ---------- Print out BYATOM 
  if (outtype_==BYATOM) {
    int idx = 0;
    if (type_ == MATRIX_DIST || type_== MATRIX_CORREL || type_ == MATRIX_IDEA) {
      // ----- DISTANCE, CORREL, IDEA ------
      if (mask2expr_==NULL) { // Half-matrix
        for (int row = 0; row < mask1tot_; ++row) {
          for (int col = 0; col < mask1tot_; ++col) {
            idx = HalfMatrixIndex( row, col );
            //mprintf("DBG:\tRow=%i Col=%i Idx=%i\n", row, col, idx);
            outfile.Printf("%6.3f ", mat_[idx]);
          }
          outfile.Printf("\n");
        }
      } else { // Full matrix
        for (int row = 0; row < mask2tot_; ++row) {
          for (int col = 0; col < mask1tot_; ++col) {
            outfile.Printf("%6.3f ", mat_[idx++]);
          }
          outfile.Printf("\n");
        }
      }
    } else if (type_ == MATRIX_COVAR || type_ == MATRIX_MWCOVAR) {
      // ----- COVAR ---------
      if (mask2expr_==NULL) { // Half-matrix
        // NOTE: Use Nelt instead of mask1tot*3?
        for (int row = 0; row < mask1tot_*3; row+=3) {
          for (int l = 0; l < 3; ++l) {
            for (int col = 0; col < mask1tot_*3; col+=3) {
              for (int ixyz = 0; ixyz < 3; ++ixyz) {
                idx = HalfMatrixIndex( row+l, col+ixyz);
                outfile.Printf("%6.3f ", mat_[idx]);
              }
            }
            outfile.Printf("\n");
          }
        }
      } else { // Full matrix
        // NOTE: Use Nelt instead of mask1tot*3?
        for (int row = 0; row < mask2tot_*3; row+=3) {
          for (int l = 0; l < 3; ++l) {
            for (int col = 0; col < mask1tot_*3; col+=3) {
              for (int ixyz = 0; ixyz < 3; ++ixyz) {
                outfile.Printf("%6.3f ", mat_[idx++]);
              }
            }
            outfile.Printf("\n");
          }
        }
      } 
    }
  } // END BYATOM

  // Close file
  outfile.Printf("\n");
  outfile.CloseFile();
}

