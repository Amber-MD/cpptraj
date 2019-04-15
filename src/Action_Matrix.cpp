#include <cmath> // sqrt
#include "Action_Matrix.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "Constants.h" // DEGRAD

// CONSTRUCTOR
Action_Matrix::Action_Matrix() :
  Mat_(0),
  matByRes_(0),
  outfile_(0),
  byMaskOut_(0),
  outtype_(BYATOM),
  debug_(0),
  order_(2),
  useMask2_(false),
  useMass_(false)
{}

void Action_Matrix::Help() const {
  mprintf("\t[out <filename>] %s\n", ActionFrameCounter::HelpText);
  mprintf("\t[name <name>] [ byatom | byres [mass] | bymask [mass] ]\n"
          "\t[ ired [order <#>] ]\n"
          "\t[ {distcovar | idea} <mask1> ]\n"
          "\t[ {dist | correl | covar | mwcovar} <mask1> [<mask2>] ]\n"
          "\t[ dihcovar dihedrals <dataset arg> ]\n"
          "  Calculate a matrix of the specified type from input coordinates.\n"
          "    dist: Distance matrix (default).\n"
          "    correl: Correlation matrix (aka dynamic cross correlation).\n"
          "    covar: Coordinate covariance matrix.\n"
          "    mwcovar: Mass-weighted coordinate covariance matrix.\n"
          "    distcovar: Distance covariance matrix.\n"
          "    idea: Isotropically Distributed Ensemble Analysis matrix.\n"
          "    ired: Isotropic Reorientational Eigenmode Dynamics matrix.\n"
          "    dihcovar: Dihedral covariance matrix.\n");
}

// Action_Matrix::Init()
Action::RetType Action_Matrix::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  trajComm_ = init.TrajComm();
# endif
  debug_ = debugIn;
  // Get Keywords
  std::string outfilename = actionArgs.GetStringKey("out");
  // Get start/stop/offset
  if (InitFrameCounter(actionArgs)) return Action::ERR;
  // Determine matrix type
  MetaData::scalarType mtype = MetaData::DIST;
  if (actionArgs.hasKey("distcovar"))
    mtype = MetaData::DISTCOVAR;
  else if (actionArgs.hasKey("mwcovar"))
    mtype = MetaData::MWCOVAR;
  else if (actionArgs.hasKey("dist"))
    mtype = MetaData::DIST;
  else if (actionArgs.hasKey("covar"))
    mtype = MetaData::COVAR;
  else if (actionArgs.hasKey("correl"))
    mtype = MetaData::CORREL;
  else if (actionArgs.hasKey("idea"))
    mtype = MetaData::IDEA;
  else if (actionArgs.hasKey("ired"))
    mtype = MetaData::IREDMAT;
  else if (actionArgs.hasKey("dihcovar"))
    mtype = MetaData::DIHCOVAR;
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
  if ( outtype_ != BYATOM && (mtype == MetaData::COVAR || 
                              mtype == MetaData::MWCOVAR || 
                              mtype == MetaData::IREDMAT ) )
  {
    mprinterr("Error: matrix: for COVAR, MWCOVAR, or IRED matrix only byatom output possible\n");
    return Action::ERR;
  }
  // Get matrix name
  std::string name = actionArgs.GetStringKey("name");
  // UseMass
  useMass_ = actionArgs.hasKey("mass");
  // NOTE: Determine matrix kind here so subsequent Actions/Analyses know about it.
  DataSet_2D::MatrixKindType mkind = DataSet_2D::HALF;
  if (mtype == MetaData::IREDMAT) { // IRED matrix
    // Setup IRED vectors and determine Legendre order
    order_ = actionArgs.getKeyInt("order",1);
    if (order_ <= 0) {
      mprinterr("Error: matrix: order parameter <= 0, ignoring command\n");
      return Action::ERR;
    }
    for ( DataSetList::const_iterator DS = init.DSL().begin(); DS != init.DSL().end(); ++DS) {
      if ( (*DS)->Type() == DataSet::VECTOR && (*DS)->Meta().ScalarType() == MetaData::IREDVEC )
        IredVectors_.push_back( (DataSet_Vector*)*DS );
    }
    if (IredVectors_.empty()) {
      mprinterr("Error: matrix: no vectors defined for IRED\n");
      return Action::ERR;
    }
  } else if (mtype == MetaData::DIHCOVAR) { // Dihedral Covariance
    // Get data set mask for dihedral covariance
    DihedralSets_.clear();
    DihedralSets_.AddTorsionSets( init.DSL().GetMultipleSets( actionArgs.GetStringKey("dihedrals") ) );
    if ( DihedralSets_.empty() ) {
      mprinterr("Error: No valid data sets found.\n");
      return Action::ERR;
    }
  } else {
    // Get masks if not IRED/DIHCOVAR
    if (mask1_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;
    std::string maskexpr = actionArgs.GetMaskNext();
    if (!maskexpr.empty()) useMask2_ = true;
    if ( useMask2_ && (mtype == MetaData::IDEA || mtype == MetaData::DISTCOVAR) )
    {
      mprinterr("Error: Mask 2 [%s] specified but not used for %s matrix\n",
                maskexpr.c_str(), MetaData::TypeString(mtype));
      useMask2_ = false;
      return Action::ERR;
    }
    if (useMask2_) {
      if (mask2_.SetMaskString( maskexpr )) return Action::ERR;
      mkind = DataSet_2D::FULL;
    }
  }

  // Create matrix BYATOM DataSet
  Mat_ = (DataSet_MatrixDbl*)init.DSL().AddSet(DataSet::MATRIX_DBL,
                                         MetaData(name, MetaData::M_MATRIX, mtype), "Mat");
  if (Mat_ == 0) return Action::ERR;
  // NOTE: Type/Kind is set here so subsequent analyses/actions know about it.
  Mat_->SetMatrixKind( mkind );
  // Set default precision for backwards compat.
  Mat_->SetupFormat().SetFormatWidthPrecision(6, 3);
  Mat_->ModifyDim(Dimension::X).SetLabel("Atom");
  // Determine what will be output.
  matByRes_ = 0;
  outfile_ = 0;
  byMaskOut_ = 0;
  if (outtype_ == BYMASK) {
    // BYMASK output - no final data set, just write to file/STDOUT.
    byMaskOut_ = init.DFL().AddCpptrajFile(outfilename, "Matrix by mask",
                                     DataFileList::TEXT, true);
    if (byMaskOut_ == 0) return Action::ERR;
  } else {
    // BYATOM or BYRESIDUE output. Set up DataFile. 
    if (outtype_ == BYRESIDUE) {
      MetaData md(Mat_->Meta().Name(), "ByRes");
      md.SetScalarMode(MetaData::M_MATRIX);
      matByRes_ = (DataSet_MatrixDbl*)init.DSL().AddSet(DataSet::MATRIX_DBL, md);
      if (matByRes_ == 0) return Action::ERR;
      matByRes_->SetupFormat().SetFormatWidthPrecision(6, 3);
      matByRes_->ModifyDim(Dimension::X).SetLabel("Res");
#     ifdef MPI
      // Since by-residue matrix allocated in Print(), no sync needed.
      matByRes_->SetNeedsSync( false );
#     endif
    } 
    outfile_ = init.DFL().AddDataFile(outfilename, "square2d noxcol noheader", actionArgs);
    if (outfile_ != 0) {
      if (outtype_ == BYATOM)
        outfile_->AddDataSet( Mat_ );
      else
        outfile_->AddDataSet( matByRes_ );
    }
  }

  mprintf("    MATRIX: Calculating %s matrix, output is", Mat_->Meta().TypeString());
  switch (outtype_) {
    case BYATOM:    mprintf(" by atom.\n"); break;
    case BYRESIDUE: mprintf(" averaged by residue.\n"); break;
    case BYMASK:    mprintf(" averaged by mask.\n"); break;
  }
  if (outtype_ != BYATOM) {
    if (useMass_)
      mprintf("\tAverages will be mass-weighted.\n");
    else
      mprintf("\tAverages will not be mass-weighted.\n");
  }
  if (mtype == MetaData::IREDMAT)
    mprintf("\t%u IRED vecs, Order of Legendre polynomials: %i\n",
            IredVectors_.size(), order_);
  else if (mtype == MetaData::DIHCOVAR)
    mprintf("\t%u data sets.\n", DihedralSets_.size());
  if (outfile_ != 0)
    mprintf("\tPrinting to file %s\n", outfile_->DataFilename().full());
  if (byMaskOut_ != 0)
    mprintf("\tAveraged by mask output to %s\n", byMaskOut_->Filename().full());
  mprintf("\tMatrix data set is '%s'\n", Mat_->legend());
  if (matByRes_ != 0)
    mprintf("\tAveraged by residue matrix data set is '%s'\n", matByRes_->legend());
  FrameCounterInfo();
  if (mtype != MetaData::IREDMAT && mtype != MetaData::DIHCOVAR) {
    mprintf("\tMask1 is '%s'\n",mask1_.MaskString());
    if (useMask2_)
      mprintf("\tMask2 is '%s'\n",mask2_.MaskString());
  }
#ifdef NEW_MATRIX_PARA
  mprintf("DEBUG: NEW COVARIANCE MATRIX PARALLELIZATION SCHEME IN USE.\n");
#endif
  return Action::OK;
}

// Action_Matrix::FillMassArray()
Action_Matrix::Darray Action_Matrix::FillMassArray(Topology const& currentParm, 
                                                   AtomMask const& mask) const
{
  Darray mass;
  mass.reserve( mask.Nselected() );
  for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom) 
    mass.push_back( currentParm[ *atom ].Mass() );
  return mass;
}

static Action::RetType PrintMask2Error() {
  mprinterr("Error: Second mask (full matrix) not supported for DISTCOVAR,\n"
            "Error:   IDEA, or IRED matrix.\n");
  return Action::ERR;
}

// Action_Matrix::MaskToMatResArray()
Action_Matrix::MatResArray Action_Matrix::MaskToMatResArray(Topology const& currentParm,
                                                            AtomMask const& mask) const
{
  MatResArray residues;
  int currentResNum = -1;
  matrix_res blank_res;
  for (int idx = 0; idx != mask.Nselected(); idx++)
  {
    int atom1 = mask[idx];
    int resNum = currentParm[atom1].ResNum();
    if (resNum != currentResNum) {
      residues.push_back( blank_res );
      residues.back().resnum_ = resNum;
      currentResNum = resNum;
    }
    residues.back().maskIdxs_.push_back( idx );
  }
  if (debug_ > 0) {
    mprintf("DEBUG: BYRES: MASK '%s'\n", mask.MaskString());
    for (MatResArray::const_iterator res = residues.begin(); res != residues.end(); ++res) {
      mprintf("\tRes %i:", res->resnum_+1);
      for (Iarray::const_iterator it = res->maskIdxs_.begin(); it != res->maskIdxs_.end(); ++it)
        mprintf(" %i (%i)", mask[*it] + 1, *it);
      mprintf("\n");
    }
  }
  return residues;
}

// Action_Matrix::Setup()
Action::RetType Action_Matrix::Setup(ActionSetup& setup) {
  size_t mask1tot = 0; // Will be # of columns
  size_t mask2tot = 0; // Will be # of rows if not symmetric matrix

  // Set up masks.
  if (Mat_->Meta().ScalarType() == MetaData::IREDMAT) {
    // IRED - matrix # cols = # of IRED vectors
    mask1tot = IredVectors_.size();
  } else if (Mat_->Meta().ScalarType() == MetaData::DIHCOVAR) {
    // Dihedral covariance - matrix # cols = # data sets
    mask1tot = DihedralSets_.size();
  } else {
    if (setup.Top().SetupIntegerMask(mask1_)) return Action::ERR;
    mask1_.MaskInfo();
    if (mask1_.None()) {
      mprintf("Warning: No atoms selected for mask1.\n");
      return Action::SKIP;
    }
    if (useMask2_) {
      if (setup.Top().SetupIntegerMask(mask2_)) return Action::ERR;
      mask2_.MaskInfo(); 
      if (mask2_.None()) {
        mprintf("Warning: No atoms selected for mask2.\n");
        return Action::SKIP;
      }
    }
    mask1tot = (size_t)mask1_.Nselected();
    mask2tot = (size_t)mask2_.Nselected();
  }
  if (mask1tot < mask2tot) {
    mprinterr("Error: # of atoms in mask1 < # of atoms in mask2\n");
    return Action::ERR;
  }
  // Store mass info for MWCOVAR matrix or when output BYRESIDUE or BYMASK.
  if ( Mat_->Meta().ScalarType() == MetaData::MWCOVAR || (outtype_ != BYATOM && useMass_) ) {
    mass1_ = FillMassArray(setup.Top(), mask1_);
    mass2_ = FillMassArray(setup.Top(), mask2_);
    if (outtype_ != BYATOM && !useMass_)
      mprintf("Warning: Using mass-weighted covar matrix, byres/bymask output will be"
              "weighted by mass.\n");
  } else if (outtype_ != BYATOM && !useMass_) {
    mass1_ = Darray(mask1_.Nselected(), 1.0);
    if (useMask2_)
      mass2_ = Darray(mask2_.Nselected(), 1.0);
  }
  // If output is BYRESIDUE, determine which mask elements correspond to 
  // which residues, as well as how many residues are selected in mask1 (cols)
  // and mask2 (rows) to determine output matrix dimensions.
  if (outtype_ == BYRESIDUE) {
    residues1_ = MaskToMatResArray( setup.Top(), mask1_ );
    residues2_ = MaskToMatResArray( setup.Top(), mask2_ );
  }
  // Determine matrix/vector dimensions. 
  size_t vectsize = 0;
  size_t ncols = 0;
  size_t nrows = 0;
  switch( Mat_->Meta().ScalarType() ) {
    case MetaData::CORREL   : // Like DIST but vectors required. 
      vectsize = (mask1tot + mask2tot) * 3;
    case MetaData::DIST     : // No vectors required.
      ncols = mask1tot;
      nrows = mask2tot;
      break;
    case MetaData::DISTCOVAR: // No Full matrix possible.
      vectsize = mask1tot * (mask1tot - 1) / 2;
      ncols = vectsize;
      if (mask2tot > 0) return PrintMask2Error();
      break;
    case MetaData::MWCOVAR  :
      // Store mass info matrix DataSet for mass-weighting eigenvectors.
      Mat_->StoreMass( mass1_ );
    case MetaData::COVAR    :
      vectsize = (mask1tot + mask2tot) * 3;
      ncols = mask1tot * 3;
      nrows = mask2tot * 3;
      break;
    case MetaData::DIHCOVAR: // Dihedral covariance
      vectsize = (mask1tot + mask2tot) * 2;
      ncols = mask1tot * 2;
      nrows = mask2tot * 2;
      break;
    case MetaData::IDEA     :
    case MetaData::IREDMAT  : // No Full matrix possible.
      vectsize = mask1tot + mask2tot;
      ncols = mask1tot;
      if (mask2tot > 0) return PrintMask2Error();
      break;
    default: return Action::ERR; // Sanity check
  }
  // Allocate vector memory.
  Mat_->AllocateVector( vectsize );
  vect2_.resize( vectsize, 0.0 );
  // Allocate matrix memory. If already allocated, do not allow sizes to change.
  if (Mat_->Size() == 0) {
    if (nrows > 0) // Full matrix - no DISTCOVAR, IDEA, or IRED possible
      Mat_->Allocate2D( ncols, nrows );
    else           // "Upper right half" matrix, including main diagonal.
      Mat_->AllocateHalf( ncols );
  } else {
    bool dimensionsHaveChanged = false;
    if (nrows > 0)
      dimensionsHaveChanged = (nrows != Mat_->Nrows() || ncols != Mat_->Ncols());
    else
      dimensionsHaveChanged = (ncols != Mat_->Ncols());
    if (dimensionsHaveChanged) {
      mprinterr("Error: Attempting to reallocate matrix with different size.\n"
                "Error:   Original # cols = %zu, new # cols = %zu.\n",
                Mat_->Ncols(), ncols);
      if (nrows > 0) mprinterr("Error:  Original # rows = %zu, new # rows = %zu.\n",
                               Mat_->Nrows(), nrows);
      mprinterr("Error:   This can occur when different #s of atoms are selected in\n"
                "Error:   different topology files.\n");
      return Action::ERR;
    }
  }
# ifdef _OPENMP
  if (
       ( Mat_->Meta().ScalarType() == MetaData::COVAR ||
         Mat_->Meta().ScalarType() == MetaData::MWCOVAR ))
  {
#   ifdef NEW_MATRIX_PARA
    // Store coordinate XYZ indices of mask 1.
    crd_indices_.clear();
    for (AtomMask::const_iterator m1 = mask1_.begin(); m1 != mask1_.end(); ++m1)
    {
      int crdidx = *m1 * 3;
      crd_indices_.push_back( crdidx );
      crd_indices_.push_back( crdidx+1 );
      crd_indices_.push_back( crdidx+2 );
    }
#   else
    if (Mat_->MatrixKind() == DataSet_2D::FULL) {
      // Store combined mask1 and mask2 for diagonal.
      crd_indices_.clear();
      crd_indices_.reserve( mask1_.Nselected() + mask2_.Nselected() );
      for (AtomMask::const_iterator at = mask1_.begin(); at != mask1_.end(); ++at)
        crd_indices_.push_back( *at * 3 );
      for (AtomMask::const_iterator at = mask2_.begin(); at != mask2_.end(); ++at)
        crd_indices_.push_back( *at * 3 );
    }
    if (debug_ > 1) {
      mprintf("DEBUG: Combined mask1+mask2 coordinate indices:\n");
      for (unsigned int i = 0; i < crd_indices_.size(); i += 2)
        mprintf("%u:\t%i %i\n", i/2, crd_indices_[i], crd_indices_[i+1]);
    }
#   endif
  }
# endif
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
void Action_Matrix::CalcIredMatrix(int frameNum) {
  v_iterator v2idx1 = vect2_.begin();
  // Store length of IRED vectors in vect2
  for (std::vector<DataSet_Vector*>::const_iterator Vtmp = IredVectors_.begin();
                                                    Vtmp != IredVectors_.end(); ++Vtmp)
    *(v2idx1++) = sqrt( (*Vtmp)->VXYZ(frameNum) * (*Vtmp)->VXYZ(frameNum) );

  // Loop over all pairs of IRED vectors.
  DataSet_MatrixDbl::iterator mat = Mat_->begin();
  v_iterator v1idx = Mat_->v1begin();
  v2idx1 = vect2_.begin();
  for (std::vector<DataSet_Vector*>::iterator Vtmp = IredVectors_.begin();
                                              Vtmp != IredVectors_.end(); ++Vtmp)
  {
    double len1 = *v2idx1;
    v_iterator v2idx2 = v2idx1;
    for (std::vector<DataSet_Vector*>::iterator Vtmp2 = Vtmp;
                                                Vtmp2 != IredVectors_.end(); ++Vtmp2)
    {
      double len2 = *(v2idx2++);
      double legendre = LegendrePoly(order_, (*Vtmp)->VXYZ(frameNum) * (*Vtmp2)->VXYZ(frameNum)  / (len1 * len2) );
      *(mat++) += legendre;
      if (Vtmp == Vtmp2)
        *(v1idx++) += legendre;
    }
    ++v2idx1;
  }
}

/** Calc Distance Matrix */
void Action_Matrix::CalcDistanceMatrix(Frame const& currentFrame) {
  DataSet_MatrixDbl::iterator mat = ((DataSet_MatrixDbl*)Mat_)->begin();
  if (!useMask2_) {
    // Upper Triangle
    for (AtomMask::const_iterator atom2 = mask1_.begin(); atom2 != mask1_.end(); ++atom2)
      for (AtomMask::const_iterator atom1 = atom2; atom1 != mask1_.end(); ++atom1)
        *(mat++) += sqrt(DIST2_NoImage(currentFrame.XYZ(*atom2), currentFrame.XYZ(*atom1)));
  } else {
    // Full matrix
    for (AtomMask::const_iterator atom2 = mask2_.begin(); atom2 != mask2_.end(); ++atom2)
      for (AtomMask::const_iterator atom1 = mask1_.begin(); atom1 != mask1_.end(); ++atom1)
        *(mat++) += sqrt(DIST2_NoImage(currentFrame.XYZ(*atom2), currentFrame.XYZ(*atom1)));
  }
}

// Action_Matrix::StoreVec()
void Action_Matrix::StoreVec(v_iterator& v1, v_iterator& v2, const double* XYZ) const 
{
  *(v1++) += XYZ[0];
  *(v2++) += (XYZ[0] * XYZ[0]);
  *(v1++) += XYZ[1];
  *(v2++) += (XYZ[1] * XYZ[1]);
  *(v1++) += XYZ[2];
  *(v2++) += (XYZ[2] * XYZ[2]);
}

/** Calc Covariance Matrix */
void Action_Matrix::CalcCovarianceMatrix(Frame const& currentFrame) {
# ifdef _OPENMP
#ifdef NEW_MATRIX_PARA /* New matrix parallelization */
  //int idx2, atomCrd2, offset2, midx, idx1, atomCrd1, offset1;
  int idx2, midx, idx1;
  double Mj;
  if (useMask2_) { // FULL MATRIX TODO
    return;
  } else { // HALF MATRIX
    DataSet_MatrixDbl& matrix = *Mat_;
    Darray& vect1 = Mat_->V1();
    //int Ncoords = mask1_.Nselected() * 3;
    int Ncoords = (int)crd_indices_.size();
//#   pragma omp parallel private(idx2, atomCrd2, offset2, midx, idx1, atomCrd1, offset1, Mj)
#   pragma omp parallel private(idx2, midx, idx1, Mj)
    {
#   pragma omp for schedule(dynamic)
    for (idx2 = 0; idx2 < Ncoords; idx2++)
    {
      //atomCrd2 = mask1_[idx2 / 3] * 3;
      //offset2 = idx2 % 3;
      //Mj = currentFrame[atomCrd2 + offset2];
      Mj = currentFrame[ crd_indices_[idx2] ];
      vect1[idx2] += Mj;
      vect2_[idx2] += (Mj * Mj);
      midx = (idx2 * (int)matrix.Ncols() - (idx2 * (idx2-1) / 2));
      for (idx1 = idx2; idx1 < Ncoords; idx1++, midx++)
      {
        //atomCrd1 = mask1_[idx1 / 3] * 3;
        //offset1 = idx1 % 3;
        //matrix[midx] += (currentFrame[atomCrd1 + offset1] * Mj);
        matrix[midx] += (currentFrame[ crd_indices_[idx1] ] * Mj);
      }
    } // END for loop
    } // END openmp pragma
  }
#else /* Original matrix parallelization */
  int m1_idx, m2_idx;
  double Vj;
  const double* XYZi;
  const double* XYZj;
  unsigned int ny;
  DataSet_MatrixDbl::iterator mat;
  v_iterator v1, v2;
  if (useMask2_) { // FULL MATRIX
    int NX = (int)Mat_->Ncols();
    int crd_max = (int)crd_indices_.size();
#   pragma omp parallel private(m1_idx, m2_idx, XYZi, XYZj, Vj, mat, v1, v2, ny)
    {
    #pragma omp for
    for (m2_idx = 0; m2_idx < mask2_.Nselected(); m2_idx++) {
      mat = Mat_->begin() + ((m2_idx*3)*NX);
      XYZj = currentFrame.XYZ( mask2_[m2_idx] );
      for (ny = 0; ny < 3; ny++) {
        Vj = XYZj[ny];
        for (m1_idx = 0; m1_idx < mask1_.Nselected(); m1_idx++) {
          XYZi = currentFrame.XYZ( mask1_[m1_idx] );
          *(mat++) += Vj * XYZi[0];
          *(mat++) += Vj * XYZi[1];
          *(mat++) += Vj * XYZi[2];
        }
      }
    }
    // Mask1/Mask2 diagonal
    #   pragma omp for
    for (m1_idx = 0; m1_idx < crd_max; m1_idx++) {
      v1 = Mat_->v1begin() + (m1_idx * 3); // Index into vect/vect2
      v2 = vect2_.begin()  + (m1_idx * 3);
      StoreVec(v1, v2, currentFrame.CRD( crd_indices_[m1_idx] ));
    }
    } // END PARALLEL BLOCK FULL
    return;
  } else {         // HALF MATRIX
    int v_idx;
    unsigned int nx;
    double d_m2_idx;
    double TwoN = (double)( Mat_->Ncols() * 2 );
#   pragma omp parallel private(m1_idx, m2_idx, d_m2_idx, v_idx, XYZi, XYZj, Vj, mat, v1, v2, ny, nx)
    {
    #pragma omp for schedule(dynamic)
    for (m2_idx = 0; m2_idx < mask1_.Nselected(); m2_idx++) {
      v_idx = m2_idx * 3;
      d_m2_idx = (double)v_idx;
      mat = Mat_->begin() + (int)(0.5*d_m2_idx*(TwoN-d_m2_idx-1.0)+d_m2_idx);
      v1 = Mat_->v1begin() + v_idx;
      v2 = vect2_.begin() + v_idx;
      XYZj = currentFrame.XYZ( mask1_[m2_idx] );
      StoreVec(v1, v2, XYZj);
      for (ny = 0; ny < 3; ny++) {
        Vj = XYZj[ny];
        // m1_idx = m2_idx, diagonal
        for (nx = ny; nx < 3; nx++)
          *(mat++) += Vj * XYZj[nx]; // Vj * i{0,1,2}, Vj * i{1,2}, Vj * i{2}
        for (m1_idx = m2_idx+1; m1_idx < mask1_.Nselected(); m1_idx++) {
          XYZi = currentFrame.XYZ( mask1_[m1_idx] );
          *(mat++) += Vj * XYZi[0];
          *(mat++) += Vj * XYZi[1];
          *(mat++) += Vj * XYZi[2];
        }
      }
    }
    } // END PARALLEL BLOCK HALF
  }
#endif
# else
  DataSet_MatrixDbl::iterator mat = Mat_->begin();
  v_iterator v1idx1 = Mat_->v1begin();
  v_iterator v2idx1 = vect2_.begin();
  if (useMask2_) {
    // Full Matrix
    // Position for mask2 halfway through vect/vect2
    v_iterator v1idx2 = Mat_->v1begin() + (mask1_.Nselected() * 3); 
    v_iterator v2idx2 = vect2_.begin()  + (mask1_.Nselected() * 3); 
    bool storeVecj = true; // Only store vecj|vecj^2 first time through inner loop
    // OUTER LOOP
    for (AtomMask::const_iterator atom2 = mask2_.begin(); atom2 != mask2_.end(); ++atom2)
    {
      const double* XYZi = currentFrame.XYZ( *atom2 );
      // Store veci and veci^2
      StoreVec(v1idx2, v2idx2, XYZi);
      // Loop over X, Y, and Z of veci
      for (int iidx = 0; iidx < 3; ++iidx) {
        double Vi = XYZi[iidx];
        // INNER LOOP
        for (AtomMask::const_iterator atom1 = mask1_.begin(); atom1 != mask1_.end(); ++atom1)
        {
          const double* XYZj = currentFrame.XYZ( *atom1 );
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
      const double* XYZi = currentFrame.XYZ( *atom2 );
      // Store veci and veci^2
      StoreVec(v1idx1, v2idx1, XYZi);
      // Loop over X, Y, and Z of veci
      for (int iidx = 0; iidx < 3; ++iidx) {
        double Vi = XYZi[iidx];
        // Diagonal
        for (int jidx = iidx; jidx < 3; jidx++)
          *(mat++) += Vi * XYZi[jidx]; // Vi * j{0,1,2}, Vi * j{1,2}, Vi * j{2}
        // INNER LOOP
        for (AtomMask::const_iterator atom1 = atom2 + 1; atom1 != mask1_.end(); ++atom1)
        {
          const double* XYZj = currentFrame.XYZ( *atom1 );
          *(mat++) += Vi * XYZj[0];
          *(mat++) += Vi * XYZj[1];
          *(mat++) += Vi * XYZj[2];
        }
      }
    }
  }
# endif
}

// Action_Matrix::StoreXY()
void Action_Matrix::StoreXY(v_iterator& v1, v_iterator& v2, const double* XY) const 
{
  *(v1++) += XY[0];
  *(v2++) += (XY[0] * XY[0]);
  *(v1++) += XY[1];
  *(v2++) += (XY[1] * XY[1]);
}

/** Dihedral covariance. */
void Action_Matrix::CalcDihedralCovariance( int frameNum ) {
  double XY2[2];
  DataSet_MatrixDbl::iterator mat = Mat_->begin();
  v_iterator v1idx1 = Mat_->v1begin();
  v_iterator v2idx1 = vect2_.begin();
    // TODO: Pre-calculate thetas
    // FIXME: Only Half Matrix for Now
    // OUTER LOOP
    for (Array1D::const_iterator ds2 = DihedralSets_.begin(); 
                                 ds2 != DihedralSets_.end(); ++ds2)
    {
      double theta2 = (*ds2)->Dval( frameNum ) * Constants::DEGRAD;
      XY2[0] = cos( theta2 );
      XY2[1] = sin( theta2 );
      // Store X and X^2
      StoreXY( v1idx1, v2idx1, XY2 );
      // Loop over X and Y of XY2
      for (int iidx = 0; iidx < 2; ++iidx) {
        double Vi = XY2[iidx];
        // Diagonal
        for (int jidx = iidx; jidx < 2; jidx++)
          *(mat++) += Vi * XY2[jidx]; // Vi * j{0,1}, Vi * j{1}
        // INNER LOOP
        for (Array1D::const_iterator ds1 = ds2 + 1; 
                                     ds1 != DihedralSets_.end(); ++ds1)
        {
          double theta1 = (*ds1)->Dval( frameNum ) * Constants::DEGRAD;
          *(mat++) += Vi * cos( theta1 );
          *(mat++) += Vi * sin( theta1 );
        }
      }
    }
}

/** Calc Isotropically distributed ensemble matrix.
  * See Proteins 2002, 46, 177; eq. 7 
  */
void Action_Matrix::CalcIdeaMatrix(Frame const& currentFrame) {
  DataSet_MatrixDbl::iterator mat = Mat_->begin();
  v_iterator v1idx1 = Mat_->v1begin();
  v_iterator v2idx1 = vect2_.begin();
  // Get COM
  Vec3 COM = currentFrame.VCenterOfMass( mask1_ );
  // Get ri, rj, and calc ri*rj
  // Matrix IDEA only uses 1 mask.
  for (AtomMask::const_iterator atomi = mask1_.begin(); atomi != mask1_.end(); ++atomi)
  {
    Vec3 ri = currentFrame.XYZ(*atomi);
    ri -= COM;
    for (AtomMask::const_iterator atomj = atomi; atomj != mask1_.end(); ++atomj)
    {
      Vec3 rj = currentFrame.XYZ(*atomj);
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
void Action_Matrix::CalcCorrelationMatrix(Frame const& currentFrame) {
  DataSet_MatrixDbl::iterator mat = Mat_->begin();
  v_iterator v1idx1 = Mat_->v1begin();
  v_iterator v2idx1 = vect2_.begin();
  if (useMask2_) {
    // Full Matrix
    // Position for mask2 halfway through vect/vect2
    v_iterator v1idx2 = Mat_->v1begin() + (mask1_.Nselected() * 3);
    v_iterator v2idx2 = vect2_.begin() + (mask1_.Nselected() * 3);
    bool storeVecj = true; // Only store vecj|vecj^2 first time through inner loop
    // OUTER LOOP
    for (AtomMask::const_iterator atom2 = mask2_.begin(); atom2 != mask2_.end(); ++atom2)
    {
      const double* XYZi = currentFrame.XYZ( *atom2 );
      // Store veci and veci^2
      StoreVec(v1idx2, v2idx2, XYZi);
      // INNER LOOP
      for (AtomMask::const_iterator atom1 = mask1_.begin(); atom1 != mask1_.end(); ++atom1)
      {
        const double* XYZj = currentFrame.XYZ( *atom1 );
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
      const double* XYZi = currentFrame.XYZ( *atom2 );
      // Store veci and veci^2
      StoreVec(v1idx1, v2idx1, XYZi);
      for (AtomMask::const_iterator atom1 = atom2; atom1 != mask1_.end(); ++atom1)
      {
        const double* XYZj = currentFrame.XYZ( *atom1 );
        *(mat++) += (XYZi[0]*XYZj[0] + XYZi[1]*XYZj[1] + XYZi[2]*XYZj[2]);
      }
    }
  }
}

/** Calculate distance covariance matrix. */
void Action_Matrix::CalcDistanceCovarianceMatrix(Frame const& currentFrame) {
  // Calculate all distance pairs for mask 1
  v_iterator pair_j = vect2_.begin();
  AtomMask::const_iterator mask1end = mask1_.end() - 1;
  for (AtomMask::const_iterator atom1 = mask1_.begin(); atom1 != mask1end; ++atom1)
    for (AtomMask::const_iterator atom2 = atom1 + 1; atom2 != mask1_.end(); ++atom2)
      *(pair_j++) = sqrt(DIST2_NoImage(currentFrame.XYZ(*atom1), currentFrame.XYZ(*atom2)));
  // Create matrix from all distance pairs
  DataSet_MatrixDbl::iterator mat = Mat_->begin();
  v_iterator v1idx = Mat_->v1begin();
  for (v_iterator pair_i = vect2_.begin(); pair_i != vect2_.end(); ++pair_i)
  {
    for (pair_j = pair_i; pair_j != vect2_.end(); ++pair_j) 
    {
      *(mat++) += ( (*pair_i) * (*pair_j) );
      if (pair_i == pair_j)
        *(v1idx++) += (*pair_i);
    }
  }
}

// Action_Matrix::DoAction()
Action::RetType Action_Matrix::DoAction(int frameNum, ActionFrame& frm) {
  // Check if this frame should be processed
  if ( CheckFrameCounter( frm.TrajoutNum() ) ) return Action::OK;
  // Increment number of snapshots
  Mat_->IncrementSnapshots();

  switch (Mat_->Meta().ScalarType()) {
    case MetaData::DIST     : CalcDistanceMatrix(frm.Frm()); break;
    case MetaData::COVAR    :
    case MetaData::MWCOVAR  : CalcCovarianceMatrix(frm.Frm()); break;
    case MetaData::CORREL   : CalcCorrelationMatrix(frm.Frm()); break;
    case MetaData::DIHCOVAR : CalcDihedralCovariance(frameNum); break;
    case MetaData::DISTCOVAR: CalcDistanceCovarianceMatrix(frm.Frm()); break;
    case MetaData::IDEA     : CalcIdeaMatrix(frm.Frm()); break;
    case MetaData::IREDMAT  : CalcIredMatrix(frameNum); break;
    default: return Action::ERR; // Sanity check
  }

  return Action::OK;
}

#ifdef MPI
int Action_Matrix::SyncAction() {
  if (!vect2_.empty()) {
    if (trajComm_.Master()) {
      Darray buf( vect2_.size() );
      trajComm_.ReduceMaster( &(buf[0]), &(vect2_[0]), vect2_.size(), MPI_DOUBLE, MPI_SUM );
      vect2_ = buf;
    } else
      trajComm_.ReduceMaster( 0,         &(vect2_[0]), vect2_.size(), MPI_DOUBLE, MPI_SUM );
  }
  return 0;
}
#endif
// -----------------------------------------------------------------------------
void Action_Matrix::Vect2MinusVect() {
  v_iterator v2 = vect2_.begin();
  for (v_iterator v1 = Mat_->v1begin(); v1 != Mat_->v1end(); ++v1)
    *(v2++) -= ( *v1 * *v1 );
}

// Action_Matrix::FinishCovariance()
void Action_Matrix::FinishCovariance(size_t element_size) {
  double Mass = 1.0;
  double mass2 = 1.0;

  DataSet_MatrixDbl::iterator mat = Mat_->begin();
  // Calc <riri> - <ri><ri>
  Vect2MinusVect(); // Is vect2 used?
  // Calc <rirj> - <ri><rj>
  if (useMask2_) {
    // Full Matrix
    M_iterator m2 = mass2_.begin();
    // Position for mask2 halfway through vect/vect2
    v_iterator v1idx2begin = Mat_->v1begin() + Mat_->Ncols();
    for (v_iterator v1idx2 = v1idx2begin; v1idx2 != Mat_->v1end(); v1idx2 += element_size)
    {
      if (Mat_->Meta().ScalarType() == MetaData::MWCOVAR)
        mass2 = *(m2++);
      for (unsigned int iidx = 0; iidx < element_size; ++iidx) {
        M_iterator m1 = mass1_.begin();
        double Vi = *(v1idx2 + iidx);
        for (v_iterator v1idx1 = Mat_->v1begin(); v1idx1 != v1idx2begin; v1idx1 += element_size)
        {
          if (Mat_->Meta().ScalarType() == MetaData::MWCOVAR)
            Mass = sqrt( mass2 * *(m1++) );
          for (unsigned int idx = 0; idx < element_size; ++idx) {
            *mat = (*mat - (Vi * *(v1idx1+idx))) * Mass;
            ++mat;
          }
        }
      }
    }
  } else {
    // Half Matrix
    M_iterator m2 = mass1_.begin();
    for (v_iterator v1idx2 = Mat_->v1begin(); v1idx2 != Mat_->v1end(); v1idx2 += element_size)
    {
      if (Mat_->Meta().ScalarType() == MetaData::MWCOVAR)
        mass2 = *m2;
      for (unsigned int iidx = 0; iidx < element_size; ++iidx) {
        M_iterator m1 = m2;
        double Vi = *(v1idx2 + iidx);
        for (v_iterator v1idx1 = v1idx2; v1idx1 != Mat_->v1end(); v1idx1 += element_size)
        {
          if (Mat_->Meta().ScalarType() == MetaData::MWCOVAR)
            Mass = sqrt( mass2 * *(m1++) );
          if ( v1idx1 == v1idx2 ) {
            for (unsigned int jidx = iidx; jidx < element_size; ++jidx) {
              *mat = (*mat - (Vi * *(v1idx1 + jidx))) * Mass;
              ++mat;
            }
          } else {
            for (unsigned int idx = 0; idx < element_size; ++idx) {
              *mat = (*mat - (Vi * *(v1idx1+idx))) * Mass;
              ++mat;
            }
          }
        }
      }
      ++m2;
    }
  }
}

// DotProdAndNorm()
void Action_Matrix::DotProdAndNorm(DataSet_MatrixDbl::iterator& mat,
                                   v_iterator& vecti, 
                                   v_iterator& vectj,
                                   v_iterator& vect2i,
                                   v_iterator& vect2j)
const {
  *(mat) -= ( *(vectj  ) * *(vecti  ) +
              *(vectj+1) * *(vecti+1) +
              *(vectj+2) * *(vecti+2)   );
  // Normalize
  *(mat++) /= sqrt( (*(vect2j) + *(vect2j+1) + *(vect2j+2)) *
                    (*(vect2i) + *(vect2i+1) + *(vect2i+2))   );
}

// Action_Matrix::FinishCorrelation()
void Action_Matrix::FinishCorrelation() {
  DataSet_MatrixDbl::iterator mat = Mat_->begin();
  // Calc <ri * ri> - <ri> * <ri>
  Vect2MinusVect();
  // Calc <ri * rj> - <ri> * <rj>
  if (useMask2_) {
    // Full Matrix
    // Position for mask2 halfway through vect/vect2
    // Vect has 3 entries per atom, but matrix only has 1 (as opposed to 
    // COVAR/MWCOVAR which has 3 for vect and matrix).
    v_iterator v1idx2begin = Mat_->v1begin() + Mat_->Ncols() * 3;
    v_iterator v2idx2      = vect2_.begin()  + Mat_->Ncols() * 3;
    for (v_iterator v1idx2 = v1idx2begin; v1idx2 != Mat_->v1end(); v1idx2 += 3)
    {
      v_iterator v2idx1 = vect2_.begin();
      for (v_iterator v1idx1 = Mat_->v1begin(); v1idx1 != v1idx2begin; v1idx1 += 3)
      {
        DotProdAndNorm( mat, v1idx2, v1idx1, v2idx2, v2idx1 );
        v2idx1 += 3;
      }
      v2idx2 += 3;
    }
  } else {
    // Half Matrix
    v_iterator v2idx2 = vect2_.begin();
    for (v_iterator v1idx2 = Mat_->v1begin(); v1idx2 != Mat_->v1end(); v1idx2 += 3)
    {
      v_iterator v2idx1 = v2idx2;
      for (v_iterator v1idx1 = v1idx2; v1idx1 != Mat_->v1end(); v1idx1 += 3)
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
  DataSet_MatrixDbl::iterator mat = Mat_->begin();
  for (v_iterator pair_i = Mat_->v1begin(); pair_i != Mat_->v1end(); ++pair_i)
    for (v_iterator pair_j = pair_i; pair_j != Mat_->v1end(); ++pair_j) 
      *(mat++) -= ((*pair_i) * (*pair_j));
}

double Action_Matrix::ByMaskAverage(unsigned int ColMax, unsigned int RowMax) const {
  double val = 0.0;
  double valnorm = 0.0;
  for (unsigned int idx2 = 0; idx2 != RowMax; idx2++) {
    for (unsigned int idx1 = 0; idx1 != ColMax; idx1++) {
      valnorm += (mass1_[idx1] * mass2_[idx2]);
      val += Mat_->GetElement(idx1, idx2);
    }
  }
  if (valnorm > 0.0) return val / valnorm;
  return 0.0;
}

// Action_Matrix::Print()
void Action_Matrix::Print() {
  if (debug_ > 1) {
    mprintf("Raw Matrix Elements:\n");
    for (unsigned int i = 0; i < Mat_->Size(); i++)
      mprintf("\t%u\t%f\n", i, (*Mat_)[i]);
    mprintf("Raw Vect1 Elements:\n");
    for (unsigned int i = 0; i < Mat_->Vect().size(); i++)
      mprintf("\t%u\t%f\n", i, Mat_->Vect()[i]);
    mprintf("Raw Vect2 Elements:\n");
    for (unsigned int i = 0; i < vect2_.size(); i++)
      mprintf("\t%u\t%f\n", i, vect2_[i]);
  }
  // ---------- Calculate average over number of sets ------
  if (Mat_->Nsnapshots() < 1) {
    mprintf("Warning: Matrix %s is empty.\n", Mat_->legend());
    return;
  }
  double norm = (double)Mat_->Nsnapshots();
  if (Mat_->Meta().ScalarType() == MetaData::IDEA) norm *= 3.0;
  norm = 1.0 / norm;
  for (v_iterator v1 = Mat_->v1begin(); v1 != Mat_->v1end(); ++v1)
    *v1 *= norm;
  for (v_iterator v2 = vect2_.begin();  v2 != vect2_.end();  ++v2)
    *v2 *= norm;
  for (DataSet_MatrixDbl::iterator m = Mat_->begin(); m != Mat_->end(); ++m)
    *m *= norm;

  switch (Mat_->Meta().ScalarType()) {
    case MetaData::COVAR    :
    case MetaData::MWCOVAR  : FinishCovariance(3); break;
    case MetaData::DIHCOVAR : FinishCovariance(2); break;
    case MetaData::CORREL   : FinishCorrelation(); break;
    case MetaData::DISTCOVAR: FinishDistanceCovariance(); break;
    default: break; 
  }

  // Perform averaging BYRESIDUE or BYMASK if necessary.
  if (outtype_ == BYRESIDUE) {
    // ---------- Print out BYRESIDUE
    // First determine which mask elements correspond to which residues,
    // as well as how many residues are selected in mask1 (cols) and
    // mask2 (rows) to determine output matrix dimensions.
    if (!useMask2_) {
      mask2_ = mask1_;
      mass2_ = mass1_;
      residues2_ = residues1_;
    }
    // Always allocate full matrix for BYRESIDUE
    matByRes_->Allocate2D( residues1_.size(), residues2_.size() ); // cols, rows
    mprintf("    MATRIX: By-residue matrix has %u rows, %u columns.\n", 
            matByRes_->Nrows(), matByRes_->Ncols());
    for (MatResArray::const_iterator resj = residues2_.begin(); resj != residues2_.end(); ++resj)
    {
      for (MatResArray::const_iterator resi = residues1_.begin(); resi != residues1_.end(); ++resi)
      {
        double val = 0.0;
        double valnorm = 0.0;
        for (Iarray::const_iterator itj = resj->maskIdxs_.begin();
                                    itj != resj->maskIdxs_.end(); ++itj)
        {
          for (Iarray::const_iterator iti = resi->maskIdxs_.begin();
                                      iti != resi->maskIdxs_.end(); ++iti)
          {
            valnorm += mass1_[*iti] * mass2_[*itj];
            val += Mat_->GetElement( *iti, *itj ); // col, row
          }
        }
        matByRes_->AddElement( val / valnorm );
      }
    }
  } else if (outtype_ == BYMASK) {
    // ---------- Print out BYMASK
    if (!useMask2_) {
      mass2_ = mass1_;
      mprintf("    MATRIX: Writing internal average over mask '%s'\n", mask1_.MaskString());
      byMaskOut_->Printf("%6.2f \n", ByMaskAverage(mask1_.Nselected(), mask1_.Nselected()));
    } else {
      mprintf("    MATRIX: Writing internal averages for mask1 '%s' and mask2 '%s':\n"
              "            mask1/mask1, mask1/mask2, mask2/mask2\n",
              mask1_.MaskString(), mask2_.MaskString());
      byMaskOut_->Printf("%6.2f %6.2f %6.2f \n",
                         ByMaskAverage(mask1_.Nselected(), mask1_.Nselected()),
                         ByMaskAverage(mask1_.Nselected(), mask2_.Nselected()),
                         ByMaskAverage(mask2_.Nselected(), mask2_.Nselected()));
    }
  }
}
