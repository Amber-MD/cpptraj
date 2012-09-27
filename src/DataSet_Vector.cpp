#include <cmath> // sqrt
#include "DataSet_Vector.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataSet_Vector::DataSet_Vector() :
  totalFrames_(-1),
  frame_(0),
  mode_(VECTOR_NOOP),
  cx_(0),
  cy_(0),
  cz_(0),
  vx_(0),
  vy_(0),
  vz_(0),
  master_(0),
  modinfo_(0),
  ibeg_(1),
  iend_(50),
  order_(2),
  npair_(-1),
  rave_(0),
  r3iave_(0),
  r6iave_(0),
  cftmp_(0),
  p2cftmp_(0),
  rcftmp_(0)
{
  width_ = 8;
  precision_ = 4;
  dType_ = VECTOR;
  SetDataSetFormat(false);
  avgcrd_[0] = 0;
  avgcrd_[1] = 0;
  avgcrd_[2] = 0;
  // DEBUG
  //debugpdb.SetupWrite("PRINCIPAL.PDB",0);
  //debugpdb.OpenFile();
  //debuginert.SetupWrite("INERT.PDB",0);
  //debuginert.OpenFile();
}

const char DataSet_Vector::ModeString[11][12] = {
  "NO_OP", "Principal X", "Principal Y", "Principal Z",
  "Dipole", "Box", "Mask", "Ired", 
  "CorrPlane", "Corr", "CorrIred"
};

// DESTRUCTOR
DataSet_Vector::~DataSet_Vector() {
  // Free shared CORRIRED mem on master
  if (mode_==VECTOR_CORRIRED) {
    if (master_==0) { 
      if (modinfo_!=0) delete modinfo_;
      if (cftmp_!=0) delete[] cftmp_;
    }
  } else {
    if (cftmp_!=0) delete[] cftmp_;
    if (p2cftmp_!=0) delete[] p2cftmp_;
    if (rcftmp_!=0) delete[] rcftmp_;
    if (cx_!=0) delete[] cx_;
    if (cy_!=0) delete[] cy_;
    if (cz_!=0) delete[] cz_;
    if (vx_!=0) delete[] vx_;
    if (vy_!=0) delete[] vy_;
    if (vz_!=0) delete[] vz_;
  }
  // DEBUG
  //debugpdb.CloseFile();
  //debuginert.CloseFile();
}

// DataSet_Vector::operator==()
/** If vector modesfile names are the same, and the beginning and ending
  * vector number is the same, then consider these two vectors the same.
  * Currently only used for CORRIRED. 
  */
bool DataSet_Vector::operator==(const DataSet_Vector& rhs) {
  if (mode_ == rhs.mode_ &&
      modesfile_ == rhs.modesfile_ &&
      ibeg_ == rhs.ibeg_ &&
      iend_ == rhs.iend_) 
    return true;
  return false; 
}

int DataSet_Vector::SetModeFromArgs( ArgList& actionArgs ) {
  // Get order for Legendre polynomial
  order_ = actionArgs.getKeyInt("order",2);
  if (order_ < 0 || order_ > 2) {
    mprintf("Warning: vector order out of bounds (<0 or >2), resetting to 2.\n");
    order_ = 2;
  }
  // Acceptable args: principal [x | y | z], dipole, box, corrplane, 
  //                  corrired, corr, ired
  if ( actionArgs.hasKey("principal") ) {
    mode_ = VECTOR_PRINCIPAL_X;
    if ( actionArgs.hasKey("x") ) mode_ = VECTOR_PRINCIPAL_X;
    if ( actionArgs.hasKey("y") ) mode_ = VECTOR_PRINCIPAL_Y;
    if ( actionArgs.hasKey("z") ) mode_ = VECTOR_PRINCIPAL_Z;
  } else if (actionArgs.hasKey("dipole"))
    mode_ = VECTOR_DIPOLE;
  else if (actionArgs.hasKey("box"))
    mode_ = VECTOR_BOX;
  else if (actionArgs.hasKey("corrplane"))
    mode_ = VECTOR_CORRPLANE;
  else if (actionArgs.hasKey("corrired"))
    mode_ = VECTOR_CORRIRED;
  else if (actionArgs.hasKey("corr"))
    mode_ = VECTOR_CORR;
  else if (actionArgs.hasKey("ired"))
    mode_ = VECTOR_IRED;
  else
    mode_ = VECTOR_MASK;

  if (mode_ == VECTOR_NOOP) {
    mprinterr("Error: No vector mode specified.\n");
    return 1;
  }
  return 0;
}

int DataSet_Vector::SetupCorrIred( DataSet_Vector* Vmaster, ArgList& actionArgs ) {
  if (mode_ == VECTOR_CORRIRED) {
    // Get Pair number
    npair_ = actionArgs.getKeyInt("npair",0);
    if (npair_ == 0) {
      mprinterr("Error: Action_Vector::Init(): No 'npair <#>' arg given, needed for 'corrired'.\n");
      return 1;
    }
    // Actual pair number needs to start from 0
    --npair_;
    modesfile_ = actionArgs.GetStringKey("modes");
    if (modesfile_.empty()) {
      mprinterr("Error: Action_Vector::Init(): No 'modes' <file> arg given, needed for 'corrired'.\n");
      return 1;
    }
    ibeg_ = actionArgs.getKeyInt("beg",1);
    iend_ = actionArgs.getKeyInt("end", 50);
    // check if modes need to be read or have been
    // previously read in by another DataSet_Vector.
    if (Vmaster != 0) {
      if (AssignMaster( Vmaster )) {
        mprinterr("Error: Could not assign vector master for CORRIRED.\n");
        return 1;
      }
    }
    // Load modes from file specified by modesfile. 
    if (modinfo_==0)
      if (ReadModesFromFile()) {
        return 1;
      }
  }
  return 0;
}

// DataSet_Vector::AssignMaster()
/** Assign this vectors ModesInfo to vecIn. vecIn is considered the
  * 'master' vector as far as ModesInfo is concerned. This means
  * among other things that this vector will not try to free
  * modinfo during destruction.
  */
int DataSet_Vector::AssignMaster( DataSet_Vector* vecIn ) {
  if (vecIn==0) return 1;
  if (vecIn->modinfo_==0) return 1;
  modinfo_ = vecIn->modinfo_;
  master_ = vecIn;
  return 0;
}

// DataSet_Vector::ReadModesFromFile()
int DataSet_Vector::ReadModesFromFile() {
  if (modesfile_.empty()) {
    mprinterr("Error: DataSet_Vector::ReadModesFromFile: filename not set.\n");
    return 1;
  }
  // NOTE: Warn here? If any other vecs linked this will mess them up.
  if (master_==0 && modinfo_!=0) {
    mprintf("Warning: Vector CORRIRED: Replacing ModesInfo.\n");
    delete modinfo_;
  }
  modinfo_ = new ModesInfo();
  return (modinfo_->ReadEvecFile( modesfile_, ibeg_, iend_ ));
}

// DataSet_Vector::Allocate()
int DataSet_Vector::Allocate(int maxFrames) {
  // Allocate data
  if (maxFrames <=0) {
    mprinterr("Error: VECTOR in cpptraj currently unsupported for unknown # of frames.\n");
    return 1;
  }

  if (mode_ == VECTOR_IRED)
    totalFrames_ = 1;
  else
    totalFrames_ = maxFrames;

  // NOTE: cftmp_ stores the complex numbers resulting from conversion of 
  //       vector coords to spherical harmonics. For each frame:
  //       [m0R][m0I][
  // CORRIRED Allocation ---------------
  if (mode_==VECTOR_CORRIRED) {
    // AssignMaster should be called if this does not have original modeinfo
    if (master_==0) {
      // This is the master with modeinfo. Allocate mem for cftmp.
      int ntotal = 2 * totalFrames_ * (2 * order_ + 1) * modinfo_->Nvect();
      cftmp_ = new double[ ntotal ];
      for (int i = 0; i < ntotal; ++i)
        cftmp_[i] = 0; 
    } else {
      // This is a child with link to master modeinfo. Also link to cftmp, which 
      // should already be allocated on master.
      if (master_->cftmp_==0) {
        mprinterr("Error: Child vector %s cannot link to master vector %s\n", Legend().c_str(),
                  master_->Legend().c_str());
        return 1;
      }
      cftmp_ = master_->cftmp_;
    }

  // CORRPLANE / CORR Allocation -------
  } else if (mode_ == VECTOR_CORRPLANE || mode_ == VECTOR_CORR) {
    // avgcrd is now static size [3]
    int ntotal = 2 * totalFrames_ * (2 * order_ + 1);
    cftmp_ = new double[ ntotal ];
    p2cftmp_ = new double[ ntotal ];
    for (int i = 0; i < ntotal; ++i) {
      cftmp_[i] = 0; 
      p2cftmp_[i] = 0;
    }
    ntotal = 2 * totalFrames_;
    rcftmp_ = new double[ ntotal ];
    for (int i = 0; i < ntotal; ++i) 
      rcftmp_[i] = 0;
    // NOTE: CORRPLANE also needs cx/cy/cz allocd to # of selected in mask

  // Everything Else -------------------
  } else {
    //C_.resize( totalFrames_ );
    //V_.resize( totalFrames_ );
    cx_ = new double[ totalFrames_ ];
    cy_ = new double[ totalFrames_ ];
    cz_ = new double[ totalFrames_ ];
    vx_ = new double[ totalFrames_ ];
    vy_ = new double[ totalFrames_ ];
    vz_ = new double[ totalFrames_ ];
  } 

  return 0; 
}

// DataSet_Vector::AllocCorrPlane()
int DataSet_Vector::AllocCorrPlane( int Nselected ) {
  if (mode_==VECTOR_CORRPLANE) {
    if (cx_!=0 || cy_!=0 || cz_!=0)
    //if (!C_.empty())
    {
      mprinterr("Error: VECTOR_CORRPLANE not yet supported for multiple topology files.\n");
      return 1;
    }
    if (Nselected < 3) {
      mprinterr("Error: < 3 atoms specified for VECTOR_CORRPLANE\n");
      return 1;
    }
    //C_.resize( mask_.Nselected() );
    cx_ = new double[ Nselected ];
    cy_ = new double[ Nselected ];
    cz_ = new double[ Nselected ];
  }
  return 0;
}

// DataSet_Vector::Info()
void DataSet_Vector::Info(const char* mask1, const char* mask2) {
  mprintf("    VECTOR: Storage to array named %s\n", Legend().c_str() );
  mprintf("            Vector Mode: %s\n", ModeString[mode_]);
  switch (mode_) {
    case VECTOR_DIPOLE:
      mprintf("\tThe dipole moment vector with respect to the center of mass\n");
      mprintf("\twill be processed for atoms in mask expression [%s]\n", mask1);
      break;

    case VECTOR_PRINCIPAL_X:
    case VECTOR_PRINCIPAL_Y:
    case VECTOR_PRINCIPAL_Z:
      mprintf("\tThe principal axis vector for");
      if (mode_==VECTOR_PRINCIPAL_X) mprintf(" (X)");
      else if (mode_==VECTOR_PRINCIPAL_Y) mprintf(" (Y)");
      else if (mode_==VECTOR_PRINCIPAL_Z) mprintf(" (Z)");
      mprintf(" will be calcd with respect to the\n");
      mprintf("\tcenter of mass of atoms in mask expression [%s]\n", mask1);
      break;

    case VECTOR_MASK:
    case VECTOR_IRED:
    case VECTOR_CORR:
    case VECTOR_CORRIRED:
      mprintf("\tCalculate the vector between the center of mass of the two atom selections\n");
      mprintf("\twhich follow (with the origin at the center of mass of the first)\n");
      mprintf("\tMask1: [%s]  Mask2: [%s]\n", mask1, mask2);
      break;

    case VECTOR_CORRPLANE:
      mprintf("\tCalculate the vector perpendicular to the least squares best plane\n");
      mprintf("\tthrough the atom selection [%s]\n",mask1);
      break;

    case VECTOR_BOX:
      mprintf("\tThe box lengths will be treated as a vector\n");
      break;

    default:
      mprinterr("Error: No vector mode defined.\n");
  }

  if (mode_ == VECTOR_CORRPLANE ||
      mode_ == VECTOR_CORR ||
      mode_ == VECTOR_CORRIRED)
    mprintf("\tThe order of Legendre polynomials is %i\n", order_);

  if (mode_ == VECTOR_CORRIRED) {
    mprintf("\tIRED modes are read from %s,\n", modesfile_.c_str());
    mprintf("\tand the pair %i is considered\n", npair_+1);
  }

}

// -----------------------------------------------------------------------------
int DataSet_Vector::Size() {
  return totalFrames_;
}

int DataSet_Vector::Xmax() {
  return totalFrames_ - 1;
}

void DataSet_Vector::WriteBuffer(CpptrajFile &cbuffer, int frameIn) {
  if (frameIn < 0 || frameIn >= frame_) {
    mprinterr("ERROR: VECTOR: Frame %i is out of range.\n",frameIn);
    return;
  }
  cbuffer.Printf(data_format_, vx_[frameIn]);
  cbuffer.Printf(data_format_, vy_[frameIn]);
  cbuffer.Printf(data_format_, vz_[frameIn]);
  cbuffer.Printf(data_format_, cx_[frameIn]);
  cbuffer.Printf(data_format_, cy_[frameIn]);
  cbuffer.Printf(data_format_, cz_[frameIn]);
  cbuffer.Printf(data_format_, vx_[frameIn]+cx_[frameIn]);
  cbuffer.Printf(data_format_, vy_[frameIn]+cy_[frameIn]);
  cbuffer.Printf(data_format_, vz_[frameIn]+cz_[frameIn]);
}

int DataSet_Vector::Width() {
  return ( ((width_+1)*9) );
}

// -----------------------------------------------------------------------------

// NOTE: Only used in 'analyze timecorr'
void DataSet_Vector::PrintAvgcrd(CpptrajFile& outfile) {
  double dnorm = 1.0 / (double)frame_;
  outfile.Printf("%10.4f %10.4f %10.4f %10.4f\n",
          rave_ * dnorm,
          sqrt(avgcrd_[0]*avgcrd_[0] + avgcrd_[1]*avgcrd_[1] + avgcrd_[2]*avgcrd_[2]) * dnorm,
          r3iave_ * dnorm,
          r6iave_ * dnorm);
}

