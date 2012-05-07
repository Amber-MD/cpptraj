#include <cmath>
#include <algorithm> // sort
#include "VectorType.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"
#include "Constants.h" // PI
#include "vectormath.h"
#include "Matrix_3x3.h"

// CONSTRUCTOR
VectorType::VectorType() :
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
  debugpdb.SetupWrite("PRINCIPAL.PDB",0);
  debugpdb.OpenFile();
  //debuginert.SetupWrite("INERT.PDB",0);
  //debuginert.OpenFile();
}

// DESTRUCTOR
VectorType::~VectorType() {
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
  debugpdb.CloseFile();
  //debuginert.CloseFile();
}

// VectorType::operator==()
bool VectorType::operator==(const VectorType& rhs) {
  if (mode_ == rhs.mode_ &&
      modesfile_ == rhs.modesfile_ &&
      ibeg_ == rhs.ibeg_ &&
      iend_ == rhs.iend_) 
    return true;
  return false; 
}

// VectorType::AssignMaster()
int VectorType::AssignMaster( VectorType* vecIn ) {
  if (vecIn==0) return 1;
  if (vecIn->modinfo_==NULL) return 1;
  modinfo_ = vecIn->modinfo_;
  master_ = vecIn;
  return 0;
}

// VectorType::ReadModesFromFile()
int VectorType::ReadModesFromFile() {
  if (modesfile_.empty()) {
    mprinterr("Error: VectorType::ReadModesFromFile: filename not set.\n");
    return 1;
  }
  // NOTE: Warn here? If any other vecs linked this will mess them up.
  if (master_==0 && modinfo_!=0) delete modinfo_;
  modinfo_ = new ModesInfo();
  return (modinfo_->ReadEvecFile( modesfile_, ibeg_, iend_ ));
}

// VectorType::init()
int VectorType::init() {
  std::string filename_ = actionArgs.GetStringKey("out");

  // Require a vector name - this behavior is consistent with ptraj
  name_ = actionArgs.GetStringNext();

  // Get order for Legendre polynomial
  order_ = actionArgs.getKeyInt("order",2);
  if (order_ < 0 || order_ > 2) {
    mprintf("Warning: vector order out of bounds (<0 or >2), resetting to 2.\n");
    order_ = 2;
  }

  // Determine vector mode.
  // Acceptable args: principal, principal x, principal y, principal z,
  //                  dipole, box, corrplane, corrired, corr, ired
  int pos = actionArgs.KeyPosition("principal");
  if (pos == actionArgs.Nargs() - 1) {
    actionArgs.MarkArg(pos);
    mode_ = VECTOR_PRINCIPAL_X;
  } else if (pos!=-1) {
    actionArgs.MarkArg(pos);
    mode_ = VECTOR_PRINCIPAL_X;
    char vecchar = actionArgs[pos+1][0];
    if (vecchar == 'x' ||
        vecchar == 'y' ||
        vecchar == 'z')
    {
      actionArgs.MarkArg(pos+1);
      switch (vecchar) {
        case 'x': mode_ = VECTOR_PRINCIPAL_X; break;
        case 'y': mode_ = VECTOR_PRINCIPAL_Y; break;
        case 'z': mode_ = VECTOR_PRINCIPAL_Z; break;
      }
    }
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

  if (mode_ == VectorType::VECTOR_NOOP) {
    mprinterr("Error: No vector mode specified.\n");
    return 1;
  }

  // VECTOR_CORRIRED
  if (mode_ == VECTOR_CORRIRED) {
    // Get Pair number
    npair_ = actionArgs.getKeyInt("npair",0);
    if (npair_ == 0) {
      mprinterr("Error: VectorType::Init(): No 'npair <#>' arg given, needed for 'corrired'.\n");
      return 1;
    }
    // Actual pair number needs to start from 0
    --npair_;
    modesfile_ = actionArgs.GetStringKey("modes");
    if (modesfile_.empty()) {
      mprinterr("Error: VectorType::Init(): No 'modes' <file> arg given, needed for 'corrired'.\n");
      return 1;
    }
    ibeg_ = actionArgs.getKeyInt("beg",1);
    iend_ = actionArgs.getKeyInt("end", 50);
    // check if modes need to be read or have been
    // previously read in by another VectorType.
    VectorType *Vtmp;
    DSL->VectorBegin();
    while ( (Vtmp = (VectorType*)DSL->NextVector()) != 0 ) {
      if ( *this == *Vtmp ) {
        if (AssignMaster( Vtmp )) {
          mprinterr("Error: Could not assign vector master for CORRIRED.\n");
          return 1;
        }
      }
    }
    // Load modes from file. modesfile name should be set by VectorType::Init.
    if (modinfo_==0)
      if (ReadModesFromFile()) {
        return 1;
      }
  }

  // Vector Mask
  char *maskexpr = actionArgs.getNextMask();
  mask_.SetMaskString( maskexpr );

  // Get second mask if necessary
  if (mode_ == VECTOR_IRED ||
      mode_ == VECTOR_CORR ||
      mode_ == VECTOR_CORRIRED ||
      mode_ == VECTOR_MASK)
  {
    maskexpr = actionArgs.getNextMask();
    if (maskexpr==NULL) {
      mprinterr("Error: vector: Specified vector mode requires a second mask.\n");
      return 1;
    }
    mask2_.SetMaskString( maskexpr );
  }

  // Allocate vector
  if (Allocate(DSL->MaxFrames())) return 1;

  // Add vector to datasetlist
  // TODO: Check for name conflicts
  DSL->AddDataSet( (DataSet*)this );
  // Since this now exists in the DataSetList and ActionList,
  // set the noDelete flag.
  SetNoDelete();

  Info();

  // Check if output is supported for the current vector mode
  if (mode_ == VECTOR_CORRPLANE ||
      mode_ == VECTOR_CORR ||
      mode_ == VECTOR_CORRIRED ||
      mode_ == VECTOR_IRED)
  {
    if (!filename_.empty()) {
      mprintf(
        "\tWarning: Output of corr, ired, corrired or corrplane vectors is not yet supported!\n");
        filename_.clear();
    }
  }

  // Add to output datafilelist
  if (!filename_.empty()) {
    mprintf("\tOutput will be dumped to a file, %s\n", filename_.c_str());
    DFL->Add( filename_.c_str(), (DataSet*)this );
  }

  return 0;
}

// VectorType::Allocate()
int VectorType::Allocate(int maxFrames) {
  // Allocate data
  if (maxFrames <=0) {
    mprinterr("Error: VECTOR in cpptraj currently unsupported for unknown # of frames.\n");
    return 1;
  }

  if (mode_ == VECTOR_IRED)
    totalFrames_ = 1;
  else
    totalFrames_ = maxFrames;

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
        mprinterr("Error: Child vector %s cannot link to master vector %s\n", c_str(),
                  master_->c_str());
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

// VectorType::Info()
void VectorType::Info() {
  mprintf("    VECTOR: Storage to array named %s\n", c_str() );
  switch (mode_) {
    case VECTOR_DIPOLE:
      mprintf("\tThe dipole moment vector with respect to the center of mass\n");
      mprintf("\twill be processed for atoms in mask expression [%s]\n", mask_.MaskString());
      break;

    case VECTOR_PRINCIPAL_X:
    case VECTOR_PRINCIPAL_Y:
    case VECTOR_PRINCIPAL_Z:
      mprintf("\tThe principal axis vector for");
      if (mode_==VECTOR_PRINCIPAL_X) mprintf(" (X)");
      else if (mode_==VECTOR_PRINCIPAL_Y) mprintf(" (Y)");
      else if (mode_==VECTOR_PRINCIPAL_Z) mprintf(" (Z)");
      mprintf(" will be calcd with respect to the\n");
      mprintf("\tcenter of mass of atoms in mask expression [%s]\n", mask_.MaskString());
      break;

    case VECTOR_MASK:
    case VECTOR_IRED:
    case VECTOR_CORR:
    case VECTOR_CORRIRED:
      mprintf("\tCalculate the vector between the center of mass of the two atom selections\n");
      mprintf("\twhich follow (with the origin at the center of mass of the first)\n");
      mprintf("\tMask1: [%s]  Mask2: [%s]\n", mask_.MaskString(), mask2_.MaskString());
      break;

    case VECTOR_CORRPLANE:
      mprintf("\tCalculate the vector perpendicular to the least squares best plane\n");
      mprintf("\tthrough the atom selection [%s]\n",mask_.MaskString());
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

// VectorType::Print()
/*void VectorType::Print() {
  if (filename_.empty()) return;

  CpptrajFile outfile;
  if (outfile.SetupWrite(filename_.c_str(), 0)) return;
  if (outfile.OpenFile()) return;

  mprintf("CPPTRAJ VECTOR: dumping vector information %s\n", c_str());
  outfile.Printf("# FORMAT: frame vx vy vz cx cy cz cx+vx cy+vy cz+vz\n");
  outfile.Printf("# FORMAT where v? is vector, c? is center of mass...\n");
  for (int i=0; i < totalFrames_; ++i) {
    outfile.Printf("%i %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
            i+1, 
            //V_[i][0], V_[i][1], V_[i][2], C_[i][0], C_[i][1], C_[i][2],
            //V_[i][0]+C_[i][0], V_[i][1]+C_[i][1], V_[i][2]+C_[i][2]
            vx_[i], vy_[i], vz_[i],
            cx_[i], cy_[i], cz_[i],
            cx_[i]+vx_[i], cy_[i]+vy_[i], cz_[i]+vz_[i]);
  }
  outfile.CloseFile();
}*/

// -----------------------------------------------------------------------------
int VectorType::Size() {
  return totalFrames_;
}

int VectorType::Xmax() {
  return totalFrames_ - 1;
}

void VectorType::WriteBuffer(CharBuffer &cbuffer, int frameIn) {
  if (frameIn < 0 || frameIn >= frame_) {
    mprinterr("ERROR: VECTOR: Frame %i is out of range.\n",frameIn);
    return;
  }
  cbuffer.WriteDouble(data_format_, vx_[frameIn]);
  cbuffer.WriteDouble(data_format_, vy_[frameIn]);
  cbuffer.WriteDouble(data_format_, vz_[frameIn]);
  cbuffer.WriteDouble(data_format_, cx_[frameIn]);
  cbuffer.WriteDouble(data_format_, cy_[frameIn]);
  cbuffer.WriteDouble(data_format_, cz_[frameIn]);
  cbuffer.WriteDouble(data_format_, vx_[frameIn]+cx_[frameIn]);
  cbuffer.WriteDouble(data_format_, vy_[frameIn]+cy_[frameIn]);
  cbuffer.WriteDouble(data_format_, vz_[frameIn]+cz_[frameIn]);
}

int VectorType::Width() {
  return ( ((width_+1)*9) );
}

// -----------------------------------------------------------------------------

// VectorType::setup()
int VectorType::setup() {
  // Setup mask 1
  if (currentParm->SetupIntegerMask(mask_)) return 1;
  mprintf("\tVector mask [%s] corresponds to %i atoms.\n",
          mask_.MaskString(), mask_.Nselected());

  // Allocate space for CORRPLANE
  if (mode_==VECTOR_CORRPLANE) {
    if (cx_!=0 || cy_!=0 || cz_!=0) 
    //if (!C_.empty())
    {
      mprinterr("Error: VECTOR_CORRPLANE not yet supported for multiple topology files.\n");
      return 1;
    }
    if (mask_.Nselected() < 3) {
      mprinterr("Error: < 3 atoms specified for VECTOR_CORRPLANE\n");
      return 1;
    }
    //C_.resize( mask_.Nselected() );
    cx_ = new double[ mask_.Nselected() ];
    cy_ = new double[ mask_.Nselected() ];
    cz_ = new double[ mask_.Nselected() ];
  }

  // Setup mask 2
  if (mask2_.MaskStringSet()) {
    if (currentParm->SetupIntegerMask(mask2_)) return 1;
    mprintf("\tVector mask [%s] corresponds to %i atoms.\n",
            mask2_.MaskString(), mask2_.Nselected());
  }
  return 0;
} 

int VectorType::action() {
  int err;

  switch (mode_) {
    case VECTOR_CORRIRED:
    case VECTOR_CORR:
    case VECTOR_CORRPLANE:
      err = Action_CORR(  );
      break;
    case VECTOR_DIPOLE:
      err = Action_DIPOLE(  );
      break;
    case VECTOR_PRINCIPAL_X:
    case VECTOR_PRINCIPAL_Y:
    case VECTOR_PRINCIPAL_Z:
      err = Action_PRINCIPAL(  );
      break;
    case VECTOR_MASK:
      err = Action_MASK(  );
      break;
    case VECTOR_IRED:
      err = Action_IRED(  );
      break;
    case VECTOR_BOX:
      err = Action_BOX(  );
      break;
    default:
      mprinterr("Error: Unhandled vector operation.\n");
      return 1;
  }
  return err;
}

int VectorType::Action_CORR() {
  double CXYZ[3], VXYZ[3];
  double Dplus[2], Dminus[2]; // 0=real, 1=imaginary
  double r3i = 0;
  int indtot;
  //Vec3 CXYZ, VXYZ;

  // Calc snapshot index
  int indsnap = (2*order_ + 1) * frame_;
  if (mode_==VECTOR_CORRIRED)
    indsnap *= modinfo_->Nvect();
  // Calc COM of masks and vector v
  //CXYZ = currentFrame->CenterOfMass(mask_);
  currentFrame->CenterOfMass( &mask_, CXYZ );
  if (mode_==VECTOR_CORR || mode_==VECTOR_CORRIRED) {
    //VXYZ = currentFrame->CenterOfMass(&mask2_);
    currentFrame->CenterOfMass( &mask2_, VXYZ );
    //VXYZ -= CXYZ;
    vector_sub( VXYZ, VXYZ, CXYZ );
  } else if (mode_==VECTOR_CORRPLANE) {
    int idx = 0;
    for (AtomMask::const_iterator atom = mask_.begin();
                                  atom != mask_.end(); ++atom) 
    {
      //VXYZ = currentFrame->GetAtomVec3( *atom );
      currentFrame->GetAtomXYZ( VXYZ, *atom );
      //C_[idx] = currentFrame->GetAtomVec3( *atom );
      //C_[idx] -= CXYZ;
      cx_[idx] = VXYZ[0] - CXYZ[0];
      cy_[idx] = VXYZ[1] - CXYZ[1];
      cz_[idx] = VXYZ[2] - CXYZ[2];
      ++idx;
    }
    leastSquaresPlane(idx, cx_, cy_, cz_, VXYZ); 
  }
  // Calc vector length
  double len = sqrt(VXYZ[0]*VXYZ[0] + VXYZ[1]*VXYZ[1] + VXYZ[2]*VXYZ[2]);
  // Update avgcrd, rave, r3iave, r6iave for VECTOR_CORR, VECTOR_CORRPLANE
  if (mode_== VECTOR_CORR || mode_ == VECTOR_CORRPLANE) {
    avgcrd_[0] += VXYZ[0];
    avgcrd_[1] += VXYZ[1];
    avgcrd_[2] += VXYZ[2];
    rave_ += len;
    double r3 = len*len*len;
    r3i = 1.0 / r3;
    r3iave_ += r3i;
    r6iave_ += r3i*r3i;
  }

  // Loop over m=0, ..., +L
  int indminus = 0;
  for (int idx = 0; idx <= order_; ++idx) {
    // Calc Spherical Harmonics
    sphericalHarmonics(order_, idx, VXYZ, len, Dplus);
    int indplus = order_ + idx; 
    if (idx > 0) {
      sphericalHarmonics(order_, -idx, VXYZ, len, Dminus);
      indminus = order_ - idx;
    }
    // CORRIRED
    // NOTE: For CORRIRED, indsnap (above), indplus, and indminus are 
    //       multiplied by nvect.
    if (mode_ == VECTOR_CORRIRED) {
      indplus *= modinfo_->Nvect();
      indminus *= modinfo_->Nvect();
      // Loop over all eigenvectors
      for (int veci = 0; veci < modinfo_->Nvect(); ++veci) {
        double Qvec = modinfo_->Evec(veci, npair_);
        indtot = 2 * (indsnap + indplus + veci);
        cftmp_[indtot  ] += (Qvec * Dplus[0]);
        cftmp_[indtot+1] += (Qvec * Dplus[1]);
        if (idx > 0) {
          indtot = 2 * (indsnap + indminus + veci);
          cftmp_[indtot  ] += Qvec * Dminus[0];
          cftmp_[indtot+1] += Qvec * Dminus[1];
        }
      } 

      // CORR / CORRPLANE
    } else if (mode_ == VECTOR_CORR || mode_ == VECTOR_CORRPLANE ) {
      indtot = 2 * (indsnap + indplus);
      cftmp_[indtot  ] += r3i * Dplus[0];
      cftmp_[indtot+1] += r3i * Dplus[1];
      p2cftmp_[indtot  ] += Dplus[0];
      p2cftmp_[indtot+1] += Dplus[1];
      if(idx > 0){
        indtot = 2 * (indsnap + indminus);
        cftmp_[indtot  ] += r3i * Dminus[0];
        cftmp_[indtot+1] += r3i * Dminus[1];
        p2cftmp_[indtot  ] += Dminus[0];
        p2cftmp_[indtot+1] += Dminus[1];
      }
      else if(idx == 0){
        indtot = 2 * frame_;
        rcftmp_[indtot  ] += r3i;
        rcftmp_[indtot+1]  = 0.0;
      }
    }
  } // END loop over m=0 ... +L
  ++frame_;
  return 0;
}

int VectorType::Action_DIPOLE()
{
  double Vec[3], CXYZ[3], VXYZ[3], XYZ[3];
  CXYZ[0] = 0;
  CXYZ[1] = 0;
  CXYZ[2] = 0;
  VXYZ[0] = 0;
  VXYZ[1] = 0;
  VXYZ[2] = 0;
  double total_mass = 0;
  for (AtomMask::const_iterator atom = mask_.begin(); 
                                atom != mask_.end(); ++atom)
  {
    currentFrame->GetAtomXYZ( XYZ, *atom );
    double mass = (*currentParm)[*atom].Mass();
    total_mass += mass;
    vector_mult_scalar( Vec, XYZ, mass );
    vector_sum( CXYZ, CXYZ, Vec );

    vector_mult_scalar( Vec, XYZ, (*currentParm)[*atom].Charge() );
    vector_sum( VXYZ, VXYZ, Vec );
  }
  ++frame_;

  return 0;
}

int VectorType::Action_PRINCIPAL( ) {
  double Inertia[9], CXYZ[3], Evec[9], Eval[3];

  currentFrame->CalculateInertia( mask_, Inertia, CXYZ );
  printMatrix_3x3("PRINCIPAL Inertia", Inertia);
  // DEBUG - Write PDB of inertia
/*  double Temp[9];
  for (int i = 0; i < 9; ++i)
    Temp[i] = 0;//Inertia[i];
  Temp[0] = 1;
  Temp[4] = 1;
  Temp[8] = 1;
  //normalize(Temp);
  //normalize(Temp+3);
  //normalize(Temp+6);
  debuginert.Printf("MODEL %i\n",frame_+1);
  PDB.pdb_write_ATOM(debuginert.IO, PDBfile::PDBATOM, 1, (char*)"Orig", (char*)"Vec", ' ', 1,
                    CXYZ[0], CXYZ[1], CXYZ[2]);
  PDB.pdb_write_ATOM(debuginert.IO, PDBfile::PDBATOM, 2, (char*)"X", (char*)"Vec", ' ', 1,
                    Temp[0]+CXYZ[0], Temp[1]+CXYZ[1], Temp[2]+CXYZ[2]);
  PDB.pdb_write_ATOM(debuginert.IO, PDBfile::PDBATOM, 3, (char*)"Y", (char*)"Vec", ' ', 1,
                    Temp[3]+CXYZ[0], Temp[4]+CXYZ[1], Temp[5]+CXYZ[2]);
  PDB.pdb_write_ATOM(debuginert.IO, PDBfile::PDBATOM, 4, (char*)"Z", (char*)"Vec", ' ', 1,
                    Temp[6]+CXYZ[0], Temp[7]+CXYZ[1], Temp[8]+CXYZ[2]);
  debuginert.Printf("ENDMDL\n");*/

  Matrix_3x3 TEMP( Inertia );
  // NOTE: Diagonalize_Sort_Chirality places sorted eigenvectors in rows.
  TEMP.Diagonalize_Sort_Chirality( Evec, Eval );
  printVector("PRINCIPAL EIGENVALUES", Eval );
  //TEMP.Print("GENERAL");
  printMatrix_3x3("PRINCIPAL EIGENVECTORS (Rows)", Evec);

  //Principal_.Diagonalize( Inertia, Eval );

  // V0x V0y V0z
  // V1x V1y V1z
  // V2x V2y V2z

  if (mode_==VECTOR_PRINCIPAL_X) {
    vx_[frame_] = Evec[0];
    vy_[frame_] = Evec[1];
    vz_[frame_] = Evec[2];
  } else if (mode_==VECTOR_PRINCIPAL_Y) {
    vx_[frame_] = Evec[3];
    vy_[frame_] = Evec[4];
    vz_[frame_] = Evec[5];
  } else if (mode_==VECTOR_PRINCIPAL_Z) {
    vx_[frame_] = Evec[6];
    vy_[frame_] = Evec[7];
    vz_[frame_] = Evec[8];
  }
  cx_[frame_] = CXYZ[0];
  cy_[frame_] = CXYZ[1];
  cz_[frame_] = CXYZ[2];

  // DEBUG - Write PDB of axes
  debugpdb.Printf("MODEL %i\n",frame_+1);
  PDB.pdb_write_ATOM(debugpdb.IO, PDBfile::PDBATOM, 1, (char*)"Orig", (char*)"Vec", ' ', 1, 
                    CXYZ[0], CXYZ[1], CXYZ[2]);
  PDB.pdb_write_ATOM(debugpdb.IO, PDBfile::PDBATOM, 2, (char*)"X", (char*)"Vec", ' ', 1,    
                    Evec[6]+CXYZ[0], Evec[7]+CXYZ[1], Evec[8]+CXYZ[2]);
  PDB.pdb_write_ATOM(debugpdb.IO, PDBfile::PDBATOM, 3, (char*)"Y", (char*)"Vec", ' ', 1,   
                    Evec[3]+CXYZ[0], Evec[4]+CXYZ[1], Evec[5]+CXYZ[2]);
  PDB.pdb_write_ATOM(debugpdb.IO, PDBfile::PDBATOM, 4, (char*)"Z", (char*)"Vec", ' ', 1,   
                    Evec[0]+CXYZ[0], Evec[1]+CXYZ[1], Evec[2]+CXYZ[2]);
  debugpdb.Printf("ENDMDL\n");

  ++frame_;

  return 0;
}

int VectorType::Action_MASK( ) {
  double CXYZ[3], VXYZ[3];

  currentFrame->CenterOfMass( &mask_, CXYZ);
  currentFrame->CenterOfMass( &mask2_, VXYZ);
  vx_[frame_] = VXYZ[0] - CXYZ[0];
  vy_[frame_] = VXYZ[1] - CXYZ[1];
  vz_[frame_] = VXYZ[2] - CXYZ[2];
  cx_[frame_] = CXYZ[0];
  cy_[frame_] = CXYZ[1];
  cz_[frame_] = CXYZ[2];
  ++frame_;
  return 0;
}

int VectorType::Action_IRED(  ) {
  double CXYZ[3], VXYZ[3];

  currentFrame->CenterOfMass( &mask_, CXYZ);
  currentFrame->CenterOfMass( &mask2_, VXYZ);
  vx_[frame_] = VXYZ[0] - CXYZ[0];
  vy_[frame_] = VXYZ[1] - CXYZ[1];
  vz_[frame_] = VXYZ[2] - CXYZ[2];
  cx_[frame_] = CXYZ[0];
  cy_[frame_] = CXYZ[1];
  cz_[frame_] = CXYZ[2];
  // IRED only ever has 1 frame, dont increment counter
  return 0;
}

int VectorType::Action_BOX(  ) {
  double XYZ[3];
  currentFrame->BoxXYZ( XYZ );
  vx_[frame_] = XYZ[0];
  vy_[frame_] = XYZ[1];
  vz_[frame_] = XYZ[2];
  cx_[frame_] = 0;
  cy_[frame_] = 0;
  cz_[frame_] = 0;
  ++frame_;
  return 0;
}
      
// -----------------------------------------------------------------------------
/** Solves a cubic equation
  * ax^3 + bx^2 + cx + d = 0
  * using "Cardan's formula"
  * (see: Bronstein, S.131f)
  */
static double solve_cubic_eq(double a, double b, double c, double d)
{

  //const double PI = 3.141592654;
  const double one3 = 1.0 / 3.0;
  const double one27 = 1.0 / 27.0;
  double droot;

  double r, s, t;
  double p, q, rho, phi;
  double D, u, v;
  std::vector<double> dtmp(3);

  /* Coeff. for normal form x^3 + rx^2 + sx + t = 0 */
  r = b / a;
  s = c / a;
  t = d / a;

  /* Coeff. for red. eq. y^3 + py + q = 0 with y = x + r/3 bzw. (x = y - r/3) */
  p = s - r * r * one3;
  q = 2.0 * r * r * r * one27 - r * s * one3 + t;

  /* Dummy variables */
  rho = sqrt(-p * p * p * one27);
  phi = acos(-q / (2.0 * rho));

  /* Discriminante(?) */
  D = pow((p * one3),3) + q * q * 0.25;

  if(D > 0){ /* x real -> one real solution */
    u = pow(-q * 0.5 + sqrt(D), one3);
    v = -p / u * one3;
    droot = (u + v) - r * one3;
  } else if(D <= 0){ 
  /* three real solutions (d < 0) | one real solution + one real double solution or 
                                                     one real triple solution (d = 0) */
    dtmp[0] = 2.0 * pow(rho, one3) * cos(phi * one3) - r * one3;
    dtmp[1] = 2.0 * pow(rho, one3) * cos((phi + 2.0 * PI) * one3) - r * one3;
    dtmp[2] = 2.0 * pow(rho, one3) * cos((phi + 4.0 * PI) * one3) - r * one3;

    sort(dtmp.begin(), dtmp.end());

    //qsort((void *) dtmp, (size_t) 3, sizeof(double), cmpdouble);
    droot = dtmp[0];
  }
  return droot;
}

/** Calcs (least-squares best) plane through a series of points
  * relative to their center of geom. (the latter has to be done outside this routine), 
  * returns (normalized) coeff. for plane eq. ax + by + cz = 0
  * following: Crystal Structure Analysis for Chem. and Biol.,
  * Glusker, Lewis, Rossi, S. 460ff
  */
void VectorType::leastSquaresPlane(int n, double* cx, double* cy, double* cz, double* XYZ) 
//void VectorType::leastSquaresPlane(std::vector<Vec3>& Cvec, Vec3& XYZ) 
{
  int i;
  double dSumXX, dSumYY, dSumZZ, dSumXY, dSumXZ, dSumYZ;
  double o, p, q, root;
  double dnorm;
  double x1, y1, z1, x2, y2, z2;

  root = 0;

  if (n == 3) {
    /*Vec3 V1 = Cvec[1];
    V1 -= Cvec[0];
    Vec3 V2 = Cvec[2];
    V2 -= Cvec[1];
    XYZ.CROSS( V1, V2 );*/

    x1 = cx[1] - cx[0];
    y1 = cy[1] - cy[0];
    z1 = cz[1] - cz[0];
    x2 = cx[2] - cx[1];
    y2 = cy[2] - cy[1];
    z2 = cz[2] - cz[1];

    XYZ[0] = y1 * z2 - z1 * y2;
    XYZ[1] = z1 * x2 - x1 * z2;
    XYZ[2] = x1 * y2 - y1 * x2;
  }
  else{
    /* Calc Var. */
    dSumXX = 0.0;
    dSumYY = 0.0;
    dSumZZ = 0.0;
    dSumXY = 0.0;
    dSumXZ = 0.0;
    dSumYZ = 0.0;

    for(i = 0; i < n; i++){
      dSumXX += cx[i] * cx[i];
      dSumYY += cy[i] * cy[i];
      dSumZZ += cz[i] * cz[i];

      dSumXY += cx[i] * cy[i];
      dSumXZ += cx[i] * cz[i];
      dSumYZ += cy[i] * cz[i];
    }

    /* Calc coeff. for -l^3 + o * l^2 + p * l + q = 0 */
    o = dSumXX + dSumYY + dSumZZ;
    p = pow(dSumXY,2) + pow(dSumXZ,2) + pow(dSumYZ,2) -
        (dSumXX * dSumYY + dSumXX * dSumZZ + dSumYY * dSumZZ);
    q = dSumXX * dSumYY * dSumZZ + 2.0 * dSumXY * dSumXZ * dSumYZ -
      (dSumXX * dSumYZ * dSumYZ + dSumYY * dSumXZ * dSumXZ + dSumZZ * dSumXY * dSumXY);

    /* Solve cubic eq. */
    root = solve_cubic_eq(-1.0, o, p, q);

    /* Calc determinantes */
    XYZ[0] = (dSumYY - root) * dSumXZ - dSumXY * dSumYZ;
    XYZ[1] = (dSumXX - root) * dSumYZ - dSumXY * dSumXZ;
    XYZ[2] =  dSumXY         * dSumXY - (dSumYY - root) * (dSumXX - root);

  }
  /* Normalize */
  dnorm = 1.0 / sqrt((XYZ[0]) * (XYZ[0]) + (XYZ[1]) * (XYZ[1]) + (XYZ[2]) * (XYZ[2]));
  XYZ[0] *= dnorm;
  XYZ[1] *= dnorm;
  XYZ[2] *= dnorm;
}

/** Calc spherical harmonics of order l=0,1,2
  * and -l<=m<=l with cartesian coordinates as input
  * (see e.g. Merzbacher, Quantum Mechanics, p. 186)
  */
void VectorType::sphericalHarmonics(int l, int m, double* XYZ, double r, double D[2])
{
  const double SH00=0.28209479;
  const double SH10=0.48860251;
  const double SH11=0.34549415;
  const double SH20=0.31539157;
  const double SH21=0.77254840;
  const double SH22=0.38627420;

  double ri;
  double x = XYZ[0];
  double y = XYZ[1];
  double z = XYZ[2];

  D[0] = 0.0;
  D[1] = 0.0;
  ri = 1.0 / r;

  if(l == 0 && m == 0){
    D[0] = SH00;
  }
  else if(l == 1){
    if(m == 0){
      D[0] = SH10 * z * ri;
    }
    else{
      D[0] = -m * SH11 * x * ri;
      D[1]  = -    SH11 * y * ri;
    }
  }
  else if(l == 2){
    if(m == 0){
      D[0] = SH20 * (2.0*z*z - x*x - y*y) * ri * ri;
    }
    else if(fabs(m) == 1){
      D[0] = -m * SH21 * x * z * ri * ri;
      D[1]  = -    SH21 * y * z * ri * ri;
    }
    else{
      D[0] = SH22 * (x*x - y*y) * ri * ri;
      D[1]  = m * SH22 * x * y * ri * ri;
    }
  }
}


