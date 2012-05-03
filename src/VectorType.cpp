#include <cmath>
#include <algorithm> // sort
#include "VectorType.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"
#include "Constants.h" // PI

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
  avgcrd_(0),
  rave_(0),
  r3iave_(0),
  r6iave_(0),
  cftmp_(0),
  p2cftmp_(0),
  rcftmp_(0)
{
  dType_ = VECTOR;
}

// DESTRUCTOR
VectorType::~VectorType() {
  // Free shared CORRIRED mem on master
  if (mode_==VECTOR_CORRIRED) {
    if (master_==0) { 
      if (modinfo_!=0) delete modinfo_;
      if (cftmp_!=0) delete cftmp_;
    }
  } else {
    if (cftmp_!=0) delete cftmp_;
    if (avgcrd_!=0) delete avgcrd_;
    if (p2cftmp_!=0) delete p2cftmp_;
    if (rcftmp_!=0) delete rcftmp_;
    if (cx_!=0) delete cx_;
    if (cy_!=0) delete cy_;
    if (cz_!=0) delete cz_;
    if (vx_!=0) delete vx_;
    if (vy_!=0) delete vy_;
    if (vz_!=0) delete vz_;
  }
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

// VectorType::Init()
int VectorType::Init(ArgList& argIn) {
  filename_ = argIn.GetStringKey("out");

  // Require a vector name - this behavior is consistent with ptraj
  name_ = argIn.GetStringNext();

  // Get order for Legendre polynomial
  order_ = argIn.getKeyInt("order",2);
  if (order_ < 0 || order_ > 2) {
    mprintf("Warning: vector order out of bounds (<0 or >2), resetting to 2.\n");
    order_ = 2;
  }

  // Determine vector mode.
  // Acceptable args: principal, principal x, principal y, principal z,
  //                  dipole, box, corrplane, corrired, corr, ired
  int pos = argIn.KeyPosition("principal");
  if (pos == argIn.Nargs() - 1) {
    argIn.MarkArg(pos);
    mode_ = VECTOR_PRINCIPAL_X;
  } else if (pos!=-1) {
    argIn.MarkArg(pos);
    mode_ = VECTOR_PRINCIPAL_X;
    char vecchar = argIn[pos+1][0];
    if (vecchar == 'x' ||
        vecchar == 'y' ||
        vecchar == 'z')
    {
      argIn.MarkArg(pos+1);
      switch (vecchar) {
        case 'x': mode_ = VECTOR_PRINCIPAL_X; break;
        case 'y': mode_ = VECTOR_PRINCIPAL_Y; break;
        case 'z': mode_ = VECTOR_PRINCIPAL_Z; break;
      }
    }
  } else if (argIn.hasKey("dipole"))
    mode_ = VECTOR_DIPOLE;
  else if (argIn.hasKey("box"))
    mode_ = VECTOR_BOX;
  else if (argIn.hasKey("corrplane"))
    mode_ = VECTOR_CORRPLANE;
  else if (argIn.hasKey("corrired"))
    mode_ = VECTOR_CORRIRED;
  else if (argIn.hasKey("corr"))
    mode_ = VECTOR_CORR;
  else if (argIn.hasKey("ired"))
    mode_ = VECTOR_IRED;
  else
    mode_ = VECTOR_MASK;

  // VECTOR_CORRIRED
  if (mode_ == VECTOR_CORRIRED) {
    // Get Pair number
    npair_ = argIn.getKeyInt("npair",0);
    if (npair_ == 0) {
      mprinterr("Error: VectorType::Init(): No 'npair <#>' arg given, needed for 'corrired'.\n");
      return 1;
    }
    modesfile_ = argIn.GetStringKey("modes");
    if (modesfile_.empty()) {
      mprinterr("Error: VectorType::Init(): No 'modes' <file> arg given, needed for 'corrired'.\n");
      return 1;
    }
    ibeg_ = argIn.getKeyInt("beg",1);
    iend_ = argIn.getKeyInt("end", 50);
  }

  // Vector Mask
  char *maskexpr = argIn.getNextMask();
  mask_.SetMaskString( maskexpr );

  // Get second mask if necessary
  if (mode_ == VECTOR_IRED ||
      mode_ == VECTOR_CORR ||
      mode_ == VECTOR_CORRIRED ||
      mode_ == VECTOR_MASK)
  {
    maskexpr = argIn.getNextMask();
    if (maskexpr==NULL) {
      mprinterr("Error: vector: Specified vector mode requires a second mask.\n");
      return 1;
    }
    mask2_.SetMaskString( maskexpr );
  }

  return 0;
}

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
    avgcrd_ = new double[ 3 ];
    avgcrd_[0] = 0;
    avgcrd_[1] = 0;
    avgcrd_[2] = 0;
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
      mprintf(" will be calcd with respect to the center of mass of atoms in\n");
      mprintf("\tmask expression [%s]\n", mask_.MaskString());
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
    mprintf("\tand the pair %i is considered\n", npair_);
  }

  if (mode_ == VECTOR_CORRPLANE ||
      mode_ == VECTOR_CORR ||
      mode_ == VECTOR_CORRIRED ||
      mode_ == VECTOR_IRED) 
  {
    if (!filename_.empty()) {
      mprintf("\tWarning: Output of corr, ired, corrired or corrplane vectors is not yet supported!\n");
      filename_.clear(); 
    }
  }

  if (!filename_.empty()) {
    mprintf("\tOutput will be dumped to a file, %s\n", filename_.c_str());
  }
}

// VectorType::Print()
void VectorType::Print() {
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
            vx_[i], vy_[i], vz_[i],
            cx_[i], cy_[i], cz_[i],
            cx_[i]+vx_[i], cy_[i]+vy_[i], cz_[i]+vz_[i]);
  }
  outfile.CloseFile();
}

int VectorType::Setup(Topology* currentParm) {
  // Setup mask 1
  if (currentParm->SetupIntegerMask(mask_)) return 1;

  // Allocate space for CORRPLANE
  if (mode_==VECTOR_CORRPLANE) {
    if (cx_!=0 || cy_!=0 || cz_!=0) {
      mprinterr("Error: VECTOR_CORRPLANE not yet supported for multiple topology files.\n");
      return 1;
    }
    cx_ = new double[ mask_.Nselected() ];
    cy_ = new double[ mask_.Nselected() ];
    cz_ = new double[ mask_.Nselected() ];
  }

  // Setup mask 2
  if (!mask2_.NoMaskString())
    if (currentParm->SetupIntegerMask(mask2_)) return 1;
  return 0;
} 

int VectorType::Action(Frame* currentFrame) {
  double CXYZ[3], VXYZ[3], len, r3, r3i;
  int indsnap, idx;

  switch (mode_) {
    case VECTOR_CORRIRED:
    case VECTOR_CORR:
    case VECTOR_CORRPLANE:
      indsnap = (2*order_ + 1) * frame_;
      // Calc COM of masks and vector v
      currentFrame->CenterOfMass(&mask_, CXYZ);
      if (mode_==VECTOR_CORR || mode_==VECTOR_CORRIRED)
        currentFrame->CenterOfMass(&mask2_, VXYZ);
      else if (mode_==VECTOR_CORRPLANE) {
        idx = 0;
        for (AtomMask::const_iterator atom = mask_.begin();
                                      atom != mask_.end(); ++atom) 
        {
          currentFrame->GetAtomXYZ( VXYZ, *atom );
          cx_[idx] = VXYZ[0] - CXYZ[0];
          cy_[idx] = VXYZ[1] - CXYZ[1];
          cz_[idx] = VXYZ[2] - CXYZ[2];
          ++idx;
        }
        lsqplane(idx, cx_, cy_, cz_, VXYZ); 
      }
      // Calc vector length
      len = sqrt(VXYZ[0]*VXYZ[0] + VXYZ[1]*VXYZ[1] + VXYZ[2]*VXYZ[2]);
      // Update avgcrd, rave, r3iave, r6iave for VECTOR_CORR, VECTOR_CORRPLANE
      if (mode_== VECTOR_CORR || mode_ == VECTOR_CORRPLANE) {
        avgcrd_[0] += VXYZ[0];
        avgcrd_[1] += VXYZ[1];
        avgcrd_[2] += VXYZ[2];
        rave_ += len;
        r3 = len*len*len;
        r3i = 1.0 / r3;
        r3iave_ += r3i;
        r6iave_ += r3i*r3i;
      }

      break;

    default:
      mprinterr("Error: Unhandled vector operation.\n");
      return 1;
  }

  return 0;
}

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
void VectorType::lsqplane(int n, double *cx, double *cy, double *cz, double *XYZ) {
  int i;
  double dSumXX, dSumYY, dSumZZ, dSumXY, dSumXZ, dSumYZ;
  double o, p, q, root;
  double dnorm;
  double x1, y1, z1, x2, y2, z2;

  root = 0;

  if(n == 3){
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


