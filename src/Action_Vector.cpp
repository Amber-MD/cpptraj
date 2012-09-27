#include "Action_Vector.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Vector::Action_Vector() :
  Vec_(0),
  totalFrames_(-1),
  frame_(0)
{ }

// Action_Vector::init()
int Action_Vector::init() {
  filename_ = actionArgs.GetStringKey("out");

  // Require a vector name - this behavior is consistent with ptraj
  std::string name = actionArgs.GetStringNext();

  // Set up vector dataset
  Vec_ = (DataSet_Vector*)DSL->AddSet(DataSet::VECTOR, name, "Vec");

  // Get order for Legendre polynomial and determine vector mode
  if (Vec_->SetModeFromArgs( actionArgs )) return 1;

  // VECTOR_CORRIRED
  if (Vec_->Mode() == DataSet_Vector::VECTOR_CORRIRED) {
    // check if modes need to be read or have been
    // previously read in by another DataSet_Vector.
    DataSet_Vector* Vmaster = 0;
    DataSet_Vector* Vtmp;
    DSL->VectorBegin();
    while ( (Vtmp = (DataSet_Vector*)DSL->NextVector()) != 0 ) {
      if ( *Vec_ == *Vtmp ) {
        Vmaster = Vtmp;
        break;
      }
    } 
    if ( Vec_->SetupCorrIred( Vmaster, actionArgs ) ) {
      mprinterr("Error: Could setup vector for CORRIRED.\n");
      return 1;
    }
  }

  // Vector Mask
  ArgList::ConstArg maskexpr = actionArgs.getNextMask();
  mask_.SetMaskString( maskexpr );

  // Get second mask if necessary
  if (Vec_->Mode() == DataSet_Vector::VECTOR_IRED ||
      Vec_->Mode() == DataSet_Vector::VECTOR_CORR ||
      Vec_->Mode() == DataSet_Vector::VECTOR_CORRIRED ||
      Vec_->Mode() == DataSet_Vector::VECTOR_MASK)
  {
    maskexpr = actionArgs.getNextMask();
    if (maskexpr==NULL) {
      mprinterr("Error: vector: Specified vector mode requires a second mask.\n");
      return 1;
    }
    mask2_.SetMaskString( maskexpr );
  }

  // Allocate vector
  if (Vec_->Allocate(DSL->MaxFrames())) return 1;

  Vec_->Info(mask_.MaskString(), mask2_.MaskString());

  // Check if output is supported for the current vector mode
  if (Vec_->Mode() == DataSet_Vector::VECTOR_CORRPLANE ||
      Vec_->Mode() == DataSet_Vector::VECTOR_CORR ||
      Vec_->Mode() == DataSet_Vector::VECTOR_CORRIRED ||
      Vec_->Mode() == DataSet_Vector::VECTOR_IRED)
  {
    if (!filename_.empty()) {
      mprintf(
        "\tWarning: Output of corr, ired, corrired or corrplane vectors is not yet supported!\n");
        filename_.clear();
    }
  }

  // Add to output datafilelist
  //if (!filename_.empty()) {
  //  mprintf("\tOutput will be dumped to a file, %s\n", filename_.c_str());
  //  DFL->Add( filename_.c_str(), (DataSet*)this );
  //}

  return 0;
}

// Action_Vector::print()
void Action_Vector::print() {
  if (filename_.empty()) return;

  CpptrajFile outfile;
  if (outfile.OpenWrite(filename_.c_str())) return;

  mprintf("CPPTRAJ VECTOR: dumping vector information %s\n", Vec_->Legend().c_str());
  outfile.Printf("# FORMAT: frame vx vy vz cx cy cz cx+vx cy+vy cz+vz\n");
  outfile.Printf("# FORMAT where v? is vector, c? is center of mass...\n");
  for (int i=0; i < totalFrames_; ++i) {
    Vec3 vxyz = Vec_->VXYZ(i);
    Vec3 cxyz = Vec_->CXYZ(i);
    Vec3 txyz  = cxyz + vxyz;
    outfile.Printf("%i %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
            i+1, 
            //V_[i][0], V_[i][1], V_[i][2], C_[i][0], C_[i][1], C_[i][2],
            //V_[i][0]+C_[i][0], V_[i][1]+C_[i][1], V_[i][2]+C_[i][2]
            vxyz[0], vxyz[1], vxyz[2],
            cxyz[0], cxyz[1], cxyz[2],
            txyz[0], txyz[1], txyz[2]);
  }
  outfile.CloseFile();
}

// Action_Vector::setup()
int Action_Vector::setup() {
  if (Vec_->Mode()==DataSet_Vector::VECTOR_BOX) {
    // Check for box info
    if (currentParm->BoxType() == Box::NOBOX) {
      mprinterr("Error: vector box: Parm %s does not have box information.\n", 
                currentParm->c_str());
      return 1;
    }
  } else {
    // Setup mask 1
    if (currentParm->SetupIntegerMask(mask_)) return 1;
    mprintf("\tVector mask [%s] corresponds to %i atoms.\n",
            mask_.MaskString(), mask_.Nselected());
  }

  // Allocate space for CORRPLANE. Does nothing if not CORRPLANE.
  Vec_->AllocCorrPlane( mask_.Nselected() );

  // Setup mask 2
  if (mask2_.MaskStringSet()) {
    if (currentParm->SetupIntegerMask(mask2_)) return 1;
    mprintf("\tVector mask [%s] corresponds to %i atoms.\n",
            mask2_.MaskString(), mask2_.Nselected());
  }
  return 0;
} 

// Action_Vector::action()
int Action_Vector::action() {
  int err;

  switch (Vec_->Mode()) {
    case DataSet_Vector::VECTOR_CORRIRED:
    case DataSet_Vector::VECTOR_CORR:
    case DataSet_Vector::VECTOR_CORRPLANE:
      err = Action_CORR(  );
      break;
    case DataSet_Vector::VECTOR_DIPOLE:
      err = Action_DIPOLE(  );
      break;
    case DataSet_Vector::VECTOR_PRINCIPAL_X:
    case DataSet_Vector::VECTOR_PRINCIPAL_Y:
    case DataSet_Vector::VECTOR_PRINCIPAL_Z:
      err = Action_PRINCIPAL(  );
      break;
    case DataSet_Vector::VECTOR_MASK:
      err = Action_MASK(  );
      break;
    case DataSet_Vector::VECTOR_IRED:
      err = Action_IRED(  );
      break;
    case DataSet_Vector::VECTOR_BOX:
      err = Action_BOX(  );
      break;
    default:
      mprinterr("Error: Unhandled vector operation.\n");
      return 1;
  }
  return err;
}

// Action_Vector::Action_CORR()
int Action_Vector::Action_CORR() {
  double Dplus[2], Dminus[2]; // 0=real, 1=imaginary
  double r3i = 0;
  int indtot;
  Vec3 CXYZ, VXYZ;

  // Calc snapshot index
  int indsnap = (2*order_ + 1) * frame_;
  if (mode_==VECTOR_CORRIRED)
    indsnap *= modinfo_->Nvect();
  // Calc COM of masks and vector v
  CXYZ = currentFrame->VCenterOfMass(mask_);
  if (mode_==VECTOR_CORR || mode_==VECTOR_CORRIRED) {
    VXYZ = currentFrame->VCenterOfMass(mask2_);
    VXYZ -= CXYZ;
  } else if (mode_==VECTOR_CORRPLANE) {
    int idx = 0;
    for (AtomMask::const_iterator atom = mask_.begin();
                                  atom != mask_.end(); ++atom) 
    {
      VXYZ = currentFrame->XYZ( *atom );
      VXYZ -= CXYZ;
      cx_[idx] = VXYZ[0];
      cy_[idx] = VXYZ[1];
      cz_[idx] = VXYZ[2];
      ++idx;
    }
    VXYZ = leastSquaresPlane(idx, cx_, cy_, cz_); 
  }
  // Calc vector length
  double len = sqrt( VXYZ.Magnitude2() );
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
  Dminus[0] = 0; // DBG
  Dminus[1] = 0; // DBG
  for (int idx = 0; idx <= order_; ++idx) {
    // Calc Spherical Harmonics
    sphericalHarmonics(order_, idx, VXYZ.Dptr(), len, Dplus);
    int indplus = order_ + idx; 
    if (idx > 0) {
      sphericalHarmonics(order_, -idx, VXYZ.Dptr(), len, Dminus);
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
        //mprintf("CDBG: Q=%.10lE\n",Qvec);
        //mprintf("CDBG: Dplus=%lE,%lE\n",Dplus[0],Dplus[1]);
        //mprintf("CDBG: Dminus=%lE,%lE\n",Dminus[0],Dminus[1]);
        indtot = 2 * (indsnap + indplus + veci);
        cftmp_[indtot  ] += (Qvec * Dplus[0]);
        cftmp_[indtot+1] += (Qvec * Dplus[1]);
        if (idx > 0) {
          indtot = 2 * (indsnap + indminus + veci);
          cftmp_[indtot  ] += (Qvec * Dminus[0]);
          cftmp_[indtot+1] += (Qvec * Dminus[1]);
        }
        //mprintf("CDBG: cftmp[%i]=%lf cftmp[%i]=%lf\n",indtot,cftmp_[indtot],indtot+1,cftmp_[indtot+1]);
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

// Action_Vector::Action_DIPOLE()
int Action_Vector::Action_DIPOLE()
{
  double CXYZ[3], VXYZ[3];
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
    const double* XYZ = currentFrame->XYZ( *atom );
    double mass = (*currentParm)[*atom].Mass();
    total_mass += mass;
    CXYZ[0] += (mass * XYZ[0]);
    CXYZ[1] += (mass * XYZ[1]);
    CXYZ[2] += (mass * XYZ[2]);

    double charge = (*currentParm)[*atom].Charge();
    VXYZ[0] += (charge * XYZ[0]);
    VXYZ[1] += (charge * XYZ[1]);
    VXYZ[2] += (charge * XYZ[2]);
  }
  vx_[frame_] = VXYZ[0];
  vy_[frame_] = VXYZ[1];
  vz_[frame_] = VXYZ[2];
  cx_[frame_] = CXYZ[0] / total_mass;
  cy_[frame_] = CXYZ[1] / total_mass;
  cz_[frame_] = CXYZ[2] / total_mass;
  ++frame_;

  return 0;
}

// Action_Vector::Action_PRINCIPAL()
int Action_Vector::Action_PRINCIPAL( ) {
  double Inertia[9], CXYZ[3], Evec[9], Eval[3];

  currentFrame->CalculateInertia( mask_, Inertia, CXYZ );
  //printMatrix_3x3("PRINCIPAL Inertia", Inertia);
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
  TEMP.Diagonalize_Sort_Chirality( Evec, Eval, debug );
  if (debug > 2) {
    printVector("PRINCIPAL EIGENVALUES", Eval );
    //TEMP.Print("GENERAL");
    printMatrix_3x3("PRINCIPAL EIGENVECTORS (Rows)", Evec);
  }

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
  /*debugpdb.Printf("MODEL %i\n",frame_+1);
  PDB.pdb_write_ATOM(debugpdb.IO, PDBfile::PDBATOM, 1, (char*)"Orig", (char*)"Vec", ' ', 1, 
                    CXYZ[0], CXYZ[1], CXYZ[2]);
  PDB.pdb_write_ATOM(debugpdb.IO, PDBfile::PDBATOM, 2, (char*)"X", (char*)"Vec", ' ', 1,    
                    Evec[6]+CXYZ[0], Evec[7]+CXYZ[1], Evec[8]+CXYZ[2]);
  PDB.pdb_write_ATOM(debugpdb.IO, PDBfile::PDBATOM, 3, (char*)"Y", (char*)"Vec", ' ', 1,   
                    Evec[3]+CXYZ[0], Evec[4]+CXYZ[1], Evec[5]+CXYZ[2]);
  PDB.pdb_write_ATOM(debugpdb.IO, PDBfile::PDBATOM, 4, (char*)"Z", (char*)"Vec", ' ', 1,   
                    Evec[0]+CXYZ[0], Evec[1]+CXYZ[1], Evec[2]+CXYZ[2]);
  debugpdb.Printf("ENDMDL\n");*/

  ++frame_;

  return 0;
}

// Action_Vector::Action_MASK()
int Action_Vector::Action_MASK( ) {
  double CXYZ[3], VXYZ[3];

  currentFrame->CenterOfMass( CXYZ, mask_ );
  currentFrame->CenterOfMass( VXYZ, mask2_ );
  vx_[frame_] = VXYZ[0] - CXYZ[0];
  vy_[frame_] = VXYZ[1] - CXYZ[1];
  vz_[frame_] = VXYZ[2] - CXYZ[2];
  cx_[frame_] = CXYZ[0];
  cy_[frame_] = CXYZ[1];
  cz_[frame_] = CXYZ[2];
  ++frame_;
  return 0;
}

// Action_Vector::Action_IRED()
int Action_Vector::Action_IRED(  ) {
  double CXYZ[3], VXYZ[3];

  currentFrame->CenterOfMass( CXYZ, mask_ );
  currentFrame->CenterOfMass( VXYZ, mask2_ );
  vx_[frame_] = VXYZ[0] - CXYZ[0];
  vy_[frame_] = VXYZ[1] - CXYZ[1];
  vz_[frame_] = VXYZ[2] - CXYZ[2];
  cx_[frame_] = CXYZ[0];
  cy_[frame_] = CXYZ[1];
  cz_[frame_] = CXYZ[2];
  // IRED only ever has 1 frame, dont increment counter
  return 0;
}

// Action_Vector::Action_BOX()
int Action_Vector::Action_BOX(  ) {
  vx_[frame_] = currentFrame->BoxX();
  vy_[frame_] = currentFrame->BoxY();
  vz_[frame_] = currentFrame->BoxZ();
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
  double droot = 0;

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
Vec3 Action_Vector::leastSquaresPlane(int n, double* cx, double* cy, double* cz) 
//void Action_Vector::leastSquaresPlane(std::vector<Vec3>& Cvec, Vec3& XYZ) 
{
  int i;
  double dSumXX, dSumYY, dSumZZ, dSumXY, dSumXZ, dSumYZ;
  double o, p, q;
  double dnorm, Xout, Yout, Zout;
  double x1, y1, z1, x2, y2, z2;

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

    Xout = y1 * z2 - z1 * y2;
    Yout = z1 * x2 - x1 * z2;
    Zout = x1 * y2 - y1 * x2;
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
    double root = solve_cubic_eq(-1.0, o, p, q);

    /* Calc determinantes */
    Xout = (dSumYY - root) * dSumXZ - dSumXY * dSumYZ;
    Yout = (dSumXX - root) * dSumYZ - dSumXY * dSumXZ;
    Zout =  dSumXY         * dSumXY - (dSumYY - root) * (dSumXX - root);

  }
  /* Normalize */
  dnorm = 1.0 / sqrt((Xout * Xout) + (Yout) * (Yout) + (Zout) * (Zout));
  Xout *= dnorm;
  Yout *= dnorm;
  Zout *= dnorm;
  return Vec3(Xout, Yout, Zout);
}

/** Calc spherical harmonics of order l=0,1,2
  * and -l<=m<=l with cartesian coordinates as input
  * (see e.g. Merzbacher, Quantum Mechanics, p. 186)
  * D[0] is the real part, D[1] is the imaginary part.
  */
void Action_Vector::sphericalHarmonics(int l, int m, const double* XYZ, double r, double D[2])
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
  //mprintf("CDBG: x=%.10lE y=%.10lE z=%.10lE\n",x,y,z);

  D[0] = 0.0;
  D[1] = 0.0;
  ri = 1.0 / r;

  if (l == 0 && m == 0) {
    D[0] = SH00; 
  } else if (l == 1) {
    if (m == 0) {
      D[0] = SH10 * z * ri;
    } else { 
      D[0] = -m * SH11 * x * ri; 
      D[1] =     -SH11 * y * ri;
    }
  } else if(l == 2) {
    if (m == 0) {
      D[0] = SH20 * (2.0*z*z - x*x - y*y) * ri * ri;
    } else if (abs(m) == 1) {
      D[0] = -m * SH21 * x * z * ri * ri;
      D[1] =     -SH21 * y * z * ri * ri;
    } else {
      D[0] = SH22 * (x*x - y*y) * ri * ri;
      D[1] = m * SH22 * x * y * ri * ri;
    }
  }
  //mprintf("CDBG: dreal=%.10lE dimg=%.10lE\n",D[0],D[1]);
}


