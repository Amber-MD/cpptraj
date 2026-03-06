#include <cmath> //sqrt
#include "Action_Vector.h"
#include "CharMask.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h" // MinImagedVec, includes Matrix_3x3 for principal
#include "DataSet_Vector.h"
#include "DataSet_3D.h"

// CONSTRUCTOR
Action_Vector::Action_Vector() :
  Vec_(0),
  Magnitude_(0),
  gridSet_(0),
  ptrajoutput_(false),
  needBoxInfo_(false),
  useMass_(true),
  dipole_in_debye_(false),
  CurrentParm_(0),
  outfile_(0),
  cmask_(0),
  debug_(0)
{}

// DESTRUCTOR
Action_Vector::~Action_Vector() {
  if (cmask_ != 0) delete cmask_;
}

// Action_Vector::Help()
void Action_Vector::Help() const {
  mprintf("\t[<name>] <Type> [out <filename> [ptrajoutput]] [<mask1>] [<mask2>]\n"
          "\t[magnitude] [geom] [ired] [gridset <grid>] [debye]\n"
          "\t<Type> = { mask     | minimage  | dipole | center   | corrplane | \n"
          "\t           box      | boxcenter | ucellx | ucelly   | ucellz    | \n"
          "\t           momentum | principal [x|y|z]  | velocity | force       }\n" 
          "  Calculate the vector of specified <Type>:\n"
          "    mask             : (Default) Vector from <mask1> to <mask2>.\n"
          "    minimage         : Store the minimum image vector between atoms in <mask1> and <mask2>.\n"
          "    dipole           : Dipole and center of mass of the atoms specified in <mask1>\n"
          "    center           : Store the center of mass of atoms in <mask1>.\n"
          "    corrplane        : Vector perpendicular to plane through the atoms in <mask1>.\n"
          "    box              : (No mask needed) Store the box lengths of the trajectory.\n"
          "    boxcenter        : (No mask needed) Store box center as vector.\n"
          "    ucell{x|y|z}     : (No mask needed) Store specified unit cell vector.\n"
          "    momentum         : Store total momentum vector of atoms in <mask1> (requires velocities).\n"
          "    principal [x|y|z]: X, Y, or Z principal axis vector for atoms in <mask1>.\n"
          "    velocity         : Store velocity of atoms in <mask1> (requires velocities).\n"
          "    force            : Store force of atoms in <mask1> (requires forces).\n"
          "  If 'magnitude' is specified, also calculate the vector magnitude.\n"
          "  If 'geom' is specified, use geometric centers, otherwise use center of mass.\n"
          "  If 'debye' is specified, report 'dipole' vector in Debye instead of e-*Ang.\n");
}

const char* Action_Vector::ModeString_[] = {
  "NO_OP", "Principal X", "Principal Y", "Principal Z",
  "Dipole", "Box", "Mask",
  "CorrPlane", "Center", "Unit cell X", "Unit cell Y", "Unit cell Z",
  "Box Center", "MinImage", "Momentum", "Velocity", "Force",
  "Bond Dipole"
};

const bool Action_Vector::NeedsOrigin_[] = {
  false, true, true, true,
  true, false, true,
  true, false, true, true, true,
  false, true, false, false, false,
  true
};

static Action::RetType WarnDeprecated() {
  mprinterr("Error: Vector: 'corrired' and 'corr' are deprecated.\n"
            "Error: 'corrired' functionality is now part of the\n"
            "Error: IRED analysis. 'corr' can now be done with a normal 2-mask\n"
            "Error: vector and TIMECORR analysis.\n");
  return Action::ERR;
}

static inline Action::RetType DeprecatedErr(const char* key) {
  mprinterr("Error: '%s' is deprecated. To write vector pseudo-traj use:\n"
            "\twritedata <filename> <vector set> [trajfmt <format>] [parmout <file>]\n", key);
  return Action::ERR;
}

// Action_Vector::Init()
Action::RetType Action_Vector::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  DataFile* df = 0;
  std::string filename = actionArgs.GetStringKey("out");
  if (actionArgs.hasKey("trajout")) return DeprecatedErr( "trajout" );
  if (actionArgs.hasKey("trajfmt")) return DeprecatedErr( "trajfmt" );
  if (actionArgs.hasKey("parmout")) return DeprecatedErr( "parmout" );
  ptrajoutput_ = actionArgs.hasKey("ptrajoutput");
  if (ptrajoutput_) {
    if (filename.empty()) {
      mprinterr("Error: 'ptrajoutput' specified but no 'out <filename>' arg given.\n");
      return Action::ERR;
    }
    outfile_ = init.DFL().AddCpptrajFile(filename, "Vector (PTRAJ)");
    if (outfile_ == 0) return Action::ERR;
  } else
    df = init.DFL().AddDataFile( filename, actionArgs );
  bool calc_magnitude = actionArgs.hasKey("magnitude");
  if (calc_magnitude && ptrajoutput_) {
    mprinterr("Error: 'ptrajoutput' and 'magnitude' are incompatible.\n");
    return Action::ERR;
  }
  dipole_in_debye_ = actionArgs.hasKey("debye");
  needBoxInfo_ = false;
  // Deprecated: corrired, corr, ired
  if ( actionArgs.hasKey("principal") ) {
    mode_ = PRINCIPAL_X;
    if ( actionArgs.hasKey("x") ) mode_ = PRINCIPAL_X;
    if ( actionArgs.hasKey("y") ) mode_ = PRINCIPAL_Y;
    if ( actionArgs.hasKey("z") ) mode_ = PRINCIPAL_Z;
  } else if (actionArgs.hasKey("center"))
    mode_ = CENTER;
  else if (actionArgs.hasKey("momentum"))
    mode_ = MOMENTUM;
  else if (actionArgs.hasKey("velocity"))
    mode_ = VELOCITY;
  else if (actionArgs.hasKey("force"))
    mode_ = FORCE;
  else if (actionArgs.hasKey("dipole"))
    mode_ = DIPOLE;
  else if (actionArgs.hasKey("bonddipole"))
    mode_ = BONDDIPOLE;
  else if (actionArgs.hasKey("box"))
    mode_ = BOX;
  else if (actionArgs.hasKey("corrplane"))
    mode_ = CORRPLANE;
  else if (actionArgs.hasKey("corrired"))
    return WarnDeprecated();
  else if (actionArgs.hasKey("corr"))
    return WarnDeprecated();
  else if (actionArgs.hasKey("mask"))
    mode_ = MASK;
  else if (actionArgs.hasKey("ucellx"))
    mode_ = BOX_X;
  else if (actionArgs.hasKey("ucelly"))
    mode_ = BOX_Y;
  else if (actionArgs.hasKey("ucellz"))
    mode_ = BOX_Z;
  else if (actionArgs.hasKey("boxcenter"))
    mode_ = BOX_CTR;
  else if (actionArgs.hasKey("minimage"))
    mode_ = MINIMAGE; 
  else
    mode_ = MASK;
  if (mode_ == BOX || mode_ == BOX_X || mode_ == BOX_Y || mode_ == BOX_Z ||
      mode_ == BOX_CTR || mode_ == MINIMAGE)
    needBoxInfo_ = true;
  gridSet_ = 0;
  if (needBoxInfo_) {
    std::string gridSetArg = actionArgs.GetStringKey("gridset");
    if (!gridSetArg.empty()) {
      DataSetList gridSetList = init.DSL().SelectGroupSets( gridSetArg, DataSet::GRID_3D );
      if (gridSetList.empty()) {
        mprinterr("Error: %s does not select any grid data set.\n", gridSetArg.c_str());
        return Action::ERR;
      }
      if (gridSetList.size() > 1) {
        mprintf("Warning: %s selects more than 1 grid data set. Only using the first set.\n", gridSetArg.c_str());
      }
      gridSet_ = (DataSet_3D*)gridSetList[0];
    }
  }
  // Check if IRED vector
  bool isIred = actionArgs.hasKey("ired"); 
  // Vector Mask
  if (mode_ != BOX && mode_ != BOX_X && mode_ != BOX_Y && mode_ != BOX_Z) {
    if (mask_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;
  }
  // Get second mask if necessary
  if (mode_ == MASK || mode_ == MINIMAGE) {
    std::string maskexpr = actionArgs.GetMaskNext();
    if (maskexpr.empty()) {
      mprinterr("Error: Specified vector mode (%s) requires a second mask.\n",
                ModeString_[ mode_ ]);
      return Action::ERR;
    }
    if (mask2_.SetMaskString( maskexpr )) return Action::ERR;
  }
  // Allocate Char mask for bond dipole
  if (mode_ == BONDDIPOLE)
    cmask_ = new CharMask();
  // Set up vector dataset and IRED status
  MetaData md(actionArgs.GetStringNext(), MetaData::M_VECTOR);
  if (isIred) md.SetScalarType( MetaData::IREDVEC );
  DataSet::DataType vtype = DataSet::VECTOR;
  //if (NeedsOrigin_[mode_])
  //  vtype = DataSet::VEC_OXYZ;
  //else
  //  vtype = DataSet::VEC_XYZ;
  Vec_ = (DataSet_Vector*)init.DSL().AddSet(vtype, md, "Vec");
  if (Vec_ == 0) return Action::ERR;
  // Add set to output file if not doing ptraj-compatible output
  if (!ptrajoutput_ && df != 0)
    df->AddDataSet( Vec_ );
  // Set up magnitude data set.
  if (calc_magnitude) {
    Magnitude_ = init.DSL().AddSet(DataSet::FLOAT, MetaData(Vec_->Meta().Name(), "Mag"));
    if (Magnitude_ == 0) return Action::ERR;
    if (df != 0) df->AddDataSet( Magnitude_ );
  }
  
  mprintf("    VECTOR: Type %s", ModeString_[ mode_ ]);
  if (calc_magnitude)
    mprintf(" (with magnitude)");
  if (isIred)
    mprintf(", IRED");
  if (mask_.MaskStringSet())
     mprintf(", mask [%s]", mask_.MaskString());
  if (mask2_.MaskStringSet())
    mprintf(", second mask [%s]", mask2_.MaskString());
  if (!filename.empty()) {
    if (ptrajoutput_)
      mprintf(", ptraj-compatible output to");
    else
      mprintf(", output to");
    mprintf(" %s", filename.c_str());
  }
  mprintf("\n");
  if (gridSet_ != 0)
    mprintf("\tExtracting box vectors from grid set '%s'\n", gridSet_->legend());
  if (mode_ == DIPOLE || mode_ == BONDDIPOLE) {
    if (mode_ == BONDDIPOLE) //TODO remove when no longer testing
      mprintf("Warning: 'bonddipole' mode is still experimental. Use with caution.\n");
    if (dipole_in_debye_)
      mprintf("\tDipole vector units will be Debye.\n");
    else
      mprintf("\tDipole vector units will be e-*Ang.\n");
  }

  return Action::OK;
}

// Action_Vector::Setup()
Action::RetType Action_Vector::Setup(ActionSetup& setup) {
  if (needBoxInfo_) {
    // Check for box info
    if (!setup.CoordInfo().TrajBox().HasBox()) {
      mprinterr("Error: vector box: No box information.\n",
                setup.Top().c_str());
      return Action::ERR;
    }
  }
  // Check for velocity/force
  if ((mode_ == MOMENTUM || mode_ == VELOCITY) && !setup.CoordInfo().HasVel()) {
    mprintf("Warning: vector %s requires velocity information. Skipping.\n", ModeString_[mode_]);
    return Action::SKIP;
  }
  if (mode_ == FORCE && !setup.CoordInfo().HasForce()) {
    mprintf("Warning: vector %s requires force information. Skipping.\n", ModeString_[mode_]);
    return Action::SKIP;
  }
  if (mask_.MaskStringSet()) {
    // Setup mask 1
    if (setup.Top().SetupIntegerMask(mask_)) return Action::ERR;
    mask_.MaskInfo();
    if (mask_.None()) {
      mprinterr("Error: First vector mask is empty.\n");
      return Action::ERR;
    }
    // Set up char mask if needed
    if (cmask_ != 0) {
      *cmask_ = CharMask( mask_.ConvertToCharMask(), mask_.Nselected() );
    }
  }

  // Setup mask 2
  if (mask2_.MaskStringSet()) {
    if (setup.Top().SetupIntegerMask(mask2_)) return Action::ERR;
    mask2_.MaskInfo();
    if (mask2_.None()) {
      mprinterr("Error: Second vector mask is empty.\n");
      return Action::ERR;
    }
  }
  // Check if center of mass is possible
  if (useMass_) {
    if (mask_.MaskStringSet() && setup.Top().MaskHasZeroMass( mask_ )) useMass_ = false;
    if (mask2_.MaskStringSet() && setup.Top().MaskHasZeroMass( mask2_ )) useMass_ = false;
  }
  CurrentParm_ = setup.TopAddress();
  return Action::OK;
}

// -----------------------------------------------------------------------------
Vec3 Action_Vector::GetVec(Frame const& currentFrame, AtomMask const& maskIn)
const
{
  if (useMass_)
    return currentFrame.VCenterOfMass(maskIn);
  else
    return currentFrame.VGeometricCenter(maskIn);
}

// Action_Vector::Mask()
void Action_Vector::Mask(Frame const& currentFrame) {
  Vec3 CXYZ = GetVec(currentFrame, mask_);
  Vec3 VXYZ = GetVec(currentFrame, mask2_);
  VXYZ -= CXYZ;
  Vec_->AddVxyzo(VXYZ, CXYZ);
}

/** Calculate dipole */
bool Action_Vector::calcBondDipole(Vec3& vBond, Vec3& vDipole, Vec3& vOrigin, double q0, Vec3 const& XYZ0, Vec3 const& XYZ1)
{
  bool hasCharge = true;
  if (q0 < 0.0) {
    vBond = XYZ0 - XYZ1;
    vDipole = vBond * (-q0);
    vOrigin = XYZ1;
  } else if (q0 > 0.0) {
    vBond = XYZ1 - XYZ0;
    vDipole = vBond * q0;
    vOrigin = XYZ0;
  } else
    hasCharge = false;
  return hasCharge;
}

/// \return true if all atoms bonded to atomIn (ignoring ignoreIdx) have only the one bond back to atomIn
static inline bool bondedAtomsHave1bond(Topology const& topIn, Atom const& atomIn, int ignoreIdx)
{
  for (int ii=0; ii < atomIn.Nbonds(); ii++) {
    int bondedAtomIdx = atomIn.Bond(ii);
    if (bondedAtomIdx != ignoreIdx) {
      if (topIn[bondedAtomIdx].Nbonds() != 1) return false;
    }
  }
  return true;
}

static inline double bondedAtomsCharge(Topology const& topIn, Atom const& atomIn, int ignoreIdx)
{
  double sumQ = 0.0;
  for (int ii=0; ii < atomIn.Nbonds(); ii++) {
    int bondedAtomIdx = atomIn.Bond(ii);
    if (bondedAtomIdx != ignoreIdx) {
      if (topIn[bondedAtomIdx].Nbonds() == 1)
        sumQ += topIn[bondedAtomIdx].Charge();
    }
  }
  return sumQ;
}

/** Calculate net dipole from individual bond dipoles. */
void Action_Vector::BondDipole_individualBonds(Frame const& currentFrame) {
  Vec3 VXYZ(0.0, 0.0, 0.0); // Negative, head
  Vec3 OXYZ(0.0, 0.0, 0.0); // Positive, tail
  Topology const& currentTop = *CurrentParm_;

  unsigned int Nbonds = 0;
  for (AtomMask::const_iterator iat = mask_.begin(); iat != mask_.end(); ++iat)
  {
    Atom const& atom0 = currentTop[*iat];
    Vec3 XYZ0( currentFrame.XYZ(*iat) );
    double q0 = atom0.Charge();
    for (Atom::bond_iterator bit = atom0.bondbegin(); bit != atom0.bondend(); ++bit)
    {
      if (*bit > *iat && cmask_->AtomInCharMask( *bit )) {
        // Both atoms in the bond are selected.
        Atom const& atom1 = currentTop[*bit];
        Vec3 XYZ1( currentFrame.XYZ(*bit) );
        double q1 = atom1.Charge();
        Vec3 vDipole(0.0, 0.0, 0.0);
        Vec3 vBond(0.0, 0.0, 0.0);
        Vec3 vOrigin(0.0, 0.0, 0.0);
        bool hasCharge = false;
        if (atom0.Element() != atom1.Element()) {
          if (atom0.Element()==Atom::HYDROGEN || atom0.Nbonds() == 1) {
            hasCharge = calcBondDipole( vBond, vDipole, vOrigin, q0, XYZ0, XYZ1 );
          } else if (atom1.Element()==Atom::HYDROGEN || atom1.Nbonds() == 1) {
            hasCharge = calcBondDipole( vBond, vDipole, vOrigin, q1, XYZ1, XYZ0 );
          } else if (atom0.Nbonds() > 1 && bondedAtomsHave1bond(currentTop, atom0, *bit)) {
            // atom0 is bonded to a single other atom, X-atom0-atom1. Sum up X and atom0 charge
            double q0X = q0 + bondedAtomsCharge(currentTop, atom0, *bit);
            if (debug_ > 1)
              mprintf("DEBUG: Total charge on %s is %f\n", currentTop.AtomMaskName(*iat).c_str(), q0X);
            hasCharge = calcBondDipole( vBond, vDipole, vOrigin, q0X, XYZ0, XYZ1 );
          } else if (atom1.Nbonds() > 1 && bondedAtomsHave1bond(currentTop, atom1, *iat)) {
            // atom1 is bonded to a single other atom, X-atom1-atom0. Sum up X and atom1 charge
            double q1X = q1 + bondedAtomsCharge(currentTop, atom1, *iat);
            if (debug_ > 1)
            mprintf("DEBUG: Total charge on %s is %f\n", currentTop.AtomMaskName(*bit).c_str(), q1X);
            hasCharge = calcBondDipole( vBond, vDipole, vOrigin, q1X, XYZ1, XYZ0 );
          } else {
            mprintf("Warning: Unhandled number of bonds: %s=%i, %s=%i\n",
                     currentTop.AtomMaskName(*iat).c_str(), atom0.Nbonds(),
                     currentTop.AtomMaskName(*bit).c_str(), atom1.Nbonds());
          }
        }

        if (hasCharge) {
          // Subtract half the dipole from the bond center to get the dipole origin
          //Vec3 vBondCtr = vBond / 2.0;
          //Vec3 vDipoleHalf = vDipole / 2.0;
          //Vec3 vOrigin = vOrigin + vBondCtr - vDipoleHalf;
          // DEBUG
          if (debug_ > 1) {
            mprintf("DEBUG: Bond %s(%f) -- %s(%f) bxyz={%f %f %f} oxyz={%f %f %f} vxyz={%f %f %f} mag=%f\n",
                    currentTop.AtomMaskName(*iat).c_str(), q0,
                    currentTop.AtomMaskName(*bit).c_str(), q1,
                    vBond[0], vBond[1], vBond[2],
                    vOrigin[0], vOrigin[1], vOrigin[2],
                    vDipole[0], vDipole[1], vDipole[2],
                    sqrt(vDipole.Magnitude2()));
          }
          // Sum
          VXYZ += vDipole;
          OXYZ += vOrigin;
          Nbonds++;
        }
      }
    } // END loop over bonded atoms
  } // END loop over mask atoms
  if (dipole_in_debye_)
    VXYZ /= Constants::DEBYE_EA;
  OXYZ /= (double)Nbonds;
  Vec_->AddVxyzo( VXYZ, OXYZ );
}

/** Calculate net dipole from atoms, origin from individual bond dipoles. */
void Action_Vector::BondDipole_net_bondOrigin(Frame const& currentFrame) {
  Vec3 VXYZ(0.0, 0.0, 0.0); // Negative, head
  Vec3 OXYZ(0.0, 0.0, 0.0); // Positive, tail
  Topology const& currentTop = *CurrentParm_;

  unsigned int Nbonds = 0;
  for (AtomMask::const_iterator iat = mask_.begin(); iat != mask_.end(); ++iat)
  {
    Atom const& atom0 = currentTop[*iat];
    Vec3 XYZ0( currentFrame.XYZ(*iat) );
    double q0 = atom0.Charge();
    double en0 = atom0.PaulingElectroNeg();
    // Note that negative charge will decrease the vector (towards the origin)
    // but by convention we want to point towards the negative charge, so at
    // the end the vector will have to be negated.
    VXYZ += ( XYZ0 * q0 );
    // Store bond origins according to electronegativity.
    for (Atom::bond_iterator bit = atom0.bondbegin(); bit != atom0.bondend(); ++bit)
    {
      if (*bit > *iat && cmask_->AtomInCharMask( *bit )) {
        // Both atoms in the bond are selected.
        Atom const& atom1 = currentTop[*bit];
        Vec3 XYZ1( currentFrame.XYZ(*bit) );
        double en1 = atom1.PaulingElectroNeg();
        Vec3 vOrigin(0.0, 0.0, 0.0);
        bool en_is_different = true;
        if (en0 > en1) {
          // atom0 is more negative
          //vBond = XYZ0 - XYZ1;
          vOrigin = XYZ1;
        } else if (en1 > en0) {
          // atom1 is more negative
          //vBond = XYZ1 - XYZ0;
          vOrigin = XYZ0;
        } else
          en_is_different = false;

        if (en_is_different) {
          // DEBUG
          if (debug_ > 1) {
            double q1 = atom1.Charge(); // DEBUG
            mprintf("DEBUG: Bond %s(%f) -- %s(%f) oxyz={%f %f %f}\n",
                    currentTop.AtomMaskName(*iat).c_str(), q0,
                    currentTop.AtomMaskName(*bit).c_str(), q1,
                    vOrigin[0], vOrigin[1], vOrigin[2]);
          }
          // Sum
          OXYZ += vOrigin;
          Nbonds++;
        }
      }
    } // END loop over bonded atoms
  } // END loop over mask atoms
  if (dipole_in_debye_)
    VXYZ /= Constants::DEBYE_EA;
  OXYZ /= (double)Nbonds;
  VXYZ.Neg();
  // Shift it so the center of VXYZ will end up at the origin.
  OXYZ -= (VXYZ / 2.0);
  Vec_->AddVxyzo( VXYZ, OXYZ );
}

// Action_Vector::Dipole()
void Action_Vector::Dipole(Frame const& currentFrame) {
  Vec3 VXYZ(0.0, 0.0, 0.0);
  Vec3 CXYZ(0.0, 0.0, 0.0);
  double total_mass = 0;
  for (AtomMask::const_iterator atom = mask_.begin();
                                atom != mask_.end(); ++atom)
  {
    double mass = (*CurrentParm_)[*atom].Mass();
    total_mass += mass;
    Vec3 XYZ = currentFrame.XYZ( *atom );
    CXYZ += ( XYZ * mass );
    double charge = (*CurrentParm_)[*atom].Charge();
    XYZ *= charge;
    VXYZ += ( XYZ );
  }
  CXYZ /= total_mass;
  if (dipole_in_debye_)
    VXYZ /= Constants::DEBYE_EA;
  Vec_->AddVxyzo( VXYZ, CXYZ );
}

// Action_Vector::Principal()
void Action_Vector::Principal(Frame const& currentFrame) {
  Matrix_3x3 Inertia;
  Vec3 Eval;

  // Origin is center of atoms in mask_ 
  Vec3 OXYZ = currentFrame.CalculateInertia( mask_, Inertia );
  // NOTE: Diagonalize_Sort_Chirality places sorted eigenvectors in rows.
  Inertia.Diagonalize_Sort_Chirality( Eval, 0 );
  //Eval.Print("PRINCIPAL EIGENVALUES");
  //Inertia.Print("PRINCIPAL EIGENVECTORS (Rows)");
  DataSet_Vector& vec = static_cast<DataSet_Vector&>( *Vec_ );
  if ( mode_ == PRINCIPAL_X ) 
    vec.AddVxyzo( Inertia.Row1(), OXYZ ); // First row = first eigenvector
  else if ( mode_ == PRINCIPAL_Y )
    vec.AddVxyzo( Inertia.Row2(), OXYZ ); // Second row = second eigenvector
  else // PRINCIPAL_Z
    vec.AddVxyzo( Inertia.Row3(), OXYZ ); // Third row = third eigenvector
}

// Action_Vector::CorrPlane()
void Action_Vector::CorrPlane(Frame const& currentFrame) {
  vcorr_.CalcLeastSquaresPlane( currentFrame, mask_, useMass_ );
  Vec_->AddVxyzo( vcorr_.Nxyz(), vcorr_.Cxyz() );
}

//  Action_Vector::UnitCell()
void Action_Vector::UnitCell(Box const& box, Vec3 const& oxyz) {
  switch ( mode_ ) {
    case BOX_X   : Vec_->AddVxyzo( box.UnitCell().Row1(), oxyz ); break;
    case BOX_Y   : Vec_->AddVxyzo( box.UnitCell().Row2(), oxyz ); break;
    case BOX_Z   : Vec_->AddVxyzo( box.UnitCell().Row3(), oxyz ); break;
    case BOX_CTR : Vec_->AddVxyz( box.UnitCell().TransposeMult(Vec3(0.5)) ); break;
    default: return;
  }
}

/** Store box vector (A, B, C) lengths as a single vector (|A|, |B|, |C|). */
void Action_Vector::BoxLengths(Box const& box) {
  Vec_->AddVxyz( box.Lengths() );
}

// Action_Vector::MinImage()
void Action_Vector::MinImage(Frame const& frm) {
  Vec3 com1 = GetVec(frm, mask_);
  Vec_->AddVxyzo( MinImagedVec(com1, GetVec(frm, mask2_), frm.BoxCrd().UnitCell(), frm.BoxCrd().FracCell()), com1 );
}

/// \return The center of selected elements in given array.
static inline Vec3 CalcCenter(const double* xyz, AtomMask const& maskIn) {
  Vec3 Coord(0.0);
  for (AtomMask::const_iterator at = maskIn.begin(); at != maskIn.end(); ++at)
  {
    int idx = *at * 3;
    Coord[0] += xyz[idx  ];
    Coord[1] += xyz[idx+1];
    Coord[2] += xyz[idx+2];
  }
  Coord /= (double)maskIn.Nselected();
  return Coord;
}

// Action_Vector::DoAction()
Action::RetType Action_Vector::DoAction(int frameNum, ActionFrame& frm) {
  switch ( mode_ ) {
    case MASK        : Mask(frm.Frm()); break;
    case CENTER      : Vec_->AddVxyz( GetVec(frm.Frm(), mask_) ); break;
    case MOMENTUM    : Vec_->AddVxyz( frm.Frm().VMomentum(mask_) ); break;
    case VELOCITY    : Vec_->AddVxyz( CalcCenter(frm.Frm().vAddress(), mask_) ); break;
    case FORCE       : Vec_->AddVxyz( CalcCenter(frm.Frm().fAddress(), mask_) ); break; 
    case DIPOLE      : Dipole(frm.Frm()); break;
//    case BONDDIPOLE  : BondDipole_individualBonds(frm.Frm()); break;
    case BONDDIPOLE  : BondDipole_net_bondOrigin(frm.Frm()); break;
    case PRINCIPAL_X :
    case PRINCIPAL_Y :
    case PRINCIPAL_Z : Principal(frm.Frm()); break;
    case CORRPLANE   : CorrPlane(frm.Frm()); break;
    case BOX         :
      if (gridSet_ != 0)
        BoxLengths( gridSet_->Bin().GridBox() );
      else
        BoxLengths( frm.Frm().BoxCrd() );
      break;
    case BOX_X       : 
    case BOX_Y       : 
    case BOX_Z       : 
    case BOX_CTR     :
      if (gridSet_ != 0)
        UnitCell( gridSet_->Bin().GridBox(), gridSet_->Bin().GridOrigin() );
      else
        UnitCell( frm.Frm().BoxCrd(), Vec3(0.0) );
      break;
    case MINIMAGE    : MinImage( frm.Frm() ); break; 
    default          : return Action::ERR; // NO_OP
  } // END switch over vectorMode
  if (Magnitude_ != 0) {
    float mag = (float)(sqrt(Vec_->Back().Magnitude2()));
    Magnitude_->Add(frameNum, &mag);
  }
  return Action::OK;
}

// Action_Vector::Print()
void Action_Vector::Print() {
  if (ptrajoutput_) {
    mprintf("    VECTOR: writing ptraj-style vector information for %s\n", Vec_->legend());
    outfile_->Printf("# FORMAT: frame vx vy vz cx cy cz cx+vx cy+vy cz+vz\n"
                   "# FORMAT where v? is vector, c? is center of mass...\n");
    int totalFrames = Vec_->Size();
    if (Vec_->HasOrigins()) {
      DataSet_Vector const& vec =
        static_cast<DataSet_Vector const&>( *Vec_ );
      for (int i=0; i < totalFrames; ++i) {
        Vec3 const& vxyz = vec[i];
        Vec3 const& cxyz = vec.OXYZ(i);
        Vec3 txyz  = cxyz + vxyz;
        outfile_->Printf("%i %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
                i+1, vxyz[0], vxyz[1], vxyz[2], cxyz[0], cxyz[1], cxyz[2],
                txyz[0], txyz[1], txyz[2]);
      }
    } else {
      // No origin
      for (int i=0; i < totalFrames; ++i) {
        Vec3 const& vxyz = (*Vec_)[i];
        const Vec3 cxyz(0.0);
        Vec3 txyz  = cxyz + vxyz;
        outfile_->Printf("%i %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
                i+1, vxyz[0], vxyz[1], vxyz[2], cxyz[0], cxyz[1], cxyz[2],
                txyz[0], txyz[1], txyz[2]);
      }
    }
  }
}
