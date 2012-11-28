// AxisType
#include <map>
#include <cstring> // memcpy
#include "AxisType.h"
#include "CpptrajStdio.h"
#ifdef NASTRUCTDEBUG
#include "StringRoutines.h" // integerToString
#endif

// ---------- NA REFERENCE BASE ATOM NAMES AND COORDS --------------------------
// TODO: Since the names get converted to NameType probably dont need extra space
const int AxisType::ADENATOM = 11;
const char AxisType::ADEnames[ADENATOM][5] = {
"C1' ","N9  ","C8  ","N7  ","C5  ","C6  ","N6  ","N1  ","C2  ","N3  ","C4  "
};
const int AxisType::ADEhbonds[ADENATOM] = {
-1,    -1,    -1,    -1,    -1,    -1,    0,     1,     -1,    -1,    -1
};
const double AxisType::ADEcoords[ADENATOM][3] = {
     {-2.479000,  5.346000,  0.000000,},
     {-1.291000,  4.498000,  0.000000,},
     { 0.024000,  4.897000,  0.000000,},
     { 0.877000,  3.902000,  0.000000,},
     { 0.071000,  2.771000,  0.000000,},
     { 0.369000,  1.398000,  0.000000,},
     { 1.611000,  0.909000,  0.000000,},
     {-0.668000,  0.532000,  0.000000,},
     {-1.912000,  1.023000,  0.000000,},
     {-2.320000,  2.290000,  0.000000,},
     {-1.267000,  3.124000,  0.000000 }
};
const int AxisType::CYTNATOM = 9;
const char AxisType::CYTnames[CYTNATOM][5] = {
"C1' ","N1  ","C2  ","O2  ","N3  ","C4  ","N4  ","C5  ","C6  "
};
const int AxisType::CYThbonds[CYTNATOM] = {
-1,    -1,    -1,    2,     1,     -1,    0,     -1,    -1
};
const double AxisType::CYTcoords[CYTNATOM][3] = {
     {-2.477000,  5.402000,  0.000000,},
     {-1.285000,  4.542000,  0.000000,},
     {-1.472000,  3.158000,  0.000000,},
     {-2.628000,  2.709000,  0.000000,},
     {-0.391000,  2.344000,  0.000000,},
     { 0.837000,  2.868000,  0.000000,},
     { 1.875000,  2.027000,  0.000000,},
     { 1.056000,  4.275000,  0.000000,},
     {-0.023000,  5.068000,  0.000000 }
};
const int AxisType::GUANATOM = 12;
const char AxisType::GUAnames[GUANATOM][5] = {
"C1' ","N9  ","C8  ","N7  ","C5  ","C6  ","O6  ","N1  ","C2  ","N2  ","N3  ","C4  "
};
const int AxisType::GUAhbonds[GUANATOM] = {
-1,    -1,    -1,    -1,    -1,    -1,    0,     1,     -1,    2,     -1,    -1
};
const double AxisType::GUAcoords[GUANATOM][3] = {
     {-2.477000,  5.399000,  0.000000,},
     {-1.289000,  4.551000,  0.000000,},
     { 0.023000,  4.962000,  0.000000,},
     { 0.870000,  3.969000,  0.000000,},
     { 0.071000,  2.833000,  0.000000,},
     { 0.424000,  1.460000,  0.000000,},
     { 1.554000,  0.955000,  0.000000,},
     {-0.700000,  0.641000,  0.000000,},
     {-1.999000,  1.087000,  0.000000,},
     {-2.949000,  0.139000, -0.001000,},
     {-2.342000,  2.364000,  0.001000,},
     {-1.265000,  3.177000,  0.000000 }
};
const int AxisType::THYNATOM = 10;
const char AxisType::THYnames[THYNATOM][5] = {
"C1' ","N1  ","C2  ","O2  ","N3  ","C4  ","O4  ","C5  ","C7  ","C6  "
};
const int AxisType::THYhbonds[THYNATOM] = {
-1,    -1,    -1,    -1,    1,     -1,    0,     -1,    -1,    -1
}; 
const double AxisType::THYcoords[THYNATOM][3] = {
     {-2.481000,  5.354000,  0.000000,},
     {-1.284000,  4.500000,  0.000000,},
     {-1.462000,  3.135000,  0.000000,},
     {-2.562000,  2.608000,  0.000000,},
     {-0.298000,  2.407000,  0.000000,},
     { 0.994000,  2.897000,  0.000000,},
     { 1.944000,  2.119000,  0.000000,},
     { 1.106000,  4.338000,  0.000000,},
     { 2.466000,  4.961000,  0.001000,},
     {-0.024000,  5.057000,  0.000000 }
};
const int AxisType::URANATOM = 9;
const char AxisType::URAnames[URANATOM][5] = {
"C1' ","N1  ","C2  ","O2  ","N3  ","C4  ","O4  ","C5  ","C6  " 
};
const int AxisType::URAhbonds[URANATOM] = {
-1,    -1,    -1,    -1,    1,     -1,    0,     -1,    -1
}; 
const double AxisType::URAcoords[URANATOM][3] = {
     {-2.481000,  5.354000,  0.000000,},
     {-1.284000,  4.500000,  0.000000,},
     {-1.462000,  3.131000,  0.000000,},
     {-2.563000,  2.608000,  0.000000,},
     {-0.302000,  2.397000,  0.000000,},
     { 0.989000,  2.884000,  0.000000,},
     { 1.935000,  2.094000, -0.001000,},
     { 1.089000,  4.311000,  0.000000,},
     {-0.024000,  5.053000,  0.000000 }
};

// UNKNOWN_BASE, ADE, CYT, GUA, THY, URA 
/// Base names corresponding to NAbaseType
const char AxisType::NAbaseName[6][4] = { "UNK", "ADE", "CYT", "GUA", "THY", "URA" };

// ------------------------- AXISTYPE FUNCTIONS -------------------------------
// CONSTRUCTOR
AxisType::AxisType() {
  X_ = 0;
  natom_ = 0;
  maxnatom_ = 0;
  Ncoord_ = 0;
  ID = UNKNOWN_BASE;
  origin[0]=0;
  origin[1]=0;
  origin[2]=0;
  R[0]=0.0; R[1]=0.0; R[2]=0.0;
  R[3]=0.0; R[4]=0.0; R[5]=0.0;
  R[6]=0.0; R[7]=0.0; R[8]=0.0;
  HbondCoord[0]=NULL;
  HbondCoord[1]=NULL;
  HbondCoord[2]=NULL;
  HbondAtom[0]=-1;
  HbondAtom[1]=-1;
  HbondAtom[2]=-1;
  residue_number=-1;
  second_resnum=-1;
  patomidx_ = -1;
  o4atomidx_ = -1;
}

// DESTRUCTOR
AxisType::~AxisType() { 
  if (X_!=0) delete[] X_;
}

// COPY CONSTRUCTOR
AxisType::AxisType(const AxisType &rhs) 
{
  natom_ = rhs.natom_;
  maxnatom_ = rhs.maxnatom_;
  Ncoord_ = rhs.Ncoord_;
  if (Ncoord_ > 0) {
    X_ = new double[ maxnatom_ * 3 ];
    memcpy(X_, rhs.X_, Ncoord_ * sizeof(double));
  } else
    X_ = 0;
  ID = rhs.ID;
  Name = rhs.Name;
  origin[0] = rhs.origin[0];
  origin[1] = rhs.origin[1];
  origin[2] = rhs.origin[2];
  memcpy(R, rhs.R, 9 * sizeof(double));
  HbondAtom[0] = rhs.HbondAtom[0];
  HbondAtom[1] = rhs.HbondAtom[1];
  HbondAtom[2] = rhs.HbondAtom[2];
  // Since HbondCoord contains memory addresses, it must be updated
  // relative to this instances X
  HbondCoord[0] = X_ + (HbondAtom[0]*3);
  HbondCoord[1] = X_ + (HbondAtom[1]*3);
  HbondCoord[2] = X_ + (HbondAtom[2]*3);
  residue_number = rhs.residue_number;
  second_resnum = rhs.second_resnum;
  patomidx_ = rhs.patomidx_;
  o4atomidx_ = rhs.o4atomidx_;
}

// Assignment Operator
// NOTE: Only allocate/deallocate when natom>maxnatom?
AxisType &AxisType::operator=(const AxisType &rhs) {
  // Check for self assignment
  if ( this == &rhs ) return *this;
  // Deallocate
  if (X_!=0) delete[] X_;
  // Allocate and copy
  natom_ = rhs.natom_;
  maxnatom_ = rhs.maxnatom_;
  Ncoord_ = rhs.Ncoord_;
  if (Ncoord_ > 0) {
    X_ = new double[ maxnatom_ * 3 ];
    memcpy(X_, rhs.X_, Ncoord_ * sizeof(double));
  } else
    X_ = 0;
  ID = rhs.ID;
  Name = rhs.Name;
  origin[0] = rhs.origin[0];
  origin[1] = rhs.origin[1];
  origin[2] = rhs.origin[2];
  memcpy(R, rhs.R, 9 * sizeof(double));
  HbondAtom[0] = rhs.HbondAtom[0];
  HbondAtom[1] = rhs.HbondAtom[1];
  HbondAtom[2] = rhs.HbondAtom[2];
  // Since HbondCoord contains memory addresses, it must be updated
  // relative to this instances X
  HbondCoord[0] = X_ + (HbondAtom[0]*3);
  HbondCoord[1] = X_ + (HbondAtom[1]*3);
  HbondCoord[2] = X_ + (HbondAtom[2]*3);
  residue_number = rhs.residue_number;
  second_resnum = rhs.second_resnum;
  patomidx_ = rhs.patomidx_;
  o4atomidx_ = rhs.o4atomidx_;
  // Return *this
  return *this;
} 

// ------------------------- VARIABLE FUNCTIONS -------------------------------
// AxisType::RX()
/** Put the X unit vector from rotation matrix into Vin.  */
void AxisType::RX(double Vin[3]) {
  Vin[0] = R[0];
  Vin[1] = R[3];
  Vin[2] = R[6];
}
// AxisType::RY()
/** Put the Y unit vector from rotation matrix into Vin.  */
void AxisType::RY(double Vin[3]) {
  Vin[0] = R[1];
  Vin[1] = R[4];
  Vin[2] = R[7];
}
// AxisType::RZ()
/** Put the Z unit vector from rotation matrix into Vin.  */
void AxisType::RZ(double Vin[3]) {
  Vin[0] = R[2];
  Vin[1] = R[5];
  Vin[2] = R[8];
}
// AxisType::OXYZ()
void AxisType::OXYZ(double Vin[3]) {
  Vin[0] = origin[0];
  Vin[1] = origin[1];
  Vin[2] = origin[2];
}
// AxisType::Origin()
/** Return pointer to the origin coords.  */
double *AxisType::Origin() {
  return origin;
}

// ------------------------- PRIVATE FUNCTIONS --------------------------------

// AxisType::ID_base()
/** Return a number indicating if this is a NA base. If not recognized, 
  * return -1. 
  * NOTE: Currently based only on amber residue names. Will not recognize
  * non-standard bases.
  */
AxisType::NAbaseType AxisType::ID_base(NameType const& resname) {
  // If residue name begins with D, assume AMBER DNA residue
  if (resname[0]=='D') {
    switch (resname[1]) {
      case 'A': return ADE;
      case 'C': return CYT;
      case 'G': return GUA;
      case 'T': return THY;
    }
  // If residue name beings with R, assume AMBER RNA residue
  } else if (resname[0]=='R') {
    switch (resname[1]) {
      case 'A': return ADE;
      case 'C': return CYT;
      case 'G': return GUA;
      case 'U': return URA;
    }
  // Look for standard 3 letter/1 letter NA residue names
  } else {
    if ( resname == "ADE " ) return ADE;
    if ( resname == "CYT " ) return CYT;
    if ( resname == "GUA " ) return GUA;
    if ( resname == "THY " ) return THY;
    if ( resname == "URA " ) return URA;
    if ( resname == "A   " ) return ADE;
    if ( resname == "C   " ) return CYT;
    if ( resname == "G   " ) return GUA;
    if ( resname == "T   " ) return THY;
    if ( resname == "U   " ) return URA;
  } 
  return UNKNOWN_BASE;
}

// AxisType::AllocAxis()
/** Allocate mem for coords and names. Replace any existing coords / names.
  */
int AxisType::AllocAxis(int natomIn) {
  if (X_!=NULL) delete[] X_;
  Name.clear(); 
  natom_ = natomIn;
  maxnatom_ = natom_;
  Ncoord_ = natom_ * 3;
  X_ = new double[ Ncoord_ ];
  Name.resize( natom_ );
  return 0;
}

// ------------------------- PUBLIC FUNCTIONS ---------------------------------
#ifdef NASTRUCTDEBUG
// AxisType::BaseName()
const char* AxisType::BaseName() {
  if (basename_num_.empty()) 
    basename_num_ = integerToString(residue_number + 1) + ":" + std::string(NAbaseName[ID]);
  return basename_num_.c_str();
}
#endif

// AxisType::ResName()
const char* AxisType::ResName() {
  return NAbaseName[ID];
}

// AxisType::AtomNameIs()
bool AxisType::AtomNameIs(int atom, char *nameIn) {
  // TODO: pass a string instead of char*
  if (atom<0 || atom>=natom_) return false;
  NameType name = nameIn;
  if (Name[atom] == name) return true;
  return false;
}

// AxisType::AtomName()
const char* AxisType::AtomName(int atom) {
  if (atom<0 || atom>=natom_) return NULL;
  return (*Name[atom]);
}

// AxisType::PrintAtomNames()
void AxisType::PrintAtomNames() {
  for (int atom = 0; atom < natom_; atom++)
    mprintf(" %s",*Name[atom]);
  mprintf("\n");
}

// AxisType::PrintAxisInfo()
/** Print origin and rotation matrix for this base.
  */
void AxisType::PrintAxisInfo(const char *title) {
  mprintf("         %s origin: %8.4lf %8.4lf %8.4lf\n",title,origin[0],origin[1],origin[2]);
  mprintf("         %s Rx vec: %8.4lf %8.4lf %8.4lf\n",title,R[0],R[3],R[6]);
  mprintf("         %s Ry vec: %8.4lf %8.4lf %8.4lf\n",title,R[1],R[4],R[7]);
  mprintf("         %s Rz vec: %8.4lf %8.4lf %8.4lf\n",title,R[2],R[5],R[8]);
}

void AxisType::SetCoordsFromFrame( Frame& frameIn ) {
  for (int i = 0; i < frameIn.size(); ++i)
    X_[i] = frameIn[i];
}

// AxisType::StoreRotMatrix()
/** Store the rotation matrix and origin coordinates associated with
  * this base.
  */
void AxisType::StoreRotMatrix(const double *RotMatrix, const double *originIn) {
  R[0] = RotMatrix[0];
  R[1] = RotMatrix[1];
  R[2] = RotMatrix[2];
  R[3] = RotMatrix[3];
  R[4] = RotMatrix[4];
  R[5] = RotMatrix[5];
  R[6] = RotMatrix[6];
  R[7] = RotMatrix[7];
  R[8] = RotMatrix[8];
  origin[0] = originIn[0];
  origin[1] = originIn[1];
  origin[2] = originIn[2];
}

// AxisType::StoreBPresnums()
void AxisType::StoreBPresnums(int r1, int r2) {
  residue_number = r1;
  second_resnum = r2;
}

// AxisType::SetRefCoord()
/** Set NA residue reference coordinates for given NA base. Ensure that
  * the atom ordering in the reference matches that in the given parm.
  * Set up a mask to select correct atoms in the given parm.
  * \return 0 on success, 1 on error.
  */
AxisType::RefReturn AxisType::SetRefCoord(Topology *currentParm, int resnum, 
                                          AtomMask &parmMask, AtomMask &fitMask,
                                          AxisType::NAbaseType customBaseType) 
{
  std::map<int,int> BaseMap;
  std::map<int,int>::iterator atom;
  std::vector<double> BaseRefCoords;
  std::vector<int> BaseHbonds;
  std::vector<NameType> BaseRefNames;
  int Nbaseatom = 0;

  // First, identify the base
  if (customBaseType==UNKNOWN_BASE) 
    ID = ID_base(currentParm->Res(resnum).Name());
  else
    ID = customBaseType;
  // If unknown exit now
  if (ID == UNKNOWN_BASE) return NA_UNKNOWN; 

  // Set up reference coord and name arrays based on the base
  switch (ID) {
    case ADE :
      if (AllocAxis(ADENATOM)) return NA_ERROR;
      Nbaseatom = ADENATOM;
      BaseRefCoords.reserve(ADENATOM);
      BaseRefNames.reserve( ADENATOM );
      BaseHbonds.reserve(ADENATOM);
      for (int i = 0; i < ADENATOM; i++) {
        BaseRefCoords.push_back( ADEcoords[i][0] );
        BaseRefCoords.push_back( ADEcoords[i][1] );
        BaseRefCoords.push_back( ADEcoords[i][2] );
        BaseRefNames.push_back( ADEnames[i] );
        BaseHbonds.push_back( ADEhbonds[i] );
      }
      break;
    case CYT :
      if (AllocAxis(CYTNATOM)) return NA_ERROR;
      Nbaseatom = CYTNATOM;
      BaseRefCoords.reserve(CYTNATOM);
      BaseRefNames.reserve( CYTNATOM );
      BaseHbonds.reserve(CYTNATOM);
      for (int i = 0; i < CYTNATOM; i++) {
        BaseRefCoords.push_back( CYTcoords[i][0] );
        BaseRefCoords.push_back( CYTcoords[i][1] );
        BaseRefCoords.push_back( CYTcoords[i][2] );
        BaseRefNames.push_back( CYTnames[i] );
        BaseHbonds.push_back( CYThbonds[i] );
      }
      break;
    case GUA :
      if (AllocAxis(GUANATOM)) return NA_ERROR;
      Nbaseatom = GUANATOM;
      BaseRefCoords.reserve(GUANATOM);
      BaseRefNames.reserve( GUANATOM );
      BaseHbonds.reserve(GUANATOM);
      for (int i = 0; i < GUANATOM; i++) {
        BaseRefCoords.push_back( GUAcoords[i][0] );
        BaseRefCoords.push_back( GUAcoords[i][1] );
        BaseRefCoords.push_back( GUAcoords[i][2] );
        BaseRefNames.push_back( GUAnames[i] );
        BaseHbonds.push_back( GUAhbonds[i] );
      }
      break;
    case THY :
      if (AllocAxis(THYNATOM)) return NA_ERROR;
      Nbaseatom = THYNATOM;
      BaseRefCoords.reserve(THYNATOM);
      BaseRefNames.reserve( THYNATOM );
      BaseHbonds.reserve(THYNATOM);
      for (int i = 0; i < THYNATOM; i++) {
        BaseRefCoords.push_back( THYcoords[i][0] );
        BaseRefCoords.push_back( THYcoords[i][1] );
        BaseRefCoords.push_back( THYcoords[i][2] );
        BaseRefNames.push_back( THYnames[i] );
        BaseHbonds.push_back( THYhbonds[i] );
      }
      break;
    case URA :
      if (AllocAxis(URANATOM)) return NA_ERROR;
      Nbaseatom = URANATOM;
      BaseRefCoords.reserve(URANATOM);
      BaseRefNames.reserve( URANATOM );
      BaseHbonds.reserve(URANATOM);
      for (int i = 0; i < URANATOM; i++) {
        BaseRefCoords.push_back( URAcoords[i][0] );
        BaseRefCoords.push_back( URAcoords[i][1] );
        BaseRefCoords.push_back( URAcoords[i][2] );
        BaseRefNames.push_back( URAnames[i] );
        BaseHbonds.push_back( URAhbonds[i] );
      }
      break;
    case UNKNOWN_BASE: // Shouldnt get here, just for completeness
      mprintf("Warning: AxisType::SetRefCoord: Missing parameters for residue %s.\n",
              currentParm->ResidueName(resnum));
      return NA_UNKNOWN;
  }

  // For each atom defined as a reference atom for this base, find the
  // corresponding atom in the parm.
  for (int ref = 0; ref < Nbaseatom; ref++) {
    int parmAtom = currentParm->FindAtomInResidue(resnum, BaseRefNames[ref]);
    // Sometimes C1' is listed as C1*; if search for C1' fails look for
    // C1*.
    if (parmAtom < 0 && BaseRefNames[ref] == "C1' ")
      parmAtom = currentParm->FindAtomInResidue(resnum, "C1* ");
    if (parmAtom<0) {
      mprinterr("Error: Ref Atom [%s] not found in NA base [%s].\n",
                *BaseRefNames[ref], currentParm->ResidueName(resnum));
      return NA_ERROR;
    } else {
      BaseMap.insert( std::pair<int,int>(parmAtom, ref) );
    }
  }
  // See if residue contains a phosphate atom
  patomidx_ = currentParm->FindAtomInResidue( resnum, "P   " );
  //if (debug_ > 0 && patomidx_ != -1)
  //  mprintf("\tPhosphorus atom found: %i\n", patomidx_+1);
  // See if residue contains an O4' atom
  o4atomidx_ = currentParm->FindAtomInResidue( resnum, "O4' ");
  //if (debug_ > 0 && o4atomidx_ != -1)
  //  mprintf("\tO4' atom found: %i\n", o4atomidx_+1);
  // Now insert ref coords in same order as parm. Also add the parm
  // atom to the mask.
  parmMask.ResetMask();
  fitMask.ResetMask();
  int coord = 0;
  int coord3 = 0;
  for (atom = BaseMap.begin(); atom != BaseMap.end(); atom++) {
    int refatom = (*atom).second;
    int refcoord = refatom * 3;
    // Check if this is an H bonding atom. If so, store the memory address
    // of the appropriate place in the coordinate array for future Hbond
    // calculations.
    int hb_index = BaseHbonds[refatom];
    if ( hb_index != -1 ) {
      HbondCoord[hb_index] = X_ + coord3;
      HbondAtom[hb_index] = coord;
    }
    X_[coord3++] = BaseRefCoords[refcoord  ]; 
    X_[coord3++] = BaseRefCoords[refcoord+1]; 
    X_[coord3++] = BaseRefCoords[refcoord+2]; 
    Name[coord] = BaseRefNames[refatom];
    parmMask.AddAtom( (*atom).first );
    // Will this atom be used for RMS fitting?
    switch (ID) {
      case ADE:
      case GUA:
        if (BaseRefNames[refatom]=="N9  " ||
            BaseRefNames[refatom]=="C8  " ||
            BaseRefNames[refatom]=="N7  " ||
            BaseRefNames[refatom]=="C5  " ||
            BaseRefNames[refatom]=="C6  " ||
            BaseRefNames[refatom]=="N1  " ||
            BaseRefNames[refatom]=="C2  " ||
            BaseRefNames[refatom]=="N3  " ||
            BaseRefNames[refatom]=="C4  "   ) 
          fitMask.AddAtom( coord );
          break;
      case CYT:
      case THY:
      case URA:
        if (BaseRefNames[refatom]=="N1  " ||
            BaseRefNames[refatom]=="C2  " ||
            BaseRefNames[refatom]=="N3  " ||
            BaseRefNames[refatom]=="C4  " ||
            BaseRefNames[refatom]=="C5  " ||
            BaseRefNames[refatom]=="C6  "   )
          fitMask.AddAtom( coord );
          break;
      case UNKNOWN_BASE: return NA_ERROR;
    }
    ++coord;
  }
  residue_number = resnum;

  return NA_OK;
}

// AxisType::FlipYZ
/** Flip the Z and Y axes. Equivalent to rotation around the X axis.
  * Done for antiparallel stranded DNA.
  */
void AxisType::FlipYZ() {
  R[1] = -R[1]; // -Yx
  R[4] = -R[4]; // -Yy
  R[7] = -R[7]; // -Yz
  R[2] = -R[2]; // -Zx
  R[5] = -R[5]; // -Zy
  R[8] = -R[8]; // -Zz
}

// AxisType::FlipXY
/** Flip the X and Y axes. Equivalent to rotation around the Z axis.
  * Done for parallel stranded DNA.
  */
void AxisType::FlipXY() {
  R[0] = -R[0]; // -Xx
  R[3] = -R[3]; // -Xy
  R[6] = -R[6]; // -Xz
  R[1] = -R[1]; // -Yx
  R[4] = -R[4]; // -Yy
  R[7] = -R[7]; // -Yz
}

