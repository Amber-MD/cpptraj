// AxisType
#include <map>
#include <cstring>
#include <cstdio> //sprintf
#include "AxisType.h"
#include "CpptrajStdio.h"
#ifdef NASTRUCTDEBUG
// DEBUG
#  include "PDBfileRoutines.h"
#endif

// ---------- NA REFERENCE BASE ATOM NAMES AND COORDS --------------------------
#define ADENATOM 11
static const NAME ADEnames[ADENATOM] = {
"C1' ","N9  ","C8  ","N7  ","C5  ","C6  ","N6  ","N1  ","C2  ","N3  ","C4  "
};
static const int ADEhbonds[ADENATOM] = {
-1,    -1,    -1,    -1,    -1,    -1,    0,     1,     -1,    -1,    -1
};
static const double ADEcoords[ADENATOM][3] = {
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
#define CYTNATOM 9
static const NAME CYTnames[CYTNATOM] = {
"C1' ","N1  ","C2  ","O2  ","N3  ","C4  ","N4  ","C5  ","C6  "
};
static const int CYThbonds[CYTNATOM] = {
-1,    -1,    -1,    2,     1,     -1,    0,     -1,    -1
};
static const double CYTcoords[CYTNATOM][3] = {
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
#define GUANATOM 12
static const NAME GUAnames[GUANATOM] = {
"C1' ","N9  ","C8  ","N7  ","C5  ","C6  ","O6  ","N1  ","C2  ","N2  ","N3  ","C4  "
};
static const int GUAhbonds[GUANATOM] = {
-1,    -1,    -1,    -1,    -1,    -1,    0,     1,     -1,    2,     -1,    -1
};
static const double GUAcoords[GUANATOM][3] = {
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
#define THYNATOM 10
static const NAME THYnames[THYNATOM] = {
"C1' ","N1  ","C2  ","O2  ","N3  ","C4  ","O4  ","C5  ","C7  ","C6  "
};
static const int THYhbonds[THYNATOM] = {
-1,    -1,    -1,    -1,    1,     -1,    0,     -1,    -1,    -1
}; 
static const double THYcoords[THYNATOM][3] = {
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
#define URANATOM 9
static const NAME URAnames[URANATOM] = {
"C1' ","N1  ","C2  ","O2  ","N3  ","C4  ","O4  ","C5  ","C6  " 
};
static const int URAhbonds[URANATOM] = {
-1,    -1,    -1,    -1,    1,     -1,    0,     -1,    -1
}; 
static const double URAcoords[URANATOM][3] = {
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
// Mostly empty - Default Frame constructor should init everything else
AxisType::AxisType() {
  ID = UNKNOWN_BASE;
  Name = NULL;
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
}

// COPY CONSTRUCTOR
// NOTE: Base class copy is part of the call
AxisType::AxisType(const AxisType &rhs) :
  Frame(rhs) 
{
  ID = rhs.ID;
  if (rhs.Name!=NULL) {
    Name = new NAME[ natom ];
    memcpy(Name, rhs.Name, natom * sizeof(NAME));
  } else
    Name=NULL;
  origin[0] = rhs.origin[0];
  origin[1] = rhs.origin[1];
  origin[2] = rhs.origin[2];
  memcpy(R, rhs.R, 9 * sizeof(double));
  HbondAtom[0] = rhs.HbondAtom[0];
  HbondAtom[1] = rhs.HbondAtom[1];
  HbondAtom[2] = rhs.HbondAtom[2];
  // Since HbondCoord contains memory addresses, it must be updated
  // relative to this instances X
  HbondCoord[0] = X + (HbondAtom[0]*3);
  HbondCoord[1] = X + (HbondAtom[1]*3);
  HbondCoord[2] = X + (HbondAtom[2]*3);
  residue_number = rhs.residue_number;
  second_resnum = rhs.second_resnum;
}

// Assignment Operator
// NOTE: Only allocate/deallocate when natom>maxnatom?
AxisType &AxisType::operator=(const AxisType &rhs) {
  // Check for self assignment
  if ( this == &rhs ) return *this;

  // Base class assignment
  Frame::operator=(rhs);

  // Deallocate
  if (Name!=NULL) delete[] Name;
  Name = NULL;

  // Allocate and copy
  ID = rhs.ID;
  if (rhs.Name!=NULL) {
    Name = new NAME[ natom ];
    memcpy(Name, rhs.Name, natom * sizeof(NAME));
  }
  origin[0] = rhs.origin[0];
  origin[1] = rhs.origin[1];
  origin[2] = rhs.origin[2];
  memcpy(R, rhs.R, 9 * sizeof(double));
  HbondAtom[0] = rhs.HbondAtom[0];
  HbondAtom[1] = rhs.HbondAtom[1];
  HbondAtom[2] = rhs.HbondAtom[2];
  // Since HbondCoord contains memory addresses, it must be updated
  // relative to this instances X
  HbondCoord[0] = X + (HbondAtom[0]*3);
  HbondCoord[1] = X + (HbondAtom[1]*3);
  HbondCoord[2] = X + (HbondAtom[2]*3);
  residue_number = rhs.residue_number;
  second_resnum = rhs.second_resnum;

  // Return *this
  return *this;
} 

// DESTRUCTOR
AxisType::~AxisType() {
  if (Name!=NULL) delete[] Name;
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
AxisType::NAbaseType AxisType::ID_base(char *resname) {
  //mprintf("DBG:\tNAresname [%s]\n",resname);
  if (resname==NULL) return UNKNOWN_BASE;
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
    if (strncmp(resname,"ADE",3)==0) return ADE;
    if (strncmp(resname,"CYT",3)==0) return CYT;
    if (strncmp(resname,"GUA",3)==0) return GUA;
    if (strncmp(resname,"THY",3)==0) return THY;
    if (strncmp(resname,"URA",3)==0) return URA;
    if (resname[0]=='A') return ADE;
    if (resname[0]=='C') return CYT;
    if (resname[0]=='G') return GUA;
    if (resname[0]=='T') return THY;
    if (resname[0]=='U') return URA;
  } 
  return UNKNOWN_BASE;
}

// AxisType::AllocAxis()
/** Allocate mem for coords and names. Replace any existing coords / names.
  */
int AxisType::AllocAxis(int natomIn) {
  if (X!=NULL) delete[] X;
  if (Name!=NULL) delete[] Name;
  natom = natomIn;
  maxnatom = natom;
  N = natom * 3;
  X = new double[ N ];
  if (X==NULL) return 1;
  Name = new NAME[ natom ]; 
  if (Name==NULL) {delete[] X; return 1;}
  return 0;
}

// ------------------------- PUBLIC FUNCTIONS ---------------------------------
// AxisType::BaseName()
char *AxisType::BaseName() {
  sprintf(basename_num,"%i:%s",residue_number+1,NAbaseName[ID]);
  return basename_num; 
}

// AxisType::ResName()
char *AxisType::ResName() {
  return (char*)NAbaseName[ID];
}

// AxisType::AtomNameIs()
bool AxisType::AtomNameIs(int atom, char *nameIn) {
  if (atom<0 || atom>=natom) return false;
  if (strcmp(Name[atom], nameIn)==0) return true;
  return false;
}

// AxisType::AtomName()
char *AxisType::AtomName(int atom) {
  if (atom<0 || atom>=natom) return NULL;
  return (Name[atom]);
}

// AxisType::PrintAtomNames()
void AxisType::PrintAtomNames() {
  for (int atom = 0; atom < natom; atom++)
    mprintf(" %s",Name[atom]);
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

// AxisType::SetAxisFromMask()
/** Set coordinates (and names if debugging) of this axis based on the
  * given axis and atom mask.
  * \param AxisIn AxisType to set from.
  * \param Mask AtomMask containing atoms to keep from AxisIn.
  */
void AxisType::SetAxisFromMask(AxisType &AxisIn, AtomMask &Mask) {
  natom = Mask.Nselected;
  N     = natom * 3;
  if (natom > maxnatom) {
    delete[] X;
    X = new double[ N ];
#   ifdef NASTRUCTDEBUG
    // Need names when debugging
    if (AxisIn.Name!=NULL) {
      delete[] Name;
      Name = new NAME[ natom ];
    }
#   endif
    maxnatom = natom;
  }
  double *newX = X;
  for (int i = 0; i < natom; i++) {
    int oldatom = Mask.Selected[i];
    int oldatom3 = oldatom * 3;
    memcpy(newX, AxisIn.X + oldatom3, COORDSIZE);
    newX += 3;
#   ifdef NASTRUCTDEBUG
    strcpy(Name[i], AxisIn.Name[oldatom]);
#   endif
  }
}

// AxisType::StoreRotMatrix()
/** Store the rotation matrix and origin coordinates associated with
  * this base.
  */
void AxisType::StoreRotMatrix(double *RotMatrix, double *originIn) {
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
AxisType::RefReturn AxisType::SetRefCoord(AmberParm *currentParm, int resnum, 
                                          AtomMask &parmMask, AtomMask &fitMask,
                                          AxisType::NAbaseType customBaseType) 
{
  std::map<int,int> BaseMap;
  std::map<int,int>::iterator atom;
  std::vector<double> BaseRefCoords;
  std::vector<int> BaseHbonds;
  NAME *BaseRefNames = NULL;
  int Nbaseatom = 0;

  // First, identify the base
  if (customBaseType==UNKNOWN_BASE) 
    ID = ID_base(currentParm->ResidueName(resnum));
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
      BaseRefNames = new NAME[ ADENATOM ];
      BaseHbonds.reserve(ADENATOM);
      for (int i = 0; i < ADENATOM; i++) {
        BaseRefCoords.push_back( ADEcoords[i][0] );
        BaseRefCoords.push_back( ADEcoords[i][1] );
        BaseRefCoords.push_back( ADEcoords[i][2] );
        strcpy(BaseRefNames[i], ADEnames[i]);
        BaseHbonds.push_back( ADEhbonds[i] );
      }
      break;
    case CYT :
      if (AllocAxis(CYTNATOM)) return NA_ERROR;
      Nbaseatom = CYTNATOM;
      BaseRefCoords.reserve(CYTNATOM);
      BaseRefNames = new NAME[ CYTNATOM ];
      BaseHbonds.reserve(CYTNATOM);
      for (int i = 0; i < CYTNATOM; i++) {
        BaseRefCoords.push_back( CYTcoords[i][0] );
        BaseRefCoords.push_back( CYTcoords[i][1] );
        BaseRefCoords.push_back( CYTcoords[i][2] );
        strcpy(BaseRefNames[i], CYTnames[i]);
        BaseHbonds.push_back( CYThbonds[i] );
      }
      break;
    case GUA :
      if (AllocAxis(GUANATOM)) return NA_ERROR;
      Nbaseatom = GUANATOM;
      BaseRefCoords.reserve(GUANATOM);
      BaseRefNames = new NAME[ GUANATOM ];
      BaseHbonds.reserve(GUANATOM);
      for (int i = 0; i < GUANATOM; i++) {
        BaseRefCoords.push_back( GUAcoords[i][0] );
        BaseRefCoords.push_back( GUAcoords[i][1] );
        BaseRefCoords.push_back( GUAcoords[i][2] );
        strcpy(BaseRefNames[i], GUAnames[i]);
        BaseHbonds.push_back( GUAhbonds[i] );
      }
      break;
    case THY :
      if (AllocAxis(THYNATOM)) return NA_ERROR;
      Nbaseatom = THYNATOM;
      BaseRefCoords.reserve(THYNATOM);
      BaseRefNames = new NAME[ THYNATOM ];
      BaseHbonds.reserve(THYNATOM);
      for (int i = 0; i < THYNATOM; i++) {
        BaseRefCoords.push_back( THYcoords[i][0] );
        BaseRefCoords.push_back( THYcoords[i][1] );
        BaseRefCoords.push_back( THYcoords[i][2] );
        strcpy(BaseRefNames[i], THYnames[i]);
        BaseHbonds.push_back( THYhbonds[i] );
      }
      break;
    case URA :
      if (AllocAxis(URANATOM)) return NA_ERROR;
      Nbaseatom = URANATOM;
      BaseRefCoords.reserve(URANATOM);
      BaseRefNames = new NAME[ URANATOM ];
      BaseHbonds.reserve(URANATOM);
      for (int i = 0; i < URANATOM; i++) {
        BaseRefCoords.push_back( URAcoords[i][0] );
        BaseRefCoords.push_back( URAcoords[i][1] );
        BaseRefCoords.push_back( URAcoords[i][2] );
        strcpy(BaseRefNames[i], URAnames[i]);
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
    int parmAtom = currentParm->FindAtomInResidue(resnum, (char*)BaseRefNames[ref]);
    if (parmAtom<0) {
      mprinterr("Error: Ref Atom [%s] not found in NA base [%s].\n",
                BaseRefNames[ref], currentParm->ResidueName(resnum));
      return NA_ERROR;
    } else {
      BaseMap.insert( std::pair<int,int>(parmAtom, ref) );
    }
  }
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
      HbondCoord[hb_index] = X + coord3;
      HbondAtom[hb_index] = coord;
    }
    X[coord3++] = BaseRefCoords[refcoord  ]; 
    X[coord3++] = BaseRefCoords[refcoord+1]; 
    X[coord3++] = BaseRefCoords[refcoord+2]; 
    strcpy(Name[coord],BaseRefNames[refatom]);
    parmMask.AddAtom( (*atom).first );
    // Will this atom be used for RMS fitting?
    switch (ID) {
      case ADE:
      case GUA:
        if (strcmp(BaseRefNames[refatom],"N9  ")==0 ||
            strcmp(BaseRefNames[refatom],"C8  ")==0 ||
            strcmp(BaseRefNames[refatom],"N7  ")==0 ||
            strcmp(BaseRefNames[refatom],"C5  ")==0 ||
            strcmp(BaseRefNames[refatom],"C6  ")==0 ||
            strcmp(BaseRefNames[refatom],"N1  ")==0 ||
            strcmp(BaseRefNames[refatom],"C2  ")==0 ||
            strcmp(BaseRefNames[refatom],"N3  ")==0 ||
            strcmp(BaseRefNames[refatom],"C4  ")==0   ) 
          fitMask.AddAtom( coord );
          break;
      case CYT:
      case THY:
      case URA:
        if (strcmp(BaseRefNames[refatom],"N1  ")==0 ||
            strcmp(BaseRefNames[refatom],"C2  ")==0 ||
            strcmp(BaseRefNames[refatom],"N3  ")==0 ||
            strcmp(BaseRefNames[refatom],"C4  ")==0 ||
            strcmp(BaseRefNames[refatom],"C5  ")==0 ||
            strcmp(BaseRefNames[refatom],"C6  ")==0   )
          fitMask.AddAtom( coord );
          break;
      case UNKNOWN_BASE: return NA_ERROR;
    }
    ++coord;
  }
  if (BaseRefNames!=NULL) delete[] BaseRefNames;
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

#ifdef NASTRUCTDEBUG
// DEBUG
// AxisType::WritePDB()
// Write coordinates to file in PDB format.
void AxisType::WritePDB(CpptrajFile &outfile, int resnum, char *resname, int *atom) {
  char buffer[82];
  int i3=0;
  for (int i=0; i<natom; i++) {
    pdb_write_ATOM(buffer,PDBATOM,(*atom)+i+1,Name[i],resname,'X',resnum+1,
                   X[i3],X[i3+1],X[i3+2],1.0,0.0,(char*)"\0",false);
    outfile.IO->Write(buffer,sizeof(char),strlen(buffer));
    i3+=3;
  }
  (*atom) += natom;
}
// AxisType::WriteAxesPDB()
void AxisType::WriteAxesPDB(CpptrajFile &outfile, int resnum, char *resname, int *atom) {
  char buffer[82];
  double dx,dy,dz;
  // Origin
  pdb_write_ATOM(buffer, PDBATOM, (*atom)+1, (char*)"Orig", resname, 'X', resnum+1,
                 origin[0],origin[1],origin[2],1.0,0.0,(char*)"\0",false);
  outfile.IO->Write(buffer,sizeof(char),strlen(buffer));
  // X vector
  dx = origin[0] + R[0];
  dy = origin[1] + R[3];
  dz = origin[2] + R[6];
  pdb_write_ATOM(buffer, PDBATOM, (*atom)+2, (char*)"X", resname, 'X', resnum+1,
                 dx,dy,dz,1.0,0.0,(char*)"\0",false);
  outfile.IO->Write(buffer,sizeof(char),strlen(buffer));
  // Y vector
  dx = origin[0] + R[1];
  dy = origin[1] + R[4];
  dz = origin[2] + R[7];
  pdb_write_ATOM(buffer, PDBATOM, (*atom)+3, (char*)"Y", resname, 'X', resnum+1,
                 dx,dy,dz,1.0,0.0,(char*)"\0",false);
  outfile.IO->Write(buffer,sizeof(char),strlen(buffer));
  // Z vector
  dx = origin[0] + R[2];
  dy = origin[1] + R[5];
  dz = origin[2] + R[8];
  pdb_write_ATOM(buffer, PDBATOM, (*atom)+4, (char*)"Z", resname, 'X', resnum+1,
                 dx,dy,dz,1.0,0.0,(char*)"\0",false);
  outfile.IO->Write(buffer,sizeof(char),strlen(buffer));

  (*atom) = (*atom) + 4;
}
#endif
