// AxisType
#include <map>
#include <cstring>
#include "AxisType.h"
#include "CpptrajStdio.h"
#ifdef NASTRUCTDEBUG
// DEBUG
#  include "PDBfileRoutines.h"
#endif

// ---------- NA REFERENCE BASE ATOM NAMES AND COORDS --------------------------
#define ADENATOM 11
const NAME AxisType::ADEnames[ADENATOM] = {
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
#define CYTNATOM 9
const NAME AxisType::CYTnames[CYTNATOM] = {
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
#define GUANATOM 12
const NAME AxisType::GUAnames[GUANATOM] = {
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
#define THYNATOM 10
const NAME AxisType::THYnames[THYNATOM] = {
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
#define URANATOM 9
const NAME AxisType::URAnames[URANATOM] = {
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

// ------------------------- AXISTYPE FUNCTIONS -------------------------------
// CONSTRUCTOR
// Mostly empty - Default Frame constructor should init everything else
AxisType::AxisType() {
  ID = UNKNOWN_BASE;
  Name = NULL;
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
}

// COPY CONSTRUCTOR
AxisType::AxisType(const AxisType &rhs) :
  Frame(rhs) 
{
  ID = rhs.ID;
  if (rhs.Name!=NULL) {
    Name = new NAME[ natom ];
    memcpy(Name, rhs.Name, natom * sizeof(NAME));
  } else
    Name=NULL;
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
}

// Assignment Operator
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

  // Return *this
  return *this;
} 

// CONSTRUCTOR - Only allocate coords
/*AxisType::AxisType(int natomIn) {
  natom = natomIn;
  maxnatom = natom;
  N = natom * 3;
  X = new double[ N ];
  ID = UNKNOWN_BASE;
  Name = NULL;
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
}*/

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
// AxisType::Origin()
/** Return pointer to the origin coords.  */
double *AxisType::Origin() {
  return (X+9);
}

// ------------------------- PRIVATE FUNCTIONS --------------------------------
// Base names corresponding to NAbaseType
// UNKNOWN_BASE, DA, DT, DG, DC, RA, RC, RG, RU
const char AxisType::NAbaseName[9][5]={"UNK","DA","DT","DG","DC","RA","RC","RG","RU" };

// AxisType::ID_base()
/** Return a number indicating if this is a NA base. If not recognized, 
  * return -1. 
  * NOTE: Currently based only on amber residue names. Will not recognize
  * non-standard bases.
  */
AxisType::NAbaseType AxisType::ID_base(char *resname) {
  char *resID;
  bool isDNA = true; // Assume DNA unless noted otherwise 
  //mprintf("        [%s]\n",resname);
  if (resname==NULL) return UNKNOWN_BASE;
  resID = resname;
  // If residue name begins with D, assume AMBER DNA residue
  if (resID[0]=='D') resID++;
  if (resID[0]=='R') {resID++; isDNA=false;}
  if (isDNA) {
    switch (resID[0]) {
      case 'A': return DA;
      case 'C': return DC;
      case 'G': return DG;
      case 'T': return DT;
    }
  } else {
    switch (resID[0]) {
      case 'A': return RA;
      case 'C': return RC;
      case 'G': return RG;
      case 'U': return RU;
    }
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
  Name = new NAME[ natom ]; // (NAME*) malloc( natom * sizeof(NAME));
  if (Name==NULL) {delete[] X; return 1;}
  return 0;
}

// ------------------------- PUBLIC FUNCTIONS ---------------------------------
// AxisType::BaseName()
char *AxisType::BaseName() {
  return (char*)NAbaseName[ID];
}

/* AxisType::AtomIndex()
 * Return atom index of given atom name in axis. 
 * Return -1 if name not found.
 */
/*int AxisType::AtomIndex(char *atomname) {
  for (int atom=0; atom < natom; atom++) 
    if (strcmp(atomname, Name[atom])==0) return atom;
  return -1;
}*/

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

// AxisType::SetFromFrame()
/** Set AxisType coords only from AxisIn. Reallocate memory only if
  * AxisIn->natom > maxnatom.
  */
void AxisType::SetFromFrame(AxisType *AxisIn) {
  natom = AxisIn->natom;
  N     = AxisIn->N;
  if (AxisIn->natom > maxnatom) {
    delete[] X;
    X = new double[ N ];
    // Dont care about Name here
#   ifdef NASTRUCTDEBUG
    // If debugging we need the names
    if (AxisIn->Name!=NULL) {
      delete[] Name;
      Name = new NAME[ natom ];
    }
#   endif
    maxnatom = natom;
  }
  memcpy(X, AxisIn->X, N * sizeof(double));
# ifdef NASTRUCTDEBUG
  if (AxisIn->Name!=NULL)
    memcpy(Name, AxisIn->Name, natom * sizeof(NAME));
# endif
}

// AxisType::StoreRotMatrix()
void AxisType::StoreRotMatrix(double *RotMatrix) {
  R[0] = RotMatrix[0];
  R[1] = RotMatrix[1];
  R[2] = RotMatrix[2];
  R[3] = RotMatrix[3];
  R[4] = RotMatrix[4];
  R[5] = RotMatrix[5];
  R[6] = RotMatrix[6];
  R[7] = RotMatrix[7];
  R[8] = RotMatrix[8];
}
    
// AxisType::SetPrincipalAxes()
/** Set this up as principal axes. Allocate memory if necessary.
  */
// NOTE: Check for mem error?
void AxisType::SetPrincipalAxes() {
  if (X==NULL) {
    AllocAxis(4);
    strcpy(Name[0],"X   ");
    strcpy(Name[1],"Y   ");
    strcpy(Name[2],"Z   ");
    strcpy(Name[3],"Orig");
  }
  X[0]=1.0; X[3]=0.0; X[6]=0.0; X[9 ]=0.0; 
  X[1]=0.0; X[4]=1.0; X[7]=0.0; X[10]=0.0; 
  X[2]=0.0; X[5]=0.0; X[8]=1.0; X[11]=0.0;
}

/* AxisType::SetRefCoord()
 * Set NA residue reference coordinates for given NA base name.
 */
/*
int AxisType::SetRefCoord(char *resname) {
  // First, identify the base
  ID = ID_base(resname);
  switch (ID) {
    case DA :
    case RA :
      if (AllocAxis(11)) return 1;
      X[0 ]=-2.479000; X[1 ]= 5.346000; X[2 ]= 0.000000; strcpy(Name[0 ],"C1' ");
      X[3 ]=-1.291000; X[4 ]= 4.498000; X[5 ]= 0.000000; strcpy(Name[1 ],"N9  ");
      X[6 ]= 0.024000; X[7 ]= 4.897000; X[8 ]= 0.000000; strcpy(Name[2 ],"C8  ");
      X[9 ]= 0.877000; X[10]= 3.902000; X[11]= 0.000000; strcpy(Name[3 ],"N7  ");
      X[12]= 0.071000; X[13]= 2.771000; X[14]= 0.000000; strcpy(Name[4 ],"C5  ");
      X[15]= 0.369000; X[16]= 1.398000; X[17]= 0.000000; strcpy(Name[5 ],"C6  ");
      X[18]= 1.611000; X[19]= 0.909000; X[20]= 0.000000; strcpy(Name[6 ],"N6  ");
      X[21]=-0.668000; X[22]= 0.532000; X[23]= 0.000000; strcpy(Name[7 ],"N1  ");
      X[24]=-1.912000; X[25]= 1.023000; X[26]= 0.000000; strcpy(Name[8 ],"C2  ");
      X[27]=-2.320000; X[28]= 2.290000; X[29]= 0.000000; strcpy(Name[9 ],"N3  ");
      X[30]=-1.267000; X[31]= 3.124000; X[32]= 0.000000; strcpy(Name[10],"C4  ");
      break;
    case DC :
    case RC :
      if (AllocAxis(9)) return 1;
      X[0 ]=-2.477000; X[1 ]= 5.402000; X[2 ]= 0.000000; strcpy(Name[0],"C1' ");
      X[3 ]=-1.285000; X[4 ]= 4.542000; X[5 ]= 0.000000; strcpy(Name[1],"N1  ");
      X[6 ]=-1.472000; X[7 ]= 3.158000; X[8 ]= 0.000000; strcpy(Name[2],"C2  ");
      X[9 ]=-2.628000; X[10]= 2.709000; X[11]= 0.000000; strcpy(Name[3],"O2  ");
      X[12]=-0.391000; X[13]= 2.344000; X[14]= 0.000000; strcpy(Name[4],"N3  ");
      X[15]= 0.837000; X[16]= 2.868000; X[17]= 0.000000; strcpy(Name[5],"C4  ");
      X[18]= 1.875000; X[19]= 2.027000; X[20]= 0.000000; strcpy(Name[6],"N4  ");
      X[21]= 1.056000; X[22]= 4.275000; X[23]= 0.000000; strcpy(Name[7],"C5  ");
      X[24]=-0.023000; X[25]= 5.068000; X[26]= 0.000000; strcpy(Name[8],"C6  ");
      break;
    case DG :
    case RG :
      if (AllocAxis(12)) return 1;
      X[0 ]=-2.477000; X[1 ]= 5.399000; X[2 ]= 0.000000; strcpy(Name[0],"C1' ");
      X[3 ]=-1.289000; X[4 ]= 4.551000; X[5 ]= 0.000000; strcpy(Name[1],"N9  ");
      X[6 ]= 0.023000; X[7 ]= 4.962000; X[8 ]= 0.000000; strcpy(Name[2],"C8  ");
      X[9 ]= 0.870000; X[10]= 3.969000; X[11]= 0.000000; strcpy(Name[3],"N7  ");
      X[12]= 0.071000; X[13]= 2.833000; X[14]= 0.000000; strcpy(Name[4],"C5  ");
      X[15]= 0.424000; X[16]= 1.460000; X[17]= 0.000000; strcpy(Name[5],"C6  ");
      X[18]= 1.554000; X[19]= 0.955000; X[20]= 0.000000; strcpy(Name[6],"O6  ");
      X[21]=-0.700000; X[22]= 0.641000; X[23]= 0.000000; strcpy(Name[7],"N1  ");
      X[24]=-1.999000; X[25]= 1.087000; X[26]= 0.000000; strcpy(Name[8],"C2  ");
      X[27]=-2.949000; X[28]= 0.139000; X[29]=-0.001000; strcpy(Name[9],"N2  ");
      X[30]=-2.342000; X[31]= 2.364000; X[32]= 0.001000; strcpy(Name[10],"N3  ");
      X[33]=-1.265000; X[34]= 3.177000; X[35]= 0.000000; strcpy(Name[11],"C4  ");
      break;
    case DT :
      if (AllocAxis(10)) return 1;
      X[0 ]=-2.481000; X[1 ]= 5.354000; X[2 ]= 0.000000; strcpy(Name[0],"C1' ");
      X[3 ]=-1.284000; X[4 ]= 4.500000; X[5 ]= 0.000000; strcpy(Name[1],"N1  ");
      X[6 ]=-1.462000; X[7 ]= 3.135000; X[8 ]= 0.000000; strcpy(Name[2],"C2  ");
      X[9 ]=-2.562000; X[10]= 2.608000; X[11]= 0.000000; strcpy(Name[3],"O2  ");
      X[12]=-0.298000; X[13]= 2.407000; X[14]= 0.000000; strcpy(Name[4],"N3  ");
      X[15]= 0.994000; X[16]= 2.897000; X[17]= 0.000000; strcpy(Name[5],"C4  ");
      X[18]= 1.944000; X[19]= 2.119000; X[20]= 0.000000; strcpy(Name[6],"O4  ");
      X[21]= 1.106000; X[22]= 4.338000; X[23]= 0.000000; strcpy(Name[7],"C5  ");
      X[24]= 2.466000; X[25]= 4.961000; X[26]= 0.001000; strcpy(Name[8],"C7  ");
      X[27]=-0.024000; X[28]= 5.057000; X[29]= 0.000000; strcpy(Name[9],"C6  ");
      break;
    case RU:
      if (AllocAxis(9)) return 1;
      X[0 ]=-2.481000; X[1 ]= 5.354000; X[2 ]= 0.000000; strcpy(Name[0],"C1' ");
      X[3 ]=-1.284000; X[4 ]= 4.500000; X[5 ]= 0.000000; strcpy(Name[1],"N1  ");
      X[6 ]=-1.462000; X[7 ]= 3.131000; X[8 ]= 0.000000; strcpy(Name[2],"C2  ");
      X[9 ]=-2.563000; X[10]= 2.608000; X[11]= 0.000000; strcpy(Name[3],"O2  ");
      X[12]=-0.302000; X[13]= 2.397000; X[14]= 0.000000; strcpy(Name[4],"N3  ");
      X[15]= 0.989000; X[16]= 2.884000; X[17]= 0.000000; strcpy(Name[5],"C4  ");
      X[18]= 1.935000; X[19]= 2.094000; X[20]=-0.001000; strcpy(Name[6],"O4  ");
      X[21]= 1.089000; X[22]= 4.311000; X[23]= 0.000000; strcpy(Name[7],"C5  ");
      X[24]=-0.024000; X[25]= 5.053000; X[26]= 0.000000; strcpy(Name[8],"C6  ");
      break;
    case UNKNOWN_BASE:
      mprintf("Warning: AxisType::SetRefCoord: Missing parameters for residue %s.\n",resname);
      return 1;
  }
  return 0;
}*/

// AxisType::SetRefCoord()
/** Set NA residue reference coordinates for given NA base. Ensure that
  * the atom ordering in the reference matches that in the given parm.
  * Set up a mask to select correct atoms in the given parm.
  * \return 0 on success, 1 on error.
  */
AxisType::RefReturn AxisType::SetRefCoord(AmberParm *currentParm, int resnum, 
                                          AtomMask &parmMask) 
{
  std::map<int,int> BaseMap;
  std::map<int,int>::iterator atom;
  std::vector<double> BaseRefCoords;
  std::vector<int> BaseHbonds;
  NAME *BaseRefNames = NULL;
  int Nbaseatom = 0;

  // First, identify the base
  ID = ID_base(currentParm->ResidueName(resnum));
  // If unknown exit now
  if (ID == UNKNOWN_BASE) return NA_UNKNOWN; 

  // Set up reference coord and name arrays based on the base
  switch (ID) {
    case DA :
    case RA :
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
    case DC :
    case RC :
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
    case DG :
    case RG :
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
    case DT :
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
    case RU:
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
    strcpy(Name[coord++],BaseRefNames[refatom]);
    parmMask.AddAtom( (*atom).first );
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
  double ox2, oy2, oz2;
    
  ox2 = X[9 ] + X[9 ];
  oy2 = X[10] + X[10];
  oz2 = X[11] + X[11];
  // Y axis  
  X[3 ] = ox2 - X[3 ];
  X[4 ] = oy2 - X[4 ];
  X[5 ] = oz2 - X[5 ];
  // Z axis  
  X[6 ] = ox2 - X[6 ];
  X[7 ] = oy2 - X[7 ];
  X[8 ] = oz2 - X[8 ]; 

  R[1] = -R[1]; // -Yx
  R[2] = -R[2]; // -Zx
  R[4] = -R[4]; // -Yy
  R[5] = -R[5]; // -Zy
  R[7] = -R[7]; // -Yz
  R[8] = -R[8]; // -Zz
}

// AxisType::FlipXY
/** Flip the X and Y axes. Equivalent to rotation around the Z axis.
  * Done for parallel stranded DNA.
  */
void AxisType::FlipXY() {
double ox2, oy2, oz2;

  ox2 = X[9 ] + X[9 ];
  oy2 = X[10] + X[10];
  oz2 = X[11] + X[11];
  // X axis
  X[0 ] = ox2 - X[0 ];
  X[1 ] = ox2 - X[1 ];
  X[2 ] = ox2 - X[2 ];
  // Y axis  
  X[3 ] = ox2 - X[3 ];
  X[4 ] = oy2 - X[4 ];
  X[5 ] = oz2 - X[5 ];

  R[0] = -R[0]; // -Xx
  R[1] = -R[1]; // -Yx
  R[3] = -R[3]; // -Xy
  R[4] = -R[4]; // -Yy
  R[6] = -R[6]; // -Xz
  R[7] = -R[7]; // -Yz
}

#ifdef NASTRUCTDEBUG
// DEBUG: AxisToPDB
void AxisType::WritePDB(CpptrajFile *outfile, int resnum, char *resname, int *atom) {
  char buffer[82];
  int i3=0;
  for (int i=0; i<natom; i++) {
    pdb_write_ATOM(buffer,PDBATOM,(*atom)+i,Name[i],resname,'X',resnum+1,
                   X[i3],X[i3+1],X[i3+2],1.0,0.0,(char*)"\0",false);
    outfile->IO->Write(buffer,sizeof(char),strlen(buffer));
    i3+=3;
  }
  (*atom) += natom;
}
#endif
