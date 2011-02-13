// AxisType
#include <cstdlib>
#include <cstring>
#include "AxisType.h"
#include "CpptrajStdio.h"
// DEBUG
#include "PDBfileRoutines.h"
/*
 * ID_base()
 * Return a number indicating if this is a NA base. If not recognized, 
 * return -1. 
 * NOTE: Currently based only on amber residue names. Will not recognize
 * non-standard bases.
 */
NAbaseType ID_base(char *resname) {
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

// ------------------------- AXISTYPE FUNCTIONS -------------------------------
// CONSTRUCTOR
// Mostly empty - Default Frame constructor should init everything
AxisType::AxisType() {
  ID = UNKNOWN_BASE;
  Name = NULL;
  maxAtom = 0;
  R[0]=0.0; R[1]=0.0; R[2]=0.0;
  R[3]=0.0; R[4]=0.0; R[5]=0.0;
  R[6]=0.0; R[7]=0.0; R[8]=0.0;
}

// CONSTRUCTOR - Only allocate coords
AxisType::AxisType(int natomIn) {
  natom = natomIn;
  maxAtom = natom;
  N = natom * 3;
  X = (double*) malloc(N * sizeof(double));
  ID = UNKNOWN_BASE;
  Name = NULL;
  R[0]=0.0; R[1]=0.0; R[2]=0.0;
  R[3]=0.0; R[4]=0.0; R[5]=0.0;
  R[6]=0.0; R[7]=0.0; R[8]=0.0;
}

// DESTRUCTOR
AxisType::~AxisType() {
  if (Name!=NULL) free(Name);
}

// ------------------------- VARIABLE FUNCTIONS -------------------------------
/*
 * AxisType::RX()
 * Put the X unit vector from rotation matrix into Vin.
 */
void AxisType::RX(double Vin[3]) {
  Vin[0] = R[0];
  Vin[1] = R[3];
  Vin[2] = R[6];
}
/*
 * AxisType::RY()
 * Put the Y unit vector from rotation matrix into Vin.
 */
void AxisType::RY(double Vin[3]) {
  Vin[0] = R[1];
  Vin[1] = R[4];
  Vin[2] = R[7];
}
/*
 * AxisType::RZ()
 * Put the Z unit vector from rotation matrix into Vin.
 */
void AxisType::RZ(double Vin[3]) {
  Vin[0] = R[2];
  Vin[1] = R[5];
  Vin[2] = R[8];
}
/*
 * AxisType::Origin()
 * Return pointer to the origin coords.
 */
double *AxisType::Origin() {
  return (X+9);
}

// ------------------------- PRIVATE FUNCTIONS --------------------------------
// Base names corresponding to NAbaseType
// UNKNOWN_BASE, DA, DT, DG, DC, RA, RC, RG, RU
const char AxisType::NAbaseName[9][5]={"UNK","DA","DT","DG","DC","RA","RC","RG","RU" };

/*
 * AxisType::AllocAxis()
 * Allocate mem for coords and names
 */
int AxisType::AllocAxis(int natomIn) {
  natom = natomIn;
  maxAtom = natom;
  N = natom * 3;
  X = (double*) malloc(N * sizeof(double));
  if (X==NULL) return 1;
  Name = (AmberParm::NAME*) malloc( natom * sizeof(AmberParm::NAME));
  if (Name==NULL) {free(X); return 1;}
  return 0;
}

// ------------------------- PUBLIC FUNCTIONS ---------------------------------
/*
 * AxisType::BaseName()
 */
char *AxisType::BaseName() {
  return (char*)NAbaseName[ID];
}

/*
 * AxisType::AtomIndex()
 * Return atom index of given atom name in axis. 
 * Return -1 if name not found.
 */
/*int AxisType::AtomIndex(char *atomname) {
  for (int atom=0; atom < natom; atom++) 
    if (strcmp(atomname, Name[atom])==0) return atom;
  return -1;
}*/
bool AxisType::AtomNameIs(int atom, char *nameIn) {
  if (atom<0 || atom>=natom) return false;
  if (strcmp(Name[atom], nameIn)==0) return true;
  return false;
}

char *AxisType::AtomName(int atom) {
  if (atom<0 || atom>=natom) return NULL;
  return (Name[atom]);
}

/*
 * AxisType::SetFromFrame()
 * Set AxisType coords only from AxisIn. Reallocate memory only if
 * AxisIn->natom > maxatom.
 */
void AxisType::SetFromFrame(AxisType *AxisIn) {
  this->natom = AxisIn->natom;
  this->N     = AxisIn->N;
  if (AxisIn->natom > this->maxAtom) { 
    X = (double*) realloc(X, this->N * sizeof(double));
    // Dont care about Name here
    this->maxAtom = this->natom;
  }
  for (int i=0; i < this->N; i++)
    this->X[i] = AxisIn->X[i];
}

/*
 * AxisType::StoreRotMatrix()
 */
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
    
/*
 * AxisType::SetPrincipalAxes()
 * Set this up as principal axes. Allocate memory if necessary.
 * NOTE: Check for mem error?
 */
void AxisType::SetPrincipalAxes() {
  if (X==NULL) {
    AllocAxis(4);
    strcpy(Name[0],"X");
    strcpy(Name[1],"Y");
    strcpy(Name[2],"Z");
    strcpy(Name[3],"Orig");
  }
  X[0]=1.0; X[3]=0.0; X[6]=0.0; X[9 ]=0.0; 
  X[1]=0.0; X[4]=1.0; X[7]=0.0; X[10]=0.0; 
  X[2]=0.0; X[5]=0.0; X[8]=1.0; X[11]=0.0;
}

/*
 * AxisType::SetRefCoord()
 * Set NA residue reference coordinates for given NA base name.
 */
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
}

/*
 * AxisType::FlipYZ
 * Flip the Z and Y axes.
 * Equiv. to rotation around the X axis.
 * Done for antiparallel stranded DNA
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

/*
 * AxisType::FlipXY
 * Flip the X and Y axes.
 * Equiv. to rotation around the Z axis.
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

/*
 * DEBUG: AxisToPDB
 */
void AxisType::WritePDB(PtrajFile *outfile, int resnum, char *resname, int *atom) {
  char buffer[82];
  int i3=0;
  for (int i=0; i<natom; i++) {
    pdb_write_ATOM(buffer,"ATOM",(*atom)+i,Name[i],resname,'X',resnum+1,
                   X[i3],X[i3+1],X[i3+2],1.0,0.0,(char*)"\0");
    outfile->IO->Write(buffer,sizeof(char),strlen(buffer));
    i3+=3;
  }
  (*atom) += natom;
}

