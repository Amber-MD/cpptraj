// NAstruct 
#include <cstdlib>
#include <cstring>
#include "Action_NAstruct.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "vectormath.h"
// DEBUG
#include "PDBfileRoutines.h"
#include <cmath>

// CONSTRUCTOR
NAstruct::NAstruct() {
  //fprintf(stderr,"NAstruct Con\n");
  BasePair=NULL;
  Nbp=0;
  Nbases=0;
  HBcut2 = 12.25; // 3.5^2
  Ocut2 = 6.25;   // 2.5^2
  Axes1=NULL;
  Axes2=NULL;
  resRange=NULL;
  outFilename=NULL;
} 

// DESTRUCTOR
NAstruct::~NAstruct() { 
  if (resRange!=NULL) delete resRange;
  ClearLists();
  if (BasePair!=NULL) free(BasePair);
  if (Axes1!=NULL) delete Axes1;
  if (Axes2!=NULL) delete Axes2;
}

// Base names corresponding to NAbaseType
// UNKNOWN_BASE, DA, DT, DG, DC, RA, RC, RG, RU
const char NAstruct::NAbaseName[9][5]={"UNK","DA","DT","DG","DC","RA","RC","RG","RU" };

// -------------------------- MEMORY FUNCTIONS --------------------------------
/*
 * NAstruct::AllocAxis()
 * Allocate memory for axis type
 */
NAstruct::AxisType *NAstruct::AllocAxis(int N) {
  AxisType *axis;
  axis = (AxisType *) malloc(sizeof(AxisType));
  if (axis==NULL) {
    mprintf("Error: NAstruct::AllocAxis: Could not allocate axis memory for %i atoms.\n",N);
  }
  axis->F= new Frame(N,NULL);
  axis->Name = (AmberParm::NAME*) malloc(N * sizeof(AmberParm::NAME));
  axis->ID=UNKNOWN_BASE;
  return axis;
}

/*  
 * NAstruct::AxisCopy()
 * Return a copy of the given axis.
 */
NAstruct::AxisType *NAstruct::AxisCopy( AxisType *axisIn ) {
  AxisType *axis;
  axis = (AxisType *) malloc(sizeof(AxisType));
  if (axis==NULL) {
    mprintf("Error: NAstruct::AxisCopy: Could not allocate axis memory.\n");
  }
  axis->F = axisIn->F->Copy();
  axis->Name = (AmberParm::NAME*) malloc(axisIn->F->natom * sizeof(AmberParm::NAME));
  for (int i=0; i<axisIn->F->natom; i++)
    strcpy(axis->Name[i], axisIn->Name[i]);
  axis->ID=axisIn->ID;
  return axis;
}

void NAstruct::AxisToFrame( AxisType *axis, Frame *frame ) {
  for (int i=0; i < 12; i++)
    frame->X[i] = axis->F->X[i];
}

/*
 * NAstruct::FreeAxis()
 * Free Memory allocated by axis type
 */
void NAstruct::FreeAxis( AxisType *axis ) {
  if (axis->F!=NULL) delete axis->F;
  if (axis->Name!=NULL) free(axis->Name);
  free(axis);
}

/*
 * NAstruct::ClearLists()
 * Clear all parm-dependent lists
 */
void NAstruct::ClearLists() {
  while (!RefCoords.empty()) {
    FreeAxis( RefCoords.back() );
    RefCoords.pop_back();
  }
  while (!BaseAxes.empty()) {
    FreeAxis( BaseAxes.back() );
    BaseAxes.pop_back();
  }
  while (!BasePairAxes.empty()) {
    FreeAxis( BasePairAxes.back() );
    BasePairAxes.pop_back();
  }
  while (!ExpMasks.empty()) {
    delete ExpMasks.back();
    ExpMasks.pop_back();
  }
  while (!ExpFrames.empty()) {
    delete ExpFrames.back();
    ExpFrames.pop_back();
  }
}

// ------------------------- PRIVATE FUNCTIONS --------------------------------
/*
 * NAstruct::ID_base()
 * Return a number indicating if this is a NA base. If not recognized, 
 * return -1. 
 * NOTE: Currently based only on amber residue names. Will not recognize
 * non-standard bases.
 */
NAstruct::NAbaseType NAstruct::ID_base(char *resname) {
  //mprintf("        [%s]\n",resname);
  if (resname==NULL) return UNKNOWN_BASE;
  if (resname[0]=='D') {
    switch (resname[1]) {
      case 'A': return DA;
      case 'C': return DC;
      case 'G': return DG;
      case 'T': return DT;
    }
  } else if (resname[0]=='R') {
    switch (resname[1]) {
      case 'A': return RA;
      case 'C': return RC;
      case 'G': return RG;
      case 'U': return RU;
    }
  }
  return UNKNOWN_BASE;
}

/* PREPROCESSOR MACRO
 * setPrincipalAxes
 * Set coords of double *X to principal Axes
 */
#define setPrincipalAxes(X) { \
  X[0]=1.0; X[3]=0.0; X[6]=0.0; X[9 ]=0.0; \
  X[1]=0.0; X[4]=1.0; X[7]=0.0; X[10]=0.0; \
  X[2]=0.0; X[5]=0.0; X[8]=1.0; X[11]=0.0; }

/*
 * NAstruct::principalAxes()
 * Set up a coordinate type with vectors of size 1.0 pointing along the x,
 * y, and z axes.
 */
NAstruct::AxisType *NAstruct::principalAxes() {
  AxisType *axis;

  axis = AllocAxis(4);
//  F->X[0]=1.0; F->X[3]=0.0; F->X[6]=0.0; F->X[9 ]=0.0;
//  F->X[1]=0.0; F->X[4]=1.0; F->X[7]=0.0; F->X[10]=0.0;
//  F->X[2]=0.0; F->X[5]=0.0; F->X[8]=1.0; F->X[11]=0.0;
  setPrincipalAxes(axis->F->X);
  strcpy(axis->Name[0],"X");
  strcpy(axis->Name[1],"Y");
  strcpy(axis->Name[2],"Z");
  strcpy(axis->Name[3],"Orig");
   
  return axis;
}


/*
 * NAstruct::getRefCoords()
 * Allocate and set the coordinates to a standard ref. for the given base type.
 * Also set target atom names. Atom Names should be 4 chars long, AMBER.
 * Coords taken from Olson et al. JMB (2001) 313, 229-237.
 */
NAstruct::AxisType *NAstruct::getRefCoords( NAbaseType btype) {
  AxisType *axis;
  Frame *AF;

  axis=NULL; 
  AF=NULL;
  switch (btype) {
    case DA :
      if ( (axis = AllocAxis(11))==NULL ) return NULL;
      axis->ID = DA;
      AF = axis->F;
      AF->X[0 ]=-2.479000; AF->X[1 ]= 5.346000; AF->X[2 ]= 0.000000; strcpy(axis->Name[0 ],"C1' ");
      AF->X[3 ]=-1.291000; AF->X[4 ]= 4.498000; AF->X[5 ]= 0.000000; strcpy(axis->Name[1 ],"N9  ");
      AF->X[6 ]= 0.024000; AF->X[7 ]= 4.897000; AF->X[8 ]= 0.000000; strcpy(axis->Name[2 ],"C8  ");
      AF->X[9 ]= 0.877000; AF->X[10]= 3.902000; AF->X[11]= 0.000000; strcpy(axis->Name[3 ],"N7  ");
      AF->X[12]= 0.071000; AF->X[13]= 2.771000; AF->X[14]= 0.000000; strcpy(axis->Name[4 ],"C5  ");
      AF->X[15]= 0.369000; AF->X[16]= 1.398000; AF->X[17]= 0.000000; strcpy(axis->Name[5 ],"C6  ");
      AF->X[18]= 1.611000; AF->X[19]= 0.909000; AF->X[20]= 0.000000; strcpy(axis->Name[6 ],"N6  ");
      AF->X[21]=-0.668000; AF->X[22]= 0.532000; AF->X[23]= 0.000000; strcpy(axis->Name[7 ],"N1  ");
      AF->X[24]=-1.912000; AF->X[25]= 1.023000; AF->X[26]= 0.000000; strcpy(axis->Name[8 ],"C2  ");
      AF->X[27]=-2.320000; AF->X[28]= 2.290000; AF->X[29]= 0.000000; strcpy(axis->Name[9 ],"N3  ");
      AF->X[30]=-1.267000; AF->X[31]= 3.124000; AF->X[32]= 0.000000; strcpy(axis->Name[10],"C4  ");
      break;
    case DC :
      if ( (axis = AllocAxis(9))==NULL ) return NULL;
      axis->ID = DC;
      AF = axis->F;
      AF->X[0 ]=-2.477000; AF->X[1 ]= 5.402000; AF->X[2 ]= 0.000000; strcpy(axis->Name[0],"C1' ");
      AF->X[3 ]=-1.285000; AF->X[4 ]= 4.542000; AF->X[5 ]= 0.000000; strcpy(axis->Name[1],"N1  ");
      AF->X[6 ]=-1.472000; AF->X[7 ]= 3.158000; AF->X[8 ]= 0.000000; strcpy(axis->Name[2],"C2  ");
      AF->X[9 ]=-2.628000; AF->X[10]= 2.709000; AF->X[11]= 0.000000; strcpy(axis->Name[3],"O2  ");
      AF->X[12]=-0.391000; AF->X[13]= 2.344000; AF->X[14]= 0.000000; strcpy(axis->Name[4],"N3  ");
      AF->X[15]= 0.837000; AF->X[16]= 2.868000; AF->X[17]= 0.000000; strcpy(axis->Name[5],"C4  ");
      AF->X[18]= 1.875000; AF->X[19]= 2.027000; AF->X[20]= 0.000000; strcpy(axis->Name[6],"N4  ");
      AF->X[21]= 1.056000; AF->X[22]= 4.275000; AF->X[23]= 0.000000; strcpy(axis->Name[7],"C5  ");
      AF->X[24]=-0.023000; AF->X[25]= 5.068000; AF->X[26]= 0.000000; strcpy(axis->Name[8],"C6  ");
      break;
    case DG :
      if ( (axis = AllocAxis(12))==NULL ) return NULL;
      axis->ID = DG;
      AF = axis->F;
      AF->X[0 ]=-2.477000; AF->X[1 ]= 5.399000; AF->X[2 ]= 0.000000; strcpy(axis->Name[0],"C1' ");
      AF->X[3 ]=-1.289000; AF->X[4 ]= 4.551000; AF->X[5 ]= 0.000000; strcpy(axis->Name[1],"N9  ");
      AF->X[6 ]= 0.023000; AF->X[7 ]= 4.962000; AF->X[8 ]= 0.000000; strcpy(axis->Name[2],"C8  ");
      AF->X[9 ]= 0.870000; AF->X[10]= 3.969000; AF->X[11]= 0.000000; strcpy(axis->Name[3],"N7  ");
      AF->X[12]= 0.071000; AF->X[13]= 2.833000; AF->X[14]= 0.000000; strcpy(axis->Name[4],"C5  ");
      AF->X[15]= 0.424000; AF->X[16]= 1.460000; AF->X[17]= 0.000000; strcpy(axis->Name[5],"C6  ");
      AF->X[18]= 1.554000; AF->X[19]= 0.955000; AF->X[20]= 0.000000; strcpy(axis->Name[6],"O6  ");
      AF->X[21]=-0.700000; AF->X[22]= 0.641000; AF->X[23]= 0.000000; strcpy(axis->Name[7],"N1  ");
      AF->X[24]=-1.999000; AF->X[25]= 1.087000; AF->X[26]= 0.000000; strcpy(axis->Name[8],"C2  ");
      AF->X[27]=-2.949000; AF->X[28]= 0.139000; AF->X[29]=-0.001000; strcpy(axis->Name[9],"N2  ");
      AF->X[30]=-2.342000; AF->X[31]= 2.364000; AF->X[32]= 0.001000; strcpy(axis->Name[10],"N3  ");
      AF->X[33]=-1.265000; AF->X[34]= 3.177000; AF->X[35]= 0.000000; strcpy(axis->Name[11],"C4  ");
      break;
    case DT :
      if ( (axis = AllocAxis(10))==NULL ) return NULL;
      axis->ID = DT;
      AF = axis->F;
      AF->X[0 ]=-2.481000; AF->X[1 ]= 5.354000; AF->X[2 ]=0.000000; strcpy(axis->Name[0],"C1' ");
      AF->X[3 ]=-1.284000; AF->X[4 ]= 4.500000; AF->X[5 ]=0.000000; strcpy(axis->Name[1],"N1  ");
      AF->X[6 ]=-1.462000; AF->X[7 ]= 3.135000; AF->X[8 ]=0.000000; strcpy(axis->Name[2],"C2  ");
      AF->X[9 ]=-2.562000; AF->X[10]= 2.608000; AF->X[11]=0.000000; strcpy(axis->Name[3],"O2  ");
      AF->X[12]=-0.298000; AF->X[13]= 2.407000; AF->X[14]=0.000000; strcpy(axis->Name[4],"N3  ");
      AF->X[15]= 0.994000; AF->X[16]= 2.897000; AF->X[17]=0.000000; strcpy(axis->Name[5],"C4  ");
      AF->X[18]= 1.944000; AF->X[19]= 2.119000; AF->X[20]=0.000000; strcpy(axis->Name[6],"O4  ");
      AF->X[21]= 1.106000; AF->X[22]= 4.338000; AF->X[23]=0.000000; strcpy(axis->Name[7],"C5  ");
      AF->X[24]= 2.466000; AF->X[25]= 4.961000; AF->X[26]=0.001000; strcpy(axis->Name[8],"C7  ");
      AF->X[27]=-0.024000; AF->X[28]= 5.057000; AF->X[29]=0.000000; strcpy(axis->Name[9],"C6  ");
      break;
    case RA:
    case RC:
    case RG:
    case RU:
    case UNKNOWN_BASE:
      mprintf("Warning: NAstruct missing parameters for residue.\n");
  }

  return axis;
}

/*
 * DEBUG: AxisToPDB
 */
void NAstruct::AxisToPDB(PtrajFile *outfile, AxisType *axis, int resnum, int *atom) {
  char buffer[82];
  int i3=0;
  for (int i=0; i<axis->F->natom; i++) {
    pdb_write_ATOM(buffer,"ATOM",(*atom)+i,axis->Name[i],P->ResidueName(resnum),'X',resnum+1,
                   axis->F->X[i3],axis->F->X[i3+1],axis->F->X[i3+2],1.0,0.0,(char*)"\0");
    outfile->IO->Write(buffer,sizeof(char),strlen(buffer));
    i3+=3;
  } 
  (*atom) += axis->F->natom; 
}

/*
 * NAstruct::GCpair()
 * Look for 3 HB based on heavy atom distances:
 * 1. G:O6 -- C:N4  6 -- 6
 * 2. G:N1 -- C:N3  7 -- 4
 * 3. G:N2 -- C:O2  9 -- 3
 * Atom positions are known in standard Ref. Multiply by 3 to get into X.
 */
bool NAstruct::GCpair(AxisType *DG, AxisType *DC) {
  int Nhbonds = 0;
  double dist2; 
  dist2 = DIST2_NoImage(DG->F->X+18, DC->F->X+18);
  if ( dist2 < HBcut2 ) {
    Nhbonds++;
    mprintf("            G:O6 -- C:N4 = %lf\n",sqrt(dist2));
  }
  dist2 = DIST2_NoImage(DG->F->X+21, DC->F->X+12);
  if ( dist2 < HBcut2 ) {
    Nhbonds++;
    mprintf("            G:N1 -- C:N3 = %lf\n",sqrt(dist2));
  }
  dist2 = DIST2_NoImage(DG->F->X+27, DC->F->X+9 );
  if ( dist2 < HBcut2 ) {
    Nhbonds++;
    mprintf("            G:N2 -- C:O2 = %lf\n",sqrt(dist2));
  }
  if (Nhbonds>0) return true;
  return false;
}

/*
 * NAstruct::ATpair()
 * Look for 2 HB based on heavy atom distances
 * 1. A:N6 -- T:O4  6 -- 6
 * 2. A:N1 -- T:N3  7 -- 4
 */
bool NAstruct::ATpair(AxisType *DA, AxisType *DT) {
  int Nhbonds = 0;
  double dist2;
  dist2 = DIST2_NoImage(DA->F->X+18, DT->F->X+18);
  if ( dist2 < HBcut2 ) {
    Nhbonds++;
    mprintf("            A:N6 -- T:O4 = %lf\n",sqrt(dist2));
  }
  dist2 = DIST2_NoImage(DA->F->X+21, DT->F->X+12);
  if ( dist2 < HBcut2 ) {
    Nhbonds++;
    mprintf("            A:N1 -- T:N3 = %lf\n",sqrt(dist2));
  }
  if (Nhbonds>0) return true;
  return false;
}

/* 
 * NAstruct::basesArePaired()
 * Given two base axes for which IDs have been given and reference coords set,
 * determine whether the bases are paired via hydrogen bonding criteria.
 */
bool NAstruct::basesArePaired(AxisType *base1, AxisType *base2) {
  // G C
  if      ( base1->ID==DG && base2->ID==DC ) return GCpair(base1,base2);
  else if ( base1->ID==DC && base2->ID==DG ) return GCpair(base2,base1);
  else if ( base1->ID==DA && base2->ID==DT ) return ATpair(base1,base2);
  else if ( base1->ID==DT && base2->ID==DA ) return ATpair(base2,base1);
//  else {
//    mprintf("Warning: NAstruct: Unrecognized pair: %s - %s\n",NAbaseName[base1->ID],
//             NAbaseName[base2->ID]);
//  }
  return false;
}

/*
 * NAstruct::determineBasePairing()
 * Determine which bases are paired from the base axes.
 */
int NAstruct::determineBasePairing() {
  double distance;
  std::vector<bool> isPaired( BaseAxes.size(), false);
  int base1,base2,numBP;

  numBP = 0;
  
  mprintf(" ==== Setup Base Pairing ==== \n");

  /* For each unpaired base, determine if it is paired with another base
   * determined by the distance between their axis origins.
   */
  for (base1=0; base1 < Nbases-1; base1++) {
    if (isPaired[base1]) continue;
    for (base2=base1+1; base2 < Nbases; base2++) {
      if (isPaired[base2]) continue;
      // First determine if origin axes coords are close enough to consider pairing
      // Origin is 4th coord
      distance = DIST2_NoImage(BaseAxes[base1]->F->X+4, BaseAxes[base2]->F->X+4);
      mprintf("  Axes distance for %i:%s -- %i:%s is %lf\n",
              base1,NAbaseName[RefCoords[base1]->ID],
              base2,NAbaseName[RefCoords[base2]->ID],distance);
      if (distance < Ocut2) {
        mprintf("    Checking %i:%s -- %i:%s\n",base1,
                NAbaseName[RefCoords[base1]->ID],base2,NAbaseName[RefCoords[base2]->ID]);
        if (basesArePaired(RefCoords[base1], RefCoords[base2])) {
          if (numBP+1 > Nbp) BasePair = (BPTYPE*) realloc(BasePair, (numBP+1)*sizeof(BPTYPE));
          BasePair[numBP][0] = base1;
          BasePair[numBP][1] = base2;
          isPaired[base1]=true;
          isPaired[base2]=true;
          numBP++;
        }
      }
    } // END Loop over base2
  } // END Loop over base1
  Nbp = numBP;

  mprintf("    NAstruct: Set up %i base pairs.\n",Nbp);
  for (base1=0; base1 < Nbp; base1++) 
    mprintf("        BP %i: Res %i:%s to %i:%s\n",base1+1,
            BasePair[base1][0]+1, NAbaseName[RefCoords[BasePair[base1][0]]->ID],
            BasePair[base1][1]+1, NAbaseName[RefCoords[BasePair[base1][1]]->ID]);

  return 0;
}

/*
 * NAstruct::flipAxesYZ
 * Given Axes, flip the Z and Y axes.
 * Equiv. to rotation around the X axis.
 */
void NAstruct::flipAxesYZ( Frame *axis ) {
  double ox2, oy2, oz2;
  
  ox2 = axis->X[9 ] + axis->X[9 ];
  oy2 = axis->X[10] + axis->X[10];
  oz2 = axis->X[11] + axis->X[11];

  axis->X[3 ] = ox2 - axis->X[3 ];
  axis->X[4 ] = oy2 - axis->X[4 ];
  axis->X[5 ] = oz2 - axis->X[5 ];
 
  axis->X[6 ] = ox2 - axis->X[6 ];
  axis->X[7 ] = oy2 - axis->X[7 ];
  axis->X[8 ] = oz2 - axis->X[8 ];
}

/*
 * NAstruct::setupBasePairAxes()
 * Given a list of base pairs and base axes, setup an 
 * Axestype structure containing reference base pair axes.
 * The axis extraction equation is based on that found in:
 *   3D game engine design: a practical approach to real-time Computer Graphics,
 *   Volume 385, By David H. Eberly, 2001, p. 16.
 */
int NAstruct::setupBasePairAxes() {
  int basepair;
  double RotMatrix[9], TransVec[6], V[3], theta;
  AxisType *Base1;
  AxisType *Base2;
  // DEBUG
  int basepairaxesatom=0;
  PtrajFile basepairaxesfile;
  basepairaxesfile.SetupFile((char*)"basepairaxes.pdb",WRITE,UNKNOWN_FORMAT,UNKNOWN_TYPE,0);
  basepairaxesfile.OpenFile();
  // END DEBUG

  mprintf(" ==== Setup Base Pair Axes ==== \n");
  // Loop over all base pairs
  for (basepair = 0; basepair < Nbp; basepair++) {
    // Set Base1 as a copy, this will become the base pair axes
    Base1 = AxisCopy( BaseAxes[ BasePair[basepair][0] ] );
    Base2 = BaseAxes[ BasePair[basepair][1] ];
    // Set frame coords for first base in the pair
    AxisToFrame( Base1, Axes1);
    // Set frame coords for second base in the pair
    AxisToFrame( Base2, Axes2);
    // Flip the axes in the second pair
    // NOTE: This flip is only correct for standard WC base-pairing, not explicitly checked
    flipAxesYZ( Axes2 );
    // RMS fit Axes of Base2 onto Base1
    Axes2->RMSD( Axes1, RotMatrix, TransVec, false);
    // Extract angle from resulting rotation matrix
    theta=matrix_to_angle(RotMatrix);
    mprintf("Base2 for pair %i will be rotated by %lf degrees.\n",basepair+1,RADDEG*theta);
    // Calc Axis of rotation
    if (axis_of_rotation(V, RotMatrix, theta)) {
      mprintf("Error: NAstruct::setupBasePairAxes(): Could not set up axis of rotation for %i.\n",
              basepair);
      FreeAxis( Base1 );
      return 1;
    }
    // Calculate new half-rotation matrix
    calcRotationMatrix(RotMatrix,V,theta/2);
    /* Rotate Base2 by half rotation towards Base1.
     * Since the new rotation axis by definition is located at
     * the origin, the coordinates of the base have to be shifted, rotated,
     * then shifted back. Use V to store the shift back.
     */
    V[0] = -TransVec[0];
    V[1] = -TransVec[1];
    V[2] = -TransVec[2];
    Base2->F->Translate( TransVec );
    Base2->F->Rotate( RotMatrix );
    Base2->F->Translate( V );
    /* Since rotation matrix for Base2 was calculated with same origin as 
     * Base1, use reverse rotation matrix (transpose) to rotate Base1. 
     * Shift, rotate, shift back. Use V to store the shift back.
     */
    V[0] = -TransVec[3];
    V[1] = -TransVec[4];
    V[2] = -TransVec[5];
    Base1->F->Translate( V );
    Base1->F->InverseRotate( RotMatrix );
    Base1->F->Translate( TransVec+3 );
    // Origin of base pair axes is midpoint between Base1 and Base2 origins
    V[0] = ( (Base1->F->X[9 ] + Base2->F->X[9 ])/2 ) - Base1->F->X[9 ];
    V[1] = ( (Base1->F->X[10] + Base2->F->X[10])/2 ) - Base1->F->X[10];
    V[2] = ( (Base1->F->X[11] + Base2->F->X[11])/2 ) - Base1->F->X[11];
    // Shift Base1 to midpoint; Base1 becomes the base pair axes
    Base1->F->Translate( V );
    // DEBUG
    V[0] = Base1->F->X[6] - Base1->F->X[9];
    V[1] = Base1->F->X[7] - Base1->F->X[10];
    V[2] = Base1->F->X[8] - Base1->F->X[11];
    //normalize(V);
    // NOTE: Axes are already normalized
    mprintf("      %i) %i:%s -- %i:%s  %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf\n",basepair,
            BasePair[basepair][0], NAbaseName[Base1->ID],
            BasePair[basepair][1], NAbaseName[Base2->ID],
            Base1->F->X[9], Base1->F->X[10], Base1->F->X[11],
            V[0], V[1], V[2]);
    // Base1 now contains absolute coords of base pair reference axes
    AxisToPDB(&basepairaxesfile, Base1, BasePair[basepair][0], &basepairaxesatom);
    BasePairAxes.push_back( Base1 );
  }
  // DEBUG
  basepairaxesfile.CloseFile();

  return 0;
}
// ----------------------------------------------------------------------------

/*
 * NAstruct::init()
 * Expected call: nastruct [resrange <range>] [out <filename>]
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int NAstruct::init() {
  char *rangeArg;
  // Get keywords
  outFilename = A->getKeyString("out",NULL);
  rangeArg = A->getKeyString("resrange",NULL); 
  resRange = A->NextArgToRange(rangeArg); 
  // Get Masks
  // Dataset
  // Add dataset to data file list

  // Temporary frames for performing fits on Axes
  Axes1 = new Frame(4,NULL);
  Axes2 = new Frame(4,NULL);

  mprintf("    NAstruct: ");
  if (resRange==NULL)
    mprintf("Scanning all NA residues");
  else
    mprintf("Scanning residues %s",rangeArg);
  if (outFilename!=NULL)
    mprintf(", output to file %s",outFilename);
  mprintf("\n");

  return 0;
}

/*
 * NAstruct::setup()
 * Set up for this parmtop. Get masks etc.
 * P is set in Action::Setup
 */
int NAstruct::setup() {
  int res, refAtom, atom;
  std::list<int>::iterator residue;
  AxisType *axis; 
  AtomMask *Mask;
  Frame *expframe;

  // Clear all lists
  ClearLists();

  // If range arg is NULL look for all NA residues.
  if (resRange==NULL) {
    resRange = new std::list<int>();
    for (res=0; res < P->nres; res++) {
      if ( ID_base(P->ResidueName(res))!=UNKNOWN_BASE )
        resRange->push_back(res);
    }

  // Otherwise, for each residue in resRange check if it is a NA
  } else {
    residue=resRange->begin();
    while (residue!=resRange->end()) {
      // User residues numbers start from 1
      (*residue) = (*residue) - 1;
      if (ID_base(P->ResidueName(*residue))==UNKNOWN_BASE) 
        residue = resRange->erase(residue);
      else 
        residue++;
    }
  }
  // Exit if no NA residues specified
  if (resRange->empty()) {
    mprintf("Error: NAstruct::setup: No NA residues found for %s\n",P->parmName);
    return 1;
  }

  // DEBUG - print all residues
  mprintf("    NAstruct: NA res:");
  for (residue=resRange->begin(); residue!=resRange->end(); residue++)
    mprintf(" %i",(*residue)+1);
  mprintf("\n");

  // Set up reference coords for each NA residue
  for (residue=resRange->begin(); residue!=resRange->end(); residue++) {
    axis = getRefCoords( ID_base(P->ResidueName(*residue)) );
    if (axis==NULL) {
      mprintf("Error: NAstruct::setup: Could not get ref coords for %i:%s\n",
              (*residue)+1, P->ResidueName(*residue));
      return 1;
    }
    RefCoords.push_back( axis );

    // Set up a mask for this NA residue in this parm. The mask will contain
    // only those atoms which are defined in the reference coords.
    Mask = new AtomMask();
    for (refAtom=0; refAtom < axis->F->natom; refAtom++) {
      res = -1; // Target atom
      //mprintf("      Ref atom: [%s]\n",axis->Name[refAtom]);
      for (atom=P->resnums[*residue]; atom < P->resnums[(*residue)+1]; atom++) {
        //mprintf("        Scanning %i [%s]\n", atom, P->names[atom]);
        if ( strcmp(axis->Name[refAtom], P->names[atom])==0 ) {
          res = atom;  
          break;
        }
      }
      if (res==-1) {
        mprintf("Error:: NAstruct::setup: Ref atom [%s] not found in residue %i:%s\n",
                 axis->Name[refAtom], (*residue)+1, P->ResidueName(*residue));
        return 1;
      }
      Mask->AddAtom(res);
    } // End Loop over reference atoms
    if (Mask->None()) {
      mprintf("Error:: NAstruct::setup: No atoms found for residue %i:%s\n",
              (*residue)+1, P->ResidueName(*residue));
      delete Mask;
      return 1;
    }
    ExpMasks.push_back( Mask );
    mprintf("      NAstruct: Res %i:%s mask atoms: ",(*residue)+1,P->ResidueName(*residue));
    Mask->PrintMaskAtoms();
    mprintf("\n");

    // Set up frame to hold input coords for this residue
    expframe = new Frame(Mask, P->mass);
    ExpFrames.push_back( expframe );

    // Set up initial axes for this NA residue.
    // NOTE: OK to overwrite axis here since it has been pushed to RefCoords already.
    axis = principalAxes();
    BaseAxes.push_back( axis );
  } // End Loop over NA residues
  Nbases = RefCoords.size(); // Also BaseAxes, ExpFrames, and ExpMasks size.
  mprintf("    NAstruct: Set up %i bases.\n",Nbases);

  return 0;  
}

/*
 * NAstruct::action()
 */
int NAstruct::action() {
  double rmsd, RotMatrix[9], TransVec[6];
  int base;
  AxisType *Origin;
  Frame *REF_TEMP=NULL;
  Frame *EXP_TEMP=NULL;
  // DEBUG
  int res = 0;
  int baseaxesatom = 0;
  int basesatom = 0;
  PtrajFile baseaxesfile;
  PtrajFile basesfile;
  baseaxesfile.SetupFile((char*)"baseaxes.pdb",WRITE,UNKNOWN_FORMAT,UNKNOWN_TYPE,0);
  baseaxesfile.OpenFile();
  basesfile.SetupFile((char*)"bases.pdb",WRITE,UNKNOWN_FORMAT,UNKNOWN_TYPE,0);
  basesfile.OpenFile();
  // END DEBUG

  // Set up Origin to be a permanent origin
  Origin = principalAxes();
  // For each axis in RefCoords, use corresponding mask in ExpMasks to set 
  // up an axis for ExpCoords.
  for (base=0; base < Nbases; base++) {
    // Reset Origin coords so it is always the origin
    setPrincipalAxes(Origin->F->X);
    // Set exp coords based on mask
    ExpFrames[base]->SetFrameCoordsFromMask( F->X, ExpMasks[base] ); 
    /* Now that we have a set of reference coords and the corresponding input
     * coords, RMS fit the reference coords to the input coords to obtain the
     * appropriate rotation and translations that will put the reference coords 
     * on top of input (experimental) coords.
     * NOTE: The RMSD routine is destructive to coords. Need copies of frames.
     */
    if (REF_TEMP!=NULL) delete REF_TEMP;
    if (EXP_TEMP!=NULL) delete EXP_TEMP;
    REF_TEMP = RefCoords[base]->F->Copy();
    EXP_TEMP = ExpFrames[base]->Copy();
    rmsd = REF_TEMP->RMSD( EXP_TEMP, RotMatrix, TransVec, false);
    mprintf("Base %i: RMS of RefCoords from ExpCoords is %lf\n",base+1,rmsd);
    //AxisToPDB(&baseaxesfile, (*baseaxis), res++, &baseaxesatom);
    // RotMatrix and TransVec now contain rotation and translation
    // that will orient refcoord to expframe.
    // Use the translation/rotation to fit principal axes in BaseAxes to experimental coords.
    BaseAxes[base]->F->Translate( TransVec );
    BaseAxes[base]->F->Rotate( RotMatrix );
    BaseAxes[base]->F->Translate( TransVec + 3);
    /* baseaxis now contains the absolute coordinates of the base reference axes.
     * Perform a second RMS fit to get the translation and rotation from
     * the absolute origin to the base reference axes.
     */
    delete REF_TEMP;
    REF_TEMP = BaseAxes[base]->F->Copy();
    Origin->F->RMSD( REF_TEMP, RotMatrix, TransVec, false);
    /* RotMatrix and TransVec now contain rotation and translation from 
     * origin (?absolute?) coords to base reference coords.
     */

    // DEBUG - Write base axis to file
    AxisToPDB(&baseaxesfile, BaseAxes[base], res, &baseaxesatom);

    // Overlap ref coords onto input coords
    // Since ref coords start at origin by default the first translation is not necessary
    RefCoords[base]->F->Rotate( RotMatrix );
    RefCoords[base]->F->Translate( TransVec + 3);
    // DEBUG - Write ref coords to file
    AxisToPDB(&basesfile, RefCoords[base], res, &basesatom);
    res++;
  }
  if (REF_TEMP!=NULL) delete REF_TEMP;
  if (EXP_TEMP!=NULL) delete EXP_TEMP;    

  // DEBUG
  baseaxesfile.CloseFile();
  basesfile.CloseFile();
  // Free up Origin
  FreeAxis( Origin );

  // Determine Base Pairing
  determineBasePairing();

  // Get base pair axes
  setupBasePairAxes();

  return 0;
} 


