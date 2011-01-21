#include <cstdlib>
#include <cstdio> // sprintf
#include <cmath>
#include <cstring> // memset
#include "Frame.h"
#include "vectormath.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR, takes number of atoms and masses
// NOTE: Probably should throw execeptions on mem errors
Frame::Frame(int natomIn, double *MassIn) {
  natom=natomIn;
  N=natom*3;
  X=(double*) malloc( N * sizeof(double));
  box[0]=0; box[1]=0; box[2]=0; box[3]=0; box[4]=0; box[5]=0;
  T=0.0;
  V=NULL;
  Mass=NULL;
  if (MassIn!=NULL) {
    Mass = (double*) malloc(natom * sizeof(double));
    for (int atom = 0; atom < natom; atom++)
      Mass[atom] = MassIn[atom];
  }
  ucell[0]=0; ucell[1]=0; ucell[2]=0;
  ucell[3]=0; ucell[4]=0; ucell[5]=0;
  ucell[6]=0; ucell[7]=0; ucell[8]=0;
  recip[0]=0; recip[1]=0; recip[2]=0;
  recip[3]=0; recip[4]=0; recip[5]=0;
  recip[6]=0; recip[7]=0; recip[8]=0;
}

// CONSTRUCTOR, Create Frame based on size of mask. Put correct masses in.
Frame::Frame(AtomMask *Mask, double *MassIn) {
  natom = Mask->Nselected;
  N = natom * 3;
  X = (double*) malloc(N * sizeof(double));
  box[0]=0; box[1]=0; box[2]=0; box[3]=0; box[4]=0; box[5]=0;
  T=0.0;
  V=NULL;
  // Copy Mass info if present
  Mass = NULL;
  if (MassIn!=NULL) {
    Mass = (double*) malloc(natom * sizeof(double));
    for (int i=0; i < Mask->Nselected; i++)
      Mass[i] = MassIn[ Mask->Selected[i] ];
  }
  ucell[0]=0; ucell[1]=0; ucell[2]=0;
  ucell[3]=0; ucell[4]=0; ucell[5]=0;
  ucell[6]=0; ucell[7]=0; ucell[8]=0;
  recip[0]=0; recip[1]=0; recip[2]=0;
  recip[3]=0; recip[4]=0; recip[5]=0;
  recip[6]=0; recip[7]=0; recip[8]=0;
}

// DESTRUCTOR
Frame::~Frame() {
  if (X!=NULL) free(X);
  if (V!=NULL) delete V;
  if (Mass!=NULL) free(Mass);
}

/*
 * Frame::Copy()
 * Return a copy of the frame
 */
Frame *Frame::Copy() {
  Frame *newFrame;
  int i;

  newFrame=new Frame(this->natom, this->Mass);
  for (i=0; i<this->N; i++)
    newFrame->X[i] = this->X[i];
  for (i=0; i<6; i++)
    newFrame->box[i] = this->box[i];
  newFrame->T = this->T;
  if (this->V!=NULL)
    newFrame->V = this->V->Copy();
  for (i=0; i<9; i++) {
    newFrame->ucell[i] = this->ucell[i];
    newFrame->recip[i] = this->recip[i];
  }

  return newFrame;
}

/*
 * Frame::printAtomCoord()
 * Print XYZ coords of given atom
 */
void Frame::printAtomCoord(int atom) {
  int natom3;

  natom3=atom*3;
  if (natom3>=N) return;
  mprintf("ATOM %i: %lf %lf %lf\n",atom,
          X[natom3],X[natom3+1],X[natom3+2]);
}

/*
 * Frame::SetCoord()
 * Set the given coordinates to those of given atom.
 */
void Frame::SetCoord(double *Coord, int atom) {
  int natom3;
  // NOTE: SHOULD CHECK FOR BOUNDARIES
  natom3 = atom * 3;
  Coord[0] = X[natom3  ];
  Coord[1] = X[natom3+1];
  Coord[2] = X[natom3+2];
}

/* ------------------- Coordinate Conversion Routines ----------------------- */
/*
 * Frame::BufferToBox()
 * Store given char buffer containing box coords with 
 * format XYZ[ABG] to box array. Each box coord has given width.
 * Return position in buffer after read. If * encountered this indicates
 * overflow in trajectory, return NULL.
 */
char *Frame::BufferToBox(char *buffer, int numBox, int width) {
  char *ptr;
  char number[64]; // Should not have to handle numbers wider than this!
  int i,atom;

//  number = (char*) malloc( (width+1) * sizeof(char));
  number[width]='\0';
  // Get box coordinates from buffer
  ptr=buffer;
  for (atom=0; atom<numBox; atom++) {
    i=0;
    while (i<width) {
      if (*ptr=='\n') ptr++;
      if (*ptr=='*') return NULL; //{free(number); return 1;}
      number[i++]=*ptr;
      ptr++;
    }
    //fprintf(stdout,"DEBUG: BOX %i: %s\n",atom,number);
    box[atom] = atof(number);
  }
//  free(number);
  return ptr;
}

/*
 * Frame::BufferToFrame()
 * Store given character buffer containing XYZ coords with
 * format X0Y0Z0X1Y1Z1...XNYNZN to the X array. Each coord has given width, 
 * newlines are skipped. buffer should be as big as N x width chars. 
 * Return position in buffer after read. If * encountered this indicates
 * overflow in trajectory, return NULL.
 */
char *Frame::BufferToFrame(char *buffer, int width) {
  char *ptr;
  char number[64]; // Should not have to handle numbers wider than this!
  int i,atom;

//  number = (char*) malloc( (width+1) * sizeof(char));
  number[width]='\0';
  ptr=buffer;
  for (atom=0; atom<N; atom++) {
    i=0;
    while (i<width) {
      if (*ptr=='\n') ptr++;
      if (*ptr=='*') return NULL;//{free(number); return 1;}
      number[i++]=*ptr;
      ptr++;
    }
    //fprintf(stdout,"DEBUG: %i: %s\n",atom/3,number);
    X[atom] = atof(number);
  }

//  free(number);
  return ptr;
}

/*
 * Frame::FrameToBuffer()
 * Given a character buffer, format, and character width corresponding
 * to format, write coords in frame to buffer. 
 * Return the position in the buffer after write. 
 */
char *Frame::FrameToBuffer(char *buffer, const char *format, int width, int numCols) {
  int coord;
  char *ptr;
  
  ptr=buffer;
  for (coord=0; coord<N; coord++) {
    sprintf(ptr,format,X[coord]);
    ptr+=width;
    if ( ((coord+1)%numCols)==0 ) {
      sprintf(ptr,"\n");
      ptr++;
    }
  }
  // If the coord record didnt end on a newline, print one
  if ( (coord%numCols)!=0 ) {
    sprintf(ptr,"\n");
    ptr++;
  }
  
  // Calculate frame size
  //coord = (int) (ptr - buffer);
  //return coord;
  return ptr;
}

/*
 * Frame::BoxToBuffer()
 * Given a character buffer, format, and character width corresponding
 * to format, write box coords to buffer.
 * Return the position in the buffer after write.
 */
char *Frame::BoxToBuffer(char *buffer, int numBox, const char *format, int width) {
//  int coord;
  char *ptr;

  ptr=buffer;
  // Box
  sprintf(ptr,format,box[0]); ptr+=width;
  sprintf(ptr,format,box[1]); ptr+=width;
  sprintf(ptr,format,box[2]); ptr+=width;
  if (numBox>3) {
    sprintf(ptr,format,box[3]); ptr+=width;
    sprintf(ptr,format,box[4]); ptr+=width;
    sprintf(ptr,format,box[5]); ptr+=width;
  }
  sprintf(ptr,"\n");
  ptr++;

  // Calculate frame size
  //coord = (int) (ptr - buffer);
  //return coord;
  return ptr;
}

/*
 * Frame::floatToFrame()
 * Convert float coords to double coords
 * NOTE: N needs to match up with size of Coord!
 */
void Frame::floatToFrame(float *Coord) {
  int i;
  for (i=0; i<N; i++)
    X[i]=(double) Coord[i];
}

/* 
 * Frame::frameToFloat()
 * Convert double coords to float coords
 */
void Frame::frameToFloat(float *Coord) {
  int i;
  for (i=0; i<N; i++)
    Coord[i]=(float) X[i];
}

/* ----------------- Coordinate Transformation Routines --------------------- */
/*
 * Frame::COM()
 * Given an AtomMask put geometric center of atoms in mask into Coord. Put 
 * center of mass instead if useMass is true.
 * Return sum of masses in Mask (useMass) or #atoms in Mask (!useMass).
 */
double Frame::COM(AtomMask *Mask, double *Coord, bool useMass) {
  int i,atom,natom3;
  double sumMass,mass;
  
  Coord[0]=0.0;
  Coord[1]=0.0;
  Coord[2]=0.0;
  sumMass=0.0;
  mass=1.0;

  //Mask->Start();
  //while ( (atom=Mask->NextAtom()) != -1 ) {
  for (i=0; i < Mask->Nselected; i++) {
      atom=Mask->Selected[i];
      natom3=atom*3;
      if (useMass) {
        mass=Mass[atom];
        sumMass+=mass;
      }
      Coord[0]+=(X[natom3]   * mass);
      Coord[1]+=(X[natom3+1] * mass);
      Coord[2]+=(X[natom3+2] * mass);
  }

  if (!useMass) sumMass=(double) Mask->Nselected;

  if (sumMass==0.0) return 0;

  Coord[0]/=sumMass;
  Coord[1]/=sumMass;
  Coord[2]/=sumMass;
  return sumMass;
}

/*
 * Frame::COM()
 * Put geometric center of all atoms between start and stop in frame into 
 * Coord. Put center of mass instead if useMass is true.
 * Return sum of masses (useMass) or #atoms (!useMass).
 * NOTE: overload to get rid of IF in loop?
 */
double Frame::COM(double *Coord, bool useMass, int startAtom, int stopAtom) {  
  int i,m;
  int startAtom3, stopAtom3;
  double sumMass,mass;
  
  Coord[0]=0.0;
  Coord[1]=0.0;
  Coord[2]=0.0;
  sumMass=0.0;
  mass=1.0;
  m=startAtom;
  startAtom3 = startAtom * 3;
  stopAtom3 = stopAtom * 3;
  
  for (i=startAtom3; i<stopAtom3; i+=3) {
    if (useMass) 
      mass=Mass[m++];
    sumMass+=mass;
    Coord[0]+=(X[i  ] * mass);
    Coord[1]+=(X[i+1] * mass);
    Coord[2]+=(X[i+2] * mass);
  }

  if (sumMass==0.0) return 0;

  Coord[0]/=sumMass;
  Coord[1]/=sumMass;
  Coord[2]/=sumMass;

  return sumMass;
}

/*
 * Frame::COM()
 * Put geometric center of all atoms in frame into Coord. Put center of mass 
 * instead if useMass is true.
 * Return sum of masses (useMass) or #atoms (!useMass).
 */
double Frame::COM(double *Coord, bool useMass) {

  return this->COM(Coord,useMass,0,natom);

}

/*
 * Frame::BoxToRecip()
 */
//void Frame::BoxToRecip(double *ucell, double *recip) {
void Frame::BoxToRecip() {
  double u12x,u12y,u12z;
  double u23x,u23y,u23z;
  double u31x,u31y,u31z;
  double volume;

  ucell[0] = box[0];
  ucell[1] = 0.0;
  ucell[2] = 0.0;
  ucell[3] = box[1]*cos(DEGRAD*box[5]);
  ucell[4] = box[1]*sin(DEGRAD*box[5]);
  ucell[5] = 0.0;
  ucell[6] = box[2]*cos(DEGRAD*box[4]);
  ucell[7] = (box[1]*box[2]*cos(DEGRAD*box[3]) - ucell[6]*ucell[3]) / ucell[4];
  ucell[8] = sqrt(box[2]*box[2] - ucell[6]*ucell[6] - ucell[7]*ucell[7]);
  
  u23x = ucell[4]*ucell[8] - ucell[5]*ucell[7];
  u23y = ucell[5]*ucell[6] - ucell[3]*ucell[8];
  u23z = ucell[3]*ucell[7] - ucell[4]*ucell[6];
  u31x = ucell[7]*ucell[2] - ucell[8]*ucell[1];
  u31y = ucell[8]*ucell[0] - ucell[6]*ucell[2];
  u31z = ucell[6]*ucell[1] - ucell[7]*ucell[0];
  u12x = ucell[1]*ucell[5] - ucell[2]*ucell[4];
  u12y = ucell[2]*ucell[3] - ucell[0]*ucell[5];
  u12z = ucell[0]*ucell[4] - ucell[1]*ucell[3];
  volume=ucell[0]*u23x + ucell[1]*u23y + ucell[2]*u23z;

  recip[0] = u23x/volume;
  recip[1] = u23y/volume;
  recip[2] = u23z/volume;
  recip[3] = u31x/volume;
  recip[4] = u31y/volume;
  recip[5] = u31z/volume;
  recip[6] = u12x/volume;
  recip[7] = u12y/volume;
  recip[8] = u12z/volume;
}

/*
 * Frame::ClosestImage()
 * Given two coordinates A and B, determine the unit XYZ vector that points 
 * towards the closest image of B to A.
 * It is assumed the coordinates are already relative to the box center.
 */
void Frame::ClosestImage(double *A, double *B, int *ixyz) {
  double halfBox[3];//, delta;
  int vectorA[3], vectorB[3], i;

  mprintf("DEBUG: CoordA     = %lf %lf %lf\n",A[0],A[1],A[2]);
  mprintf("DEBUG: CoordB     = %lf %lf %lf\n",B[0],B[1],B[2]);

  halfBox[0] = box[0] / 6.0; 
  halfBox[1] = box[1] / 6.0; 
  halfBox[2] = box[2] / 6.0;
  mprintf("DEBUG: Half box = %lf %lf %lf\n",halfBox[0],halfBox[1],halfBox[2]);

  // Vector A
  vectorA[0] = 0; vectorA[1] = 0; vectorA[2] = 0;
  for (i=0; i<3; i++) {
    if (A[i] < -halfBox[i]) vectorA[i] = -1;
    if (A[i] >  halfBox[i]) vectorA[i] =  1;
/*    delta = A[i] - boxCenter[i];
    if (delta > 0.0) vectorA[i] = 1;
    if (delta < 0.0) vectorA[i] = -1;
*/
  }
  mprintf("DEBUG:  VectorA = %2i %2i %2i\n",vectorA[0],vectorA[1],vectorA[2]);

  // NOT Vector B
  vectorB[0] = 0; vectorB[1] = 0; vectorB[2] = 0;
  for (i=0; i<3; i++) {
    if (B[i] < -halfBox[i]) vectorB[i] =  1; // NOT
    if (B[i] >  halfBox[i]) vectorB[i] = -1; // NOT
/*    delta = B[i] - boxCenter[i];
    if (delta > 0.0) vectorB[i] = -1; // NOT
    if (delta < 0.0) vectorB[i] = 1;  // NOT
*/
  }
  mprintf("DEBUG: !VectorB = %2i %2i %2i\n",vectorB[0],vectorB[1],vectorB[2]);

  // A & !B
  ixyz[0]=0; ixyz[1]=0; ixyz[2]=0;
  for (i=0; i<3; i++) {
    if (vectorA[i] == vectorB[i]) ixyz[i] = vectorA[i];
    //ixyz[i] = vectorA[i] & vectorB[i];
  }
}

/*
 * Frame::MinImageNonOrtho2()
 * Given two sets of coordinates and reciprocal space information based on
 * the current non-orthorhombic box, return the shortest imaged distance^2
 * between the coordinates.
 * The integer coefficients describing the closest reflection in reciprocal
 * space will be placed in ixyz.
 */
double Frame::MinImageNonOrtho2(double *Coord1, double *Coord2, bool origin, int *ixyz) {
  double min, f[3], f2[3];

  min = 100.0 * (box[0]*box[0]+box[1]*box[1]+box[2]*box[2]);

  //if (prnlev > 6) {
  //  fprintf(stdout, "ATOM      0  XXX A1      1     %7.3f %7.3f %7.3f\n",
  //          x1, y1, z1);
  //  fprintf(stdout, "ATOM      1  XXX A2      1     %7.3f %7.3f %7.3f\n",
  //          x2, y2, z2);
  //}

  f[0] = Coord1[0]*recip[0] + Coord1[1]*recip[1] + Coord1[2]*recip[2];
  f[1] = Coord1[0]*recip[3] + Coord1[1]*recip[4] + Coord1[2]*recip[5];
  f[2] = Coord1[0]*recip[6] + Coord1[1]*recip[7] + Coord1[2]*recip[8];

  f2[0] = Coord2[0]*recip[0] + Coord2[1]*recip[1] + Coord2[2]*recip[2];
  f2[1] = Coord2[0]*recip[3] + Coord2[1]*recip[4] + Coord2[2]*recip[5];
  f2[2] = Coord2[0]*recip[6] + Coord2[1]*recip[7] + Coord2[2]*recip[8];

  if (origin) {
    f[0] += 0.5;
    f[1] += 0.5;
    f[2] += 0.5;
    f2[0] += 0.5;
    f2[1] += 0.5;
    f2[2] += 0.5;
  }

  min = this->DIST2_ImageNonOrtho(f, f2, min, ixyz);

  return min;
}

/*
 * Frame::DIST2_ImageNonOrtho()
 * Given two coordinates and reciprocal space information based on 
 * the current non-orthorhombic box, return the shortest imaged distance^2 
 * between the coordinates.
 */
double Frame::DIST2_ImageNonOrtho(double *a1, double *a2) { // double closest2
  double f[3], f2[3];
  int ixyz[3];

  f[0] = a2[0]*recip[0] + a2[1]*recip[1] + a2[2]*recip[2];
  f[1] = a2[0]*recip[3] + a2[1]*recip[4] + a2[2]*recip[5];
  f[2] = a2[0]*recip[6] + a2[1]*recip[7] + a2[2]*recip[8];

  f2[0] = a1[0]*recip[0] + a1[1]*recip[1] + a1[2]*recip[2];
  f2[1] = a1[0]*recip[3] + a1[1]*recip[4] + a1[2]*recip[5];
  f2[2] = a1[0]*recip[6] + a1[1]*recip[7] + a1[2]*recip[8];

  return this->DIST2_ImageNonOrtho(f, f2, -1.0, ixyz);
}

/*
 * Frame::DIST2_ImageNonOrtho()
 * Given two coordinate sets in reciprocal space, return the minimum imaged
 * distance^2 between them.
 * If minIn is > 0.0 it is considered a possible minimum distance.
 * The integer coefficients describing the closest reflection in reciprocal
 * space will be placed in ixyz.
 */
double Frame::DIST2_ImageNonOrtho(double *f, double *f2, double minIn, int *ixyz) { // double closest2
  double fx, fy, fz, f2x, f2y, f2z;
  double x,y,z,D,min;
  int ix,iy,iz;
  // DEBUG

  /*
   *  NON-ORTHORHOMBIC CASE: find shortest distance in periodic reference
   *  This is a brute force check requiring up to 26 distance evaluations.
   *  It has been adapted to be smarter by returning the first distance that
   *  is shorter than the minimum possible distance between images.
   */

  fx = f[0] - floor(f[0]);
  fy = f[1] - floor(f[1]);
  fz = f[2] - floor(f[2]); 
  
  f2x = f2[0] - floor(f2[0]);
  f2y = f2[1] - floor(f2[1]);
  f2z = f2[2] - floor(f2[2]);

  // Calc ix iy iz = 0 case
  x = (fx*ucell[0] + fy*ucell[3] + fz*ucell[6]) - (f2x*ucell[0] + f2y*ucell[3] + f2z*ucell[6]);
  y = (fx*ucell[1] + fy*ucell[4] + fz*ucell[7]) - (f2x*ucell[1] + f2y*ucell[4] + f2z*ucell[7]);
  z = (fx*ucell[2] + fy*ucell[5] + fz*ucell[8]) - (f2x*ucell[2] + f2y*ucell[5] + f2z*ucell[8]);
  // DEBUG
  //fprintf(stdout,"DEBUG: a2: fx  fy  fz  = %lf %lf %lf\n",fx,fy,fz);
  //fprintf(stdout,"DEBUG: a1: f2x f2y f2z = %lf %lf %lf\n",f2x,f2y,f2z);
  min = (x*x) + (y*y) + (z*z);

  if (minIn > 0.0 && minIn < min) min = minIn;

  ixyz[0] = 0;
  ixyz[1] = 0;
  ixyz[2] = 0;

  //if (closest2 != 0.0 && min < closest2) return (min);
//  this->ClosestImage(a1, a2, ixyz);
//  fprintf(stdout,"DEBUG: Predict  = %2i %2i %2i\n",ixyz[0],ixyz[1],ixyz[2]);

//  ix = ixyz[0];
//  iy = ixyz[1];
//  iz = ixyz[2];

  for (ix = -1; ix <= 1; ix++) {
    for (iy = -1; iy <= 1; iy++) {
      for (iz = -1; iz <= 1; iz++) {

        if (! (ix == 0 && iy == 0 && iz == 0) ) {
          x = ((fx+ix)*ucell[0] + (fy+iy)*ucell[3] + (fz+iz)*ucell[6]) - 
              (f2x*ucell[0] + f2y*ucell[3] + f2z*ucell[6]);
          y = ((fx+ix)*ucell[1] + (fy+iy)*ucell[4] + (fz+iz)*ucell[7]) - 
              (f2x*ucell[1] + f2y*ucell[4] + f2z*ucell[7]);
          z = ((fx+ix)*ucell[2] + (fy+iy)*ucell[5] + (fz+iz)*ucell[8]) - 
              (f2x*ucell[2] + f2y*ucell[5] + f2z*ucell[8]);
          D = (x*x) + (y*y) + (z*z);

          //if (debug > 3) 
          //  printf("DISTANCE + %2i*X %2i*Y %2i*Z unit cells is %8.3f\n", ix, iy, iz, sqrt(D));
          
          if (D < min) {
            min = D;
            ixyz[0] = ix;
            ixyz[1] = iy;
            ixyz[2] = iz;
            // DEBUG
            //if (closest2 != 0.0 && min < closest2)
            //  return(min);
          }
        }

      }
    }
  }

  //D = sqrt(min);
//  fprintf(stdout,"DEBUG: MinDist  = %2i %2i %2i = %8.3f\n", ixmin, iymin, izmin, D);
//  printf("---------------------------------------------------------------\n");
  return(min);
}

/*
 * Frame::DIST2_ImageOrtho()
 * Return the minimum orthorhombic imaged distance^2 between coordinates a1 
 * and a2.
 */
double Frame::DIST2_ImageOrtho(double *a1, double *a2) {
  double x,y,z,D;

  x = a1[0] - a2[0];
  y = a1[1] - a2[1];
  z = a1[2] - a2[2];

  // Get rid of sign info
  if (x<0) x=-x;
  if (y<0) y=-y;
  if (z<0) z=-z;

  // Get rid of multiples of box lengths 
  while (x > box[0]) x = x - box[0];
  while (y > box[1]) y = y - box[1];
  while (z > box[2]) z = z - box[2];

  // Find shortest distance in periodic reference
  D = box[0] - x;
  if (D < x) x = D;
  D = box[1] - y;
  if (D < y) y = D;  
  D = box[0] - z;
  if (D < z) z = D;

  x = x * x;
  y = y * y;
  z = z * z;
 
  //D = sqrt(x + y + z);
  D = x + y + z;

  return D;
}

/*
 * Frame::DIST2_NoImage()
 * Return distance^2 between coordinates in a1 and a2.
 */
double Frame::DIST2_NoImage(double *a1, double *a2) {
  double x,y,z,D;

  x = a1[0] - a2[0];
  y = a1[1] - a2[1];
  z = a1[2] - a2[2];

  x=x*x;
  y=y*y;
  z=z*z;

  //D=sqrt(x + y + z);
  D = x + y + z;

  //fprintf(stdout,"Mask1=%8.3lf %8.3lf %8.3lf Mask2=%8.3lf %8.3lf %8.3lf D=%8.3lf\n",
  //        a1[0],a1[1],a1[2],a2[0],a2[1],a2[2],D);

  return D;
}

/*
 * Frame::DIST2()
 * Call the appropriate distance calc for atoms in Mask1 and Mask2 based on
 * given box type.
 *   0 = None
 *   1 = Orthorhombic
 *   2 = Non-orthorhombic
 * Based on useMassIn, calculate geometric center (false) or center of mass 
 * (true) of the atoms in each mask.
 */
double Frame::DIST2(AtomMask *Mask1, AtomMask *Mask2, bool useMassIn, int ifbox) {
  //double ucell[9], recip[9];
  double a1[3], a2[3];

  COM(Mask1, a1, useMassIn);
  COM(Mask2, a2, useMassIn);

  if (ifbox == 0) 
    return this->DIST2_NoImage(a1, a2);
  else if (ifbox == 1) 
    return this->DIST2_ImageOrtho(a1, a2);
  else if (ifbox == 2) {
    //this->BoxToRecip(ucell, recip);
    return this->DIST2_ImageNonOrtho(a1, a2);
  }

  mprintf("    Error: Frame::DIST: Unrecognized box type (%i)\n.", ifbox);

  return (-1.0);
}

/*
 * Frame::DIST2()
 * Return the distance between atoms A1 and A2 with optional imaging.
 */
double Frame::DIST2(int A1, int A2, int ifbox) {
  int atom3;
  //double ucell[9], recip[9];
  double a1[3], a2[3];

  atom3 = A1 * 3;
  a1[0] = X[atom3  ];
  a1[1] = X[atom3+1];
  a1[2] = X[atom3+2];
  atom3 = A2 * 3;
  a2[0] = X[atom3  ];
  a2[1] = X[atom3+1];
  a2[2] = X[atom3+2];

  if (ifbox == 0)
    return this->DIST2_NoImage(a1, a2);
  else if (ifbox == 1)
    return this->DIST2_ImageOrtho(a1, a2);
  else if (ifbox == 2) {
    //this->BoxToRecip(ucell, recip);
    return this->DIST2_ImageNonOrtho(a1, a2);
  }

  mprintf("    Error: Frame::DIST: Unrecognized box type (%i)\n.", ifbox);

  return (-1.0);
}

/*
 * Frame::DIST()
 * Return the distance between atoms A1 and A2, no imaging.
 */
double Frame::DIST(int A1, int A2) {
  int i, j; // Actual indices into X
  double x,y,z,D;

  i = A1 * 3;
  j = A2 * 3;

  x = X[i  ] - X[j  ];
  y = X[i+1] - X[j+1];
  z = X[i+2] - X[j+2];

  x=x*x;
  y=y*y;
  z=z*z;

  D=sqrt(x + y + z);
  return D;
}

/* 
 * Frame::ANGLE()
 * Return the angle between atoms in M1, M2, M3.
 * Adapted from PTRAJ.
 */
double Frame::ANGLE(AtomMask *M1, AtomMask *M2, AtomMask *M3) {
  double a1[3],a2[3],a3[3];
  double angle, xij, yij, zij, xkj, ykj, zkj, rij, rkj;
  
  COM(M1,a1,false);
  COM(M2,a2,false);
  COM(M3,a3,false);

  xij = a1[0] - a2[0];
  yij = a1[1] - a2[1];
  zij = a1[2] - a2[2];

  xkj = a3[0] - a2[0];
  ykj = a3[1] - a2[1]; 
  zkj = a3[2] - a2[2];

  rij = xij*xij + yij*yij + zij*zij;
  rkj = xkj*xkj + ykj*ykj + zkj*zkj;

  if (rij > SMALL && rkj > SMALL) {
    angle = (xij*xkj + yij*ykj + zij*zkj) / sqrt(rij*rkj);
    if (angle > 1.0)
      angle = 1.0;
    else if (angle < -1.0)
      angle = -1.0;
    angle = acos(angle) * RADDEG;
  } else
    angle = 0.0;

  return angle;
}

/* 
 * Frame::ANGLE()
 * Return the angle between atoms specified by A1, A2, A3.
 * Adapted from PTRAJ.
 */
double Frame::ANGLE(int A1, int A2, int A3) {
  double angle, xij, yij, zij, xkj, ykj, zkj, rij, rkj;
  int a1, a2, a3;

  a1 = A1 * 3;
  a2 = A2 * 3;
  a3 = A3 * 3;
  xij = X[a1  ] - X[a2  ];
  yij = X[a1+1] - X[a2+1];
  zij = X[a1+2] - X[a2+2];

  xkj = X[a3  ] - X[a2  ];
  ykj = X[a3+1] - X[a2+1];
  zkj = X[a3+2] - X[a2+2];

  rij = xij*xij + yij*yij + zij*zij;
  rkj = xkj*xkj + ykj*ykj + zkj*zkj;

  if (rij > SMALL && rkj > SMALL) {
    angle = (xij*xkj + yij*ykj + zij*zkj) / sqrt(rij*rkj);
    if (angle > 1.0)
      angle = 1.0;
    else if (angle < -1.0)
      angle = -1.0;
    angle = acos(angle) * RADDEG;
  } else
    angle = 0.0;

  return angle;
}

/*
 * Frame::DIHEDRAL()
 * Return dihedral angle between COM of atoms in M1-M4
 * Adapted from PTRAJ
 */
double Frame::DIHEDRAL(AtomMask *M1, AtomMask *M2, AtomMask *M3, AtomMask *M4) {
  double a1[3],a2[3],a3[3],a4[3];
  double Lx, Ly, Lz, Lnorm;
  double Rx, Ry, Rz, Rnorm;
  double Sx, Sy, Sz;
  double angle;

  COM(M1,a1,false);
  COM(M2,a2,false);
  COM(M3,a3,false);
  COM(M4,a4,false);

  CROSS_PRODUCT(     Lx,      Ly,      Lz,
                (a2[0]-a1[0]), (a2[1]-a1[1]), (a2[2]-a1[2]),
                (a3[0]-a2[0]), (a3[1]-a2[1]), (a3[2]-a2[2]));

  CROSS_PRODUCT(     Rx,      Ry,      Rz,
                (a4[0]-a3[0]), (a4[1]-a3[1]), (a4[2]-a3[2]),
                (a2[0]-a3[0]), (a2[1]-a3[1]), (a2[2]-a3[2]));

  Lnorm = sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
  Rnorm = sqrt(Rx*Rx + Ry*Ry + Rz*Rz);

  CROSS_PRODUCT(Sx, Sy, Sz,
                Lx, Ly, Lz,
                Rx, Ry, Rz);

  angle = (Lx*Rx + Ly*Ry + Lz*Rz) / (Lnorm * Rnorm);

  if ( angle > 1.0 ) angle = 1.0;
  if ( angle < -1.0 ) angle = -1.0;

  angle = acos( angle );
  angle = angle * RADDEG;

  if ( (Sx * (a3[0]-a2[0]) + Sy * (a3[1]-a2[1]) + Sz * (a3[2]-a2[2])) < 0 )
    angle = -angle;

  return angle;
}

/*
 * Frame::RADGYR()
 * Return the radius of gyration of atoms in mask. Also set the maximum 
 * distance from center. Use center of mass if useMassIn is true.
 */
double Frame::RADGYR(AtomMask *Mask, bool useMassIn, double *max) {
  double mid[3], Coord[3];
  double currentMass, total_mass, maxMass, dist2, sumDist2;
  int i,atom,natom3;

  total_mass=0.0;
  maxMass=1.0;
  currentMass=1.0;
  sumDist2=0.0;
  *max=0.0;

  this->COM(Mask,mid,useMassIn);

  for (i=0; i < Mask->Nselected; i++) {
      atom=Mask->Selected[i];
      natom3=atom*3;
      if (useMassIn) {
        currentMass=Mass[atom];
        total_mass+=currentMass;
      }
      Coord[0] = X[natom3  ] - mid[0];
      Coord[1] = X[natom3+1] - mid[1];
      Coord[2] = X[natom3+2] - mid[2];
      dist2 = (Coord[0]*Coord[0]) + (Coord[1]*Coord[1]) + (Coord[2]*Coord[2]);
      dist2 *= currentMass;
      if (dist2 > *max) {
        *max = dist2;
        maxMass = currentMass;
      }
      sumDist2 += dist2;
  }

  if (!useMassIn) total_mass=(double) Mask->Nselected;

  if (total_mass==0.0) return 0;  

  currentMass = sqrt(sumDist2 / total_mass); // Radius of Gyration
  *max = sqrt(*max / maxMass);

  return currentMass;
}

/*
 * Frame::SetFrameFromMask()
 * Given an existing Frame and an AtomMask (with Nselected atoms and atom 
 * numbers corresponding to FrameIn), set this Frame to be a copy of FrameIn
 * according to Mask. If the number of atoms in the Mask is greater than
 * this frame has been set up for, resize this Frame accordingly. 
 * For this to work properly AtomMask needs to have been setup based on 
 * FrameIn, although this is not explicitly checked for.
 */
void Frame::SetFrameFromMask(Frame *FrameIn, AtomMask *Mask) {
  int i,oldatom3,newatom3;

  // Check if number of atoms in mask is greater than what is allocd for frame.
  // Realloc X if necessary, set natom and N.
  if (Mask->Nselected > natom) {
    natom = Mask->Nselected;
    N = natom * 3;
    X = (double*) realloc(X, N * sizeof(double));

  // Otherwise just set the new natom and N.
  } else {
    natom = Mask->Nselected;
    N = natom * 3;
  }

  newatom3 = 0;
  for (i=0; i < Mask->Nselected; i++) {
    oldatom3 = Mask->Selected[i] * 3;
    this->X[newatom3  ] = FrameIn->X[oldatom3  ];
    this->X[newatom3+1] = FrameIn->X[oldatom3+1];
    this->X[newatom3+2] = FrameIn->X[oldatom3+2];
    newatom3 += 3;
  }

  // Copy box/T as well
  for (i=0; i<6; i++)
    this->box[i] = FrameIn->box[i];
  this->T = FrameIn->T;

  // Copy Mass info if present
  if (FrameIn->Mass!=NULL) {
    if (Mask->Nselected>natom || Mass==NULL)
      Mass = (double*) realloc(Mass, natom * sizeof(double));
    for (i=0; i < Mask->Nselected; i++) 
      this->Mass[i] = FrameIn->Mass[Mask->Selected[i]];
  }

  // Copy velocities if present
  if (FrameIn->V!=NULL) {
    if (this->V==NULL)
      this->V = new Frame(FrameIn->V->natom, NULL);
    this->V->SetFrameFromMask(FrameIn->V, Mask);
  }
}

/*
 * Frame::SetFrameCoordsFromMask()
 * Like SetFrameFromMask except only copy the coordinates. Unlike
 * SetFrameFromMask the #atoms in this Frame must match the 
 * #selected atoms in Mask.
 */
int Frame::SetFrameCoordsFromMask(double *Xin, AtomMask *Mask) {
  int i,oldatom3,newatom3;

  if (Mask->Nselected != natom) {
    mprintf("Internal Error: Frame::SetFrameCoordsFromMask: Frame set for %i atoms\n",natom);
    mprintf("                but attempting to set %i atoms.\n",Mask->Nselected);
    return 1;
  }

  newatom3 = 0;
  for (i=0; i < Mask->Nselected; i++) {
    oldatom3 = Mask->Selected[i] * 3;
    this->X[newatom3  ] = Xin[oldatom3  ];
    this->X[newatom3+1] = Xin[oldatom3+1];
    this->X[newatom3+2] = Xin[oldatom3+2];
    newatom3 += 3;
  }
  return 0;
}

/*
 * Frame::Translate()
 * Translate X by V. 
 */
void Frame::Translate(double *V) {
  int i;

  for (i=0; i<N; i+=3) {
    X[i  ]+=V[0];
    X[i+1]+=V[1];
    X[i+2]+=V[2];
  }
}

/*
 * Frame::Translate()
 * Translate atoms between start and end in X by V.
 * NOTE: SHOULD CHECK BOUNDS! 
 */
void Frame::Translate(double *V, int startAtom, int endAtom) {
  int atom, startAtom3, endAtom3;

  startAtom3 = startAtom * 3;
  endAtom3 = endAtom * 3;

  for (atom=startAtom3; atom<endAtom3; atom+=3) {
    X[atom  ] += V[0];
    X[atom+1] += V[1];
    X[atom+2] += V[2];
  }
}

/*  
 * Frame::Rotate()
 * Multiply natomx3 matrix X by 3x3 matrix T. If T is a rotation matrix
 * this rotates the coords in X. 
 */
void Frame::Rotate(double *T) {
  int i;
  double x,y,z;
  
  for (i=0; i<N; i+=3) {
    x=X[i]; y=X[i+1]; z=X[i+2];

    X[i  ]=(x*T[0]) + (y*T[1]) + (z*T[2]);
    X[i+1]=(x*T[3]) + (y*T[4]) + (z*T[5]);
    X[i+2]=(x*T[6]) + (y*T[7]) + (z*T[8]);
  }
} 

/* 
 * Frame::Center()
 * Center coordinates in Mask to the coordinates in box[0-2]. When called from
 * Action_Center box will be either 0.0 or box center. Use geometric center if 
 * mass is NULL, otherwise center of mass will be used.
 */
void Frame::Center(AtomMask *Mask, double *box, bool useMassIn) {
  double center[3];

  this->COM(Mask, center, useMassIn);
  //fprintf(stderr,"  FRAME CENTER: %lf %lf %lf\n",center[0],center[1],center[2]); //DEBUG

  // Shift to whatever is in box (origin or center of box in Action_Center) 
  center[0] = box[0] - center[0]; 
  center[1] = box[1] - center[1]; 
  center[2] = box[2] - center[2];
  this->Translate(center);
}

/*
 * Frame::ShiftToCenter()
 * Shift Frame and Ref to common COM
 */
void Frame::ShiftToCenter( Frame *Ref ) {
  double frameCOM[3], refCOM[3];

  // Rotation will occur around geometric center.
  // NOTE could pass in Mass to make true Center of Mass rotation
  this->COM(frameCOM, false);
  Ref->COM(refCOM, false);
  //fprintf(stderr,"  FRAME COM: %lf %lf %lf\n",frameCOM[0],frameCOM[1],frameCOM[2]); //DEBUG
  //fprintf(stderr,"  REF   COM: %lf %lf %lf\n",refCOM[0],refCOM[1],refCOM[2]); //DEBUG
  
  // Shift to common COM
  frameCOM[0]=-frameCOM[0]; frameCOM[1]=-frameCOM[1]; frameCOM[2]=-frameCOM[2];
  refCOM[0]  =-refCOM[0];   refCOM[1]  =-refCOM[1];   refCOM[2]  =-refCOM[2];
  this->Translate(frameCOM);
  Ref->Translate(refCOM);
}

/* 
 * Frame::RMSD()
 * Get the RMSD of this Frame to Ref Frame. Ref frame must contain the same
 * number of atoms as this Frame - should be checked for before this routine
 * is called. Put the best-fit rotation matrix in U and the COM translation 
 * vectors in T. The translation is composed of two XYZ vectors; the first is 
 * the shift of the XYZ coords to origin, and the second is the shift to Ref 
 * origin. To reproduce the fit perform the first translation (Trans[0...2]), 
 * then rotate (U), then the second translation (Trans[3...5]).
 * Adapted from PTRAJ.
 */
double Frame::RMSD( Frame *Ref, double *U, double *Trans, bool useMassIn) {
  double frameCOM[3], refCOM[3], rms_return, total_mass;
  double mwss, rot[9], rtr[9];
  double xt,yt,zt,xr,yr,zr;
  double *Evector[3], Eigenvalue[3], Emat[9];
  double b[9];
  double cp[3], sig3;
  int i, m;
  int i3,k,k3,j;

  memset(U,     0, 9*sizeof(double));
  memset(Trans, 0, 6+sizeof(double));
 
  // Rotation will occur around geometric center.
  // NOTE could pass in Mass to make true Center of Mass rotation
  total_mass = this->COM(frameCOM,useMassIn);
  Ref->COM(refCOM,useMassIn);
  //fprintf(stderr,"  FRAME COM: %lf %lf %lf\n",frameCOM[0],frameCOM[1],frameCOM[2]); //DEBUG
  //fprintf(stderr,"  REF   COM: %lf %lf %lf\n",refCOM[0],refCOM[1],refCOM[2]); //DEBUG

  // Shift to common COM
  frameCOM[0]=-frameCOM[0]; frameCOM[1]=-frameCOM[1]; frameCOM[2]=-frameCOM[2];
  refCOM[0]  =-refCOM[0];   refCOM[1]  =-refCOM[1];   refCOM[2]  =-refCOM[2];
  this->Translate(frameCOM);
  Ref->Translate(refCOM);
  //fprintf(stderr,"  SHIFTED FRAME 0: %lf %lf %lf\n",X[0],X[1],X[2]); //DEBUG
  //fprintf(stderr,"  SHIFTED REF 0  : %lf %lf %lf\n",Ref->X[0],Ref->X[1],Ref->X[2]); //DEBUG

  // Use Kabsch algorithm to calculate optimum rotation matrix.
  // U = [(RtR)^.5][R^-1]
  mwss=0.0;
  memset(rot, 0, 9*sizeof(double));
  memset(rtr, 0, 9*sizeof(double));
  // Calculate covariance matrix of Coords and Reference (R = Xt * Ref)
  m=0;
  for (i=0; i<N; i+=3) {
    xt = X[i];
    yt = X[i+1];
    zt = X[i+2];
    xr = Ref->X[i];
    yr = Ref->X[i+1];
    zr = Ref->X[i+2];

    // Use rms_return to hold mass for this atom if specified
    rms_return = 1.0;
    if (useMassIn) 
      rms_return = Mass[m++];
    //total_mass+=rms_return;

    mwss += rms_return * ( (xt*xt)+(yt*yt)+(zt*zt)+(xr*xr)+(yr*yr)+(zr*zr) );

    // Calculate the Kabsch matrix: R = (rij) = Sum(yni*xnj)
    rot[0] += rms_return*xt*xr;
    rot[1] += rms_return*xt*yr;
    rot[2] += rms_return*xt*zr;

    rot[3] += rms_return*yt*xr;
    rot[4] += rms_return*yt*yr;
    rot[5] += rms_return*yt*zr;

    rot[6] += rms_return*zt*xr;
    rot[7] += rms_return*zt*yr;
    rot[8] += rms_return*zt*zr;
  }
  mwss *= 0.5;    // E0 = 0.5*Sum(xn^2+yn^2) 

  //DEBUG
  //fprintf(stderr,"ROT:\n%lf %lf %lf\n%lf %lf %lf\n%lf %lf %lf\n",
  //        rot[0],rot[1],rot[2],rot[3],rot[4],rot[5],rot[6],rot[7],rot[8]);
  //fprintf(stderr,"MWSS: %lf\n",mwss);

  // calculate Kabsch matrix multiplied by its transpose: RtR 
  rtr[0] = rot[0]*rot[0] + rot[1]*rot[1] + rot[2]*rot[2];
  rtr[1] = rot[0]*rot[3] + rot[1]*rot[4] + rot[2]*rot[5];
  rtr[2] = rot[0]*rot[6] + rot[1]*rot[7] + rot[2]*rot[8];
  rtr[3] = rot[3]*rot[0] + rot[4]*rot[1] + rot[5]*rot[2];
  rtr[4] = rot[3]*rot[3] + rot[4]*rot[4] + rot[5]*rot[5];
  rtr[5] = rot[3]*rot[6] + rot[4]*rot[7] + rot[5]*rot[8];
  rtr[6] = rot[6]*rot[0] + rot[7]*rot[1] + rot[8]*rot[2];
  rtr[7] = rot[6]*rot[3] + rot[7]*rot[4] + rot[8]*rot[5];
  rtr[8] = rot[6]*rot[6] + rot[7]*rot[7] + rot[8]*rot[8];

  // Diagonalize?
  if (!diagEsort(rtr, Emat, Evector, Eigenvalue))
    return(0);

  // a3 = a1 x a2 
  CROSS_PRODUCT(Evector[2][0], Evector[2][1], Evector[2][2],
                Evector[0][0], Evector[0][1], Evector[0][2],
                Evector[1][0], Evector[1][1], Evector[1][2]);

  // Evector dot transpose rot: b = R . ak 
  b[0] = Evector[0][0]*rot[0] + Evector[0][1]*rot[3] + Evector[0][2]*rot[6];
  b[1] = Evector[0][0]*rot[1] + Evector[0][1]*rot[4] + Evector[0][2]*rot[7];
  b[2] = Evector[0][0]*rot[2] + Evector[0][1]*rot[5] + Evector[0][2]*rot[8];
  normalize(&b[0]);
  b[3] = Evector[1][0]*rot[0] + Evector[1][1]*rot[3] + Evector[1][2]*rot[6];
  b[4] = Evector[1][0]*rot[1] + Evector[1][1]*rot[4] + Evector[1][2]*rot[7];
  b[5] = Evector[1][0]*rot[2] + Evector[1][1]*rot[5] + Evector[1][2]*rot[8];
  normalize(&b[3]);
  b[6] = Evector[2][0]*rot[0] + Evector[2][1]*rot[3] + Evector[2][2]*rot[6];
  b[7] = Evector[2][0]*rot[1] + Evector[2][1]*rot[4] + Evector[2][2]*rot[7];
  b[8] = Evector[2][0]*rot[2] + Evector[2][1]*rot[5] + Evector[2][2]*rot[8];
  normalize(&b[6]);

 /* b3 = b1 x b2 */
  CROSS_PRODUCT(cp[0], cp[1], cp[2],
                 b[0],  b[1],  b[2],
                 b[3],  b[4],  b[5]);

  if ( (cp[0]*b[6] + cp[1]*b[7] + cp[2]*b[8]) < 0.0 )
    sig3 = -1.0;
  else
    sig3 = 1.0;

  b[6] = cp[0];
  b[7] = cp[1];
  b[8] = cp[2];

  // U has the best rotation 
  for (k=k3=0; k<3; k++,k3+=3)
    for (i=i3=0; i<3; i++,i3+=3)
      for (j=0; j<3; j++)
        U[i3+j] += Evector[k][j] * b[k3+i]; 

  // E=E0-sqrt(mu1)-sqrt(mu2)-sig3*sqrt(mu3) 
  rms_return = mwss
               - sqrt(fabs(Eigenvalue[0]))
               - sqrt(fabs(Eigenvalue[1]))
               - sig3*sqrt(fabs(Eigenvalue[2]));

  if (rms_return<0) {
    //fprintf(stderr,"RMS returned is <0 before sqrt, setting to 0 (%lf)\n",rms_return);
    rms_return=0.0;
  }
  else
    rms_return = sqrt((2.0*rms_return)/total_mass);

  /* Translation matrix - coords are shifted to common CoM first (origin), then
   * to original reference location.
   * Remember frameCOM and refCOM were negated above to facilitate translation to COM.
   */
  Trans[0] = frameCOM[0];
  Trans[1] = frameCOM[1];
  Trans[2] = frameCOM[2];
  Trans[3] = -refCOM[0];
  Trans[4] = -refCOM[1];
  Trans[5] = -refCOM[2];

  //DEBUG
  //printRotTransInfo(U,Trans);
  //fprintf(stdout,"RMS is %lf\n",rms_return);

  return rms_return;
}


/*
 * Frame::RMSD()
 * Calculate RMSD of Frame to Ref with no fitting. Frames must contain
 * same # atoms.
 */
double Frame::RMSD( Frame *Ref, bool useMass ) {
  double rms_return, total_mass, xx,yy,zz, currentMass;
  int i, m;

  rms_return = 0.0;
  total_mass = 0.0;
  currentMass = 1.0;
  m=0;

  for (i=0; i < N; i+=3) {
    if (useMass) 
      currentMass = Mass[m++];
    total_mass += currentMass;
    xx = Ref->X[i]   - X[i];
    yy = Ref->X[i+1] - X[i+1];
    zz = Ref->X[i+2] - X[i+2];
    rms_return += currentMass * (xx*xx + yy*yy + zz*zz);
  }
  rms_return = sqrt(rms_return / total_mass);

  return rms_return;
}
