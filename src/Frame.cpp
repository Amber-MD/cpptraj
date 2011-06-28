#include <cstdlib>
#include <cstdio> // sprintf
#include <cmath>
#include <cstring> // memset
#include "Frame.h"
#include "vectormath.h"
#include "DistRoutines.h"
#include "TorsionRoutines.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Frame::Frame() {
  natom=0;
  N=0;
  X=NULL;
  box[0]=0; box[1]=0; box[2]=0; box[3]=0; box[4]=0; box[5]=0;
  T=0.0;
  V=NULL;
  Mass=NULL;
}

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
}

// DESTRUCTOR
Frame::~Frame() {
  if (X!=NULL) free(X);
  if (V!=NULL) delete V;
  if (Mass!=NULL) free(Mass);
}

/* Frame::ZeroCoords()
 * Set all coords to 0.0
 */
void Frame::ZeroCoords() {
  for (int coord=0; coord < N; coord++)
    X[coord]=0.0;
}

/* Frame::AddCoord()
 * Add the coord values from the input frame to the coord values of 
 * this frame.
 */
void Frame::AddCoord(Frame *FrameIn) {
  if (FrameIn->N != this->N) {
    mprintf("Error: Frame::AddCoord: Attempting to add %i coords to %i coords.\n",
            FrameIn->N,this->N);
  } else {
    for (int coord=0; coord < N; coord++)
      this->X[coord] += FrameIn->X[coord];
  }
}

/* Frame::Divide()
 * Divide all coord values by input. Dont do it if the number is too small.
 */
void Frame::Divide(double divisor) {
  if (divisor < SMALL) return;
  for (int coord=0; coord < N; coord++)
    X[coord] /= divisor;
}

/* Frame::Copy()
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
 * Frame::GetCoord()
 * Get coordinates of specified atom and put into Coord.
 */
void Frame::GetCoord(double *Coord, int atom) {
  int natom3;
  // NOTE: SHOULD CHECK FOR BOUNDARIES
  if (Coord==NULL) return;
  natom3 = atom * 3;
  Coord[0] = X[natom3  ];
  Coord[1] = X[natom3+1];
  Coord[2] = X[natom3+2];
}

/*
 * Frame::SetCoord()
 * Set coordinates of specified atom to those of Coord.
 */
void Frame::SetCoord(int atom, double *Coord) {
  int natom3;
  // NOTE: SHOULD CHECK FOR BOUNDARIES
  if (Coord==NULL) return;
  //mprintf("      Setting coord %i/%i: %lf %lf %lf\n",atom,natom,Coord[0],Coord[1],Coord[2]);
  natom3 = atom * 3;
  X[natom3  ] = Coord[0];
  X[natom3+1] = Coord[1];
  X[natom3+2] = Coord[2];
}

/*
 * Frame::Coord()
 * Return double pointer to XYZ coord of given atom.
 */
double *Frame::Coord(int atom) {
  if (atom<0 || atom>=natom) return NULL;
  return ( X+(atom*3) );
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

  // NOTE: Not using == since it is unreliable for floating point numbers.
  // Should NEVER have a mass smaller than SMALL (vectormath.h)
  if (sumMass < SMALL) return 0;

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

  // NOTE: Not using == since it is unreliable for floating point numbers.
  // Should NEVER have a mass smaller than SMALL (vectormath.h)
  if (sumMass < SMALL) return 0;

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
void Frame::BoxToRecip(double *ucell, double *recip) {
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
 * Frame::DIST2()
 * Call the appropriate distance calc for atoms in Mask1 and Mask2 based on
 * given box type.
 *   0 = None
 *   1 = Orthorhombic
 *   2 = Non-orthorhombic
 * Based on useMassIn, calculate geometric center (false) or center of mass 
 * (true) of the atoms in each mask.
 */
double Frame::DIST2(AtomMask *Mask1, AtomMask *Mask2, bool useMassIn, int boxType,
                    double *ucell, double *recip) {
  double a1[3], a2[3];

  COM(Mask1, a1, useMassIn);
  COM(Mask2, a2, useMassIn);

  if (boxType == 0) 
    return DIST2_NoImage(a1, a2);
  else if (boxType == 1) 
    return DIST2_ImageOrtho(a1, a2, this->box);
  else if (boxType == 2) 
    return DIST2_ImageNonOrtho(a1, a2, ucell, recip);

  mprintf("    Error: Frame::DIST: Unrecognized box type (%i)\n.", boxType);

  return (-1.0);
}

/*
 * Frame::DIST2()
 * Return the distance between atoms A1 and A2 with optional imaging.
 *   0 = None
 *   1 = Orthorhombic
 *   2 = Non-orthorhombic
 */
double Frame::DIST2(int A1, int A2, int boxType, double *ucell, double *recip) {
  int atom3;
  double a1[3], a2[3];

  atom3 = A1 * 3;
  a1[0] = X[atom3  ];
  a1[1] = X[atom3+1];
  a1[2] = X[atom3+2];
  atom3 = A2 * 3;
  a2[0] = X[atom3  ];
  a2[1] = X[atom3+1];
  a2[2] = X[atom3+2];

  if (boxType == 0)
    return DIST2_NoImage(a1, a2);
  else if (boxType == 1)
    return DIST2_ImageOrtho(a1, a2, this->box);
  else if (boxType == 2) 
    return DIST2_ImageNonOrtho(a1, a2, ucell, recip);

  mprintf("    Error: Frame::DIST: Unrecognized box type (%i)\n.", boxType);

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

  COM(M1,a1,false);
  COM(M2,a2,false);
  COM(M3,a3,false);
  COM(M4,a4,false);

  return Torsion(a1,a2,a3,a4);
}

/*
 * Frame::PUCKER()
 * Return the pseudorotation between atoms in masks M1-M5 for the given
 * puckerMethod:
 *   0: Use Altona & Sundaralingam method/conventions
 *   1: Use Cremer & Pople method
 * If amplitude is true, return amplitude instead of pseudorotation.
 */
double Frame::PUCKER(AtomMask *M1, AtomMask *M2, AtomMask *M3, AtomMask *M4, AtomMask *M5,
                     int puckerMethod, bool amplitude, bool useMassIn) {
  double a1[3],a2[3],a3[3],a4[3],a5[3]; 
  double angle, amp;

  COM(M1,a1,useMassIn);
  COM(M2,a2,useMassIn);
  COM(M3,a3,useMassIn);
  COM(M4,a4,useMassIn);
  COM(M5,a5,useMassIn);

  angle = 0.0;
  amp = 0.0;
  switch (puckerMethod) {
    case 0 : angle = Pucker_AS(a1,a2,a3,a4,a5,&amp); break;
    case 1 : angle = Pucker_CP(a1,a2,a3,a4,a5,&amp); break;
  }

  if (amplitude) return amp;
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

  // NOTE: Not using == since it is unreliable for floating point numbers.
  // Should NEVER have a mass smaller than SMALL (vectormath.h)
  if (total_mass < SMALL) return 0;

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
 * Translate atom in X by V.
 * NOTE: SHOULD CHECK BOUNDS! 
 */
void Frame::Translate(double *V, int Atom) {
  int atom3;

  atom3 = Atom * 3;

  X[atom3  ] += V[0];
  X[atom3+1] += V[1];
  X[atom3+2] += V[2];
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
 * Frame::InverseRotate()
 * Multiply natomx3 matrix X by transpose of 3x3 matrix T. If T is a rotation
 * matrix this rotates the coords in X in the opposite direction.
 */
void Frame::InverseRotate(double *T) {
  int i;
  double x,y,z;

  for (i=0; i<N; i+=3) {
    x=X[i]; y=X[i+1]; z=X[i+2];

    X[i  ]=(x*T[0]) + (y*T[3]) + (z*T[6]);
    X[i+1]=(x*T[1]) + (y*T[4]) + (z*T[7]);
    X[i+2]=(x*T[2]) + (y*T[5]) + (z*T[8]);
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
  memset(Trans, 0, 6*sizeof(double));
 
  // Rotation will occur around geometric center/center of mass
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
               - (sig3*sqrt(fabs(Eigenvalue[2])));

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
