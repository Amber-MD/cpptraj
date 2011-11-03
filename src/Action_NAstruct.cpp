/* Action_NAstruct
 * Calculate nucleic acid base/base pair structural parameters.
 * Algorithms for calculation of base/base pair structural parameters
 * adapted from:
 *   Babcock MS, Pednault EPD, Olson WK, "Nucleic Acid Structure Analysis: 
 *   Mathematics for Local Cartesian and Helical Structure Parameters That
 *   Are Truly Comparable Between Structures", J. Mol. Biol. (1994) 237,
 *   125-156.
 * NA base reference frame coordinates taken from:
 *   Olson WK, Bansal M, Burley SK, Dickerson RE, Gerstein M, Harvey SC,
 *   Heinemann U, Lu XJ, Neidle S, Shekked Z, Sklenar H, Suzuki M, Tung CS,
 *   Westhof E, Wolberger C, Berman H, "A Standard Reference Frame for the 
 *   Description of Nucleic Acid Base-pair Geometry", J. Mol. Biol. (2001)
 *   313, 229-237.
 */
#include "Action_NAstruct.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "Constants.h" // RADDEG
#include "vectormath.h"
#include <cmath>
#include <cstdio> // sprintf

// CONSTRUCTOR
NAstruct::NAstruct() {
  //fprintf(stderr,"NAstruct Con\n");
  Nbp=0;
  Nbases=0;
  HBcut2=12.25; // 3.5^2
  Ocut2=6.25;   // 2.5^2
  Nframe=0;
  outFilename=NULL;
  naoutFilename=NULL;
  noheader = false;
} 

// DESTRUCTOR
NAstruct::~NAstruct() { 
  ClearLists();
  while (!BasePairAxes.empty()) {
    delete BasePairAxes.back();
    BasePairAxes.pop_back();
  }
}

/* NAstruct::ClearLists()
 * Clear all parm-dependent lists
 */
void NAstruct::ClearLists() {
  while (!RefCoords.empty()) {
    delete RefCoords.back();
    RefCoords.pop_back();
  }
  while (!BaseAxes.empty()) {
    delete BaseAxes.back();
    BaseAxes.pop_back();
  }
// NOTE: Since BasePairAxes are set up to correspond with SHEAR etc dont
// free in this routine - should only be freed at the very end.
//  while (!BasePairAxes.empty()) {
//    delete BasePairAxes.back();
//    BasePairAxes.pop_back();
//  }
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

/* NAstruct::GCpair()
 * Look for 3 HB based on heavy atom distances:
 * 1. G:O6 -- C:N4  6 -- 6
 * 2. G:N1 -- C:N3  7 -- 4
 * 3. G:N2 -- C:O2  9 -- 3
 * Atom positions are known in standard Ref. Multiply by 3 to get into X.
 */
bool NAstruct::GCpair(AxisType *DG, AxisType *DC) {
  int Nhbonds = 0;
  double dist2; 
  dist2 = DIST2_NoImage(DG->X+18, DC->X+18);
  if ( dist2 < HBcut2 ) {
    Nhbonds++;
    //mprintf("            G:O6 -- C:N4 = %lf\n",sqrt(dist2));
  }
  dist2 = DIST2_NoImage(DG->X+21, DC->X+12);
  if ( dist2 < HBcut2 ) {
    Nhbonds++;
    //mprintf("            G:N1 -- C:N3 = %lf\n",sqrt(dist2));
  }
  dist2 = DIST2_NoImage(DG->X+27, DC->X+9 );
  if ( dist2 < HBcut2 ) {
    Nhbonds++;
    //mprintf("            G:N2 -- C:O2 = %lf\n",sqrt(dist2));
  }
  if (Nhbonds>0) return true;
  return false;
}

/* NAstruct::ATpair()
 * Look for 2 HB based on heavy atom distances
 * 1. A:N6 -- T:O4  6 -- 6
 * 2. A:N1 -- T:N3  7 -- 4
 */
bool NAstruct::ATpair(AxisType *DA, AxisType *DT) {
  int Nhbonds = 0;
  double dist2;
  dist2 = DIST2_NoImage(DA->X+18, DT->X+18);
  if ( dist2 < HBcut2 ) {
    Nhbonds++;
    //mprintf("            A:N6 -- T:O4 = %lf\n",sqrt(dist2));
  }
  dist2 = DIST2_NoImage(DA->X+21, DT->X+12);
  if ( dist2 < HBcut2 ) {
    Nhbonds++;
    //mprintf("            A:N1 -- T:N3 = %lf\n",sqrt(dist2));
  }
  if (Nhbonds>0) return true;
  return false;
}

/* NAstruct::basesArePaired()
 * Given two base axes for which IDs have been given and reference coords set,
 * determine whether the bases are paired via hydrogen bonding criteria.
 * NOTE: Currently only set up for Antiparallel WC detection
 */
bool NAstruct::basesArePaired(AxisType *base1, AxisType *base2) {
  // G C
  if      ( base1->ID==DG && base2->ID==DC ) return GCpair(base1,base2);
  else if ( base1->ID==DC && base2->ID==DG ) return GCpair(base2,base1);
  else if ( base1->ID==RG && base2->ID==RC ) return GCpair(base1,base2);
  else if ( base1->ID==RC && base2->ID==RG ) return GCpair(base2,base1);
  // A T
  else if ( base1->ID==DA && base2->ID==DT ) return ATpair(base1,base2);
  else if ( base1->ID==DT && base2->ID==DA ) return ATpair(base2,base1);
  // A U
  else if ( base1->ID==RA && base2->ID==RU ) return ATpair(base1,base2);
  else if ( base1->ID==RU && base2->ID==RA ) return ATpair(base2,base1);
//  else {
//    mprintf("Warning: NAstruct: Unrecognized pair: %s - %s\n",NAbaseName[base1->ID],
//             NAbaseName[base2->ID]);
//  }
  return false;
}

/* NAstruct::determineBasePairing()
 * Determine which bases are paired from the individual base axes.
 */
int NAstruct::determineBasePairing() {
  double distance;
  std::vector<bool> isPaired( BaseAxes.size(), false);
  int base1,base2;
  double Z1[3], Z2[3];
  bool AntiParallel = false;

  Nbp = 0;
  BasePair.clear();
  
  //mprintf(" ==== Setup Base Pairing ==== \n");

  /* For each unpaired base, determine if it is paired with another base
   * determined by the distance between their axis origins.
   */
  for (base1=0; base1 < Nbases-1; base1++) {
    if (isPaired[base1]) continue;
    for (base2=base1+1; base2 < Nbases; base2++) {
      if (isPaired[base2]) continue;
      // First determine if origin axes coords are close enough to consider pairing
      distance = DIST2_NoImage(BaseAxes[base1]->Origin(), BaseAxes[base2]->Origin());
      /*mprintf("  Axes distance for %i:%s -- %i:%s is %lf\n",
              base1,RefCoords[base1]->BaseName(),
              base2,RefCoords[base2]->BaseName(),sqrt(distance));*/
      if (distance < Ocut2) {
        //mprintf("    Checking %i:%s -- %i:%s\n",base1,RefCoords[base1]->BaseName(),
        //        base2,RefCoords[base2]->BaseName());
        // Figure out if z vectors point in same (<90 deg) or opposite (>90 deg) direction
        BaseAxes[base1]->RZ(Z1);
        BaseAxes[base2]->RZ(Z2);
        //printVector("Base1Z",Z1); printVector("Base2Z",Z2);
        distance = dot_product_angle(Z1, Z2);
        //mprintf("    Dot product of Z vectors: %lf\n",distance);
        if (distance > (PIOVER2)) { // If theta(Z) > 90 deg.
          //mprintf("      Base2 %i is anti-parallel to Base1 %i\n",base2,base1);
          AntiParallel = true;
        } else {
          //mprintf("      Base2 %i is parallel to Base1 %i\n",base2,base1);
          AntiParallel = false;
        }
        if (basesArePaired(RefCoords[base1], RefCoords[base2])) {
          BasePair.push_back(base1);
          BasePair.push_back(base2);
          if (AntiParallel) 
            BasePair.push_back(1);
          else
            BasePair.push_back(0);
          isPaired[base1]=true;
          isPaired[base2]=true;
          Nbp++;
        }
      }
    } // END Loop over base2
  } // END Loop over base1

  if (debug>0) mprintf("    NAstruct: Set up %i base pairs.\n",Nbp);
  base2=1;
  //mprintf("DEBUG: BasePair.size = %i\n",(int)BasePair.size());
  //mprintf("DEBUG: SHEAR.size = %i\n",(int)SHEAR.size());
  //mprintf("DEBUG: BasePairAxes.size = %i\n",(int)BasePairAxes.size());
  DataSet *na_dataset = NULL;
  for (base1=0; base1 < (int)BasePair.size(); base1+=3) {
    // For each base pair, set up a dataset for each structural parameter
    // if one is not already set up.
    if ( base2 > SHEAR.Size() ) {
      na_dataset = SHEAR.AddMultiN(DOUBLE, "SHR", "BP", base2);
      DFL->Add(outFilename, na_dataset);
      na_dataset = STRETCH.AddMultiN(DOUBLE, "STR", "BP", base2);
      DFL->Add(outFilename, na_dataset);
      na_dataset = STAGGER.AddMultiN(DOUBLE, "STA", "BP", base2);
      DFL->Add(outFilename, na_dataset);
      na_dataset = BUCKLE.AddMultiN(DOUBLE, "BCK", "BP", base2);
      DFL->Add(outFilename, na_dataset);
      na_dataset = PROPELLER.AddMultiN(DOUBLE, "PRP", "BP", base2);
      DFL->Add(outFilename, na_dataset);
      na_dataset = OPENING.AddMultiN(DOUBLE, "OPN", "BP", base2);
      DFL->Add(outFilename, na_dataset);
      // Also set up a place to hold the base pair axes
      BasePairAxes.push_back( new AxisType() );
      BasePairAxes[base2-1]->SetPrincipalAxes();
    } 
    // Print base pair info
    if (debug>1) {
      mprintf("        BP %i: Res %i:%s to %i:%s",base2,
              BasePair[base1  ]+1, RefCoords[ BasePair[base1  ] ]->BaseName(),
              BasePair[base1+1]+1, RefCoords[ BasePair[base1+1] ]->BaseName());
      if ( BasePair[base1+2] )
        mprintf(" AntiParallel.\n");
      else
        mprintf(" Parallel.\n");
    }
    base2++;
  }
  // Also set up base pair step data. One less than # base pairs
  //mprintf("DEBUG: SHIFT.size() = %i\n",(int)SHIFT.size());
  for (base1=0; base1 < Nbp-1; base1++) {
    if ( base1+1 > SHIFT.Size() ) {
      na_dataset = SHIFT.AddMultiN(DOUBLE, "SHF", "BS", base1 + 1);
      DFL->Add(outFilename, na_dataset);
      na_dataset = SLIDE.AddMultiN(DOUBLE, "SLD", "BS", base1 + 1);
      DFL->Add(outFilename, na_dataset);
      na_dataset = RISE.AddMultiN(DOUBLE, "RIS", "BS", base1 + 1);
      DFL->Add(outFilename, na_dataset);
      na_dataset = TILT.AddMultiN(DOUBLE, "TLT", "BS", base1 + 1);
      DFL->Add(outFilename, na_dataset);
      na_dataset = ROLL.AddMultiN(DOUBLE, "RLL", "BS", base1 + 1);
      DFL->Add(outFilename, na_dataset);
      na_dataset = TWIST.AddMultiN(DOUBLE, "TWS", "BS", base1 + 1);
      DFL->Add(outFilename, na_dataset);
    }
    // Print base pair step info
    if (debug>1) mprintf("        BP step %i: \n",base1+1);
  }

  return 0;
}

/* CURRENTLY NOT USED
 * NAstruct::setupBasePairAxes()
 * Given a list of base pairs and base axes, setup an 
 * Axestype structure containing reference base pair axes.
 */
/*
int NAstruct::setupBasePairAxes() {
  int basepair, BP;
  double RotMatrix[9], TransVec[6], V[3], theta;
  AxisType *Base1;
  AxisType *Base2;
  // DEBUG
  int basepairaxesatom=0;
  CpptrajFile basepairaxesfile;
  basepairaxesfile.SetupFile((char*)"basepairaxes.pdb",WRITE,UNKNOWN_FORMAT,UNKNOWN_TYPE,0);
  basepairaxesfile.OpenFile();
  // END DEBUG

  mprintf(" ==== Setup Base Pair Axes ==== \n");
  // Loop over all base pairs
  BP = 0;
  for (basepair = 0; basepair < (int)BasePair.size(); basepair+=3) {
    // Set Base1 as a copy, this will become the base pair axes
    Base1 = new AxisType();
    Base1->SetPrincipalAxes();
    Base1->SetFromFrame( BaseAxes[ BasePair[basepair] ] );
    Base2 = BaseAxes[ BasePair[basepair+1] ];
    // Set frame coords for first base in the pair
    RefFrame.SetFromFrame( Base1 );
    // Set frame coords for second base in the pair
    ExpFrame.SetFromFrame( Base2 );
    // Flip the axes in the second pair
    // NOTE: This flip is only correct for standard WC base-pairing, not explicitly checked
    ExpFrame.FlipYZ();
    // RMS fit Axes of Base2 onto Base1
    ExpFrame.RMSD( &RefFrame, RotMatrix, TransVec, false);
    // DEBUG
    //mprintf("  BP: %i Rotation matrix/Translation vector:\n",BP+1);
    //printRotTransInfo(RotMatrix, TransVec);
    // Extract angle from resulting rotation matrix
    theta=matrix_to_angle(RotMatrix);
    mprintf("Base2 for pair %i will be rotated by %lf degrees.\n",BP+1,RADDEG*theta);
    // Calc Axis of rotation
    if (axis_of_rotation(V, RotMatrix, theta)) {
      mprintf("Error: NAstruct::setupBasePairAxes(): Could not set up axis of rotation for %i.\n",
              basepair);
      delete Base1;
      return 1;
    }
    // Calculate new half-rotation matrix
    calcRotationMatrix(RotMatrix,V,theta/2);
    printMatrix("Rhalf",RotMatrix);
    // Rotate Base2 by half rotation towards Base1.
    // Since the new rotation axis by definition is located at
    // the origin, the coordinates of the base have to be shifted, rotated,
    // then shifted back. Use V to store the shift back.
    V[0] = -TransVec[0];
    V[1] = -TransVec[1];
    V[2] = -TransVec[2];
    Base2->Translate( TransVec );
    Base2->Rotate( RotMatrix );
    Base2->Translate( V );
    // Since rotation matrix for Base2 was calculated with same origin as 
    // Base1, use reverse rotation matrix (transpose) to rotate Base1. 
    // Shift, rotate, shift back. Use V to store the shift back.
    V[0] = -TransVec[3];
    V[1] = -TransVec[4];
    V[2] = -TransVec[5];
    Base1->Translate( V );
    Base1->InverseRotate( RotMatrix );
    Base1->Translate( TransVec+3 );
    // Origin of base pair axes is midpoint between Base1 and Base2 origins
    V[0] = ( (Base1->X[9 ] + Base2->X[9 ])/2 ) - Base1->X[9 ];
    V[1] = ( (Base1->X[10] + Base2->X[10])/2 ) - Base1->X[10];
    V[2] = ( (Base1->X[11] + Base2->X[11])/2 ) - Base1->X[11];
    // Shift Base1 to midpoint; Base1 becomes the base pair axes
    Base1->Translate( V );
    // X, Y, and Z unit vectors compose the base pair rotation matrix
    Base1->R[0] = Base1->X[0] - Base1->X[9];
    Base1->R[3] = Base1->X[1] - Base1->X[10];
    Base1->R[6] = Base1->X[2] - Base1->X[11];
    Base1->R[1] = Base1->X[3] - Base1->X[9];
    Base1->R[4] = Base1->X[4] - Base1->X[10];
    Base1->R[7] = Base1->X[5] - Base1->X[11]; 
    Base1->R[2] = Base1->X[6] - Base1->X[9];
    Base1->R[5] = Base1->X[7] - Base1->X[10];
    Base1->R[8] = Base1->X[8] - Base1->X[11];
    // NOTE: Axes are by definition already normalized
    mprintf("      %i) %i:%s -- %i:%s  %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf %8.2lf\n",BP,
            BasePair[basepair  ], Base1->BaseName(),
            BasePair[basepair+1], Base2->BaseName(),
            Base1->X[9], Base1->X[10], Base1->X[11],
            Base1->R[2], Base1->R[5], Base1->R[8]);
    printMatrix("BP R",Base1->R);
    // Base1 contains absolute coords and rotation matrix of base pair reference axes
    Base1->WritePDB(&basepairaxesfile, BasePair[basepair], P->ResidueName(BP),&basepairaxesatom);
    BasePairAxes.push_back( Base1 );
    BP++;
  }
  // DEBUG
  basepairaxesfile.CloseFile();

  return 0;
}
*/

/* NAstruct::setupBaseAxes()
 * For each residue defined in reference coords, get the corresponding input
 * coords and fit the reference coords (and reference axes) on top of input 
 * coords. This sets up the reference axes for each base.
 */
int NAstruct::setupBaseAxes(Frame *InputFrame) {
  double rmsd, RotMatrix[9], TransVec[6];
  int base;
  // DEBUG
  //int res = 0;
  //int baseaxesatom = 0;
  //int basesatom = 0;
  //CpptrajFile baseaxesfile;
  //CpptrajFile basesfile;
  //baseaxesfile.SetupFile((char*)"baseaxes.pdb",WRITE,UNKNOWN_FORMAT,UNKNOWN_TYPE,0);
  //baseaxesfile.OpenFile();
  //basesfile.SetupFile((char*)"bases.pdb",WRITE,UNKNOWN_FORMAT,UNKNOWN_TYPE,0);
  //basesfile.OpenFile();
  // END DEBUG

  // For each axis in RefCoords, use corresponding mask in ExpMasks to set 
  // up an axis for ExpCoords.
  for (base=0; base < Nbases; base++) {
    // Set exp coords based on previously set-up mask
    ExpFrames[base]->SetFrameCoordsFromMask( InputFrame->X, ExpMasks[base] ); 
    /* Now that we have a set of reference coords and the corresponding input
     * coords, RMS fit the reference coords to the input coords to obtain the
     * appropriate rotation and translations that will put the reference coords 
     * on top of input (experimental) coords.
     * NOTE: The RMSD routine is destructive to coords. Need copies of frames.
     */
    RefFrame.SetFromFrame( RefCoords[base] );
    ExpFrame.SetFromFrame( ExpFrames[base] );
    rmsd = RefFrame.RMSD( &ExpFrame, RotMatrix, TransVec, false);
    if (debug>0) 
      mprintf("Base %i: RMS of RefCoords from ExpCoords is %lf\n",base+1,rmsd);
    // BaseAxes start at origin
    BaseAxes[base]->SetPrincipalAxes();
    // Store the Rotation matrix.
    BaseAxes[base]->StoreRotMatrix( RotMatrix );
    // DEBUG
    //mprintf("         Rotation matrix/Translation vector:\n");
    //printRotTransInfo(RotMatrix, TransVec);
    //AxisToPDB(&baseaxesfile, (*baseaxis), res++, &baseaxesatom);
    /* RotMatrix and TransVec now contain rotation and translation
     * that will orient refcoord to expframe. The first translation is that of
     * the reference frame to the absolute origin, the second translation is
     * that of the reference frame to the exp. coords after rotation.
     * The rotation matrix essentially contains the absolute coordinates of the
     * X, Y, and Z unit vectors of the base reference coordinates.
     */

    // Use the translation/rotation to fit principal axes in BaseAxes to experimental coords.
    BaseAxes[base]->Trans_Rot_Trans(TransVec, RotMatrix);
    // This BaseAxis now contains the absolute coordinates of the base reference axes.
    
    // DEBUG
    if (debug>0) {
      mprintf("         BaseAxes origin: %8.4lf %8.4lf %8.4lf\n",
              BaseAxes[base]->X[9],BaseAxes[base]->X[10],BaseAxes[base]->X[11]);
      mprintf("         BaseAxes X vec.: %8.4lf %8.4lf %8.4lf\n",
              BaseAxes[base]->X[0 ]-BaseAxes[base]->X[9 ],
              BaseAxes[base]->X[1 ]-BaseAxes[base]->X[10],
              BaseAxes[base]->X[2 ]-BaseAxes[base]->X[11]);
      mprintf("         BaseAxes Y vec.: %8.4lf %8.4lf %8.4lf\n",
              BaseAxes[base]->X[3 ]-BaseAxes[base]->X[9 ],
              BaseAxes[base]->X[4 ]-BaseAxes[base]->X[10],
              BaseAxes[base]->X[5 ]-BaseAxes[base]->X[11]);
      mprintf("         BaseAxes Z vec.: %8.4lf %8.4lf %8.4lf\n",
              BaseAxes[base]->X[6 ]-BaseAxes[base]->X[9 ],
              BaseAxes[base]->X[7 ]-BaseAxes[base]->X[10],
              BaseAxes[base]->X[8 ]-BaseAxes[base]->X[11]);
    }
    // DEBUG - Write base axis to file
    //BaseAxes[base]->WritePDB(&baseaxesfile, res, P->ResidueName(res), &baseaxesatom);

    // Overlap ref coords onto input coords. Rotate, then translate to baseaxes origin
    //RefFrame.SetFromFrame( RefCoords[base] );
    //RefFrame.Rotate( RotMatrix );
    //RefFrame.Translate( BaseAxes[base]->Origin() );
    // DEBUG - Write ref coords to file
    //RefCoords[base]->WritePDB(&basesfile, res, P->ResidueName(res), &basesatom);
    //res++;
  }

  // DEBUG
  //baseaxesfile.CloseFile();
  //basesfile.CloseFile();

  return 0;
}

/* FLIP_V
 * Multiply each element of vector V by corresponding element of flip vector.
 */
#define FLIP_V( V, F ) { \
  V[0]*=F[0]; \
  V[1]*=F[1]; \
  V[2]*=F[2]; } 
/* NAstruct::determineBaseParameters()
 * For each base in a base pair, get the values of buckle, propeller twist,
 * opening, shear, stretch, and stagger. Also determine the origin and 
 * rotation matrix for each base pair reference frame.
 */
int NAstruct::determineBaseParameters() {
  int base1, base2, BP, nbasepair;
  AxisType *Base1, *Base2;
  double V1[3], V2[3], Flip[3];
  double X1[3], X2[3], Y1[3], Y2[3], Z1[3], Z2[3];
  double Phi, dpX, dpY, dpZ;
  double numerator, denominator;
  double Kappa, Omega, Sigma, Sign, absK, absO, absS;
  double y1z2, z1y2, z1x2, x1z2, x1y2, y1x2;
  double O21[3], R1V1[3], Vec[3], FV2[3];
  double Rhalf[9], R2t[9], Rb[9];
  double Shear, Stretch, Stagger;
  // DEBUG
  //int basepairaxesatom=0;
  //CpptrajFile basepairaxesfile;
  //basepairaxesfile.SetupFile((char*)"basepairaxes.pdb",WRITE,UNKNOWN_FORMAT,UNKNOWN_TYPE,0);
  //basepairaxesfile.OpenFile();
  // END DEBUG

  // Default pivot points
  V1[0]=0.0; V1[1]=1.808; V1[2]=0.0;
  V2[0]=0.0; V2[1]=1.808; V2[2]=0.0;

  nbasepair=0;
  for (BP=0; BP < (int)BasePair.size(); BP+=3) {
    base1 = BasePair[BP  ];
    base2 = BasePair[BP+1];
    //mprintf("\n*** Determining base parameters for pair %i -- %i\n",base1+1,base2+1);
    Base1 = BaseAxes[base1];
    Base2 = BaseAxes[base2];
    // Fill X, Y, and Z axis unit vectors from Rotation matrices
    Base1->RX(X1);
    Base1->RY(Y1);
    Base1->RZ(Z1);
    Base2->RX(X2);
    Base2->RY(Y2);
    Base2->RZ(Z2);
    // Calculate dot products common to both parallel and anti-parallel bases  
    dpX = dot_product(X1,X2);
    dpY = dot_product(Y1,Y2);
    dpZ = dot_product(Z1,Z2);
    // Dot Products for sign determination
    y1z2 = dot_product(Y1, Z2);
    z1y2 = dot_product(Z1, Y2);
    z1x2 = dot_product(Z1, X2);
    x1z2 = dot_product(X1, Z2);
    x1y2 = dot_product(X1, Y2);
    y1x2 = dot_product(Y1, X2);

    // -------------------- Anti-parallel --------------------
    if (BasePair[BP+2]) {
      Phi = acos( 0.5 * ( dpX - dpY - dpZ - 1.0 ) );
      //mprintf("PHI= %lf\n",Phi*RADDEG);
      denominator = ( dpX - dpY - dpZ - 3.0 );
      // Kappa (Buckle)
      numerator = ( - dpX - dpY - dpZ - 1.0 );
      Kappa = Phi * sqrt( fabs( numerator / denominator ) );
      absK = fabs(Kappa);
      // Omega (Propeller Twist)
      numerator = (   dpX + dpY - dpZ - 1.0 );
      Omega = Phi * sqrt( fabs( numerator / denominator ) );
      absO = fabs(Omega);
      // Sigma (Opening)
      numerator = (   dpX - dpY + dpZ - 1.0 ); 
      Sigma = Phi * sqrt( fabs( numerator / denominator ) );
      absS = fabs(Sigma);
      // Debug
      //mprintf("DEBUG: K/O/S = %lf %lf %lf\n",Kappa,Omega,Sigma);
      // Kappa Sign determination
      Sign=-1.0;
      if ( absK>=absO && absK>=absS && y1z2<=z1y2 ) Sign=1.0;
      else if   ( absO>=absK && absO>=absS ) {
          if      ( z1x2>=(-x1z2) && (x1y2-y1x2)<=0 ) Sign=1.0;
          else if ( z1x2 <(-x1z2) && (x1y2-y1x2)> 0 ) Sign=1.0;
      } else if ( absS>=absK && absS>=absO ) {
          if      ( x1y2<=(-y1x2) && (x1z2-z1x2)<=0 ) Sign=1.0;
          else if ( x1y2 >(-y1x2) && (x1z2-z1x2)> 0 ) Sign=1.0;
      }
      Kappa *= Sign;
      //mprintf("        Kappa = %8.2lf Sign=%lf\n",Kappa*RADDEG,Sign);
      // Omega Sign determination
      // NOTE: In Table 2 of Babcock et. al (1994) the final sign determination
      //       rule is listed as x1y2>-y1x2 and y1z2+z1x2>0. The final term
      //       appears to contain a typo and has been changed to :
      //          y1z2+z1y2>0 
      Sign=-1.0;
      if        ( absK>=absO && absK>=absS ) {
          if      ( y1z2<=z1y2 && (x1y2-y1x2)<=0 ) Sign=1.0;
          else if ( y1z2> z1y2 && (x1y2-y1x2)> 0 ) Sign=1.0;
      } else if ( absO>=absK && absO>=absS && z1x2>=(-x1z2) ) Sign=1.0;
      else if   ( absS>=absK && absS>=absO ) {
          if      ( x1y2<=(-y1x2) && (y1z2+z1y2)<=0 ) Sign=1.0;
          else if ( x1y2> (-y1x2) && (y1z2+z1y2)> 0 ) Sign=1.0;
      }
      Omega *= Sign;
      //mprintf("        Omega = %8.2lf Sign=%lf\n",Omega*RADDEG,Sign);
      // Sigma Sign determination
      Sign=-1.0;
      if        ( absK>=absO && absK>=absS ) {
          if      ( y1z2<=z1y2 && (x1z2-z1x2)<=0 ) Sign=1.0;
          else if ( y1z2> z1y2 && (x1z2-z1x2)> 0 ) Sign=1.0;
      } else if ( absO>=absK && absO>=absS ) {
          if      ( z1x2>=(-x1z2) && (y1z2+z1y2)<=0 ) Sign=1.0;
          else if ( z1x2< (-x1z2) && (y1z2+z1y2)> 0 ) Sign=1.0;
      } else if ( absS>=absK && absS>=absO && x1y2<=(-y1x2) ) Sign=1.0;
      Sigma *= Sign;
      //mprintf("        Sigma = %8.2lf Sign=%lf\n",Sigma*RADDEG,Sign);
      Flip[0]=1.0; Flip[1]=-1.0; Flip[2]=-1.0;

    // -------------------- Parallel --------------------
    } else {
      Phi = acos( -0.5 * ( dpX + dpY - dpZ + 1.0 ) );
      //mprintf("PHI= %lf\n",Phi*RADDEG);
      denominator = ( dpX + dpY - dpZ + 3.0 );
      // Kappa (Buckle)
      numerator = ( - dpX + dpY - dpZ + 1.0 );
      Kappa = Phi * sqrt( fabs( numerator / denominator ) );
      absK = fabs(Kappa);
      // Omega (Propeller Twist)
      numerator = (   dpX - dpY - dpZ + 1.0 );
      Omega = Phi * sqrt( fabs( numerator / denominator ) );
      absO = fabs(Omega);
      // Sigma (Opening)
      numerator = (   dpX + dpY + dpZ + 1.0 ); 
      Sigma = Phi * sqrt( fabs( numerator / denominator ) );
      absS = fabs(Sigma);
      // Debug
      //mprintf("DEBUG: K/O/S = %lf %lf %lf\n",Kappa,Omega,Sigma);
      // Kappa Sign determination
      Sign=-1.0;
      if ( absK>=absO && absK>=absS && y1z2>=(-z1y2) ) Sign=1.0;
      else if   ( absO>=absK && absO>=absS ) {
          if      ( z1x2<=(-x1z2) && (x1y2+y1x2)<=0 ) Sign=1.0;
          else if ( z1x2 >(-x1z2) && (x1y2+y1x2)> 0 ) Sign=1.0;
      } else if ( absS>=absK && absS>=absO ) {
          if      ( x1y2<=y1x2 && (x1z2-z1x2)>=0 ) Sign=1.0;
          else if ( x1y2 >y1x2 && (x1z2-z1x2)< 0 ) Sign=1.0;
      }
      Kappa *= Sign;
      //mprintf("        Kappa = %8.2lf Sign=%lf\n",Kappa*RADDEG,Sign);
      // Omega Sign determination
      Sign=-1.0;
      if        ( absK>=absO && absK>=absS ) {
          if      ( y1z2>=(-z1y2) && (x1y2+y1x2)<=0 ) Sign=1.0;
          else if ( y1z2< (-z1y2) && (x1y2+y1x2)> 0 ) Sign=1.0;
      } else if ( absO>=absK && absO>=absS && z1x2<=(-x1z2) ) Sign=1.0;
      else if   ( absS>=absK && absS>=absO ) {
          if      ( x1y2<=y1x2 && (y1z2-z1y2)>=0 ) Sign=1.0;
          else if ( x1y2> y1x2 && (y1z2-z1y2)< 0 ) Sign=1.0;
      }
      Omega *= Sign;
      //mprintf("        Omega = %8.2lf Sign=%lf\n",Omega*RADDEG,Sign);
      // Sigma Sign determination
      Sign=-1.0;
      if        ( absK>=absO && absK>=absS ) {
          if      ( y1z2>=(-z1y2) && (x1z2-z1x2)>=0 ) Sign=1.0;
          else if ( y1z2< (-z1y2) && (x1z2-z1x2)< 0 ) Sign=1.0;
      } else if ( absO>=absK && absO>=absS ) {
          if      ( z1x2<=(-x1z2) && (y1z2-z1y2)>=0 ) Sign=1.0;
          else if ( z1x2> (-x1z2) && (y1z2-z1y2)< 0 ) Sign=1.0;
      } else if ( absS>=absK && absS>=absO && x1y2<=y1x2 ) Sign=1.0;
      Sigma *= Sign;
      //mprintf("        Sigma = %8.2lf Sign=%lf\n",Sigma*RADDEG,Sign);
      Flip[0]=-1.0; Flip[1]=-1.0; Flip[2]=1.0;
    }

    // Shear, stretch, and stagger
    // Half-rotation matrix
    calcRotationMatrix(Rhalf,-Kappa/2.0,-Omega/2.0,-Sigma/2.0);
    // Calc O21 = R2t(o1 - o2)
    matrix_transpose(R2t, Base2->R);
    vector_sub(O21, Base1->Origin(), Base2->Origin());
    matrix_times_vector(O21, R2t, O21);
    // Flip O21
    FLIP_V( O21, Flip);
    // Calc FV2
    FV2[0]=V2[0]; FV2[1]=V2[1]; FV2[2]=V2[2];
    FLIP_V( FV2, Flip);
    // Calc R1V1
    matrix_times_vector(R1V1, Base1->R, V1);
    // Calc R2tR1V1
    matrix_times_vector(Vec, R2t, R1V1);
    // Flip R2tR1V1
    FLIP_V( Vec, Flip);
    // Calc F021 + FR2tR1V1
    vector_sum(Vec, O21, Vec);
    // Calc F021+FR2tR1V1 - FV2
    vector_sub(Vec, Vec, FV2);
    // Calc Rhalf[F021+FR2tR1V1-FV2]
    matrix_times_vector(Vec, Rhalf, Vec);
    // Calc Rhalf[F021+FR2tR1V1-FV2] + FV2
    vector_sum(Vec, Vec, FV2);
    // Calc Rhalf[F021+FR2tR1V1-FV2]+FV2 - V1
    vector_sub(Vec, Vec, V1);
    Shear   = Vec[0]; 
    Stretch = Vec[1]; 
    Stagger = Vec[2];

    // Base pair rotation matrix, R1Rhalf
    matrix_multiply(Rb, Base1->R, Rhalf);

    // Base pair origin
    // Calc V1 + FV2
    vector_sum(Vec, V1, FV2);
    // Calc Rb[V1+FV2]
    matrix_times_vector(Vec, Rb, Vec);
    // Calc R2V2
    matrix_times_vector(O21, Base2->R, V2);
    // Calc R2V2 - Rb[V1+FV2]
    vector_sub(Vec, O21, Vec);
    // Calc R1V1 + R2V2-Rb[V1+FV2]
    vector_sum(Vec, R1V1, Vec);
    // Calc o2 + R1V1+R2V2-Rb[V1+FV2]
    vector_sum(Vec, Base2->Origin(), Vec);
    // Calc o1 + o2+R1V1+R2V2-Rb[V1+FV2]
    vector_sum(Vec, Base1->Origin(), Vec);
    // 0.5 * (o1+o2+R1V1+R2V2-Rb[V1+FV2])
    Vec[0] *= 0.5;
    Vec[1] *= 0.5;
    Vec[2] *= 0.5;

    if (debug>1) {
      mprintf("Buckle= %8.2lf  Propeller= %8.2lf  Opening= %8.2lf\n",
              Kappa*RADDEG, Omega*RADDEG, Sigma*RADDEG);
      mprintf(" Shear= %8.2lf    Stretch= %8.2lf  Stagger= %8.2lf\n",Shear,Stretch,Stagger);
      mprintf("    Ox= %8.2lf         Oy= %8.2lf       Oz= %8.2lf\n",Vec[0],Vec[1],Vec[2]);
      mprintf("    Nx= %8.2lf         Ny= %8.2lf       Nz= %8.2lf\n",Rb[2],Rb[5],Rb[8]);
    }

    // Store data
    Kappa*=RADDEG; Omega*=RADDEG; Sigma*=RADDEG;
    SHEAR.AddData(frameNum, &Shear, nbasepair);
    STRETCH.AddData(frameNum, &Stretch, nbasepair);
    STAGGER.AddData(frameNum, &Stagger, nbasepair);
    BUCKLE.AddData(frameNum, &Kappa, nbasepair);
    PROPELLER.AddData(frameNum, &Omega, nbasepair);
    OPENING.AddData(frameNum, &Sigma, nbasepair);

    // Store BP axes
    BasePairAxes[nbasepair]->X[9 ] = Vec[0];
    BasePairAxes[nbasepair]->X[10] = Vec[1];
    BasePairAxes[nbasepair]->X[11] = Vec[2];
    BasePairAxes[nbasepair]->X[0] = Vec[0] + Rb[0];
    BasePairAxes[nbasepair]->X[1] = Vec[1] + Rb[3];
    BasePairAxes[nbasepair]->X[2] = Vec[2] + Rb[6];
    BasePairAxes[nbasepair]->X[3] = Vec[0] + Rb[1];
    BasePairAxes[nbasepair]->X[4] = Vec[1] + Rb[4];
    BasePairAxes[nbasepair]->X[5] = Vec[2] + Rb[7];
    BasePairAxes[nbasepair]->X[6] = Vec[0] + Rb[2];
    BasePairAxes[nbasepair]->X[7] = Vec[1] + Rb[5];
    BasePairAxes[nbasepair]->X[8] = Vec[2] + Rb[8];
    BasePairAxes[nbasepair]->StoreRotMatrix(Rb);
    // DEBUG - write base pair axes
    //BPaxes->WritePDB(&basepairaxesfile, base1, P->ResidueName(base1), &basepairaxesatom);
    //BasePairAxes.push_back( BPaxes );

    nbasepair++; // Actual base pair count; BP is nbasepair*3
  }
  // DEBUG
  //basepairaxesfile.CloseFile(); 

  return 0;
}

/* NAstruct::determineBasepairParameters() 
 * For each base pair step, determine values of Tilt, Roll, Twist, Shift,
 * Slide, and Rise.
 */
int NAstruct::determineBasepairParameters() {
  double Xi[3], Yi[3], Zi[3], Xj[3], Yj[3], Zj[3];
  double Vi[3], Vj[3];
  int bpi, bpj;
  AxisType *BPi, *BPj;
  double Omegah, dpX, dpY, dpZ;
  double numerator, denominator;
  double Tau, Rho, Omega, Sign, absT, absR, absO;
  double Shift, Slide, Rise;
  double yizj, ziyj, zixj, xizj, xiyj, yixj;
  double Oij[3], Vec[3];
  double Rti[9], Rhalf[9];

  // Default pivot points
  Vi[0]=0.0; Vi[1]=0.0; Vi[2]=0.0;
  Vj[0]=0.0; Vj[1]=0.0; Vj[2]=0.0;

  for (bpi = 0; bpi < ((int)BasePairAxes.size()) - 1; bpi++) {
    bpj = bpi+1;
    //mprintf("\n*** Determining BP parameters for step %i to %i\n",bpi+1,bpj+1);
    BPi = BasePairAxes[bpi];
    BPj = BasePairAxes[bpj];
    // Fill X Y and Z unit vectors
    BPi->RX(Xi);
    BPi->RY(Yi);
    BPi->RZ(Zi);
    BPj->RX(Xj);
    BPj->RY(Yj);
    BPj->RZ(Zj);
    // Calculate dot products 
    dpX = dot_product(Xi,Xj);
    dpY = dot_product(Yi,Yj);
    dpZ = dot_product(Zi,Zj);
    // Dot Products for sign determination
    yizj = dot_product(Yi, Zj);
    ziyj = dot_product(Zi, Yj);
    zixj = dot_product(Zi, Xj);
    xizj = dot_product(Xi, Zj);
    xiyj = dot_product(Xi, Yj);
    yixj = dot_product(Yi, Xj);
    // Calculate Omegah
    Omegah = acos( 0.5 * (dpX + dpY + dpZ - 1.0) );
    // Common denominator
    denominator = ( dpX + dpY + dpZ - 3.0 );
    // Tau (Tilt)
    numerator = ( -dpX + dpY + dpZ - 1.0 );
    Tau   = Omegah * sqrt( fabs( numerator / denominator ) );
    absT = fabs(Tau);
    // Rho (Roll)
    numerator = (  dpX - dpY + dpZ - 1.0 );
    Rho   = Omegah * sqrt( fabs( numerator / denominator ) );
    absR = fabs(Rho);
    // Omega (Twist)
    numerator = (  dpX + dpY - dpZ - 1.0 );
    Omega = Omegah * sqrt( fabs( numerator / denominator ) );
    absO = fabs(Omega);
    // Tau sign determination
    Sign=-1.0;
    if        ( absT>=absR && absT>=absO && ziyj>=yizj ) Sign=1.0;
    else if   ( absR>=absT && absR>=absO ) {
        if      ( xizj>=zixj && (yixj+xiyj)>=0 ) Sign=1.0;
        else if ( xizj< zixj && (yixj+xiyj)< 0 ) Sign=1.0;
    } else if ( absO>=absT && absO>=absR ) {
        if      ( yixj>=xiyj && (zixj+xizj)>=0 ) Sign=1.0;
        else if ( yixj< xiyj && (zixj+xizj)< 0 ) Sign=1.0;
    }
    Tau *= Sign;
    // Rho sign determination
    Sign=-1.0;
    if        ( absT>=absR && absT>=absO ) {
        if      ( ziyj>=yizj && (yixj+xiyj)>=0 ) Sign=1.0;
        else if ( ziyj< yizj && (yixj+xiyj)< 0 ) Sign=1.0;
    } else if ( absR>=absT && absR>=absO && xizj>=zixj ) Sign=1.0;
    else if   ( absO>=absT && absO>=absR ) {
        if      ( yixj>=xiyj && (ziyj+yizj)>=0 ) Sign=1.0;
        else if ( yixj< xiyj && (ziyj+yizj)< 0 ) Sign=1.0;
    }
    Rho *= Sign;
    // Omega sign determination
    Sign=-1.0;
    if        ( absT>=absR && absT>=absO ) {
        if      ( ziyj>=yizj && (zixj+xizj)>=0 ) Sign=1.0;
        else if ( ziyj< yizj && (zixj+xizj)< 0 ) Sign=1.0;
    } else if ( absR>=absT && absR>=absO ) {
        if      ( xizj>=zixj && (ziyj+yizj)>=0 ) Sign=1.0;
        else if ( xizj< zixj && (ziyj+yizj)< 0 ) Sign=1.0;
    } else if ( absO>=absT && absO>=absR && yixj>=xiyj ) Sign=1.0;
    Omega *= Sign;

    // Calculate half rotation matrix
    calcRotationMatrix(Rhalf, -Tau/2.0, -Rho/2.0, -Omega/2.0);
    // Calc Oij = Rti(oj - oi)
    matrix_transpose(Rti, BPi->R);
    vector_sub(Oij, BPj->Origin(), BPi->Origin());
    matrix_times_vector(Oij, Rti, Oij);
    // Calc RtiRjvj
    matrix_times_vector(Vec, BPj->R, Vj);
    matrix_times_vector(Vec, Rti, Vec);
    // Calc Oij + RtiRjvj - vi
    vector_sub(Vec, Vec, Vi);
    vector_sum(Vec, Oij, Vec);
    // Calc Rhalf[Oij + RtiRjvj - vi]
    matrix_times_vector(Vec, Rhalf, Vec);
    // Calc Rhalf[Oij + RtiRjvj - vi] + vi - vj
    vector_sum(Vec, Vec, Vi);
    vector_sub(Vec, Vec, Vj);
    Shift = Vec[0];
    Slide = Vec[1];
    Rise  = Vec[2];

    if (debug>1) {
      mprintf(" Tilt= %8.2lf   Roll= %8.2lf  Twist= %8.2lf\n", Tau*RADDEG, Rho*RADDEG, Omega*RADDEG);
      mprintf("Shift= %8.2lf  Slide= %8.2lf   Rise= %8.2lf\n", Shift, Slide, Rise);
    }

    // Store data
    Tau *= RADDEG;
    Rho *= RADDEG;
    Omega *= RADDEG;
    SHIFT.AddData(frameNum, &Shift, bpi);
    SLIDE.AddData(frameNum, &Slide, bpi);
    RISE.AddData(frameNum, &Rise, bpi);
    TILT.AddData(frameNum, &Tau, bpi);
    ROLL.AddData(frameNum, &Rho, bpi);
    TWIST.AddData(frameNum, &Omega, bpi);
  }

  return 0;
}
// ----------------------------------------------------------------------------

/* NAstruct::init()
 * Expected call: nastruct [resrange <range>] [out <filename>] [naout <nafilename>]
 *                         [noheader]
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int NAstruct::init() {
  // Get keywords
  outFilename = actionArgs.getKeyString("out",NULL);
  naoutFilename = actionArgs.getKeyString("naout",NULL);
  resRange.SetRange( actionArgs.getKeyString("resrange",NULL) );
  noheader = actionArgs.hasKey("noheader");
  // Get Masks
  // Dataset
  // Add dataset to data file list

  mprintf("    NAstruct: ");
  if (resRange.Empty())
    mprintf("Scanning all NA residues");
  else
    mprintf("Scanning residues %s",resRange.RangeArg());
  if (naoutFilename!=NULL) {
      mprintf(", formatted output to file %s",naoutFilename);
    if (noheader) mprintf(", no header");
  }
  mprintf(".\n");

  return 0;
}

/* NAstruct::setup()
 * Set up for this parmtop. Get masks etc.
 * P is set in Action::Setup
 * Determine the number of NA bases that will be analyzed, along with the masks
 * that correspond to the reference frame atoms.
 */
int NAstruct::setup() {
  int res, refAtom, atom;
  int residue;
  AxisType *axis; 
  AtomMask *Mask;
  Range actualRange;

  // Clear all lists
  ClearLists();

  // If range is empty (i.e. no resrange arg given) look for all NA residues.
  if (resRange.Empty()) {
    for (res=0; res < P->nres; res++) {
      if ( ID_base(P->ResidueName(res))!=UNKNOWN_BASE )
        actualRange.AddToRange(res);
    }

  // Otherwise, for each residue in resRange check if it is a NA
  } else {
    resRange.Begin();
    while (resRange.NextInRange(&res)) {
      // User residues numbers start from 1
      res--;
      if (ID_base(P->ResidueName(res))!=UNKNOWN_BASE) 
        actualRange.AddToRange(res);
    }
  }
  // Exit if no NA residues specified
  if (actualRange.Empty()) {
    mprintf("Error: NAstruct::setup: No NA residues found for %s\n",P->parmName);
    return 1;
  }

  // DEBUG - print all residues
  if (debug>0)
    actualRange.PrintRange("    NAstruct: NA res:",1);

  // Set up frame to hold reference coords for each NA residue
  actualRange.Begin();
  while (actualRange.NextInRange(&residue)) {
    axis = new AxisType();
    if ( axis->SetRefCoord( P->ResidueName(residue) ) ) {
      mprintf("Error: NAstruct::setup: Could not get ref coords for %i:%s\n",
              residue+1, P->ResidueName(residue));
      return 1;
    }
    RefCoords.push_back( axis );

    // Set up a mask for this NA residue in this parm. The mask will contain
    // only those atoms which are defined in the reference coords.
    Mask = new AtomMask();
    for (refAtom=0; refAtom < axis->natom; refAtom++) {
      res = -1; // Target atom
      //mprintf("      Ref atom: [%s]\n",axis->Name[refAtom]);
      for (atom=P->resnums[residue]; atom < P->resnums[residue+1]; atom++) {
        //mprintf("        Scanning %i [%s]\n", atom, P->names[atom]);
        if ( axis->AtomNameIs(refAtom, P->names[atom]) ) {
          res = atom;  
          break;
        }
      }
      if (res==-1) {
        mprintf("Error:: NAstruct::setup: Ref atom [%s] not found in residue %i:%s\n",
                 axis->AtomName(refAtom), residue+1, P->ResidueName(residue));
        return 1;
      }
      Mask->AddAtom(res);
    } // End Loop over reference atoms
    if (Mask->None()) {
      mprintf("Error:: NAstruct::setup: No atoms found for residue %i:%s\n",
              residue+1, P->ResidueName(residue));
      delete Mask;
      return 1;
    }
    ExpMasks.push_back( Mask );
    if (debug>1) {
      mprintf("      NAstruct: Res %i:%s mask atoms: ",residue+1,P->ResidueName(residue));
      Mask->PrintMaskAtoms();
      mprintf("\n");
    }

    // Set up empty frame to hold input coords for this residue
    axis = new AxisType(Mask->Nselected);
    ExpFrames.push_back( axis );

    // Set up initial axes for this NA residue.
    axis = new AxisType(); 
    axis->SetPrincipalAxes();
    BaseAxes.push_back( axis );
  } // End Loop over NA residues

  Nbases = (int)RefCoords.size(); // Also BaseAxes, ExpFrames, and ExpMasks size.
  mprintf("    NAstruct: Set up %i bases.\n",Nbases);

  return 0;  
}

/* NAstruct::action()
 */
int NAstruct::action() {

  // Set up base axes
  if ( setupBaseAxes(F) ) return 1;

  // Determine Base Pairing
  if ( determineBasePairing() ) return 1;

  // Determine base parameters
  determineBaseParameters();

  // Determine base pair parameters
  determineBasepairParameters();

  Nframe++;

  // Get base pair axes
  //setupBasePairAxes();

  return 0;
} 

/* NAstruct::print()
 */
void NAstruct::print() {
  CpptrajFile outfile;
  CharBuffer buffer;
  int frame, nbasepair;
  char tempName[64]; // NOTE: Replce with string?
  size_t dataFileSize;
  DataSet *na_dataset = NULL;

  // Set precision of all datasets
/*  SHEAR.SetPrecisionOfDatasets(8,3);
  STRETCH.SetPrecisionOfDatasets(8,3);
  STAGGER.SetPrecisionOfDatasets(8,3);
  BUCKLE.SetPrecisionOfDatasets(8,3);
  PROPELLER.SetPrecisionOfDatasets(8,3);
  OPENING.SetPrecisionOfDatasets(8,3);
  SHIFT.SetPrecisionOfDatasets(8,3);
  SLIDE.SetPrecisionOfDatasets(8,3);
  RISE.SetPrecisionOfDatasets(8,3);
  TILT.SetPrecisionOfDatasets(8,3);
  ROLL.SetPrecisionOfDatasets(8,3);
  TWIST.SetPrecisionOfDatasets(8,3);*/

  if (naoutFilename==NULL) return;

  // ---------- Base pair parameters ----------
  sprintf(tempName,"BP.%s",naoutFilename);
  if ( outfile.SetupFile(tempName, WRITE, UNKNOWN_FORMAT, UNKNOWN_TYPE, debug) ) {
    mprinterr("Error: NAstruct::print(): Could not set up naout file %s\n",tempName);
    return;
  }
  if ( outfile.OpenFile() ) return;
  // Calculate dataFileSize
  // First 2 cols 8, rest are 12 because all double datasets being used.
  // 1) Header: [8*2] + [12*6] + 8
  if (!noheader)
    dataFileSize = 96;
  else 
    dataFileSize = 0;
  // 2) Data: 
  //   ((header size * nbasepair)+1) * nframes
  dataFileSize += (((96 * SHEAR.Size())+1) * Nframe);
  // Allocate buffer
  buffer.Allocate( dataFileSize );
  // Write File header
  if (!noheader)
    buffer.Sprintf("%-8s %8s %12s %12s %12s %12s %12s %12s\n","#Frame","BasePair",
                   "Shear","Stretch","Stagger","Buckle","Propeller","Opening");
  // Write Base pair data for each frame
  for (frame=0; frame < Nframe; frame++) {
    // Base-pair parameters
    for (nbasepair=0; nbasepair < SHEAR.Size(); nbasepair++) {
      // Frame and base pair #
      buffer.Sprintf("%8i %8i",frame,nbasepair);
      na_dataset = SHEAR.GetDataSetN(nbasepair);
      na_dataset->WriteBuffer(buffer,frame);
      na_dataset = STRETCH.GetDataSetN(nbasepair);
      na_dataset->WriteBuffer(buffer,frame);
      na_dataset = STAGGER.GetDataSetN(nbasepair);
      na_dataset->WriteBuffer(buffer,frame);
      na_dataset = BUCKLE.GetDataSetN(nbasepair);
      na_dataset->WriteBuffer(buffer,frame);
      na_dataset = PROPELLER.GetDataSetN(nbasepair);
      na_dataset->WriteBuffer(buffer,frame);
      na_dataset = OPENING.GetDataSetN(nbasepair);
      na_dataset->WriteBuffer(buffer,frame);
      buffer.NewLine();
    }
    buffer.NewLine(); 
  }
  outfile.IO->Write(buffer.Buffer(), sizeof(char), buffer.CurrentSize());
  outfile.CloseFile();
  
  // ---------- Base pair step parameters ----------
  sprintf(tempName,"BPstep.%s",naoutFilename);
  if ( outfile.SetupFile(tempName, WRITE, UNKNOWN_FORMAT, UNKNOWN_TYPE, debug) ) {
    mprinterr("Error: NAstruct::print(): Could not set up naout file %s\n",tempName);
    return;
  }
  if ( outfile.OpenFile() ) return;
  // Calculate dataFileSize
  // First 2 cols 8, rest are 12 because all double datasets being used.
  // 1) Header: [8*2] + [12*6] + 8
  if (!noheader)
    dataFileSize = 96;
  else
    dataFileSize = 0;
  // 2) Data: 
  //   ((header size * nbpstep)+1) * nframes
  dataFileSize += (((96 * SHIFT.Size())+1) * Nframe);
  // Allocate buffer
  buffer.Allocate( dataFileSize );
  // Write File header
  if (!noheader)
    buffer.Sprintf("%-8s %8s %12s %12s %12s %12s %12s %12s\n","#Frame","BPstep",
                   "Shift","Slide","Rise","Tilt","Roll","Twist");
  // Write base pair step data for each frame
  for (frame=0; frame < Nframe; frame++) {
    // Base-pair step parameters
    for (nbasepair=0; nbasepair < SHIFT.Size(); nbasepair++) {
      // Frame and base pair #
      buffer.Sprintf("%8i %8i",frame,nbasepair);
      na_dataset = SHIFT.GetDataSetN(nbasepair);
      na_dataset->WriteBuffer(buffer,frame);
      na_dataset = SLIDE.GetDataSetN(nbasepair);
      na_dataset->WriteBuffer(buffer,frame);
      na_dataset = RISE.GetDataSetN(nbasepair);
      na_dataset->WriteBuffer(buffer,frame);
      na_dataset = TILT.GetDataSetN(nbasepair);
      na_dataset->WriteBuffer(buffer,frame);
      na_dataset = ROLL.GetDataSetN(nbasepair);
      na_dataset->WriteBuffer(buffer,frame);
      na_dataset = TWIST.GetDataSetN(nbasepair);
      na_dataset->WriteBuffer(buffer,frame);
      buffer.NewLine();
    }
    buffer.NewLine();
  }
  outfile.IO->Write(buffer.Buffer(), sizeof(char), buffer.CurrentSize());
  outfile.CloseFile();
}

