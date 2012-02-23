#include <cmath>
#include <cstdio> // sprintf
#include "Action_NAstruct.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "Constants.h" // RADDEG
#include "vectormath.h"

// CONSTRUCTOR
NAstruct::NAstruct() {
  //fprintf(stderr,"NAstruct Con\n");
  Nbp=0;
  Nbases=0;
  HBdistCut2=9.61;  // Hydrogen Bond distance cutoff^2: 3.1^2
  // NOTE: Currently not used
  HBangleCut2=2.53; // Hydrogen Bond angle cutoff (in radians, ~145 degs)
  // NOTE: Is this too big?
  originCut2=6.25;  // Origin cutoff^2 for base-pairing: 2.5^2
  Nframe=0;
  outFilename=NULL;
  naoutFilename=NULL;
  noheader = false;
# ifdef NASTRUCTDEBUG
  calcparam = true;
# endif
} 

// DESTRUCTOR
NAstruct::~NAstruct() { 
  ClearLists();
  // NOTE: Since BasePairAxes are set up to correspond with SHEAR etc dont
  // free in this routine - should only be freed at the very end.
  BasePairAxes.clear();
}

// NAstruct::ClearLists()
/** Clear all parm-dependent lists
  */
void NAstruct::ClearLists() {
  RefCoords.clear();
  BaseAxes.clear();
  //ExpFrames.clear();
  ExpMasks.clear();
  FitMasks.clear();
}

// ------------------------- PRIVATE FUNCTIONS --------------------------------
// NAstruct::setupBaseAxes()
/** For each residue defined in reference coords, get the corresponding input
  * coords and fit the reference coords (and reference axes) on top of input 
  * coords. This sets up the reference axes for each base.
  */
int NAstruct::setupBaseAxes(Frame *InputFrame) {
  double rmsd, RotMatrix[9], TransVec[6];
  AxisType refFrame; // Hold copy of base reference coords
  AxisType expFrame; // Hold copy of input base coords
# ifdef NASTRUCTDEBUG
  AxisPDBwriter baseaxesfile;
  baseaxesfile.Open("baseaxes.pdb");
  AxisPDBwriter basesfile;
  basesfile.Open("bases.pdb");
  mprintf("\n=================== Setup Base Axes ===================\n");
# endif

  // For each axis in RefCoords, use corresponding mask in ExpMasks to set 
  // up an axis for ExpCoords.
  for (int base=0; base < Nbases; base++) {
    // Set exp coords based on previously set-up mask
    BaseAxes[base].SetFrameCoordsFromMask( InputFrame->X, &ExpMasks[base] );
#   ifdef NASTRUCTDEBUG
    int expbasenum = BaseAxes[base].BaseNum();
    mprintf("Base REF %i:%4s   EXP %i:%4s\n",
            RefCoords[base].BaseNum()+1,RefCoords[base].ResName(),
            expbasenum+1,currentParm->ResidueName(expbasenum));
    ExpMasks[base].PrintMaskAtoms("ExpMask");
    FitMasks[base].PrintMaskAtoms("FitMask");
    mprintf("#  %4s %8s %8s %8s   %4s %8s %8s %8s\n","Ref","Rx","Ry","Rz","Exp","Ex","Ey","Ez");
    for (int i = 0; i < RefCoords[base].natom; i++) {
      int j = i * 3;
      mprintf("%-2i %4s %8.3lf %8.3lf %8.3lf   %4s %8.3lf %8.3lf %8.3lf\n",i,
              RefCoords[base].AtomName(i),
              RefCoords[base].X[j],RefCoords[base].X[j+1],RefCoords[base].X[j+2],
              BaseAxes[base].AtomName(i),
              BaseAxes[base].X[j],BaseAxes[base].X[j+1],BaseAxes[base].X[j+2]);
    }
#   endif 
    /* Now that we have a set of reference coords and the corresponding input
     * coords, RMS fit the reference coords to the input coords to obtain the
     * appropriate rotation and translations that will put the reference coords 
     * on top of input (experimental) coords. Per 3DNA procedure, not all 
     * reference atoms are used in the RMS fit; only ring atoms are used. 
     */
    refFrame.SetAxisFromMask( RefCoords[base], FitMasks[base] );
    expFrame.SetAxisFromMask( BaseAxes[base], FitMasks[base] );
    rmsd = refFrame.RMSD( &expFrame, RotMatrix, TransVec, false);
    /* RotMatrix and TransVec now contain rotation and translation
     * that will orient refcoord to expframe. The first translation is that of
     * the reference frame to the absolute origin, the second translation is
     * that of the reference frame to the exp. coords after rotation.
     * The rotation matrix contains the coordinates of the X, Y, and Z unit 
     * vectors of the base axes.
     */
    // Store the Rotation matrix and the rotated and translated origin.
    double Vec[3];
    Vec[0]=(TransVec[0]*RotMatrix[0]+TransVec[1]*RotMatrix[1]+TransVec[2]*RotMatrix[2])+TransVec[3];
    Vec[1]=(TransVec[0]*RotMatrix[3]+TransVec[1]*RotMatrix[4]+TransVec[2]*RotMatrix[5])+TransVec[4];
    Vec[2]=(TransVec[0]*RotMatrix[6]+TransVec[1]*RotMatrix[7]+TransVec[2]*RotMatrix[8])+TransVec[5];
    BaseAxes[base].StoreRotMatrix( RotMatrix, Vec );
    if (debug>0) { 
      mprintf("Base %i: RMS of RefCoords from ExpCoords is %lf\n",base+1,rmsd);
      //printMatrix_3x3("Rotation matrix:",RotMatrix);
      //printRotTransInfo(RotMatrix,TransVec);
      BaseAxes[base].PrintAxisInfo("BaseAxes");
    }
#   ifdef NASTRUCTDEBUG
    // DEBUG - Write base axis to file
    baseaxesfile.WriteAxes(BaseAxes[base], base, RefCoords[base].ResName());
    // Overlap ref coords onto input coords. 
    refFrame.SetFromFrame( &RefCoords[base] );
    refFrame.Trans_Rot_Trans(TransVec,RotMatrix);
    basesfile.Write(refFrame, base, RefCoords[base].ResName());
#   endif
  } // END loop over bases

  return 0;
}

// NAstruct::GCpair()
/** Look for 3 HB based on heavy atom distances:
  * 1. G:O6 -- C:N4  6 -- 6
  * 2. G:N1 -- C:N3  7 -- 4
  * 3. G:N2 -- C:O2  9 -- 3
  * Atom positions are known in standard Ref. Multiply by 3 to get into X.
  */
bool NAstruct::GCpair(AxisType *DG, AxisType *DC) {
  int Nhbonds = 0;
  double dist2;
  for (int hb = 0; hb < 3; hb++) {
    dist2 = DIST2_NoImage(DG->HbondCoord[hb], DC->HbondCoord[hb]);
    if ( dist2 < HBdistCut2 ) {
      ++Nhbonds;
#     ifdef NASTRUCTDEBUG
      int dg_hbatom = DG->HbondAtom[hb];
      int dc_hbatom = DC->HbondAtom[hb];
      mprintf("\t\t%s:%s -- %s:%s = %lf\n",
              DG->BaseName(),DG->AtomName(dg_hbatom),
              DC->BaseName(),DC->AtomName(dc_hbatom),sqrt(dist2));
#     endif
    }
  }
  if (Nhbonds>0) return true;
  return false;
}

// NAstruct::ATpair()
/** Look for 2 HB based on heavy atom distances
  * 1. A:N6 -- T:O4  6 -- 6
  * 2. A:N1 -- T:N3  7 -- 4
  */
bool NAstruct::ATpair(AxisType *DA, AxisType *DT) {
  int Nhbonds = 0;
  double dist2;
  for (int hb = 0; hb < 2; hb++) {
    dist2 = DIST2_NoImage(DA->HbondCoord[hb], DT->HbondCoord[hb]);
    if ( dist2 < HBdistCut2 ) {
      ++Nhbonds;
#     ifdef NASTRUCTDEBUG
      int da_hbatom = DA->HbondAtom[hb];
      int dt_hbatom = DT->HbondAtom[hb];
      mprintf("\t\t%s:%s -- %s:%s = %lf\n",
              DA->BaseName(),DA->AtomName(da_hbatom),
              DT->BaseName(),DT->AtomName(dt_hbatom),sqrt(dist2));
#     endif
    }
  }
  if (Nhbonds>0) return true;
  return false;
}

// NAstruct::basesArePaired()
/** Given two base axes for which IDs have been given and reference coords set,
  * determine whether the bases are paired via hydrogen bonding criteria.
  * NOTE: Currently only set up for Antiparallel WC detection
  */
bool NAstruct::basesArePaired(AxisType *base1, AxisType *base2) {
  // G C
  if      ( base1->ID==AxisType::GUA && base2->ID==AxisType::CYT ) return GCpair(base1,base2);
  else if ( base1->ID==AxisType::CYT && base2->ID==AxisType::GUA ) return GCpair(base2,base1);
  // A T
  else if ( base1->ID==AxisType::ADE && base2->ID==AxisType::THY ) return ATpair(base1,base2);
  else if ( base1->ID==AxisType::THY && base2->ID==AxisType::ADE ) return ATpair(base2,base1);
  // A U
  else if ( base1->ID==AxisType::ADE && base2->ID==AxisType::URA ) return ATpair(base1,base2);
  else if ( base1->ID==AxisType::URA && base2->ID==AxisType::ADE ) return ATpair(base2,base1);
//  else {
//    mprintf("Warning: NAstruct: Unrecognized pair: %s - %s\n",NAbaseName[base1->ID],
//             NAbaseName[base2->ID]);
//  }
  return false;
}

// NAstruct::determineBasePairing()
/** Determine which bases are paired from the individual base axes.
  */
int NAstruct::determineBasePairing() {
  double distance;
  std::vector<bool> isPaired( BaseAxes.size(), false);
  int base1,base2;
  double Z1[3], Z2[3];
  bool AntiParallel = false;
  double minDistance;
  int minBaseNum;

  Nbp = 0;
  BasePair.clear();
# ifdef NASTRUCTDEBUG  
  mprintf("\n=================== Setup Base Pairing ===================\n");
# endif

  /* For each unpaired base, determine if it is paired with another base
   * determined by the distance between their axis origins.
   */
  for (base1=0; base1 < Nbases-1; base1++) {
    if (isPaired[base1]) continue;
    minBaseNum = -1;
    minDistance = 0;
    for (base2=base1+1; base2 < Nbases; base2++) {
      if (isPaired[base2]) continue;
      // First determine which origin axes coords are close enough to 
      // consider pairing.
      distance = DIST2_NoImage(BaseAxes[base1].Origin(), BaseAxes[base2].Origin());
      if (distance < originCut2) {
#       ifdef NASTRUCTDEBUG
        mprintf("  Axes distance for %s -- %s is %lf\n",
                BaseAxes[base1].BaseName(),BaseAxes[base2].BaseName(),sqrt(distance));
#       endif
        if (minBaseNum==-1) {
          minDistance = distance;
          minBaseNum = base2;
        } else {
          if (distance < minDistance) {
            minDistance = distance;
            minBaseNum = base2;
          }
        }
      }
    }
/*
        // Figure out angle between y vectors, determines whether bases
        // are able to hydrogen bond based on angle cutoff.
        BaseAxes[base1].RY(Z1);
        BaseAxes[base2].RY(Z2);
        distance = dot_product_angle(Z1, Z2);
#       ifdef NASTRUCTDEBUG
        mprintf("      Angle between Y vectors is %lf deg. (%lf)\n",distance * RADDEG,distance);
#       endif
        if (distance > HBangleCut2) {*/
    if (minBaseNum!=-1) {
          base2 = minBaseNum;
#         ifdef NASTRUCTDEBUG
          mprintf("    Closest base is %i, %lf Ang.\n",base2+1,sqrt(minDistance));
#         endif
          // Figure out if z vectors point in same (<90 deg) or opposite (>90 deg) direction
          BaseAxes[base1].RZ(Z1);
          BaseAxes[base2].RZ(Z2);
          //printVector("Base1Z",Z1); printVector("Base2Z",Z2);
          distance = dot_product_angle(Z1, Z2);
          //mprintf("    Dot product of Z vectors: %lf\n",distance);
          if (distance > (PIOVER2)) { // If theta(Z) > 90 deg.
#           ifdef NASTRUCTDEBUG
            mprintf("      %s is anti-parallel to %s\n",BaseAxes[base1].BaseName(),
                    BaseAxes[base2].BaseName());
#           endif
            AntiParallel = true;
          } else {
#           ifdef NASTRUCTDEBUG
            mprintf("      %s is parallel to %s\n",BaseAxes[base1].BaseName(),
                    BaseAxes[base2].BaseName());
#           endif
            AntiParallel = false;
          }
          if (basesArePaired(&BaseAxes[base1], &BaseAxes[base2])) {
            BasePair.push_back(base1);
            BasePair.push_back(base2);
            if (AntiParallel) 
              BasePair.push_back(1);
            else
              BasePair.push_back(0);
            isPaired[base1]=true;
            isPaired[base2]=true;
            ++Nbp;
          }
//        } // END if distance > HBangleCut2
//      } // END if distance < originCut2
    } // END if minBaseNum!=-1
  } // END Loop over base1

  if (debug>0) mprintf("    NAstruct: Set up %i base pairs.\n",Nbp);
  base2=1;
  //mprintf("DEBUG: BasePair.size = %i\n",(int)BasePair.size());
  //mprintf("DEBUG: SHEAR.size = %i\n",(int)SHEAR.size());
  //mprintf("DEBUG: BasePairAxes.size = %i\n",(int)BasePairAxes.size());
  AxisType bpaxis;
  //bpaxis.SetPrincipalAxes();
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
      BasePairAxes.push_back( bpaxis );
    } 
    // Print base pair info
    if (debug>1) {
      int bp_1 = BasePair[base1  ];
      int bp_2 = BasePair[base1+1];
      mprintf("        BP %i: Res %s to %s",base2,
              RefCoords[bp_1].BaseName(), RefCoords[bp_2].BaseName());
      if ( BasePair[base1+2] )
        mprintf(" AntiParallel.\n");
      else
        mprintf(" Parallel.\n");
    }
    ++base2;
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

// NAstruct::calculateParameters()
/** Given two axes, calculate translational and rotational parameters
  * between them.
  */
int NAstruct::calculateParameters(AxisType &BaseAxis1, AxisType &BaseAxis2, 
                                  AxisType *BPaxis, double *Param) 
{
  double hingeAxis[3],X1[3],Y1[3],Z1[3],X2[3],Y2[3],Z2[3],O1[3],O2[3];
  double R[9], Rinv[9], RotatedR2[9], RotatedR1[9];
  double r2, OM[3], O21[3], Vec[3];
  AxisType Axis1, Axis2;
# ifdef NASTRUCTDEBUG
  AxisType tempAxis;
  AxisPDBwriter paramfile;
  if (calcparam)
    paramfile.Open("Param.pdb");
# endif

  Axis1 = BaseAxis1;
  Axis2 = BaseAxis2;

  Axis1.RX(X1);
  Axis1.RY(Y1);
  Axis1.RZ(Z1);
  Axis1.OXYZ(O1);

  Axis2.RX(X2);
  Axis2.RY(Y2);
  Axis2.RZ(Z2);
  Axis2.OXYZ(O2);
# ifdef NASTRUCTDEBUG
  printVector("O1",O1);
  printMatrix_3x3("R1",Axis1.R);
  printVector("O2",O2);
  printMatrix_3x3("R2",Axis2.R);
# endif

  // Hinge axis is cross product between Z1 and Z2
  cross_product(hingeAxis, Z1, Z2);
# ifdef NASTRUCTDEBUG
  printVector("hinge",hingeAxis);
# endif
  // Normalize hinge axis
  vector_norm( hingeAxis, &r2 );
  //double rolltilt = asin( vector_norm( hingeAxis, &r2 ));
# ifdef NASTRUCTDEBUG
  printVector("norm(hinge)",hingeAxis);
# endif

  // Angle between Z1 and Z2, Z1 dot Z2 
  double rolltilt = dot_product_angle(Z1, Z2);
# ifdef NASTRUCTDEBUG
  mprintf("\tAngle between Z1 and Z2= %lf\n",rolltilt*RADDEG);
# endif

  // Calculate forward and backwards half rolltilt rotation around
  // hinge axis.
  calcRotationMatrix(R, hingeAxis, -0.5*rolltilt);
  //calcRotationMatrix(R, hingeAxis, 0.5*rolltilt);
  matrix_transpose(Rinv, R);
# ifdef NASTRUCTDEBUG
  printMatrix_3x3("Rhalf",R);
# endif

  // Rotate R2 by -0.5 * rolltilt degrees around the hinge
  matrix_multiply_3x3(RotatedR2, R,    Axis2.R);
  // Rotate R1 by 0.5 * rolltilt degrees around the hinge (inverse rotation)
  matrix_multiply_3x3(RotatedR1, Rinv, Axis1.R);

# ifdef NASTRUCTDEBUG
  // Print rotated R1 and R2
  printMatrix_3x3("Rotated R1",RotatedR1);
  printMatrix_3x3("Rotated R2",RotatedR2);
  if (calcparam) {
    tempAxis.StoreRotMatrix(RotatedR1,O1);
    paramfile.WriteAxes(tempAxis,0,(char*)"R1'");
    tempAxis.StoreRotMatrix(RotatedR2,O2);
    paramfile.WriteAxes(tempAxis,1,(char*)"R2'");
  }
# endif

  // Average R1 and R2 to get the middle frame
  for (int i = 0; i < 9; i++)
    R[i] = (RotatedR1[i] + RotatedR2[i]) / 2;
  // Normalize X, Y and Z vectors
  r2 = sqrt( R[0]*R[0] + R[3]*R[3] + R[6]*R[6] );
  R[0] /= r2;
  R[3] /= r2;
  R[6] /= r2;
  r2 = sqrt( R[1]*R[1] + R[4]*R[4] + R[7]*R[7] );
  R[1] /= r2;
  R[4] /= r2;
  R[7] /= r2;
  r2 = sqrt( R[2]*R[2] + R[5]*R[5] + R[8]*R[8] );
  R[2] /= r2;
  R[5] /= r2;
  R[8] /= r2;

  // Take average of origins
  OM[0] = (O1[0] + O2[0]) / 2;
  OM[1] = (O1[1] + O2[1]) / 2;
  OM[2] = (O1[2] + O2[2]) / 2;

# ifdef NASTRUCTDEBUG
  printVector("Origin Mean",OM);
  // Print Rm and hinge axis
  printMatrix_3x3("Rm",R);
  if (calcparam) {
    // Use Rinv to store hinge axis in Z
    for (int i=0; i<9; i++) Rinv[i]=0;
    Rinv[2]=hingeAxis[0]; Rinv[5]=hingeAxis[1]; Rinv[8]=hingeAxis[2];
    tempAxis.StoreRotMatrix(Rinv, OM);
    paramfile.WriteAxes(tempAxis,2,(char*)"Hng");
    // Store middle frame
    tempAxis.StoreRotMatrix(R,OM);
    paramfile.WriteAxes(tempAxis,3,(char*)"Rm");
  }
# endif

  // If BPaxis is not NULL, store Rm and OM as BP axis.
  if (BPaxis!=NULL) 
    BPaxis->StoreRotMatrix(R, OM);

  // Shift Slide Rise / Shear Stretch Stagger
  vector_sub(O21, O2, O1);
  // Since this is really vector times matrix, use matrix transpose times vec
  matrixT_times_vector(Vec, R, O21);
# ifdef NASTRUCTDEBUG
  printVector("O21",O21);
  printVector("Vec",Vec);
# endif
  Param[0] = Vec[0];
  Param[1] = Vec[1];
  Param[2] = Vec[2];

  // Set Z1 to Z from middle frame
  Z1[0] = R[2];
  Z1[1] = R[5];
  Z1[2] = R[8];

  // Twist / Opening
  // Angle between rotated Y1 and rotated Y2
  Y1[0] = RotatedR1[0];
  Y1[1] = RotatedR1[3];
  Y1[2] = RotatedR1[6];
  Y2[0] = RotatedR2[0];
  Y2[1] = RotatedR2[3];
  Y2[2] = RotatedR2[6];
  double twistopen = dot_product_angle(Y1, Y2);
  // Sign of twistopen related to (Y1'xY2') dot Z of middle frame
  cross_product(Vec,Y1,Y2);
  double sign = dot_product(Vec,Z1);
# ifdef NASTRUCTDEBUG
  mprintf("\tAngle between Y1' and Y2' is %lf\n",twistopen*RADDEG);
  mprintf("\tDot product of Y1'xY2' with Zm is %lf radians\n",sign);
# endif
  if (sign<0) 
    sign=-1.0; 
  else 
    sign=1.0;
  twistopen*=sign;
# ifdef NASTRUCTDEBUG
  mprintf("\tFinal Twist/Opening is %10.4lf\n",twistopen*RADDEG);
# endif
  Param[3] = twistopen;

  // Phase angle
  // Angle between hinge axis and middle frame Y axis
  Y1[0] = R[1];
  Y1[1] = R[4];
  Y1[2] = R[7];
  double phi = dot_product_angle(hingeAxis, Y1);
  // Sign of phi related to (hingeAxis x Ym) dot Z of middle frame
  cross_product(Vec,hingeAxis,Y1);
  sign = dot_product(Vec,Z1);
# ifdef NASTRUCTDEBUG
  mprintf("\tAngle between hinge axis and Ym is %lf\n",phi*RADDEG);
  mprintf("\tDot product of (hingeAxis x Ym) with Zm is %lf radians\n",sign);
# endif
  if (sign<0) 
    sign=-1.0;
  else 
    sign=1.0;
  phi *= sign;
  double sinphi = sin( phi );
  double cosphi = cos( phi );
# ifdef NASTRUCTDEBUG
  mprintf("\tPhase angle is %lf, sinphi is %lf, cosphi is %lf\n",phi*RADDEG,sinphi,cosphi);
# endif

  // Roll / Propeller
  double rollprop = rolltilt * cosphi;
  Param[4] = rollprop;

  // Tilt / Buckle
  double tiltbuck = rolltilt * sinphi;
  Param[5] = tiltbuck;

# ifdef NASTRUCTDEBUG
  mprintf("\tRoll/Propeller %10.4lf\n",rollprop*RADDEG);
  mprintf("\tTilt/Buckle %10.4lf\n",tiltbuck*RADDEG);
  if (calcparam) calcparam=false;
# endif
  return 0;
}

// NAstruct::determineBaseParameters()
/** For each base in a base pair, get the values of buckle, propeller twist,
  * opening, shear, stretch, and stagger. Also determine the origin and 
  * rotation matrix for each base pair reference frame.
  */
int NAstruct::determineBaseParameters() {
  double Param[6];
# ifdef NASTRUCTDEBUG
  AxisPDBwriter basepairaxesfile;
  basepairaxesfile.Open("basepairaxes.pdb");
  mprintf("\n=================== Determine BP Parameters ===================\n");
# endif

  int nbasepair=0;
  for (unsigned int BP=0; BP < BasePair.size(); BP+=3) {
    int base1 = BasePair[BP  ];
    int base2 = BasePair[BP+1];
#   ifdef NASTRUCTDEBUG
    mprintf("BasePair %s to %s",BaseAxes[base1].BaseName(),BaseAxes[base2].BaseName());
    if (BasePair[BP+2])
      mprintf(" Anti-parallel.\n");
    else
      mprintf(" Parallel.\n");
#   endif
    // Check Antiparallel / Parallel
    // Flip YZ (rotate around X) for antiparallel
    // Flip XY (rotate around Z) for parallel
    if (BasePair[BP+2])
      BaseAxes[base2].FlipYZ();
    else
      BaseAxes[base2].FlipXY();
    // Calc BP parameters, set up basepair axes
    //calculateParameters(BaseAxes[base1],BaseAxes[base2],&BasePairAxes[nbasepair],Param);
    calculateParameters(BaseAxes[base2],BaseAxes[base1],&BasePairAxes[nbasepair],Param);
    // Store data
    Param[3] *= RADDEG;
    Param[4] *= RADDEG;
    Param[5] *= RADDEG;
    SHEAR.AddData(frameNum, Param, nbasepair);
    STRETCH.AddData(frameNum, Param+1, nbasepair);
    STAGGER.AddData(frameNum, Param+2, nbasepair);
    OPENING.AddData(frameNum, Param+3, nbasepair);
    PROPELLER.AddData(frameNum, Param+4, nbasepair);
    BUCKLE.AddData(frameNum, Param+5, nbasepair);

#   ifdef NASTRUCTDEBUG
    // DEBUG - write base pair axes
    basepairaxesfile.WriteAxes(BasePairAxes[nbasepair], base1, RefCoords[base1].ResName());
#   endif

    ++nbasepair; // Actual base pair count; BP is nbasepair*3
  }

  return 0;
}

// NAstruct::determineBasepairParameters() 
/** For each base pair step, determine values of Tilt, Roll, Twist, Shift,
  * Slide, and Rise.
  */
int NAstruct::determineBasepairParameters() {
  double Param[6];
# ifdef NASTRUCTDEBUG
  mprintf("\n=================== Determine BPstep Parameters ===================\n");
# endif
  for (int bpi = 0; bpi < (int)BasePairAxes.size() - 1; bpi++) {
    int bpj = bpi+1;
#   ifdef NASTRUCTDEBUG
    mprintf("BasePair step %i to %i\n",bpi+1,bpj+1);
#   endif
    // Calc step parameters
    calculateParameters(BasePairAxes[bpi], BasePairAxes[bpj], NULL, Param);
    // Store data
    Param[3] *= RADDEG;
    Param[4] *= RADDEG;
    Param[5] *= RADDEG;
    SHIFT.AddData(frameNum, Param, bpi);
    SLIDE.AddData(frameNum, Param+1, bpi);
    RISE.AddData(frameNum, Param+2, bpi);
    TWIST.AddData(frameNum, Param+3, bpi);
    ROLL.AddData(frameNum, Param+4, bpi);
    TILT.AddData(frameNum, Param+5, bpi);
  }

  return 0;
}
// ----------------------------------------------------------------------------

// NAstruct::init()
/** Expected call: nastruct [resrange <range>] [out <filename>] [naout <nafilename>]
  *                         [noheader]
  */
// Dataset name will be the last arg checked for. Check order is:
//    1) Keywords
//    2) Masks
//    3) Dataset name
int NAstruct::init() {
  char *resrange_arg, *maparg;
  ArgList maplist;
  AxisType::NAbaseType mapbase;
  std::string resname;
  NAME mapresname;
  mapresname[4]='\0';
  // Get keywords
  outFilename = actionArgs.getKeyString("out",NULL);
  naoutFilename = actionArgs.getKeyString("naout",NULL);
  resrange_arg = actionArgs.getKeyString("resrange",NULL);
  if (resrange_arg != NULL)
    if (resRange.SetRange( resrange_arg )) return 1;
  noheader = actionArgs.hasKey("noheader");
  // Get custom residue maps
  while ( (maparg = actionArgs.getKeyString("resmap",NULL))!=NULL ) {
    // Split maparg at ':'
    maplist.SetList(maparg,":");
    // Expect only 2 args
    if (maplist.Nargs()!=2) {
      mprinterr("Error: nastruct: resmap format should be '<ResName>:{A,C,G,T,U}' (%s)\n",maparg);
      return 1;
    }
    // Check that second arg is A,C,G,T,or U
    if (maplist.ArgIs(1,"A")) mapbase=AxisType::ADE;
    else if (maplist.ArgIs(1,"C")) mapbase=AxisType::CYT;
    else if (maplist.ArgIs(1,"G")) mapbase=AxisType::GUA;
    else if (maplist.ArgIs(1,"T")) mapbase=AxisType::THY;
    else if (maplist.ArgIs(1,"U")) mapbase=AxisType::URA;
    else {
      mprinterr("Error: nastruct: resmap format should be '<ResName>:{A,C,G,T,U}' (%s)\n",maparg);
      return 1;
    }
    // Check that residue name is <= 4 chars
    resname = maplist[0]; 
    if (resname.size() > 4) {
      mprinterr("Error: nastruct: resmap resname > 4 chars (%s)\n",maparg);
      return 1;
    }
    // Format residue name
    int i = 0;
    for (std::string::iterator it = resname.begin(); it != resname.end(); it++)
      mapresname[i++]=*it;
    mapresname[i]='\0';
    PadWithSpaces(mapresname);
    resname.assign(mapresname);
    mprintf("\tCustom Map: [%s]\n",resname.c_str());
    //maplist.PrintList();
    // Add to CustomMap
    customRes = CustomMap.find(resname);
    if (customRes!=CustomMap.end()) {
      mprintf("Warning: nastruct: resmap: %s already mapped.\n",resname.c_str());
    } else {
      CustomMap.insert( std::pair<std::string,AxisType::NAbaseType>(resname,mapbase) );
    }
  }
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

// NAstruct::setup()
/** Determine the number of NA bases that will be analyzed, along with 
  * the masks that correspond to the reference frame atoms.
  */
int NAstruct::setup() {
  int resnum;
  AxisType axis; 
  AtomMask Mask;
  AtomMask fitMask;
  Range actualRange;
  AxisType::RefReturn refreturn;
  std::string resname; // Used to check for custom mapped residues
  AxisType::NAbaseType customBaseType;

  // Clear all lists
  ClearLists();

  // If range is empty (i.e. no resrange arg given) look through all 
  // solute residues.
  if (resRange.Empty()) 
    actualRange.SetRange(0, currentParm->FinalSoluteRes());
  // If user range specified, create new range shifted by -1 since internal
  // resnums start from 0.
  else {
    actualRange = resRange;
    actualRange.ShiftBy(-1);
  }

  // Exit if no residues specified
  if (actualRange.Empty()) {
    mprinterr("Error: NAstruct::setup: No residues specified for %s\n",currentParm->parmName);
    return 1;
  }

  // DEBUG - print all residues
  //if (debug>0)
  //  actualRange.PrintRange("    NAstruct: NA res:",1);

  // Set up frame to hold reference coords for each NA residue
  actualRange.Begin();
  while (actualRange.NextInRange(&resnum)) {
    customBaseType = AxisType::UNKNOWN_BASE;
    // Check if the residue at resnum matches any of the custom maps
    if (!CustomMap.empty()) {
      resname.assign( currentParm->ResidueName(resnum) );
      customRes = CustomMap.find( resname );
      if (customRes!=CustomMap.end()) {
        mprintf("\tCustom map found: %i [%s]\n",resnum+1,(*customRes).first.c_str());
        customBaseType = (*customRes).second;
      }
    }
    // Set up ref coords in correct order, along with corresponding 
    // parm mask for this residue. SetRefCoord should overwrite all
    // previously set up information in axis and Mask.
    refreturn = axis.SetRefCoord( currentParm, resnum, Mask, fitMask, customBaseType );
    // If not recognized as a NA residue, continue to next.
    // Print a warning if the user specified this range.
    if ( refreturn == AxisType::NA_UNKNOWN ) {
      if (!resRange.Empty()) {
        mprintf("Warning: Residue %i:%s not recognized as NA residue.\n",
                resnum+1, currentParm->ResidueName(resnum));
      }
      continue;
    } else if ( refreturn == AxisType::NA_ERROR  ) {
      mprinterr("Error: NAstruct::setup: Could not get ref coords for %i:%s\n",
                resnum+1, currentParm->ResidueName(resnum));
      return 1;
    }
    if (Mask.None()) {
      mprintf("Error:: NAstruct::setup: No atoms found for residue %i:%s\n",
              resnum+1, currentParm->ResidueName(resnum));
      return 1;
    }
    if (fitMask.None()) {
      mprintf("Error:: NAstruct::setup: No fit atoms found for residue %i:%s\n",
              resnum+1, currentParm->ResidueName(resnum));
      return 1;
    }
    RefCoords.push_back( axis );
    ExpMasks.push_back( Mask );
    FitMasks.push_back( fitMask );
    if (debug>1) {
      mprintf("\tNAstruct: Res %i:%s ",resnum+1,currentParm->ResidueName(resnum));
      Mask.PrintMaskAtoms("NAmask");
      mprintf("\t          Ref %i:%s ",resnum+1,axis.BaseName());
      axis.PrintAtomNames();
    }

    // Set up empty frame to hold input coords for this residue.
    // Initially set to be a copy of the reference frame. The coordinates
    // will be overwritten later.
    BaseAxes.push_back( axis );
  } // End Loop over NA residues

  Nbases = (int)RefCoords.size(); // Also BaseAxes, ExpFrames, and ExpMasks size.
  mprintf("    NAstruct: Set up %i bases.\n",Nbases);

  return 0;  
}

// NAstruct::action()
int NAstruct::action() {

  // Set up base axes
  if ( setupBaseAxes(currentFrame) ) return 1;

  // Determine Base Pairing
  if ( determineBasePairing() ) return 1;

  // Determine base parameters
  determineBaseParameters();

  // Determine base pair parameters
  determineBasepairParameters();

  ++Nframe;

  return 0;
} 

// NAstruct::print()
void NAstruct::print() {
  CpptrajFile outfile;
  CharBuffer buffer;
  int frame, nbasepair;
  char tempName[64]; // NOTE: Replce with string?
  size_t dataFileSize;
  DataSet *na_dataset = NULL;

  // Set precision of all datasets
/*  int dsw = 12;
  int dsp = 2;
  SHEAR.SetPrecisionOfDatasets(dsw,dsp);
  STRETCH.SetPrecisionOfDatasets(dsw,dsp);
  STAGGER.SetPrecisionOfDatasets(dsw,dsp);
  BUCKLE.SetPrecisionOfDatasets(dsw,dsp);
  PROPELLER.SetPrecisionOfDatasets(dsw,dsp);
  OPENING.SetPrecisionOfDatasets(dsw,dsp);
  SHIFT.SetPrecisionOfDatasets(dsw,dsp);
  SLIDE.SetPrecisionOfDatasets(dsw,dsp);
  RISE.SetPrecisionOfDatasets(dsw,dsp);
  TILT.SetPrecisionOfDatasets(dsw,dsp);
  ROLL.SetPrecisionOfDatasets(dsw,dsp);
  TWIST.SetPrecisionOfDatasets(dsw,dsp);*/

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
      buffer.Sprintf("%8i %8i",frame+OUTPUTFRAMESHIFT,nbasepair+1);
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
      buffer.Sprintf("%8i %8i",frame+OUTPUTFRAMESHIFT,nbasepair+1);
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

