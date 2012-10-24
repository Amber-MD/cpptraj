#include <cmath>
#include "Action_NAstruct.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "Constants.h" // RADDEG
#include "vectormath.h"

// CONSTRUCTOR
Action_NAstruct::Action_NAstruct() :
  HBdistCut2_(12.25),  // Hydrogen Bond distance cutoff^2: 3.5^2
  //HBangleCut2_(2.53),  // Hydrogen Bond angle cutoff (in radians, ~145 degs)
  // NOTE: Is this too big?
  originCut2_(6.25),   // Origin cutoff^2 for base-pairing: 2.5^2
  maxResSize_(0),
  debug_(0),
  printheader_(true),
  useReference_(false),
  masterDSL_(0)
# ifdef NASTRUCTDEBUG
  ,calcparam(true)
# endif
{ 
  //fprintf(stderr,"NAstruct Con\n");
}

void Action_NAstruct::Help() {
  mprintf("nastruct [resrange <range>] [naout <nafilename>]\n");
  mprintf("         [noheader] [resmap <ResName>:{A,C,G,T,U} ...]\n");
  mprintf("         [hbcut <hbcut>] [origincut <origincut>]\n");
  mprintf("         [ reference | refindex <#> | ref <REF> ]\n");
}

// DESTRUCTOR
Action_NAstruct::~Action_NAstruct() { 
  ClearLists();
  // NOTE: Since BasePairAxes are set up to correspond with SHEAR etc they 
  //       are only freed at the very end.
  BasePairAxes.clear();
}

// Output Format Strings
static const char BP_OUTPUT_FMT[66] = "%8i %8i %8i %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %2i\n";
static const char NA_OUTPUT_FMT[73] = "%8i %4i-%-4i %4i-%-4i %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf\n";

// ------------------------- PRIVATE FUNCTIONS --------------------------------
// Action_NAstruct::ClearLists()
/** Clear all parm-dependent lists */
void Action_NAstruct::ClearLists() {
  RefCoords.clear();
  BaseAxes.clear();
  ExpMasks.clear();
  FitMasks.clear();
}

// Action_NAstruct::setupBaseAxes()
/** For each residue defined in reference coords, get the corresponding input
  * coords and fit the reference coords (and reference axes) on top of input 
  * coords. This sets up the reference axes for each base.
  */
int Action_NAstruct::setupBaseAxes(Frame *InputFrame) {
  double rmsd, RotMatrix[9], TransVec[6];
  Frame refFrame(maxResSize_); // Hold copy of base reference coords for RMS fit
  Frame inpFrame(maxResSize_); // Hold copy of input base coords for RMS fit
# ifdef NASTRUCTDEBUG
  AxisPDBwriter baseaxesfile;
  baseaxesfile.Open("baseaxes.pdb");
  AxisPDBwriter basesfile;
  basesfile.Open("bases.pdb");
  mprintf("\n=================== Setup Base Axes ===================\n");
# endif

  // For each axis in RefCoords, use corresponding mask in ExpMasks to set 
  // up an axis for ExpCoords.
  int Nbases = (int)RefCoords.size();
  for (int base=0; base < Nbases; base++) {
    // Set exp coords based on previously set-up mask
    Frame expCoords( *InputFrame, ExpMasks[base] );
    BaseAxes[base].SetCoordsFromFrame( expCoords );
    // If P atom defined store phosphorus coords
    if (BaseAxes[base].HasPatom()) 
      BaseAxes[base].SetPcrd( InputFrame->XYZ( BaseAxes[base].Pidx() ) );
    // If O4' atom defined store phosphorus coords
    if (BaseAxes[base].HasO4atom())
      BaseAxes[base].SetO4crd( InputFrame->XYZ( BaseAxes[base].O4idx() ) );
#   ifdef NASTRUCTDEBUG
    int expbasenum = BaseAxes[base].ResNum();
    mprintf("Base REF %i:%4s   EXP %i:%4s\n",
            RefCoords[base].ResNum()+1,RefCoords[base].ResName(),
            expbasenum+1,currentParm->ResidueName(expbasenum));
    ExpMasks[base].PrintMaskAtoms("ExpMask");
    FitMasks[base].PrintMaskAtoms("FitMask");
    mprintf("#  %4s %8s %8s %8s   %4s %8s %8s %8s\n","Ref","Rx","Ry","Rz","Exp","Ex","Ey","Ez");
    for (int i = 0; i < RefCoords[base].Natom(); i++) {
      int j = i * 3;
      mprintf("%-2i %4s %8.3lf %8.3lf %8.3lf",i,
              RefCoords[base].AtomName(i),
              RefCoords[base][j],RefCoords[base][j+1],RefCoords[base][j+2]);
      mprintf("   %4s %8.3lf %8.3lf %8.3lf\n",
              BaseAxes[base].AtomName(i),
              BaseAxes[base][j],BaseAxes[base][j+1],BaseAxes[base][j+2]);
    }
#   endif 
    /* Now that we have a set of reference coords and the corresponding input
     * coords, RMS fit the reference coords to the input coords to obtain the
     * appropriate rotation and translations that will put the reference coords 
     * on top of input (experimental) coords. Per 3DNA procedure, not all 
     * reference atoms are used in the RMS fit; only ring atoms are used. 
     */
    refFrame.SetCoordinatesByMask( RefCoords[base].xAddress(), FitMasks[base] );
    inpFrame.SetCoordinatesByMask( BaseAxes[base].xAddress(), FitMasks[base] );
    rmsd = refFrame.RMSD( inpFrame, RotMatrix, TransVec, false);
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
    if (debug_>0) { 
      mprintf("Base %i: RMS of RefCoords from ExpCoords is %lf\n",base+1,rmsd);
      //printMatrix_3x3("Rotation matrix:",RotMatrix);
      //printRotTransInfo(RotMatrix,TransVec);
      BaseAxes[base].PrintAxisInfo("BaseAxes");
    }
#   ifdef NASTRUCTDEBUG
    // DEBUG - Write base axis to file
    baseaxesfile.WriteAxes(BaseAxes[base], base, BaseAxes[base].ResName());
    // Overlap ref coords onto input coords.
    Frame reftemp(RefCoords[base].Natom(), RefCoords[base].xAddress() ); 
    reftemp.Trans_Rot_Trans(TransVec,RotMatrix);
    basesfile.Write(RefCoords[base], reftemp.xAddress(), base, RefCoords[base].ResName());
#   endif
  } // END loop over bases

  return 0;
}

// Action_NAstruct::GCpair()
/** Look for 3 HB based on heavy atom distances:
  * 1. G:O6 -- C:N4  6 -- 6
  * 2. G:N1 -- C:N3  7 -- 4
  * 3. G:N2 -- C:O2  9 -- 3
  */
int Action_NAstruct::GCpair(AxisType *DG, AxisType *DC) {
  int Nhbonds = 0;
  double dist2;
  for (int hb = 0; hb < 3; hb++) {
    dist2 = DIST2_NoImage(DG->HbondCoord[hb], DC->HbondCoord[hb]);
    if ( dist2 < HBdistCut2_ ) {
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
  return Nhbonds;
}

// Action_NAstruct::ATpair()
/** Look for 2 HB based on heavy atom distances
  * 1. A:N6 -- T:O4  6 -- 6
  * 2. A:N1 -- T:N3  7 -- 4
  */
int Action_NAstruct::ATpair(AxisType *DA, AxisType *DT) {
  int Nhbonds = 0;
  double dist2;
  for (int hb = 0; hb < 2; hb++) {
    dist2 = DIST2_NoImage(DA->HbondCoord[hb], DT->HbondCoord[hb]);
    if ( dist2 < HBdistCut2_ ) {
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
  return Nhbonds;
}

// Action_NAstruct::basesArePaired()
/** Given two base axes for which IDs have been given and reference coords set,
  * determine whether the bases are paired via hydrogen bonding criteria.
  * NOTE: Currently only set up for WC detection
  */
int Action_NAstruct::basesArePaired(AxisType *base1, AxisType *base2) {
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
  return 0;
}

// Action_NAstruct::determineBasePairing()
/** Determine which bases are paired from the individual base axes. Also 
  * sets up BP and BP step parameter DataSets.
  */
int Action_NAstruct::determineBasePairing() {
  double distance;
  std::vector<bool> isPaired( BaseAxes.size(), false);
  int base1,base2;
  double Z1[3], Z2[3];
  bool AntiParallel = false;
  double minDistance;
  int minBaseNum;

  BasePair.clear();
  NumberOfHbonds_.clear();
# ifdef NASTRUCTDEBUG  
  mprintf("\n=================== Setup Base Pairing ===================\n");
# endif

  /* For each unpaired base, find the closest potential pairing base 
   * determined by the distance between their axis origins.
   */
  int Nbases = (int)RefCoords.size();
  int Nbases1 = Nbases - 1;
  for (base1=0; base1 < Nbases1; base1++) {
    if (isPaired[base1]) continue;
    minBaseNum = -1;
    minDistance = 0;
    for (base2=base1+1; base2 < Nbases; base2++) {
      if (isPaired[base2]) continue;
      distance = DIST2_NoImage(BaseAxes[base1].Origin(), BaseAxes[base2].Origin());
      if (distance < originCut2_) {
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
        if (distance > HBangleCut2_) {*/
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
          int NHB = basesArePaired(&BaseAxes[base1], &BaseAxes[base2]);
          if (NHB > 0) {
            BasePair.push_back(base1);
            BasePair.push_back(base2);
            if (AntiParallel) 
              BasePair.push_back(1);
            else
              BasePair.push_back(0);
            isPaired[base1]=true;
            isPaired[base2]=true;
          }
          NumberOfHbonds_.push_back( NHB );
//        } // END if distance > HBangleCut2
//      } // END if distance < originCut2
    } // END if minBaseNum!=-1
  } // END Loop over base1

  unsigned int Nbp = BasePair.size() / 3;
  if (debug_>0) mprintf("    NAstruct: Set up %i base pairs.\n",Nbp);
  // Resize BasePairAxes
  BasePairAxes.resize( Nbp );
  // Print Base Pair info
  if (debug_>1) {
    base2=1; //  Base pair #
    for (base1 = 0; base1 < (int)BasePair.size(); base1+= 3) {
      int bp_1 = BasePair[base1  ];
      int bp_2 = BasePair[base1+1];
      mprintf("        BP %i: Res %i:%s to %i:%s",base2++,
              BaseAxes[bp_1].ResNum()+1, BaseAxes[bp_1].ResName(), 
              BaseAxes[bp_2].ResNum()+1, BaseAxes[bp_2].ResName());
      if ( BasePair[base1+2] )
        mprintf(" AntiParallel.\n");
      else
        mprintf(" Parallel.\n");
    }
  }
  // For each BP, set up a dataset for each structural parameter if
  // one is not already set up.
  base2 = 1; // Base pair # and DataSet idx
  for (base1 = 0; base1 < (int)BasePair.size(); base1+= 3) {
      if ( base2 > (int)SHEAR_.size() ) {
        // Create legend
        int bp_1 = BasePair[base1  ];
        int bp_2 = BasePair[base1+1];
        std::string bpname = integerToString( BaseAxes[bp_1].ResNum()+1 ) +
                             BaseAxes[bp_1].ResName() +
                             integerToString( BaseAxes[bp_2].ResNum()+1 ) + 
                             BaseAxes[bp_2].ResName();
        // Create sets
        SHEAR_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"shear",bpname) );
        STRETCH_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"stretch",bpname));
        STAGGER_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"stagger",bpname));
        BUCKLE_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"buckle",bpname) );
        PROPELLER_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"prop",bpname) );
        OPENING_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"open",bpname) );
        BPHBONDS_.push_back( masterDSL_->AddSetIdxAspect(DataSet::INT,dataname_,base2,"hb",bpname) );
        MAJOR_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"major",bpname) );
        MINOR_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"minor",bpname) );
      }
    ++base2;
  }
  // For each BP step, set up a dataset for each structural parameter 
  // if one is not already set up. One less than total # BP.
  base2 = 1; // Base pair step # and DataSet idx
  int NBPstep = (int)(BasePairAxes.size() - 1);
  for (base1 = 0; base1 < NBPstep; ++base1) {
    if ( base2 > (int)SHIFT_.size() ) {
      // Create legend
      int bpidx = base1 * 3;
      int bp_1 = BasePair[bpidx  ];
      int bp_2 = BasePair[bpidx+1];
      int bp_3 = BasePair[bpidx+3];
      int bp_4 = BasePair[bpidx+4];
      std::string sname = integerToString( BaseAxes[bp_1].ResNum()+1 ) +
                          BaseAxes[bp_1].ResName()[0] +
                          integerToString( BaseAxes[bp_2].ResNum()+1 ) +
                          BaseAxes[bp_2].ResName()[0] + "-" +
                          integerToString( BaseAxes[bp_3].ResNum()+1 ) +
                          BaseAxes[bp_3].ResName()[0] +
                          integerToString( BaseAxes[bp_4].ResNum()+1 ) +
                          BaseAxes[bp_4].ResName()[0];
      // Create Sets
      SHIFT_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"shift",sname) );
      SLIDE_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"slide",sname) );
      RISE_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"rise",sname) );
      TILT_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"tilt",sname) );
      ROLL_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"roll",sname) );
      TWIST_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"twist",sname) );
      XDISP_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"xdisp",sname) );
      YDISP_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"ydisp",sname) );
      HRISE_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"hrise",sname) );
      INCL_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"incl",sname) );
      TIP_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"tip",sname) );
      HTWIST_.push_back( masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,base2,"htwist",sname) );
    }
    ++base2;
  }
  
/*
  //mprintf("DEBUG: BasePair.size = %i\n",(int)BasePair.size());
  //mprintf("DEBUG: SHEAR.size = %i\n",(int)SHEAR.size());
  //mprintf("DEBUG: BasePairAxes.size = %i\n",(int)BasePairAxes.size());
*/

  return 0;
}

// AverageMatrices()
static void AverageMatrices(double *R, double *RotatedR1, double *RotatedR2) {
  double r2;
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
}

// Action_NAstruct::calculateParameters()
/** Given two axes, calculate translational and rotational parameters
  * between them.
  */
int Action_NAstruct::calculateParameters(AxisType &Axis1, AxisType &Axis2, 
                                  AxisType *BPaxis, double *Param) 
{
  double hingeAxis[3],Y1[3],Z1[3],Y2[3],Z2[3],O1[3],O2[3];
  double R[9], Rinv[9], RotatedR2[9], RotatedR1[9];
  double r2, OM[3], O21[3], Vec[3];
# ifdef NASTRUCTDEBUG
  AxisType tempAxis;
  AxisPDBwriter paramfile;
  if (calcparam)
    paramfile.Open("Param.pdb");
  printVector("O1",Axis1.Origin());
  printMatrix_3x3("R1",Axis1.R);
  printVector("O2",Axis2.Origin());
  printMatrix_3x3("R2",Axis2.R);
# endif

  // Hinge axis is cross product between Z1 and Z2
  Axis1.RZ(Z1);
  Axis2.RZ(Z2);
  cross_product(hingeAxis, Z1, Z2);
# ifdef NASTRUCTDEBUG
  printVector("hinge",hingeAxis);
# endif
  // Normalize hinge axis
  vector_norm( hingeAxis, &r2 );
# ifdef NASTRUCTDEBUG
  printVector("norm(hinge)",hingeAxis);
# endif

  // Roll/Tilt is Angle between Z1 and Z2
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
    tempAxis.StoreRotMatrix(RotatedR1,Axis1.Origin());
    paramfile.WriteAxes(tempAxis,0,(char*)"R1'");
    tempAxis.StoreRotMatrix(RotatedR2,Axis2.Origin());
    paramfile.WriteAxes(tempAxis,1,(char*)"R2'");
  }
# endif

  // Average R1 and R2 to get the middle frame
  AverageMatrices(R, RotatedR1, RotatedR2);

  // Take average of origins
  Axis1.OXYZ(O1);
  Axis2.OXYZ(O2);
  vector_sum(OM, O1, O2);
  OM[0] /= 2;
  OM[1] /= 2;
  OM[2] /= 2;
  //OM[0] = (O1[0] + O2[0]) / 2;
  //OM[1] = (O1[1] + O2[1]) / 2;
  //OM[2] = (O1[2] + O2[2]) / 2;

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
  if (BPaxis!=NULL) { 
    BPaxis->StoreRotMatrix(R, OM);
    BPaxis->StoreBPresnums(Axis2.ResNum(), Axis1.ResNum());
  }

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
  // Sign of twistopen related to (Y1'xY2') dot Z of middle frame
  Y1[0] = RotatedR1[1];
  Y1[1] = RotatedR1[4];
  Y1[2] = RotatedR1[7];
  Y2[0] = RotatedR2[1];
  Y2[1] = RotatedR2[4];
  Y2[2] = RotatedR2[7];
  double twistopen = dot_product_sign(Y1, Y2, Z1);
# ifdef NASTRUCTDEBUG
  mprintf("\tFinal Twist/Opening is %10.4lf\n",twistopen*RADDEG);
# endif
  Param[3] = twistopen;

  // Phase angle
  // Angle between hinge axis and middle frame Y axis
  // Sign of phi related to (hingeAxis x Ym) dot Z of middle frame
  Y1[0] = R[1];
  Y1[1] = R[4];
  Y1[2] = R[7];
  double phi = dot_product_sign(hingeAxis,Y1,Z1);
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

// Action_NAstruct::helicalParameters()
int Action_NAstruct::helicalParameters(AxisType &Axis1, AxisType &Axis2, double *Param) {
  double X1[3],X2[3],Y1[3],Y2[3],Z1[3],Z2[3],O1[3],O2[3],helicalAxis[3];
  double hingeAxis[3], R[9], RotatedR1[9], RotatedR2[9], Vec[3], r2;
  // NOTE: Just use Vec for hingeAxis?
  // X2 - X1
  Axis1.RX(X1);
  Axis2.RX(X2);
  vector_sub(O1, X2, X1);
  // Y2 - Y1
  Axis1.RY(Y1);
  Axis2.RY(Y2);
  vector_sub(O2, Y2, Y1);
  // Local helical axis: (X2-X1) x (Y2-Y1)
  cross_product( helicalAxis, O1, O2 );
# ifdef NASTRUCTDEBUG
  printVector("X2 - X1",O1);
  printVector("Y2 - Y1",O2);
  printVector("(X2-X1) x (Y2-Y1)",helicalAxis);
# endif
  normalize( helicalAxis );
# ifdef NASTRUCTDEBUG
  printVector("NORM[(X2-X1)x(Y2-Y1)]",helicalAxis);
# endif

  // Tip/inclination is angle between helical axis and z1
  Axis1.RZ(Z1);
  double tipinc = dot_product_angle(helicalAxis, Z1);
  // Hinge axis is normalized cross product of helical axis to z1
  cross_product(hingeAxis, helicalAxis, Z1);
  normalize( hingeAxis );
  // Rotate R1 around hinge axis by -tipinc
  calcRotationMatrix(R, hingeAxis, -tipinc);
  matrix_multiply_3x3(RotatedR1, R, Axis1.R);
# ifdef NASTRUCTDEBUG
  mprintf("\tTip/Inclination: %lf\n",tipinc*RADDEG);
  printVector("Hinge axis 1",hingeAxis);
  printMatrix_3x3("Rotated R1", RotatedR1);
# endif

  // Tip/inclination should be same for z2
  Axis2.RZ(Z2);
  //mprintf("\tTipCheck= %lf\n",dot_product_angle(helicalAxis, Z2)*RADDEG);
  // Hinge axis is normalized cross product from h to z2
  cross_product(Vec, helicalAxis, Z2);
  normalize( Vec );
  // Rotate R2 around hinge axis by -tipinc
  calcRotationMatrix(R, Vec, -tipinc); 
  matrix_multiply_3x3(RotatedR2, R, Axis2.R);
# ifdef NASTRUCTDEBUG
  printVector("Hinge axis 2",Vec);
  printMatrix_3x3("Rotated R2",RotatedR2);
# endif

  // Average Rotated R1 and R2 to get middle helical frame
  AverageMatrices(R, RotatedR1, RotatedR2);

  // Helical twist is angle from Rotated Y1 to Rotated Y2
  // Sign is given by (Y1'xY2' dot helicalAxis)
  Y1[0] = RotatedR1[1];
  Y1[1] = RotatedR1[4];
  Y1[2] = RotatedR1[7];
  Y2[0] = RotatedR2[1];
  Y2[1] = RotatedR2[4];
  Y2[2] = RotatedR2[7];
  double Twist = dot_product_sign(Y1,Y2,helicalAxis);
  Param[5] = Twist;

  // Calc O2 - O1
  vector_sub(Vec, Axis2.Origin(), Axis1.Origin());
  // Project (O2-O1) onto helical axis
  double Rise = dot_product(Vec,helicalAxis);
  Param[2] = Rise;
# ifdef NASTRUCTDEBUG
  printMatrix_3x3("Hm",R);
  mprintf("\tTwist is %lf\n",Twist*RADDEG);
  mprintf("\tRise is %lf\n",Rise);
# endif

  // Phase angle is angle from hinge Axis 1 to RotatedR1 Y
  // Sign is given by (hingeAxis x Y1') dot helicalAxis
  double phase = dot_product_sign(hingeAxis, Y1, helicalAxis);

  // Tip is tipinc * cos( phase )
  double Tip = tipinc * cos( phase );
  Param[4] = Tip;
  // Inclination is tipinc * sin( phase )
  double Inc = tipinc * sin( phase );
  Param[3] = Inc;
# ifdef NASTRUCTDEBUG
  mprintf("\tPhase angle is %lf\n",phase*RADDEG);
  mprintf("\tTip is %lf\n",Tip*RADDEG);
  mprintf("\tInclination is %lf\n",Inc*RADDEG);
# endif

  // Calc vector AB (store in X1)
  // Vec contains O2-O1
  Z1[0] = helicalAxis[0] * Rise; 
  Z1[1] = helicalAxis[1] * Rise; 
  Z1[2] = helicalAxis[2] * Rise;
  vector_sub(X1, Vec, Z1);

  // Calc vector AD (store in X2)
  double AD_angle = PIOVER2 - (0.5 * Twist);
  // rotation of AD_angle around helicalAxis
  // NOTE: Assuming we dont need RotatedR2 anymore
  calcRotationMatrix(RotatedR2, helicalAxis, AD_angle);
  // rotate AB
  matrix_times_vector(X2,RotatedR2,X1);
  normalize( X2 );
# ifdef NASTRUCTDEBUG
  printVector("AB",X1);
  mprintf("\tAD_angle is %lf\n",AD_angle*RADDEG);
  printVector("AD",X2);
# endif

  // Calc magnitude of AD; 0.5 * |AB| / sin( 0.5 * Twist )
  double AB_mag = vector_norm( X1, &r2 );
  double AD_mag = (0.5 * AB_mag) / sin( 0.5 * Twist );

  // Calc origin of local helical frame for BP 1
  // Origin1 + (AD_mag * AD)
  X2[0] *= AD_mag;
  X2[1] *= AD_mag;
  X2[2] *= AD_mag;
  vector_sum(O1, Axis1.Origin(), X2);

  // Calc origin of local helical frame for BP 2
  // O1 + (Rise * helicalAxis)
  // Z1 contains helicalAxis * Rise
  vector_sum(O2, O1, Z1);

  // Calculate origin of middle helical frame
  vector_sum(Z2, O2, O1);
  Z2[0] /= 2; 
  Z2[1] /= 2; 
  Z2[2] /= 2;
# ifdef NASTRUCTDEBUG
  mprintf("\t|AD| = %lf\n",AD_mag);
  printVector("o1_h",O1);
  printVector("o2_h",O2);
  printVector("Om_h",Z2);
# endif

  // Calc vector from O1 to Origin1
  vector_sub(Vec, Axis1.Origin(), O1);

  // X-disp is projection of vector from O1 to Origin1 onto 
  // X axis of RotatedR1.
  X1[0] = RotatedR1[0];
  X1[1] = RotatedR1[3];
  X1[2] = RotatedR1[6];
  double X_disp = dot_product(Vec, X1);
  Param[0] = X_disp;

  // Y-disp is projection of vector from O1 to Origin1 onto 
  // Y axis of RotatedR1.
  Y1[0] = RotatedR1[1];
  Y1[1] = RotatedR1[4];
  Y1[2] = RotatedR1[7];
  double Y_disp = dot_product(Vec, Y1);
  Param[1] = Y_disp;
# ifdef NASTRUCTDEBUG
  mprintf("\tX-displacement= %lf\n",X_disp);
  mprintf("\tY-displacement= %lf\n",Y_disp);
# endif

  return 0;
}

// Action_NAstruct::determineBaseParameters()
/** For each base in a base pair, get the values of buckle, propeller twist,
  * opening, shear, stretch, and stagger. Also determine the origin and 
  * rotation matrix for each base pair reference frame.
  */
int Action_NAstruct::determineBaseParameters(int frameNum) {
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
    // TEST - calc P--P distance
    float dPtoP = 0.0;
    //mprintf("\tDEBUG: %i %i:", BaseAxes[base1].Pidx(), BaseAxes[base2].Pidx() );
    if ( BaseAxes[base1].HasPatom() && BaseAxes[base2].HasPatom() ) {
      double DP = DIST2_NoImage( BaseAxes[base1].Pcrd(), BaseAxes[base2].Pcrd() );
      //mprintf(" %i to %i P--P D= %f", BaseAxes[base1].ResNum()+1, BaseAxes[base2].ResNum()+1,
      //        sqrt(dPtoP) );
      DP = sqrt(DP);
      dPtoP = (float)DP;
    }
    //mprintf("\n");
    float dOtoO = 0.0;
    //mprintf("\tDEBUG: %i %i:", BaseAxes[base1].O4idx(), BaseAxes[base2].O4idx() );
    if ( BaseAxes[base1].HasO4atom() && BaseAxes[base2].HasO4atom() ) {
      double DO4 = DIST2_NoImage( BaseAxes[base1].O4crd(), BaseAxes[base2].O4crd() );
      //mprintf(" %i to %i O4'--O4' D= %f", BaseAxes[base1].ResNum()+1, BaseAxes[base2].ResNum()+1,
      //        sqrt(dOtoO) );
      DO4 = sqrt(DO4);
      dOtoO = (float)DO4;
    }
    //mprintf("\n");
    // Calc BP parameters, set up basepair axes
    //calculateParameters(BaseAxes[base1],BaseAxes[base2],&BasePairAxes[nbasepair],Param);
    calculateParameters(BaseAxes[base2],BaseAxes[base1],&BasePairAxes[nbasepair],Param);
    // Store data
    Param[3] *= RADDEG;
    Param[4] *= RADDEG;
    Param[5] *= RADDEG;
    //mprintf("DBG: BP %i # hbonds = %i\n", nbasepair+1, NumberOfHbonds_[nbasepair]);
    // Convert everything to float to save space
    float shear = (float)Param[0];
    float stretch = (float)Param[1];
    float stagger = (float)Param[2];
    float opening = (float)Param[3];
    float prop = (float)Param[4];
    float buckle = (float)Param[5];
    int n_of_hb = NumberOfHbonds_[nbasepair];
    // Add to DataSets
    SHEAR_[nbasepair]->Add(frameNum, &shear);
    STRETCH_[nbasepair]->Add(frameNum, &stretch);
    STAGGER_[nbasepair]->Add(frameNum, &stagger);
    OPENING_[nbasepair]->Add(frameNum, &opening);
    PROPELLER_[nbasepair]->Add(frameNum, &prop);
    BUCKLE_[nbasepair]->Add(frameNum, &buckle);
    BPHBONDS_[nbasepair]->Add(frameNum, &n_of_hb);
    MAJOR_[nbasepair]->Add(frameNum, &dPtoP);
    MINOR_[nbasepair]->Add(frameNum, &dOtoO);
#   ifdef NASTRUCTDEBUG
    // DEBUG - write base pair axes
    basepairaxesfile.WriteAxes(BasePairAxes[nbasepair], base1, BaseAxes[base1].ResName());
#   endif
    ++nbasepair; // Actual base pair count; BP is nbasepair*3
  }

  return 0;
}

// Action_NAstruct::determineBasepairParameters() 
/** For each base pair step, determine values of Tilt, Roll, Twist, Shift,
  * Slide, and Rise.
  */
int Action_NAstruct::determineBasepairParameters(int frameNum) {
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
    // Convert everything to float to save space
    float shift = (float)Param[0];
    float slide = (float)Param[1];
    float rise = (float)Param[2];
    float twist = (float)Param[3];
    float roll = (float)Param[4];
    float tilt = (float)Param[5];
    SHIFT_[bpi]->Add(frameNum, &shift);
    SLIDE_[bpi]->Add(frameNum, &slide);
    RISE_[bpi]->Add(frameNum, &rise);
    TWIST_[bpi]->Add(frameNum, &twist);
    ROLL_[bpi]->Add(frameNum, &roll);
    TILT_[bpi]->Add(frameNum, &tilt);
    // Calc helical parameters
    helicalParameters(BasePairAxes[bpi], BasePairAxes[bpj], Param);
    Param[3] *= RADDEG;
    Param[4] *= RADDEG;
    Param[5] *= RADDEG;
    // Convert to float
    float xdisp = (float)Param[0];
    float ydisp = (float)Param[1];
    float hrise = (float)Param[2];
    float incl = (float)Param[3];
    float tip = (float)Param[4];
    float htwist = (float)Param[5];
    XDISP_[bpi]->Add(frameNum, &xdisp);
    YDISP_[bpi]->Add(frameNum, &ydisp);
    HRISE_[bpi]->Add(frameNum, &hrise);
    INCL_[bpi]->Add(frameNum, &incl);
    TIP_[bpi]->Add(frameNum, &tip);
    HTWIST_[bpi]->Add(frameNum, &htwist);
  }

  return 0;
}
// ----------------------------------------------------------------------------

// Action_NAstruct::init()
// Dataset name will be the last arg checked for. Check order is:
//    1) Keywords
//    2) Masks
//    3) Dataset name
Action::RetType Action_NAstruct::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  Frame* refframe = NULL;
  Topology* refparm = NULL;

  // Get keywords
  outputsuffix_ = actionArgs.GetStringKey("naout");
  double hbcut = actionArgs.getKeyDouble("hbcut", -1);
  if (hbcut > 0) 
    HBdistCut2_ = hbcut * hbcut;
  double origincut = actionArgs.getKeyDouble("origincut", -1);
  if (origincut > 0)
    originCut2_ = origincut * origincut;
  ArgList::ConstArg resrange_arg = actionArgs.getKeyString("resrange");
  if (resrange_arg != NULL)
    if (resRange.SetRange( resrange_arg )) return Action::ERR;
  printheader_ = !actionArgs.hasKey("noheader");
  // Reference for setting up basepairs
  int refindex = actionArgs.getKeyInt("refindex", -1);
  if (actionArgs.hasKey("reference")) refindex = 0;
  std::string refname = actionArgs.GetStringKey("ref");
  if (refindex!=-1 || !refname.empty()) {
    useReference_ = true;
    // Reference by name/tag
    if (!refname.empty())
      refindex = FL->FindName( refname );
    // Get reference by index
    refframe = FL->GetFrame( refindex );
    if (refframe==NULL) {
      mprinterr("Error: nastruct: Could not get ref frame, index=%i\n",refindex);
      return Action::ERR;
    }
    // Get parm for reference
    refparm = FL->GetFrameParm( refindex );
    if (refparm == NULL) {
      mprinterr("Error: nastruct: Could not get parm for frame %s\n", FL->FrameName(refindex));
      return Action::ERR;
    }
  }

  // Get custom residue maps
  ArgList::ConstArg maparg;
  ArgList maplist;
  AxisType::NAbaseType mapbase;
  while ( (maparg = actionArgs.getKeyString("resmap"))!=NULL ) {
    // Split maparg at ':'
    maplist.SetList(maparg,":");
    // Expect only 2 args
    if (maplist.Nargs()!=2) {
      mprinterr("Error: nastruct: resmap format should be '<ResName>:{A,C,G,T,U}' (%s)\n",maparg);
      return Action::ERR;
    }
    // Check that second arg is A,C,G,T,or U
    if      (maplist[1] == "A") mapbase=AxisType::ADE;
    else if (maplist[1] == "C") mapbase=AxisType::CYT;
    else if (maplist[1] == "G") mapbase=AxisType::GUA;
    else if (maplist[1] == "T") mapbase=AxisType::THY;
    else if (maplist[1] == "U") mapbase=AxisType::URA;
    else {
      mprinterr("Error: nastruct: resmap format should be '<ResName>:{A,C,G,T,U}' (%s)\n",maparg);
      return Action::ERR;
    }
    // Check that residue name is <= 4 chars
    std::string resname = maplist[0]; 
    if (resname.size() > 4) {
      mprinterr("Error: nastruct: resmap resname > 4 chars (%s)\n",maparg);
      return Action::ERR;
    }
    // Format residue name
    // TODO: Use NameType in map
    NameType mapresname = resname;
    resname.assign( *mapresname );
    mprintf("\tCustom Map: [%s]\n",resname.c_str());
    //maplist.PrintList();
    // Add to CustomMap
    ResMapType::iterator customRes = CustomMap.find(resname);
    if (customRes!=CustomMap.end()) {
      mprintf("Warning: nastruct: resmap: %s already mapped.\n",resname.c_str());
    } else {
      CustomMap.insert( std::pair<std::string,AxisType::NAbaseType>(resname,mapbase) );
    }
  }
  // Get Masks
  // Dataset
  dataname_ = actionArgs.GetStringNext();
  if (dataname_.empty())
    dataname_ = DSL->GenerateDefaultName("NA");
  // DataSets are added to data file list in print()

  mprintf("    NAstruct: ");
  if (resRange.Empty())
    mprintf("Scanning all NA residues");
  else
    mprintf("Scanning residues %s",resRange.RangeArg());
  if (!outputsuffix_.empty()) {
      mprintf(", formatted output using file suffix %s",outputsuffix_.c_str());
    if (!printheader_) mprintf(", no header");
  }
  mprintf(".\n");
  mprintf("\tHydrogen bond cutoff for determining base pairs is %.2lf Angstroms.\n",
          sqrt( HBdistCut2_ ) );
  mprintf("\tBase reference axes origin cutoff for determining base pairs is %.2lf Angstroms.\n",
          sqrt( originCut2_ ) );

  // Use reference to determine base pairing
  if (useReference_) {
    mprintf("\tUsing reference %s to determine base-pairing.\n",FL->FrameName(refindex));
    if (Setup(refparm, 0)) return Action::ERR;
    // Set up base axes
    if ( setupBaseAxes(refframe) ) return Action::ERR;
    // Determine Base Pairing
    if ( determineBasePairing() ) return Action::ERR;
    mprintf("\tSet up %zu base pairs.\n", BasePairAxes.size() ); 
  } else {
    mprintf("\tUsing first frame to determine base pairing.\n");
  }
  masterDSL_ = DSL;
  return Action::OK;
}

// Action_NAstruct::setup()
/** Determine the number of NA bases that will be analyzed, along with 
  * the masks that correspond to the reference frame atoms.
  */
Action::RetType Action_NAstruct::Setup(Topology* currentParm, Topology** parmAddress) {
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
    mprinterr("Error: NAstruct::setup: No residues specified for %s\n",currentParm->c_str());
    return Action::ERR;
  }

  // DEBUG - print all residues
  //if (debug>0)
  //  actualRange.PrintRange("    NAstruct: NA res:",1);

  // Set up frame to hold reference coords for each NA residue
  for (Range::const_iterator resnum = actualRange.begin();
                             resnum != actualRange.end(); ++resnum)
  {
    customBaseType = AxisType::UNKNOWN_BASE;
    // Check if the residue at resnum matches any of the custom maps
    if (!CustomMap.empty()) {
      resname.assign( currentParm->ResidueName(*resnum) );
      ResMapType::iterator customRes = CustomMap.find( resname );
      if (customRes!=CustomMap.end()) {
        mprintf("\tCustom map found: %i [%s]\n",*resnum+1,(*customRes).first.c_str());
        customBaseType = (*customRes).second;
      }
    }
    // Set up ref coords in correct order, along with corresponding 
    // parm mask for this residue. SetRefCoord should overwrite all
    // previously set up information in axis and Mask.
    refreturn = axis.SetRefCoord( currentParm, *resnum, Mask, fitMask, customBaseType );
    // If not recognized as a NA residue, continue to next.
    // Print a warning if the user specified this range.
    if ( refreturn == AxisType::NA_UNKNOWN ) {
      if (!resRange.Empty()) {
        mprintf("Warning: Residue %i:%s not recognized as NA residue.\n",
                *resnum+1, currentParm->ResidueName(*resnum));
      }
      continue;
    } else if ( refreturn == AxisType::NA_ERROR  ) {
      mprinterr("Error: NAstruct::setup: Could not get ref coords for %i:%s\n",
                *resnum+1, currentParm->ResidueName(*resnum));
      return Action::ERR;
    }
    if (Mask.None()) {
      mprintf("Error:: NAstruct::setup: No atoms found for residue %i:%s\n",
              *resnum+1, currentParm->ResidueName(*resnum));
      return Action::ERR;
    }
    if (fitMask.None()) {
      mprintf("Error:: NAstruct::setup: No fit atoms found for residue %i:%s\n",
              *resnum+1, currentParm->ResidueName(*resnum));
      return Action::ERR;
    }
    RefCoords.push_back( axis );
    ExpMasks.push_back( Mask );
    // Determine the largest residue for setting up frames for RMS fit later.
    // ExpMask is larger than or equal to FitMask so use that.
    if (Mask.Nselected() > maxResSize_)
      maxResSize_ = Mask.Nselected();
    FitMasks.push_back( fitMask );
    if (debug_>1) {
      mprintf("\tNAstruct: Res %i:%s ",*resnum+1,currentParm->ResidueName(*resnum));
      Mask.PrintMaskAtoms("NAmask");
      mprintf("\t          Ref %i:%s ",*resnum+1,axis.ResName());
      axis.PrintAtomNames();
    }

    // Set up empty frame to hold input coords for this residue.
    // Initially set to be a copy of the reference frame. The coordinates
    // will be overwritten later.
    BaseAxes.push_back( axis );
  } // End Loop over NA residues

  // NOTE: RefCoords size is also BaseAxes, ExpFrames, and ExpMasks size.
  mprintf("\tSet up %zu bases.\n", RefCoords.size());

  return Action::OK;  
}

// Action_NAstruct::action()
Action::RetType Action_NAstruct::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {

  // Set up base axes
  if ( setupBaseAxes(currentFrame) ) return Action::ERR;

  if (!useReference_) {
    // Determine Base Pairing based on first frame
    if ( determineBasePairing() ) return Action::ERR;
    useReference_ = true;
  } else {
    // Base pairing determined from ref. Just calc # hbonds for each pair.
    int bp = 0;
    for (int bp3 = 0; bp3 < (int)BasePair.size(); bp3 += 3) 
      NumberOfHbonds_[bp++] = basesArePaired(&BaseAxes[BasePair[bp3]], &BaseAxes[BasePair[bp3+1]]); 
  }

  // Determine base parameters
  determineBaseParameters(frameNum);

  // Determine base pair parameters
  determineBasepairParameters(frameNum);

  return Action::OK;
} 

// Action_NAstruct::print()
void Action_NAstruct::Print() {
  CpptrajFile outfile;
  int nframes;
  if (outputsuffix_.empty()) return;

  // ---------- Base pair parameters ----------
  std::string outfilename = "BP." + outputsuffix_;
  // Check that there is actually data
  if ( SHEAR_.empty() || SHEAR_[0]->Empty() )
    mprinterr("Error: nastruct: Could not write BP file %s: No BP data.\n",outfilename.c_str()); 
  else {
    if (outfile.OpenWrite( outfilename ) == 0) {
      // Determine number of frames from SHEAR[0] DataSet
      nframes = SHEAR_[0]->Size();
      mprintf("\tBase pair output file %s; %i frames, %zu base pairs.\n", 
              outfilename.c_str(), nframes, BasePair.size() / 3);
      //  File header
      if (printheader_)
        outfile.Printf("%-8s %8s %8s %10s %10s %10s %10s %10s %10s %2s\n","#Frame","Base1","Base2",
                       "Shear","Stretch","Stagger","Buckle","Propeller","Opening", "HB");
      // Loop over all frames
      for (int frame = 0; frame < nframes; ++frame) {
        int nbp = 0;
        for (unsigned int bpidx = 0; bpidx < BasePair.size(); bpidx+= 3) {
          int bp_1 = BasePair[bpidx  ];
          int bp_2 = BasePair[bpidx+1];
          // FIXME: Hack for integer
          int n_of_hb = (int)BPHBONDS_[nbp]->Dval(frame);
          outfile.Printf(BP_OUTPUT_FMT, frame+OUTPUTFRAMESHIFT, 
                         BaseAxes[bp_1].ResNum()+1, BaseAxes[bp_2].ResNum()+1,
                         SHEAR_[nbp]->Dval(frame), STRETCH_[nbp]->Dval(frame),
                         STAGGER_[nbp]->Dval(frame), BUCKLE_[nbp]->Dval(frame),
                         PROPELLER_[nbp]->Dval(frame), OPENING_[nbp]->Dval( frame),
                         n_of_hb);
          ++nbp;
        }
        outfile.Printf("\n");
      }
      outfile.CloseFile();
    } else {
      mprinterr("Error: nastruct: Could not open %s for writing.\n", outfilename.c_str());
    }
  }

  // ---------- Base pair step parameters ----------
  CpptrajFile outfile2;
  outfilename = "BPstep." + outputsuffix_;
  std::string outfilename2 = "Helix." + outputsuffix_;
  // Check that there is actually data
  // TODO: Check helix data as well
  if ( SHIFT_.empty() || SHIFT_[0]->Empty() )
    mprinterr("Error: nastruct: Could not write BPstep / helix files: No data.\n"); 
  else {
    int err = 0;
    err += outfile.OpenWrite( outfilename );
    err += outfile2.OpenWrite( outfilename2 );
    if (err == 0) {
      // Determine number of frames from SHIFT[0] DataSet. Should be same as SHEAR.
      nframes = SHIFT_[0]->Size();
      mprintf("\tBase pair step output file %s;",outfilename.c_str());
      mprintf("Helix output file %s; %i frames, %zu base pair steps.\n", outfilename2.c_str(),
              nframes, BasePairAxes.size() - 1);
      //  File headers
      if (printheader_) {
        outfile.Printf("%-8s %-9s %-9s %10s %10s %10s %10s %10s %10s\n","#Frame","BP1","BP2",
                       "Shift","Slide","Rise","Tilt","Roll","Twist");
        outfile2.Printf("%-8s %-9s %-9s %10s %10s %10s %10s %10s %10s\n","#Frame","BP1","BP2",
                        "X-disp","Y-disp","Rise","Incl.","Tip","Twist");
      }
      // Loop over all frames
      for (int frame = 0; frame < nframes; ++frame) {
        int nstep = 0;
        for (unsigned int bpi = 0; bpi < BasePairAxes.size() - 1; bpi++) {
          unsigned int bpj = bpi+1;
          // BPstep write
          outfile.Printf(NA_OUTPUT_FMT, frame+OUTPUTFRAMESHIFT, 
                         BasePairAxes[bpi].ResNum()+1, BasePairAxes[bpi].ResNum2()+1,
                         BasePairAxes[bpj].ResNum()+1, BasePairAxes[bpj].ResNum2()+1,
                         SHIFT_[nstep]->Dval(frame), SLIDE_[nstep]->Dval(frame),
                         RISE_[nstep]->Dval(frame), TILT_[nstep]->Dval(frame),
                         ROLL_[nstep]->Dval(frame), TWIST_[nstep]->Dval(frame));
          // Helix write
          outfile2.Printf(NA_OUTPUT_FMT, frame+OUTPUTFRAMESHIFT,
                          BasePairAxes[bpi].ResNum()+1, BasePairAxes[bpi].ResNum2()+1,
                          BasePairAxes[bpj].ResNum()+1, BasePairAxes[bpj].ResNum2()+1,
                          XDISP_[nstep]->Dval(frame), YDISP_[nstep]->Dval(frame),
                          HRISE_[nstep]->Dval(frame), INCL_[nstep]->Dval(frame),
                          TIP_[nstep]->Dval(frame), HTWIST_[nstep]->Dval(frame));
          ++nstep;
        }
        outfile.Printf("\n");
        outfile2.Printf("\n");
      }
      outfile.CloseFile();
      outfile2.CloseFile();
    } else {
      mprinterr("Error: nastruct: Could not open BPstep files for writing.\n");
    }
  }
}

