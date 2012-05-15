// ptraj_analyze.c
// ORIGINALLY FROM PTRAJ

/*  ________________________________________________________________________
 */


/*
 *  The code defined herein is used to analyze and data accumulated by the
 *  various actions during the trajectory processing.
 *
 *  Among others, the following routines (along with supplemental routines
 *  as necessary) are defined:
 *
 *  analyze matrix  --- diagonalization of matrices generated with the 
 *                      "matrix" command (Holger Gohlke, Scripps)
 *
 *  analyze modes   --- vibrational analysis of modes generated with the
 *                      "analyze matrix" command (Holger Gohlke, Scripps)
 *
 */

// ---------- CSTDLIB includes -------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// ---------- Defines ----------------------------------------------------------
#define ANALYZE_MODULE

#ifdef CLINK_PLAIN
#  define dsaupd_ dsaupd
#  define dseupd_ dseupd
#  define dspev_  dspev
#  define thermo_ thermo
#  define cffti_  cffti
#  define cfftf_  cfftf
#  define cfftb_  cfftb
#endif
#ifdef CLINK_CAPS
#  define dsaupd_ DSAUPD
#  define dseupd_ DSEUPD
#  define dspev_  DSPEV
#  define thermo_ THERMO
#  define cffti_  CFFTI
#  define cfftf_  CFFTF
#  define cfftb_  CFFTB
#endif

// ---------- PTRAJ includes ---------------------------------------------------
#include "ptraj_analyze.h"
#include "ptraj_scalar.h"
#include "ptraj_arg.h"
#include "ptraj_common.h"

// ---------- STATIC variables -------------------------------------------------
static const char pucker_ss[10][9] = { 
  "C3'-endo", "C4'-exo ", "O4'-endo", "C1'-exo ", "C2'-endo", "C3'-exo ", 
  "C4'-endo", "O4'-exo ", "C1'-endo", "C2'-exo " 
};
static const char torsion_ss_2D[6][6][6] = {
  {"g+ g+", "g+ a+", "g+ t", "g+ a-", "g+ g-", "g+ c"},
  {"a+ g+", "a+ a+", "a+ t", "a+ a-", "a+ g-", "a+ c"},
  {"t  g+", "t  a+", "t  t", "t  a-", "t  g-", "t  c"},
  {"a- g+", "a- a+", "a- t", "a- a-", "a- g-", "a- c" },
  {"g- g+", "g- a+", "g- t", "g- a-", "g- g-", "g- c" },
  {"c  g+", "c  a+", "c  t", "c  a-", "c  g-", "c  c" }
};
static const char torsion_ss[6][8] = {
  "g+     ", "a+     ", "t      ", "a-     ", "g-     ", "c      "
};
static const double torsion_offset[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 180.0 };
static const char distance_ss_2D[6][6][9] = {
  {"< 2, < 2", "< 2, 2-3", "< 2, 3-4", "< 2, 4-5", "< 2, 5-6", "< 2, > 6" },
  {"2-3, < 2", "2-3, 2-3", "2-3, 3-4", "2-3, 4-5", "2-3, 5-6", "2-3, > 6" },
  {"3-4, < 2", "3-4, 2-3", "3-4, 3-4", "3-4, 4-5", "3-4, 5-6", "3-4, > 6" },
  {"4-5, < 2", "4-5, 2-3", "4-5, 3-4", "4-5, 4-5", "4-5, 5-6", "4-5, > 6" },
  {"5-6, < 2", "5-6, 2-3", "5-6, 3-4", "5-6, 4-5", "5-6, 5-6", "5-6, > 6" },
  {"> 6, < 2", "> 6, 2-3", "> 6, 3-4", "> 6, 4-5", "> 6, 5-6", "> 6, > 6" }
};
static const char distance_ss[6][8] = {
  " < 2.5 ", "2.5-3.5", "3.5-4.5", "4.5-5.5", "5.5-6.5", " > 6.5 "
};
// =============================================================================
   int
analyzeTest(analyzeInformation *analyze, stackType *sp, int mode)
{
  argStackType **argumentStackPointer;

  if (mode == PTRAJ_SETUP) {

    /*
     *  ANALYZE: PTRAJ_SETUP
     *
     *  Parse arguments off the stack and fill in any necessary
     *  information in the analyzeInformation structure.
     *
     *  This mode is invoked by ptrajSetupAnalyze().  The current
     *  argumentStack is passed in as action->carg1.  The
     *  "analyze" structure is already partially setup via
     *  the following code in ptrajSetupAnalyze():
     *
     *     statep = ptrajCurrentState();
     *     state = *statep;
     *
     *     analyze = (analyzeInformation *)
     *        safe_malloc(sizeof(analyzeInformation));
     *     INITIALIZE_analyzeInformation(analyze);
     *
     *     analyze->state = ptrajCopyState(ptrajCurrentState());
     *
     *     analyze->type = TRANSFORM_TEST;
     *     analyze->fxn  = (analyzeFunction) analyzeTest;
     */

    argumentStackPointer = (argStackType **) analyze->carg1;
    analyze->carg1 = NULL;

    //if (prnlev > 4) {
    //  printStack(argumentStackPointer, printString, NULL);
    //}

    /*
     *  See the comments on argument processing in dispatch.c or in the similar
     *  routine actionTest() in actions.c
     *
     *  Note that after argument processing, assuming successful return 
     *  (i.e. return 0 is performed by this routine), the analyze structure is
     *  placed on the transformAnalyzeStack by ptrajSetupAnalyze() with
     *  the following code:
     *
     *     pushBottomStack( &transformAnalyzeStack, (void *) action );
     *
     *  If an error was encountered during argument processing, return -1
     *  and this action will be freed and not placed on the 
     *  transformActionStack
     */

    printf("In analyzeTest, PTRAJ_SETUP mode\n");

    return 0;
  }


  if (mode == PTRAJ_STATUS) {

    /*
     *  Print out a summary of information about this action
     */
    printf("In analyzeTest, PTRAJ_STATUS mode\n");
  }

  /*
   *  Other possible modes:
   *
   *  PTRAJ_PRINT   -- dump information to output files
   *  PTRAJ_CLEANUP -- free any allocated memory
   *  PTRAJ_NOOP    -- NULL mode
   */

  if (mode != PTRAJ_ACTION) return 0;

  /*
   *  Perform action on coordinates, etc.
   */

  printf("In analyzeTest, PTRAJ_ACTION mode\n");
  return 1;


}

   int
analyzeCrankshaft(analyzeInformation *analyze, stackType *sp, int mode)
{
  /*
   *  usage:
   *
   *  analyze crankshaft {angle | distance} <scalar-name1> <scalar-name2> info <string>
   *
   *  argument usage:
   *
   *    iarg1: = 0 angle, =1 distance
   *    iarg2: start
   *    iarg3: stop
   *    iarg4: offset
   *    carg1: scalar1
   *    carg2: scalar2
   *    carg3: info string
   */

  argStackType **argumentStackPointer;
  char *buffer, *info;
  scalarInfo *scalar1, *scalar2;
  stackType *stack;

  double v1_avg[6][6], v2_avg[6][6], v1_sd[6][6], v2_sd[6][6];
  double v1, v2, initial_v1, initial_v2, final_v1, final_v2;
  int visits[6][6], transitions[6][6];
  int previous_i1, previous_i2;
  int initial_i1, initial_i2;
  int final_i1, final_i2;

  int i,j,k,i1,i2;
  char *filename;
  FILE *outfile;

  int start, stop, offset, totalFrames;

  if (mode == PTRAJ_SETUP) {

    /*
     *  PTRAJ_SETUP:
     */


    argumentStackPointer = (argStackType **) analyze->carg1;
    analyze->carg1 = NULL;
    
    info = argumentStackKeyToString(argumentStackPointer, "info", NULL);
    analyze->carg3 = info;

    if (argumentStackContains(argumentStackPointer, "angle"))
      analyze->iarg1 = 0;
    if (argumentStackContains(argumentStackPointer, "distance"))
      analyze->iarg1 = 1;

    filename = argumentStackKeyToString( argumentStackPointer, "out", NULL );
    analyze->carg4 = (void *) filename;

    analyze->iarg2 = argumentStackKeyToInteger(argumentStackPointer, "start", 1);
    analyze->iarg3 = argumentStackKeyToInteger(argumentStackPointer, "stop", -1);
    analyze->iarg4 = argumentStackKeyToInteger(argumentStackPointer, "offset", 1);


    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL) {
      fprintf(stdout, "ptraj(), analyzeCrankshaft: No name specified for angle 1\n");
      return -1;
    }

    scalar1 = scalarStackGetName(&sp, buffer);
    if (scalar1 == NULL) {
      fprintf(stdout, "ptraj(), analyzeCrankshaft: Name (%s) not found, ignoring\n",
	      buffer);
      return -1;
    }
    analyze->carg1 = (void *) scalar1;

    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL) {
      fprintf(stdout, "ptraj(), analyzeCrankshaft: No name specified for angle 2\n");
      return -1;
    }

    scalar2 = scalarStackGetName(&sp, buffer);
    if (scalar2 == NULL) {
      fprintf(stdout, "ptraj(), analyzeCrankshaft: Name (%s) not found, ignoring\n",
	      buffer);
      return -1;
    }

    analyze->carg2 = (void *) scalar2;
    return 0;  
  }

  scalar1 = (scalarInfo *) analyze->carg1;
  scalar2 = (scalarInfo *) analyze->carg2;
  info = (char *) analyze->carg3;
  filename = (char *) analyze->carg4;

  if (mode == PTRAJ_STATUS) {

    /*
     *  PTRAJ_STATUS:
     */
    fprintf(stdout, "  ANALYZE CRANKSHAFT: %s ", info);
    fprintf(stdout, "%s named %s and %s\n", (analyze->iarg1 ? "distances" : "angles"),
	    scalar1->name, scalar2->name);
    return 0;
  }

  /*
   *  Other possible modes:
   *
   *  PTRAJ_CLEANUP -- free any allocated memory
   *  PTRAJ_NOOP    -- NULL mode
   */




  if (mode != PTRAJ_ACTION) return 0;

  /*
   *  PTRAJ_ACTION:
   */

  start = analyze->iarg2-1;
  stop = analyze->iarg3;
  if (stop < 0)
    stop = scalar1->totalFrames;
  offset = analyze->iarg4;
  totalFrames = (stop - start) / offset;

  if (start >= scalar1->totalFrames) return 0;

  for (i=start; i < stop; i+=offset) {

    if (i==start) {

      /*
       *  do initialization
       */

      for (j=0;j<6;j++) {
	for (k=0;k<6;k++) {
	  v1_avg[j][k] = 0.0;
	  v2_avg[j][k] = 0.0;
	  v1_sd[j][k] = 0.0;
	  v2_sd[j][k] = 0.0;
	
	  visits[j][k] = 0;
	  transitions[j][k] = 0;
	}
      }

      if (filename != NULL) {
	outfile = safe_fopen(filename, "w");
	if (outfile == NULL) {
	  fprintf(stderr, "WARNING in ptraj(), analyze crankshaft: error on opening %s for output\n",
		  filename);
	  safe_free(filename);
	  filename = NULL;
	  analyze->carg4 = NULL;
	}
      }
    }
      
    v1 = scalar1->value[i];
    v2 = scalar2->value[i];

    if (analyze->iarg1) {  
      /*
       *  this is a DISTANCE
       */

      i1 = (v1 - 1.0) / 1;     /*  the current algorithm is a test and aims to bin things   */
      if (i1 > 5) i1 = 5;      /*  from 0->5 starting from values < 2A, in increments of 1  */
      i2 = (v2 - 1.0) / 1;     /*  angstrom to > 6 A.  i.e. value-1.0/1                     */
      if (i2 > 5) i2 = 5;

    } else { 
      /*
       *  this is an ANGLE
       *    -- bin from 0->5
       *    -- subtract 30 from value, such that 0->60 = g+, 60->120 = a+, 120->180 = t
       */

      v1 -= 30.0;
      v2 -= 30.0;

      if (v1 < 0)    v1 += 360.0;
      if (v1 > 360)  v1 -= 360.0;
      if (v2 < 0)    v2 += 360.0;
      if (v2 > 360)  v2 -= 360.0;

      i1 = v1 / 60;
      i2 = v2 / 60;

    }

    /*
     *  Store initial and final bins/values
     */
    if (i == start) {
      initial_i1 = i1;
      initial_i2 = i2;
      initial_v1 = scalar1->value[i];
      initial_v2 = scalar2->value[i];
    }

    final_i1 = i1;
    final_i2 = i2;
    final_v1 = scalar1->value[i];
    final_v2 = scalar2->value[i];

    /*
     *  debugging info
     */
    if (prnlev > 5) {
      fprintf(stdout, "Binning %s values %6.2f %6.2f into %i x %i\n",
	      (analyze->iarg1 ? "distance" : "angle"), 
	      scalar1->value[i], scalar2->value[i], i1, i2);
    }

    /*
     *  update bin counter and averages/standard deviations
     */
    visits[i1][i2] += 1;

    v1 = scalar1->value[i];
    v2 = scalar2->value[i]; 
    v1 += torsion_offset[i1];
    v2 += torsion_offset[i2];
    v1_avg[i1][i2] += v1;
    v2_avg[i1][i2] += v2;

    v1_sd[i1][i2]  = v1_sd[i1][i2] + (v1*v1);
    v2_sd[i1][i2]  = v2_sd[i1][i2] + (v2*v2);

    if (filename != NULL) {
      /*
       *  hack to map substate numbers 0->36 to a-z0-9
       *
       *  j = i1*6 + i2 + 97;
       *  if (j > 122)
       *  j -= 75;
       *
       *       g+    a+    t    a-    g-    c
       *  g+   0     1     2    3     4     5
       *  a+   6     7     8    9    10    11
       *   
       */
      fprintf(outfile, "%7i %i\n", i, i1*6 + i2);
    }

    /*
     *  check for transitions from one bin to another
     */
    if (i > start) {
      if (i1 != previous_i1 || i2 != previous_i2) {

	transitions[previous_i1][previous_i2] += 1;

	if ( ! (i1 != previous_i1 && i2 != previous_i2) &&
	     (i1 != previous_i1 && (abs(i1 - previous_i1) == 1)) ||
	     (i2 != previous_i2 && (abs(i2 - previous_i2) == 1)) ||
	     (i1 != previous_i1 && (abs(i1 - previous_i1) == 5)) ||
	     (i2 != previous_i2 && (abs(i2 - previous_i2) == 5)) ) {
	  
	  if (prnlev > 3) {
	    if (analyze->iarg1) {
	      fprintf(stdout, "SMALL TRANSITION frame %6i (%s,%s): substate (%s) to (%s)\n",
		      i, scalar1->name, scalar2->name, 
		      distance_ss_2D[previous_i1][previous_i2], distance_ss_2D[i1][i2]);
	    } else {
	      fprintf(stdout, "SMALL TRANSITION frame %6i (%s,%s): substate (%s) to (%s)\n",
		      i, scalar1->name, scalar2->name, 
		      torsion_ss_2D[previous_i1][previous_i2], torsion_ss_2D[i1][i2]);
	    }
	  }
	} else {

	  if (prnlev > 2) {
	    if (analyze->iarg1) {
	      fprintf(stdout, "LARGE TRANSITION frame %6i (%s,%s): substate (%s) to (%s)\n",
		      i, scalar1->name, scalar2->name, 
		      distance_ss_2D[previous_i1][previous_i2], distance_ss_2D[i1][i2]);
	    } else {
	      fprintf(stdout, "LARGE TRANSITION frame %6i (%s,%s): substate (%s) to (%s)\n",
		      i, scalar1->name, scalar2->name, 
		      torsion_ss_2D[previous_i1][previous_i2], torsion_ss_2D[i1][i2]);
	    }

	  }
	}
      }
    }

    previous_i1 = i1;
    previous_i2 = i2;
  }

  /*
   *  after processing all frames, compute the averages and standard deviations
   */
  for (j=0;j<6;j++) {
    for (k=0;k<6;k++) {

      if (visits[j][k]) {
	v1_avg[j][k] = v1_avg[j][k]/visits[j][k];
	v2_avg[j][k] = v2_avg[j][k]/visits[j][k];
      }
      if (visits[j][k] > 1) {
	v1_sd[j][k] = v1_sd[j][k]/visits[j][k];
	v2_sd[j][k] = v2_sd[j][k]/visits[j][k];
	v1_sd[j][k] = sqrt( v1_sd[j][k] - v1_avg[j][k]*v1_avg[j][k] );
	v2_sd[j][k] = sqrt( v2_sd[j][k] - v2_avg[j][k]*v2_avg[j][k] );
      } else {
	v1_sd[j][k] = 0.0;
	v2_sd[j][k] = 0.0;
      }
      if (visits[j][k]) {
	v1_avg[j][k] -= torsion_offset[j];
	v2_avg[j][k] -= torsion_offset[k];
      }
    }
  }

  if (filename != NULL) {
    safe_fclose(outfile);
    safe_free(filename);
    analyze->carg4 = NULL;
  }

  /*
   *  print out results
   */
  fprintf(stdout, "\n\nCRANKSHAFT: %s.\n", info);
  fprintf(stdout, "  start at frame %i, stop after frame %i, offset between frames is %i.\n", 
	  start+1, stop, offset);

  fprintf(stdout, "  total number of frames is %i.  Table values are\n", totalFrames);
  fprintf(stdout, "  %%occupied, #transitions to another substate, average angles and stddev\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "  INITIAL VALUE: %s (%6.1f, %6.1f)\n",
	  (analyze->iarg1 ? distance_ss_2D[initial_i1][initial_i2] : torsion_ss_2D[initial_i1][initial_i2]),
	  initial_v1, initial_v2);
  fprintf(stdout, "  FINAL VALUE:   %s (%6.1f, %6.1f)\n\n",
	  (analyze->iarg1 ? distance_ss_2D[final_i1][final_i2] : torsion_ss_2D[final_i1][final_i2]),
	  final_v1, final_v2);

    /*
     *  supplementary information based on type of crankshaft
     */
  if (scalar1->type == SCALAR_TYPE_EPSILON && scalar2->type == SCALAR_TYPE_ZETA) {

      /*
       *  epsilon/zeta in nucleic acids!
       */
    i1 = 0;
    i2 = 0;
    for (i=start; i < stop; i+=offset) {
      v1 = scalar1->value[i];
      v2 = scalar2->value[i];
      if (v1 < 0) v1 += 360.0;
      if (v2 < 0) v2 += 360.0;
      v1 = v1 - v2;
      if (v1 > 60 && v1 < 120) i2++;
      if (v1 < -30 && v1 > -120) i1++;
    } 
    fprintf(stdout, "    EPSILON/ZETA crankshaft\n");
    fprintf(stdout, "      BI  = (t, g-) or eps-zeta ~ -90 [currently = %.1f%%]\n",
	    i1*100.0 / totalFrames);
    fprintf(stdout, "      BII = (g-, t) or eps-zeta ~ +90 [currently = %.1f%%]\n\n",
	    i2*100.0 / totalFrames);


  } else if (scalar1->type == SCALAR_TYPE_ALPHA && scalar2->type == SCALAR_TYPE_GAMMA) {

      /*
       *  alpha/gamma in nucleic acids!
       */
    fprintf(stdout, "    ALPHA/GAMMA crankshaft\n");
    fprintf(stdout, "      canonical is (g-, g+) [currently at %.1f%%]\n",
	    visits[4][0]*100.0/totalFrames);
    fprintf(stdout, "      other possible states are (t, t) {%.1f%%} or (t, g-) {%.1f%%}\n",
	    visits[2][2]*100.0/totalFrames,
	    visits[2][4]*100.0/totalFrames);
    fprintf(stdout, "      (g+, t) is found < 5%% in protein/DNA complexes {%.1f%%}\n",
	    visits[0][2]*100.0/totalFrames);
    if ( (visits[0][2] + visits[1][2])*100.0/totalFrames > 10.0 ) 
      fprintf(stdout, "    *** > 10%% population in (g+, t) / (a+, t) states!!!\n");
    fprintf(stdout, "\n");
  }
    

  fprintf(stdout, "                %s         %s         %s         %s         %s          %s\n", 
	  torsion_ss[0],torsion_ss[1],torsion_ss[2],
	  torsion_ss[3],torsion_ss[4],torsion_ss[5]);
  fprintf(stdout, "        -------------------------------------------------------------------------------------------------\n");
  for (i=0; i < 6; i++) {

    fprintf(stdout, "        |");
    for (j=0; j < 6; j++) {
      if (visits[i][j] == 0)
	fprintf(stdout, "               |");
      else
	fprintf(stdout, "   %7.1f%%    |", 100.0 * visits[i][j]/totalFrames);
    }
    fprintf(stdout, "\n");

    fprintf(stdout, " %s|", torsion_ss[i]);
    for (j=0; j < 6; j++) {
      if (visits[i][j] == 0)
	fprintf(stdout, "               |");
      else
	fprintf(stdout, " %8i      |", transitions[i][j]);
    }
    fprintf(stdout, "\n");
    
    fprintf(stdout, "        |");
    for (j=0; j < 6; j++) {
      if (visits[i][j] == 0)
	fprintf(stdout, "               |");
      else
	fprintf(stdout, " %6.1f %6.1f |", v1_avg[i][j], v2_avg[i][j]);
    }
    fprintf(stdout, "\n");

    fprintf(stdout, "        |");
    for (j=0; j < 6; j++) {
      if (visits[i][j] < 2)
	fprintf(stdout, "               |");
      else
	fprintf(stdout, " %6.1f %6.1f |", v1_sd[i][j], v2_sd[i][j]);
    }
    fprintf(stdout, "\n");
    
    fprintf(stdout, "        |-----------------------------------------------------------------------------------------------|\n");
  }
  
  return 1;

}


/** ANALYZE ROUTINE *************************************************************
 *
 *  analyzeModes() --- vibrational analysis
 *
 *  Supplementary routines:
 *    readEvecFile (in evec.h/.c)
 *    freeAnalyzeModesMemory (below)
 *
 ******************************************************************************/
typedef enum _analyzeModesType {
  ANALYZEMODES_UNKNOWN,
  ANALYZEMODES_FLUCT,
  ANALYZEMODES_DISPL,
  ANALYZEMODES_CORR
} analyzeModesType;

   void
freeAnalyzeModesMemory(analyzeInformation *analyze){

  modesInfo *modinfo;

  modinfo = (modesInfo *) analyze->carg1;
  if(modinfo != NULL && analyze->darg2 > 0.0){
    /*
     * source == MS_FILE -> free modinfo 
     * (source == MS_STACK will be handled by analyzeMatrix)
     */

    if(modinfo->name != NULL)
      safe_free(modinfo->name);
    if(modinfo->avg != NULL)
      safe_free(modinfo->avg);
    if(modinfo->freq != NULL)
      safe_free(modinfo->freq);
    if(modinfo->evec != NULL)
      safe_free(modinfo->evec);

    INITIALIZE_modesInfo(modinfo);
    safe_free(modinfo);
  }

  if(analyze->carg2 != NULL)
    safe_free(analyze->carg2);
  if(analyze->carg3 != NULL)
    safe_free(analyze->carg3);

  if(analyze->carg4 != NULL){
    /* extern void clearStack( stackType ** ); */
    clearStack( (stackType **) &(analyze->carg4) );
  }

}

   int
analyzeModes(analyzeInformation *analyze, stackType *sp, int mode)
{
  /* Constant definitions */
  const double TKBC2 = 0.46105E-34;
  const double AVO   = 6.023E23;
  const double CNST  = TKBC2 * AVO;
  const double CMTOA = 1.000E8;
  const double TWOPI = 6.2832;
  const double CONT  = CMTOA / TWOPI;
  const double CONSQ = 2.39805E-3; /* = hc/2kT in cm, with T=300K; use for quantum Bose statistics) */

  /* Variable declarations */
  argStackType **argumentStackPointer;
  char *buffer;

  modesInfo *modinfo;
  modesSource source;
  analyzeModesType type;
  FILE *fp;
  stackType *atompairStack, *apsp;
  ptrajState *state;
  int i, j, k, ind, indi1, indi2, indj1, indj2, ncnt;
  int nvect, nvectelem, ibose, ncoords, ibeg, iend, natoms;
  int *mask1, *mask2, *mp, nm1, nm2, at1, at2, npair;
  double factor;
  double f, fi, sumx, sumy, sumz, distx, disty, distz, argq, qcorr;
  double *freq, *evec, *avg, *results;
  double sqrtcnst;
  double dnorm, e[3], del[3][3], val;
  char *outfile;

  /*
   *  USAGE:
   *
   *  analyze modes fluct|displ|corr
   *                            stack <stackname> | file <filename> 
   *                            [beg <beg>] [end <end>] 
   *                            [bose] [factor <factor>]
   *                            [out <outfile>] [maskp <mask1> <mask2> [...]]
   *  
   *  - fluct: rms fluctations from normal modes
   *  - displ: displacement of cartesian coordinates along normal mode directions
   *
   *  analyze argument usage:
   *    iarg1:
   *      analysis type
   *    iarg2:
   *      start value for inclusion of modes (beg)
   *    iarg3:
   *      stop value for inclusion of modes (end)
   *    iarg4:
   *      0 -- Boltzmann statistics
   *      1 -- Bose statistics
   *
   *    darg1:
   *      factor
   *    darg2:
   *      > 0.0: modes read from file,   modesInfo will     be freed in freeAnalyzeModesMemory 
   *      < 0.0: modes taken from stack, modesInfo will not be freed in freeAnalyzeModesMemory
   *                                       (but in freeAnalyzeMatrixMemory)
   *
   *    carg1:
   *      pointer to modesInfo
   *    carg2:
   *      pointer to outfile name
   *    carg3:
   *      pointer to results vector
   *    carg4:
   *      pointer to mask pair stack
   *
   *  results vector usage:
   *    fluct:
   *      [rmsx(atom1), rmsy(atom1), rmsz(atom1), rms(atom1), ..., rmsx(atomN), ..., rms(atomN)]
   *    displ:
   *      [displx(atom1), disply(atom1), displz(atom1), ..., displx(atomN), ..., displz(atomN)]
   *    corr:
   *      [corr(pair1, vec1), ..., corr(pair1, vecN), ..., corr(pairM, vec1), ..., corr(pairM, vecN)]
   */

  if (mode == PTRAJ_SETUP) {

    /*
     *  -------- ANALYZE: PTRAJ_SETUP
     */
    argumentStackPointer = (argStackType **) analyze->carg1;
    analyze->carg1 = NULL;

    /*
     *  Process analysis type
     */ 
    if(argumentStackContains(argumentStackPointer, "fluct")){
      analyze->iarg1 = ANALYZEMODES_FLUCT;
    }
    else if(argumentStackContains(argumentStackPointer, "displ")){
      analyze->iarg1 = ANALYZEMODES_DISPL;
    }
    else if(argumentStackContains(argumentStackPointer, "corr")){
      analyze->iarg1 = ANALYZEMODES_CORR;
    }
    else{
      analyze->iarg1 = ANALYZEMODES_UNKNOWN;
      fprintf(stderr, "ptraj(), analyzeModes: no analysis type given, returning\n");
      freeAnalyzeModesMemory(analyze);
      return -1;
    }

    /*
     *  Get beg, end, factor, bose
     */
    analyze->iarg2 = ibeg = argumentStackKeyToInteger(argumentStackPointer, "beg", 7);
    analyze->iarg3 = iend = argumentStackKeyToInteger(argumentStackPointer, "end", 50);
    analyze->iarg4 = 0;
    if( argumentStackContains( argumentStackPointer, "bose" ) )
      analyze->iarg4 = 1;
    analyze->darg1 = (double) argumentStackKeyToFloat(argumentStackPointer, "factor", 1.0);

    /*
     *  Get modes 
     */
    if((buffer = argumentStackKeyToString(argumentStackPointer, "file", NULL)) != NULL){
      /*
       * Set flag that modinfo is freed in freeAnalyzeModesMemory
       */
      analyze->darg2 = 1.0;

      /*
       *  Allocate modesInfo structure
       */
      modinfo = (modesInfo *) safe_malloc(sizeof(modesInfo));
      analyze->carg1 = (void *) modinfo;
      INITIALIZE_modesInfo(modinfo);
      modinfo->name = buffer;
      modinfo->type = MT_UNKNOWN;
      modinfo->source = MS_FILE;

      /*
       *  Read evec file
       */
      fp = safe_fopen(buffer, "r");
      if(fp == NULL){
        fprintf(stderr,
                "WARNING in ptraj(), analyzeModes: file %s not opened, ignoring command\n", buffer);
        freeAnalyzeModesMemory(analyze);
        return -1;
      }
      if(readEvecFile(fp, analyze->iarg2, analyze->iarg3, modinfo)){
        fprintf(stderr,
                "WARNING in ptraj(), analyze modes: error while reading %s, ignoring command\n", buffer);
        freeAnalyzeModesMemory(analyze);
        return -1;
      }
      if(modinfo->nvect != (iend - ibeg + 1)){
        fprintf(stderr,
                "FYI: Number of read evecs is %i, number of requested evecs is %i\n", 
                modinfo->nvect, iend - ibeg + 1);
      }
      safe_fclose(fp);
    }        
    else if((buffer = argumentStackKeyToString(argumentStackPointer, "stack", NULL)) != NULL){
      /*
       * Set flag that modinfo is not freed in freeAnalyzeModesMemory
       */
      analyze->darg2 = -1.0;

      /*
       * Get modes from stack 
       */
      modinfo = modesInfoStackGetName(&modesStack, buffer); 
      if (modinfo == NULL) {
        /*
         *  try searching for a match without a leading $ just in case the name
         *  was padded
         */
        modinfo = modesInfoStackGetName(&modesStack, buffer+1); 
      }
      if (modinfo == NULL) {
        fprintf(stderr, "ptraj(), analyzeModes: cannot find a match in the\n");
        fprintf(stderr, "modesInfoStack for name (%s), returning.\n", buffer);
        freeAnalyzeModesMemory(analyze);
        return -1;
      }

      analyze->carg1 = (void *) modinfo;
      safe_free(buffer);
    }
    else{
      analyze->iarg1 = ANALYZEMODES_UNKNOWN;
      fprintf(stderr, "ptraj(), analyzeModes: no stack/file information given, returning\n");
      freeAnalyzeModesMemory(analyze);
      return -1;
    }

    if(modinfo->type != MT_COVAR && modinfo->type != MT_MWCOVAR){
      fprintf(stderr,
              "WARNING in ptraj(), analyzeModes: evecs not of type COVAR or MWCOVAR, ignoring command\n");
      freeAnalyzeModesMemory(analyze);
      return -1;
    }

    /*
     *  Get outfile name
     */
    buffer = argumentStackKeyToString(argumentStackPointer, "out", NULL);
    analyze->carg2 = (void *) buffer;

    /*
     *  Get mask pair info for ANALYZEMODES_CORR option and build the atom pair stack
     */
    analyze->carg4 = atompairStack = NULL;
    state = analyze->state; //ptrajCopyState(ptrajCurrentState());
    npair = 0;
    if(analyze->iarg1 == ANALYZEMODES_CORR){
      while(argumentStringContains(argumentStackPointer, "maskp")){
	/*
	 *  Next two arguments should be one-atom masks
	 */
	buffer = getArgumentString(argumentStackPointer, NULL);
	mask1 = NULL;
	if(buffer != NULL){
	  mask1 = processAtomMask(buffer, state);
	  safe_free(buffer);
	}
	buffer = getArgumentString(argumentStackPointer, NULL);
	mask2 = NULL;
	if(buffer != NULL){
	  mask2 = processAtomMask(buffer, state);
	  safe_free(buffer);
	}

	if(mask1 != NULL && mask2 != NULL){
	  /*
	   *  Check to see if each mask only represents a single atom or not
	   */
	  nm1 = nm2 = 0;
	  for(i=0; i < state->atoms; i++){
	    if(mask1[i] == 1){
	      nm1++;
	      at1 = i+1;
	    }
	    if(mask2[i] == 1){
	      nm2++;
	      at2 = i+1;
	    }
	  }
	  safe_free(mask1);
	  safe_free(mask2);
	
	  if(nm1 == 1 && nm2 == 1){
	    /*
	     *  Store atom pair
	     */
	    mp = safe_malloc(sizeof(int) * 2);
	    mp[0] = at1;
	    mp[1] = at2;
	    pushBottomStack(&atompairStack, (void *) mp);
	    npair++;
	  }
	  else{
	    fprintf(stderr,
		    "WARNING in ptraj(), analyzeModes: masks should only contain one atom, ignoring command\n");
	    freeAnalyzeModesMemory(analyze);
	    return -1;
	  }
	}
	else{
	  fprintf(stderr,
		  "WARNING in ptraj(), analyzeModes: mask pair not valid, ignoring command\n");
	  if(mask1 != NULL)
	    safe_free(mask1);
	  if(mask2 != NULL)
	    safe_free(mask2);
	  freeAnalyzeModesMemory(analyze);
	  return -1;
	}
      }

      if(npair == 0){
	fprintf(stderr,
		"WARNING in ptraj(), analyzeModes: no atom pairs found, ignoring command\n");
	freeAnalyzeModesMemory(analyze);
	return -1;
      }
    }
    analyze->carg4 = (void *) atompairStack;

    /*
     *  Allocate memory for results vector
     */
    if(analyze->iarg1 == ANALYZEMODES_FLUCT){
      analyze->carg3 = (void *) safe_malloc(sizeof(double) * modinfo->navgelem * 4 / 3);
    }
    else if(analyze->iarg1 == ANALYZEMODES_DISPL){
      analyze->carg3 = (void *) safe_malloc(sizeof(double) * modinfo->navgelem);
    }
    else if(analyze->iarg1 == ANALYZEMODES_CORR){
      analyze->carg3 = (void *) safe_malloc(sizeof(double) * npair * (iend - ibeg + 1));
    }

    return 0;

  } else if (mode == PTRAJ_STATUS) {

    /*
     *  -------- ANALYZE: PTRAJ_STATUS
     */

    fprintf(stdout, "  ANALYZE MODES: Calculating %s using modes %i to %i from %s\n",
	    (analyze->iarg1 == ANALYZEMODES_FLUCT ? "rms fluctuations" : 
	     (analyze->iarg1 == ANALYZEMODES_DISPL ? "displacements" :
	      (analyze->iarg1 == ANALYZEMODES_CORR ? "correlation functions" : "????"
	       )
	      )
	     ),
            analyze->iarg2,
            analyze->iarg3,
	    ((modesInfo *) analyze->carg1)->name);
    fprintf(stdout, "      Results are written to %s\n", analyze->carg2 != NULL ? 
                                                                  (char *) analyze->carg2 : "STDOUT");
    fprintf(stdout, "      Statistic used: %s\n", analyze->iarg4 ? "Bose" : "Boltzmann");
    if(analyze->iarg1 == ANALYZEMODES_DISPL)
      fprintf(stdout, "      Factor for displacement: %f\n", analyze->darg1);
    if(analyze->iarg1 == ANALYZEMODES_CORR){
      fprintf(stdout, "      Using the following atom pairs: ");
      for(apsp = ((stackType*) analyze->carg4); apsp != NULL; apsp = apsp->next)
	fprintf(stdout, "(%i,%i)  ", ((int *) apsp->entry)[0], ((int *) apsp->entry)[1]);
      fprintf(stdout, "\n");
    }
  } else if (mode == PTRAJ_CLEANUP) {
    /*
     *  -------- ANALYZE: PTRAJ_CLEANUP
     */
    freeAnalyzeModesMemory(analyze);
  }

  /*
   *  Other possible modes:
   *
   *  PTRAJ_PRINT   -- dump information to output files
   *  PTRAJ_NOOP    -- NULL mode
   */

  if (mode != PTRAJ_ACTION) return 0;

  /*
   *  -------- ANALYZE: PTRAJ_ACTION
   */

  type    = (analyzeModesType) analyze->iarg1;
  ibeg    = analyze->iarg2;
  iend    = analyze->iarg3;
  ibose   = analyze->iarg4;
  factor  = analyze->darg1;
  modinfo = (modesInfo *) analyze->carg1;
  outfile = (char *) analyze->carg2;
  results = (double *) analyze->carg3;
  atompairStack = (stackType *) analyze->carg4;

  source    = modinfo->source;
  ncoords   = modinfo->navgelem;
  nvect     = modinfo->nvect;
  nvectelem = modinfo->nvectelem;
  freq      = modinfo->freq;
  evec      = modinfo->evec;
  avg       = modinfo->avg;

  natoms    = ncoords / 3;

  if(type == ANALYZEMODES_FLUCT){
    /*
     * Calc rms atomic fluctuations
     */
    for(i = 0; i < natoms; i++){
      sumx = sumy = sumz = 0.0;
      for(j = 0; j < nvect; j++){
        if(source == MS_FILE ||
           source == MS_STACK && j+1 >= ibeg && j+1 <= iend){
          f = freq[j];
          if(f >= 0.5){
            /* 
             * Don't use eigenvectors associated with 
             * zero or negative eigenvalues 
             */
            ind = nvectelem*j + 3*i;
            distx = evec[ind  ] * evec[ind  ];
            disty = evec[ind+1] * evec[ind+1];
            distz = evec[ind+2] * evec[ind+2];
            fi = 1.0 / (f*f);
            if(ibose){
              argq = CONSQ * f;
              fi = fi * argq / tanh(argq);
            }
            sumx += distx * fi;
            sumy += disty * fi;
            sumz += distz * fi;
          }
        }
      }

      sumx *= CNST;
      sumy *= CNST;
      sumz *= CNST;

      ind = 4*i;
      results[ind  ] = sqrt(sumx) * CONT;
      results[ind+1] = sqrt(sumy) * CONT;
      results[ind+2] = sqrt(sumz) * CONT;
      results[ind+3] = sqrt(sumx + sumy + sumz) * CONT;
    }
  }
  else if(type == ANALYZEMODES_DISPL){
    /*
     * Calc displacement of coordinates along normal mode directions
     */
    sqrtcnst = sqrt(CNST);

    for(i = 0; i < nvect; i++){
      if(source == MS_FILE ||
         source == MS_STACK && i+1 >= ibeg && i+1 <= iend){
        f = freq[i];
        if(f >= 0.5){
          /* 
           * Don't use eigenvectors associated with 
           * zero or negative eigenvalues 
           */
          fi = 1.0 / f;
          if(ibose){
            argq = CONSQ * f;
            fi = fi * fi * argq / tanh(argq);
            fi = sqrt(fi);
          }
          fi *= sqrtcnst * CONT * factor;

          for(j = 0; j < natoms; j++){
            if(i == 0){
              ind = 3 * j;
              results[ind] = evec[nvectelem * i + ind] * fi;
              ind = 3 * j + 1;
              results[ind] = evec[nvectelem * i + ind] * fi;
              ind = 3 * j + 2;
              results[ind] = evec[nvectelem * i + ind] * fi;
            }
            else{
              ind = 3 * j;
              results[ind] += evec[nvectelem * i + ind] * fi;
              ind = 3 * j + 1;
              results[ind] += evec[nvectelem * i + ind] * fi;
              ind = 3 * j + 2;
              results[ind] += evec[nvectelem * i + ind] * fi;
            }
          }
        }
      }      
    }
  }
  else if(type == ANALYZEMODES_CORR){
    /*
     *  Calc dipole-dipole correlation functions
     */

    ncnt = 0;
    for(apsp = atompairStack; apsp != NULL; apsp = apsp->next){
      at1 = ((int *) apsp->entry)[0];
      at2 = ((int *) apsp->entry)[1];

      /*
       *  Calc unit vector along at2->at1 bond
       */
      dnorm = 0.0;
      for(i = 0; i < 3; i++){
	e[i] = avg[3 * (at1 - 1) + i] - avg[3 * (at2 - 1) + i];
	dnorm += e[i] * e[i];
      }
      dnorm = sqrt(dnorm);
      for(i = 0; i < 3; i++){
	e[i] /= dnorm;
      }

      /*
       *  Loop over desired modes
       */
      for(i = 0; i < nvect; i++){
	if(source == MS_FILE ||
	   source == MS_STACK && i+1 >= ibeg && i+1 <= iend){
	  if(freq[i] >= 0.5){
	    /* 
	     * Don't use eigenvectors associated with 
	     * zero or negative eigenvalues 
	     */
	    f = freq[i]*freq[i] / 11791.79;
	    if(ibose){
	      argq = CONSQ * freq[i];
	      qcorr = argq / tanh(argq);
	    }
	    else{
	      qcorr = 1.0;
	    }

	    /*
	     *  Calc the correlation matrix for delta
	     *    as in eq. 7.16 of lamm and szabo, J Chem Phys 1986, 85, 7334.
	     *  Note that the rhs of this eq. should be multiplied by kT
	     */
	    for(j = 0; j < 3; j++){
	      indi1 = 3 * (at1 - 1) + j;
	      indi2 = 3 * (at2 - 1) + j;
	      for(k = 0; k < 3; k++){
		indj1 = 3 * (at1 - 1) + k;
		indj2 = 3 * (at2 - 1) + k;
		del[j][k] = 0.6 * (qcorr / f) * (evec[nvectelem * i + indi1] - evec[nvectelem * i + indi2]) * 
		                                (evec[nvectelem * i + indj1] - evec[nvectelem * i + indj2]);
	      }
	    }

	    /*
	     *  Correlation in length, eq. 10.2 of lamm and szabo
	     */
	    /*****
	    rtr0 = 0.0;
	    for(j = 0; j < 3; j++)
	      for(k = 0; k < 3; k++)
		rtr0 += e[j] * e[k] * del[j][k];
	    *****/
	     
	    /*
	     *  Librational correlation function, using eq. 7.12 of lamm and szabo
	     *    (w/o beta on the lhs)
	     */
	    val = 0.0;
	    for(j = 0; j < 3; j++){
	      val -= del[j][j];
	      for(k = 0; k < 3; k++){
		val += e[j] * e[k] * del[j][k];
	      }
	    }
	    val *= (3.0 / (dnorm * dnorm));
	    
            results[ncnt] = val;
	    ncnt++;
          }
        }
      }      
    }
  }

  /*
   * Output of results
   */
  if(outfile != NULL){
    fp = safe_fopen(outfile, "w");
    if(fp == NULL){
      fprintf(stderr, "WARNING in ptraj(), analyzeModes: error on opening %s for output\n",
	      outfile);
      return 0;
    }
  }
  else{
    fp = stdout;
  }

  fprintf(fp,"Analysis of modes: %s\n", (type == ANALYZEMODES_FLUCT ? "RMS FLUCTUATIONS" :
                                         (type == ANALYZEMODES_DISPL ? "DISPLACEMENT" :
					  (type == ANALYZEMODES_CORR ? "CORRELATION FUNCTIONS" : "UNKNOWN"
					   )
					  )
					 ));  

  if(type == ANALYZEMODES_FLUCT){
    fprintf(fp, "%10s %10s %10s %10s %10s\n", "Atom no.", "rmsX", "rmsY", "rmsZ", "rms");
    for(i = 0; i < natoms; i++)
      fprintf(fp, "%10i %10.3f %10.3f %10.3f %10.3f\n", i+1, 
                                                        results[4*i  ],
                                                        results[4*i+1],
                                                        results[4*i+2],
                                                        results[4*i+3]);
  }
  else if(type == ANALYZEMODES_DISPL){
    fprintf(fp, "%10s %10s %10s %10s\n", "Atom no.", "displX", "displY", "displZ");
    for(i = 0; i < natoms; i++)
      fprintf(fp, "%10i %10.3f %10.3f %10.3f\n", i+1, 
                                                 results[3*i  ],
                                                 results[3*i+1],
                                                 results[3*i+2]);
  }
  else if(type == ANALYZEMODES_CORR){
    fprintf(fp, "%10s %10s %10s %10s %10s %10s\n", "Atom1", "Atom2", "Mode", "Freq", "1-S^2", "P2(cum)");
    ncnt = 0;
    for(apsp = atompairStack; apsp != NULL; apsp = apsp->next){
      at1 = ((int *) apsp->entry)[0];
      at2 = ((int *) apsp->entry)[1];
      fprintf(fp, "%10i %10i\n", at1, at2);

      val = 1.0;
      for(i = 0; i < nvect; i++){
	if(source == MS_FILE ||
	   source == MS_STACK && i+1 >= ibeg && i+1 <= iend){
	  if(freq[i] >= 0.5){
	    val += results[ncnt];
	    fprintf(fp, "%10s %10s %10i %10.5f %10.5f %10.5f\n", 
		        "",   "",  i,   freq[i], results[ncnt], val);
	    ncnt++;
	  }
	}
      }
    }
  }

  if(outfile != NULL)
    safe_fclose(fp);

  return 1;
}

/** ANALYZE ROUTINE *************************************************************
 *
 *  CorrelationCoefficient() 
 *
 *  Supplementary routines:
 *
 ******************************************************************************/

typedef struct _correlationCoefficientResults {
  int N;
  double average;
  double a2;
  double stddev;
  double coeff;
  double a;
  double b;
  double significance;
} correlationCoefficientResults;

#define INITIALIZE_correlationCoefficientResults(_p_) \
  _p_->N = 0; \
  _p_->average = 0.0; \
  _p_->a2 = 0.0; \
  _p_->stddev = 0.0; \
  _p_->coeff = 0.0; \
  _p_->a = 0.0; \
  _p_->b = 0.0; \
  _p_->significance = 0.0

   int
analyzeCorrelationCoefficient(analyzeInformation *analyze, stackType *sp, int mode)
{
  argStackType **argumentStackPointer;
  char *buffer;

  /*
   *  usage:
   *
   *    analyze correlationcoefficient {name | ALL} {name | ALL}
   *
   *  argument usage:
   *
   *    iarg1:
   *      0 -- a specific scalarInfo entry for first argument
   *      1 -- ALL scalarInfo entries for first argument
   *    iarg2:
   *      0 -- a specific scalarInfo entry for the second argument
   *      1 -- ALL scalarInfo entries for the second argument
   *    carg1: if iarg1 == 0, the scalarInfo entry
   *    carg2: if iarg2 == 0, the scalarInfo entry
   */

  scalarInfo *scalar1, *scalar2;
  correlationCoefficientResults Results1, Results2, *results1, *results2;
  stackType *stack1;
  int i, total;

  results1 = &Results1;
  results2 = &Results2;

  if (mode == PTRAJ_SETUP) {

    /*
     *  PTRAJ_SETUP:
     */

    argumentStackPointer = (argStackType **) analyze->carg1;
    analyze->carg1 = NULL;

       /*
        *  process argument 1
        */
    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL) {
      fprintf(stderr, "ptraj(), analyzeCorrelationCoefficient: missing the first name on\n");
      fprintf(stderr, "scalarStack to process...  Returning.\n");
      return -1;
    }
    if (strcmp(buffer, "ALL") == 0 ||
	strcmp(buffer, "all") == 0) {
      analyze->iarg1 = 1;
    } else {
      scalar1 = scalarStackGetName(&sp, buffer);
      if (scalar1 == NULL) {
	/*
	 *  try searching for a match without a leading $ just in case the name
	 *  was padded
	 */
	scalar1 = scalarStackGetName(&sp, buffer+1);
      }
      if (scalar1 == NULL) {
	fprintf(stderr, "ptraj(), analyzeCorrelationCoefficient: cannot find a match in the\n");
	fprintf(stderr, "scalarStack for name (%s), returning.\n", buffer);
	return -1;
      }
      analyze->carg1 = (void *) scalar1;
    }
    safe_free(buffer);

       /*
        *  process argument 2
        */
    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL) {
      fprintf(stderr, "ptraj(), analyzeCorrelationCoefficient: missing the second name on\n");
      fprintf(stderr, "scalarStack to process...  Returning.\n");
      return -1;
    }
    if (strcmp(buffer, "ALL") == 0 ||
	strcmp(buffer, "all") == 0) {
      analyze->iarg2 = 1;
    } else {
      scalar2 = scalarStackGetName(&sp, buffer);
      if (scalar2 == NULL) {
	/*
	 *  try searching for a match without a leading $ just in case the name
	 *  was padded
	 */
	scalar2 = scalarStackGetName(&sp, buffer+1);
      }
      if (scalar2 == NULL) {
	fprintf(stderr, "ptraj(), analyzeCorrelationCoefficient: cannot find a match in the\n");
	fprintf(stderr, "scalarStack for name (%s), returning.\n", buffer);
	return -1;
      }
      analyze->carg2 = (void *) scalar2;
    }
    safe_free(buffer);


  } else if (mode == PTRAJ_STATUS) {

    /*
     *  PTRAJ_STATUS:
     */

    scalar1 = (scalarInfo *) analyze->carg1;
    scalar2 = (scalarInfo *) analyze->carg2;
    fprintf(stdout, "  ANALYZE CORRELATIONCOEFFICIENT: comparing %s to %s\n",
	    (analyze->iarg1 ? "ALL" : scalar1->name),
	    (analyze->iarg2 ? "ALL" : scalar2->name));
  }

  /*
   *  Other possible modes:
   *
   *  PTRAJ_PRINT   -- dump information to output files
   *  PTRAJ_CLEANUP -- clean up
   *  PTRAJ_NOOP    -- NULL mode
   */

  if (mode != PTRAJ_ACTION) return 0;

  /*
   *  PTRAJ_ACTION:
   */

  scalar1 = (scalarInfo *) analyze->carg1;
  scalar2 = (scalarInfo *) analyze->carg2;

     /*
      *  If ALL is set, count the number of entries on the stack
      */
  total = 0;
  if (analyze->iarg1 == 1 || analyze->iarg2 == 1) {
    for (stack1 = sp; stack1 != NULL; stack1 = stack1->next)
      total++;
  }

     /*
      *  Handle the all-to-all case first
      */
  if (analyze->iarg1 == 1 && analyze->iarg2 == 1) {




  } else if (analyze->iarg1 == 1 || analyze->iarg2 == 1) {
    if (analyze->iarg2 == 1) scalar1 = scalar2;
       /*
        *  If ALL is specified (once): no sorting is done
        */

    INITIALIZE_correlationCoefficientResults(results1);
    INITIALIZE_correlationCoefficientResults(results2);

    results1->N = scalar1->totalFrames;
    for (i=0; i < results1->N; i++) {
      results1->average += scalar1->value[i];
      results1->a2 += (scalar1->value[i]*scalar1->value[i]);
    }
    results1->average /= results1->N;
    results1->a2 /= results1->N;
    results1->stddev = sqrt(results1->a2 - results1->average*results1->average);
      
    for (stack1 = sp; stack1 != NULL; stack1 = stack1->next) {
      scalar2 = (scalarInfo *) stack1->entry;

      if (scalar1 != scalar2) {
	results2->N = scalar2->totalFrames;
	if (results2->N != results1->N) {
	  fprintf(stderr, "Attempting to calculate a correlation coefficient for two data\n");
	  fprintf(stderr, "sets (%s) and (%s) that have differing numbers of values\n",
		  scalar1->name, scalar2->name);
	  fprintf(stderr, "(%i) and (%i); ignoring this calculation...\n", 
		  results1->N, results2->N);
	} else {
	
	  for (i=0; i < results2->N; i++) {
	    results2->average += scalar2->value[i];
	    results2->a2 += (scalar2->value[i]*scalar2->value[i]);
	  }
	  results2->average /= results2->N;
	  results2->a2 /= results2->N;
	  results2->stddev = sqrt(results2->a2 - results2->average*results2->average);

	  for (i=0; i < results1->N; i++)
	    results1->coeff += scalar1->value[i]*scalar2->value[i];
	  results1->coeff /= results1->N;

	  results1->a = results1->coeff - results1->average - results2->average /
	    (results1->a2 - results2->a2);
	  results1->b = results2->average - results1->a * results1->average;

	  results1->coeff = results1->coeff - 
	    (results1->average*results2->average) / (results1->stddev*results2->stddev);

	  if (results1->coeff != 1.0)
	    results1->significance = results1->coeff * 
	      sqrt( (results1->N-2)/(1-results1->coeff*results1->coeff));
	  else
	    results1->significance = 0.0;

	  fprintf(stdout, "  CORRELATION COEFFICIENT %6s to %6s IS %10.4f\n",
		  scalar1->name, scalar2->name, results1->coeff);
	  /*
	  fprintf(stdout, "     A = %8.4f, B = %8.4f, significance = %8.4f\n",
		  results1->a, results1->b, results1->significance);
	  */
	  INITIALIZE_correlationCoefficientResults(results2);
	}
      }
    }
  } else {


    INITIALIZE_correlationCoefficientResults(results1);
    INITIALIZE_correlationCoefficientResults(results2);

    results1->N = scalar1->totalFrames;
    for (i=0; i < results1->N; i++) {
      results1->average += scalar1->value[i];
      results1->a2 += (scalar1->value[i]*scalar1->value[i]);
    }
    results1->average /= results1->N;
    results1->a2 /= results1->N;
    results1->stddev = sqrt(results1->a2 - results1->average*results1->average);
      
    results2->N = scalar2->totalFrames;
    if (results2->N != results1->N) {
      fprintf(stderr, "Attempting to calculate a correlation coefficient for two data\n");
      fprintf(stderr, "sets (%s) and (%s) that have differing numbers of values\n",
	      scalar1->name, scalar2->name);
      fprintf(stderr, "(%i) and (%i); ignoring this calculation...\n", 
	      results1->N, results2->N);
    } else {
	
      for (i=0; i < results2->N; i++) {
	results2->average += scalar2->value[i];
	results2->a2 += (scalar2->value[i]*scalar2->value[i]);
      }
      results2->average /= results2->N;
      results2->a2 /= results2->N;
      results2->stddev = sqrt(results2->a2 - results2->average*results2->average);

      for (i=0; i < results1->N; i++)
	results1->coeff += scalar1->value[i]*scalar2->value[i];
      results1->coeff /= results1->N;

      results1->a = (results1->coeff - results1->average - results2->average) /
	(results1->a2 - results2->a2);
      results1->b = results2->average - results1->a * results1->average;

      results1->coeff = (results1->coeff - results1->average*results2->average) / 
	(results1->stddev*results2->stddev);

      if (results1->coeff != 1.0)
	results1->significance = results1->coeff * 
	  sqrt( (results1->N-2)/(1-results1->coeff*results1->coeff));
      else
	results1->significance = 0.0;

      fprintf(stdout, "  CORRELATION COEFFICIENT %6s to %6s IS %10.4f\n",
	      scalar1->name, scalar2->name, results1->coeff);
      /*
      fprintf(stdout, "     A = %8.4f, B = %8.4f, significance = %8.4f\n",
	      results1->a, results1->b, results1->significance);
      */
      INITIALIZE_correlationCoefficientResults(results2);
    }
  }


  return 1;

}

/* DISABLED FOR CPPTRAJ
   int
analyzeHBond(analyzeInformation *analyze, stackType *sp, int mode)
{
  argStackType **argumentStackPointer;
  char *buffer;

  //
  //  Usage:
  //
  //  analyze hbond <name>
  //          [lifetimes [window <value>]]
  //          [maxoccupied]
  //          [scalar <name>]
  //          [donor { solvent | resname atomname | mask <mask> }]
  //          [acceptor { solvent | res atom res atom | mask <mask> <maskH>}
  //
  //  argument usage:
  //
  //    carg1 -- the transformHBondInfo structure pointer
  //    carg2 -- the scalarStack entry
  //    iarg1 -- if set do lifetime analysis, value is window
  //    iarg2 -- if set do maxoccupied analysis, value is window
  //    iarg3 -- if non-zero, a single solvent donor selection
  //    iarg4 -- if non-zero, a single solvent acceptor selection

  transformHBondInfo *hbondInfo, *hbondInfoSelect;
  scalarInfo *scalar;
  int *maskD, *maskA, *maskAH1, *maskAH2, *maskAH3;
  int i, j, k, idonor, iacceptor, iacceptorH, donors, acceptors;
  int *sortindex, *maxoccupied;
  float *lifetimes, a, d;
  int lifetime, lifetimecounter, index, ia, id;

  if (mode == PTRAJ_SETUP) {

    argumentStackPointer = (argStackType **) analyze->carg1;
    analyze->carg1 = NULL;

        // name
    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL) {
      fprintf(stdout, "WARNING in ptraj(), analyze: hbond name not specified...\n");
      return -1;
    }

    hbondInfo = hbondInfoStackGetName(&hbondStack, buffer);
    if (hbondInfo == NULL) {
      fprintf(stdout, "WARNING in ptraj(), analyze: hbond name %s not found...\n",
	      buffer);
      safe_free(buffer);
      return -1;
    }
    safe_free(buffer);

    analyze->carg1 = (void *) hbondInfo;

       //  lifetimes, maxoccupied, window
    analyze->iarg1 = argumentStackContains(argumentStackPointer, "lifetimes");
    analyze->iarg2 = argumentStackContains(argumentStackPointer, "maxoccupied");
    i = argumentStackKeyToInteger(argumentStackPointer, "window", 0);
    if (i) {
      if (analyze->iarg1) analyze->iarg1 = i;
      if (analyze->iarg2) analyze->iarg2 = i;
    }

       //  scalar <name>
       //
       //  Place values on the scalar stack using "name" as the key.  Note
       //  that this requires selecting an individual h-bond

    scalar = NULL;
    if (argumentStackContains(argumentStackPointer, "scalar")) {

      scalar = (scalarInfo *) safe_malloc(sizeof(scalarInfo));
      INITIALIZE_scalarInfo(scalar);
      buffer = argumentStackKeyToString(argumentStackPointer, "distance", NULL);
      if (buffer != NULL) {
	scalar->mode = SCALAR_DISTANCE;
      } else {
	buffer = argumentStackKeyToString(argumentStackPointer, "angle", NULL);
	if (buffer == NULL) {	
	  fprintf(stdout, "WARNING in ptraj, analyze hbond: scalar command used without\n");
	  fprintf(stdout, "specifying distance or angle\n");
	  safe_free(scalar);
	  scalar = NULL;
	  return -1;
	}
	scalar->mode = SCALAR_ANGLE;
      }

      if ( scalarStackGetName(&scalarStack, buffer) != NULL ) {
	fprintf(stdout, "WARNING: ptraj(), analyze hbond: The chosen name (%s) has already\n",
		buffer);
	fprintf(stdout, "been used.  Ignoring scalar command\n");
	safe_free(scalar);
	return -1;
      }
      scalar->name = buffer;
      scalar->totalFrames = -1;
    }
    analyze->carg2 = (void *) scalar;

       //  SELECTION

    hbondInfoSelect = NULL;
    maskD = NULL;
    maskA = NULL;
    maskAH1 = NULL;
    maskAH2 = NULL;
    maskAH3 = NULL;
    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL) 
      return 0;
    else if (strcmp(buffer, "donor") != 0) {
      fprintf(stdout, "WARNING in ptraj, analyze hbond: Error on donor specification in\n");
      fprintf(stdout, "select specification...  Ignoring command...\n");
      return -1;
    }
    safe_free(buffer);

    buffer = (char *) popStack(argumentStackPointer);
    if ( strstr(buffer, "solvent") != NULL) {

      if (strcmp(buffer, "solvent") == 0) {
	analyze->iarg3 = -1;
      } else if (strncmp(buffer, "solvent", 7) == 0) {
	sscanf(buffer+7, "%i", &analyze->iarg3);
	if (analyze->iarg3 > hbondInfo->solventNeighbor) {
	  fprintf(stdout, "WARNING in ptraj, analyze hbond: Number solvent selection out of\n");
	  fprintf(stdout, "range (%s), max is %i.  Ignoring...\n", 
		  buffer, hbondInfo->solventNeighbor);
	  return -1;
	}
      }
      safe_free(buffer);
    } else {
      pushStack(argumentStackPointer, buffer);

      maskD = (int *) safe_malloc(sizeof(int) * hbondInfo->state->atoms);
      for (i=0; i < hbondInfo->state->atoms; i++)
	maskD[i] = 0;
      parseHBondDonor(argumentStackPointer, hbondInfo->state, maskD);
    }

    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL)
      return 0;
    else if (strcmp(buffer, "acceptor") != 0) {
      fprintf(stdout, "WARNING in ptraj, analyze hbond: Error on acceptor specification in\n");
      fprintf(stdout, "select specification...  Ignoring command...\n");
      return -1;
    }
    safe_free(buffer);

    buffer = (char *) popStack(argumentStackPointer);
    if ( strstr(buffer, "solvent") != NULL) {

      if (strcmp(buffer, "solvent") == 0) {
	analyze->iarg4 = -1;
      } else if (strncmp(buffer, "solvent", 7) == 0) {
	sscanf(buffer+7, "%i", &analyze->iarg4);
	if (analyze->iarg4 > hbondInfo->solventNeighbor) {
	  fprintf(stdout, "WARNING in ptraj, analyze hbond: Number solvent selection out of\n");
	  fprintf(stdout, "range (%s), max is %i.  Ignoring...\n",
		  buffer, hbondInfo->solventNeighbor);
	  return -1;
	}
      }
      safe_free(buffer);
    } else {
      pushStack(argumentStackPointer, buffer);
      
      maskA   = (int *) safe_malloc(sizeof(int) * hbondInfo->state->atoms);
      maskAH1 = (int *) safe_malloc(sizeof(int) * hbondInfo->state->atoms);
      maskAH2 = (int *) safe_malloc(sizeof(int) * hbondInfo->state->atoms);
      maskAH3 = (int *) safe_malloc(sizeof(int) * hbondInfo->state->atoms);
      for (i=0; i < hbondInfo->state->atoms; i++) {
	maskA[i] = 0;
	maskAH1[i] = -1;
	maskAH2[i] = -1;
	maskAH3[i] = -1;
      }

      parseHBondAcceptor(argumentStackPointer, hbondInfo->state,
			 maskA, maskAH1, maskAH2, maskAH3);
    }

    if (maskD || maskA) {
      hbondInfoSelect = (transformHBondInfo *) safe_malloc(sizeof(transformHBondInfo));
      INITIALIZE_transformHBondInfo(hbondInfoSelect);

      if (maskD)
	atomMaskIsActive(maskD, hbondInfo->state, &hbondInfoSelect->numdonor, &i);
      if (maskA) {
	k = 0;
	for (i = 0; i < hbondInfo->state->atoms; i++) {
	  if (maskA[i] == 1) {
	    k++;
	    if (maskAH2[i] >= 0) k++;
	    if (maskAH3[i] >= 0) k++;
	  }
	}
	hbondInfoSelect->numacceptor = k;
      }

      if (maskD) hbondInfoSelect->donor     = (int *) 
		   safe_malloc(sizeof(int) * hbondInfoSelect->numdonor);
      if (maskA) {
	hbondInfoSelect->acceptor  = (int *) 
	  safe_malloc(sizeof(int) * hbondInfoSelect->numacceptor);
	hbondInfoSelect->acceptorH = (int *) 
	  safe_malloc(sizeof(int) * hbondInfoSelect->numacceptor);
      }
      
      idonor = 0;
      iacceptor = 0;
      for (i=0; i < hbondInfo->state->atoms; i++) {

	if (maskD && maskD[i] == 1)
	  hbondInfoSelect->donor[idonor++] = i;
	if (maskA && maskA[i] == 1) {
	  hbondInfoSelect->acceptor[iacceptor] = i;
	  hbondInfoSelect->acceptorH[iacceptor] = maskAH1[i];
	  iacceptor++;
	  if (maskAH2[i] >= 0) {
	    hbondInfoSelect->acceptor[iacceptor] = i;
	    hbondInfoSelect->acceptorH[iacceptor] = maskAH2[i];
	    iacceptor++;
	  }
	  if (maskAH3[i] >= 0) {
	    hbondInfoSelect->acceptor[iacceptor] = i;
	    hbondInfoSelect->acceptorH[iacceptor] = maskAH3[i];
	    iacceptor++;
	  }
	}
      }
    }


    
    analyze->carg3 = (void *) hbondInfoSelect;

    if (scalar != NULL) {
      if (!( (analyze->iarg3 != 0 || (maskD && hbondInfoSelect->numdonor == 1)) &&
	     (analyze->iarg4 != 0 || (maskA && hbondInfoSelect->numacceptor == 1)) ) ) {
	fprintf(stderr, "WARNING in ptraj, analyze hbond: scalar option.  More than a single\n");
	fprintf(stderr, "hbond selected...  This is not supported.\n");
	safe_free(scalar);
	return -1;
      }
    }

    return 0;
  }


  hbondInfo = (transformHBondInfo *) analyze->carg1;
  scalar = (scalarInfo *) analyze->carg2;
  hbondInfoSelect = (transformHBondInfo *) analyze->carg3;

  if (mode == PTRAJ_STATUS) {

    fprintf(stdout, "  ANALYZE HBOND: %s -- ", hbondInfo->name);
    if (analyze->iarg1)
      fprintf(stdout, "lifetimes, ");
    if (analyze->iarg2)
      fprintf(stdout, "max occupied, ");
    if (analyze->iarg1 > 1 || analyze->iarg2 > 1)
      fprintf(stdout, "window %i, ", (analyze->iarg1 > 1 ? analyze->iarg1 : analyze->iarg2));
    fprintf(stdout, "\n");

    if (scalar != NULL) {
      fprintf(stdout, "      Saving hbond %s data to scalar value %s\n",
	      (scalar->mode == SCALAR_DISTANCE ? "distance" : "angle"),
	      scalar->name);
    }

    if ( (analyze->iarg3 != 0 || hbondInfoSelect->numdonor > 0) &&
	 (analyze->iarg4 != 0 || hbondInfoSelect->numacceptor > 0) ) {
      fprintf(stdout, "      hydrogen bond selection:\n");
      if (analyze->iarg3 != 0) {
	if (analyze->iarg3 < 0) 
	  fprintf(stdout, "      solvent donor (avg)  ");
	else
	  fprintf(stdout, "      solvent donor %i  ", analyze->iarg3);
      } else
	fprintf(stdout, "      donors: %i  ",
	      hbondInfoSelect->numdonor);

      if (analyze->iarg4 != 0) {
	if (analyze->iarg4 < 0)
	  fprintf(stdout, "      solvent acceptor (avg)\n");
	else
	  fprintf(stdout, "      solvent acceptor %i\n", analyze->iarg4);
      } else
	fprintf(stdout, "acceptors %i\n", hbondInfoSelect->numacceptor);

      if (prnlev > 2 && hbondInfoSelect != NULL)
	printHBondInfo(hbondInfoSelect->numdonor, hbondInfoSelect->donor, 
		       hbondInfoSelect->numacceptor, hbondInfoSelect->acceptor, 
		       hbondInfoSelect->acceptorH, hbondInfo->state);

    }
  }

  //  Other possible modes:
  //
  //  PTRAJ_PRINT   -- dump information to output files
  //  PTRAJ_CLEANUP -- free any allocated memory
  //  PTRAJ_NOOP    -- NULL mode

  if (mode != PTRAJ_ACTION) return 0;

  //  Perform action on coordinates, etc.


  donors = hbondInfo->numdonor;
  if (hbondInfo->numSolventAcceptor)
    donors += hbondInfo->solventNeighbor;

  acceptors = hbondInfo->numacceptor;
  if (hbondInfo->numSolventDonor)
    acceptors += hbondInfo->solventNeighbor;

  //  Don't do hbond selection for now...
  if (analyze->iarg1) {
    lifetimes = (float *) safe_malloc(sizeof(float) * donors * acceptors);
    for (i=0; i < donors*acceptors; i++) 
      lifetimes[i] = 0.0;
  }

  if (analyze->iarg2) {
    maxoccupied = (int *) safe_malloc(sizeof(int) * donors * acceptors);
    for (i=0; i < donors*acceptors; i++) 
      maxoccupied[i] = 0;
  }


  if (analyze->iarg1) {
    sortindex = (int *) safe_malloc(sizeof(int) * donors * acceptors);
    for (k = 0; k < donors*acceptors; k++)
      sortindex[k] = k;

    for (i=0; i < hbondInfo->numdonor; i++) {
      idonor = hbondInfo->donor[i];
      for (j = 0; j < hbondInfo->numacceptor; j++) {

	index = i*acceptors + j;
	iacceptor = hbondInfo->acceptor[j];
	iacceptorH = hbondInfo->acceptorH[j];

	lifetime++;
	lifetimecounter = 0;
	for (k = 0; k < hbondInfo->visit; k++) {

	
	  d = hbondInfo->seriesDistance[k*donors*acceptors + index];
	  a = 180.0 - hbondInfo->seriesAngle[k*donors*acceptors + index];

	  if (prnlev > 4 ) {
	    printf("series %i, index %i, distance is %.3f, angle is %.3f\n",
		   k, index, d, a);
	  }

	 // if (hbondInfo->seriesDistance[k*donors*acceptors+index] < hbondInfo->distanceCutoff && 
	 //     (180.0-hbondInfo->seriesAngle[k*donors*acceptors+index])>hbondInfo->angleCutoff) {
	  if (d < hbondInfo->distanceCutoff &&
	      a > hbondInfo->angleCutoff) {
	    lifetimecounter++;
	  } else {
	    lifetimes[index] += lifetimecounter * hbondInfo->timeinterval;
	    lifetimecounter = 0;
	    lifetime++;
	  }
	}
	if (lifetimecounter == hbondInfo->visit) {
	  lifetimes[index] = lifetimecounter * hbondInfo->timeinterval;
	} else if (lifetime > 0) {
	  lifetimes[index] /= lifetime;
	}
      }
    }

    sortIndexFloat(lifetimes, sortindex, acceptors*donors);

    fprintf(stdout, "  HBOND LIFETIMES SUMMARY:\n");
    fprintf(stdout, "        DONOR                   ACCEPTORH            ACCEPTOR\n");
    fprintf(stdout, "atom# res#  atom res -- atom# res#  atom res atom# res#  atom res");
    fprintf(stdout, " lifetimes\n");
    for (k = acceptors*donors-1; k >= 0; k--) {
      j = sortindex[k];
      if (lifetimes[j] > 0.0) {

	id = j / acceptors;
	ia = j % acceptors;

	idonor     = hbondInfo->donor[     id ];
	iacceptor  = hbondInfo->acceptor[  ia ];
	iacceptorH = hbondInfo->acceptorH[ ia ];

	if ( id <= hbondInfo->numdonor && ia <= hbondInfo->numacceptor &&
	     ! (id == hbondInfo->numdonor && ia == hbondInfo->numacceptor) ) {
	  if (id == hbondInfo->numdonor) 
	    fprintf(stdout, "       solvent donor ");
	  else {
	    fprintf(stdout, "%5i ", idonor+1);
	    printAtom(stdout, idonor, hbondInfo->state);
	  }
	  
	  if (ia == hbondInfo->numacceptor) 
	    fprintf(stdout, "--  solvent acceptor                         ");
	  else {
	    fprintf(stdout, "-- %5i ", iacceptorH+1);
	    printAtom(stdout, iacceptorH, hbondInfo->state);
	    fprintf(stdout, "%5i ", iacceptor+1);
	    printAtom(stdout, iacceptor, hbondInfo->state);
	  }
	  
	  fprintf(stdout," %6.3f\n", lifetimes[j]);
	}
      }
    }

    safe_free(sortindex);
  }

  return 1;


}
*/






   int
analyzeSet(analyzeInformation *analyze, stackType *sp, int mode)
{
  argStackType **argumentStackPointer;

  if (mode == PTRAJ_SETUP) {

    /*
     *  PTRAJ_SETUP:
     */

    argumentStackPointer = (argStackType **) analyze->carg1;
    analyze->carg1 = NULL;

    //if (prnlev > 4) {
    //  printStack(argumentStackPointer, printString, NULL);
    //}

  } else if (mode == PTRAJ_STATUS) {

    /*
     *  PTRAJ_STATUS:
     */

  }

  /*
   *  Other possible modes:
   *
   *  PTRAJ_PRINT   -- dump information to output files
   *  PTRAJ_CLEANUP -- free any allocated memory
   *  PTRAJ_NOOP    -- NULL mode
   */

  if (mode != PTRAJ_ACTION) return 0;

  /*
   *  PTRAJ_ACTION:
   */

  return 1;


}



   int
analyzeStatistics(analyzeInformation *analyze, stackType *sp, int mode)
{
  /*
   *  usage:
   *
   *  analyze statistics {name | ALL} [shift value]
   *
   *  argument usage:
   *
   *    iarg1: 1 if ALL
   *    darg1: shift value
   *    carg1: the scalarInfo * if iarg1 != 1
   */

  argStackType **argumentStackPointer;
  char *buffer;
  scalarInfo *scalar;
  stackType *stack;
  double average, stddev, value;
  int i, j, k, l, periodic, curbin, prevbin;

  int pucker_visits[10];
  int pucker_transitions[10][10];
  double pucker_avg[10];
  double pucker_sd[10];

  int torsion_visits[6];
  int torsion_transitions[6][6];
  double torsion_avg[10];
  double torsion_sd[10];

  int distance_visits[6];
  int distance_transitions[6][6];
  double distance_avg[6];
  double distance_sd[6];

  if (mode == PTRAJ_SETUP) {

    /*
     *  PTRAJ_SETUP:
     */

    argumentStackPointer = (argStackType **) analyze->carg1;
    analyze->carg1 = NULL;

    if (argumentStackContains(argumentStackPointer, "all"))
      analyze->iarg1 = 1;
    else
      analyze->iarg1 = 0;
    analyze->darg1 = argumentStackKeyToDouble(argumentStackPointer, "shift", 0.0);


    if (analyze->iarg1 != 1) {
      buffer = getArgumentString(argumentStackPointer, NULL);
      if (buffer == NULL) {
	fprintf(stdout, "ptraj(), analyzeStatistics: No name specified\n");
	return -1;
      }
      scalar = scalarStackGetName(&sp, buffer);
      if (scalar == NULL) {
	fprintf(stdout, "ptraj(), analyzeStatistics: Name (%s) not found, ignoring\n",
		buffer);
	return -1;
      }
      analyze->carg1 = (void *) scalar;
    }
    return 0;
  }


  scalar = (scalarInfo *) analyze->carg1;


  if (mode == PTRAJ_STATUS) {

    /*
     *  PTRAJ_STATUS:
     */
    fprintf(stdout, "  ANALYZE STATISTICS: ");
    if (analyze->iarg1 == 1)
      fprintf(stdout, "ALL accumulated values ");
    else
      fprintf(stdout, "name %s ", scalar->name);
    if (analyze->darg1 != 0.0) {
      fprintf(stdout, "shift (about %5.2f) is being applied",
	      analyze->darg1);
    }
    fprintf(stdout, "\n");

  }

  /*
   *  Other possible modes:
   *
   *  PTRAJ_PRINT   -- dump information to output files
   *  PTRAJ_CLEANUP -- free any allocated memory
   *  PTRAJ_NOOP    -- NULL mode
   */

  if (mode != PTRAJ_ACTION) return 0;

  /*
   *  PTRAJ_ACTION:
   */

  for (stack = sp; stack != NULL; stack=stack->next) {
    if (analyze->iarg1 == 1) 
      scalar = (scalarInfo *) stack->entry;

    average = 0.0;
    stddev = 0.0;

    periodic = 0;
    if (scalar->mode == SCALAR_ANGLE ||
	scalar->mode == SCALAR_TORSION ||
	scalar->mode == SCALAR_PUCKER)
      periodic = 1;
	
    for (i=0; i < scalar->totalFrames; i++) {
      value = scalar->value[i] - analyze->darg1;
      if (periodic) {
	if (value > 180.0) 
	  value -= 360.0;
	else if (value < -180.0)
	  value += 360.0;
      }
      average += value;
      stddev += value*value;
    }
    average /= scalar->totalFrames;
    stddev /= scalar->totalFrames;
    stddev = stddev - (average*average);
    if (stddev > 0)
      stddev = sqrt(stddev);
    else
      stddev = 0.0;

    average += analyze->darg1;

    fprintf(stdout, "__________________________________________________________________\n\n");
    fprintf(stdout, "STATISTICS %6s\n", scalar->name);
    if (scalar->mode == SCALAR_DISTANCE) {
      if (scalar->mask1 == NULL && scalar->mask2 == NULL) {
	fprintf(stdout, "              DISTANCE between atoms ");
	printAtomCompact2(stdout, scalar->atom1, scalar->state);
	fprintf(stdout, " & ");
	printAtomCompact2(stdout, scalar->atom2, scalar->state);
	fprintf(stdout, "\n");
      } else {
	fprintf(stdout, "              DISTANCE between masks ");
	if (scalar->mask1)
	  printAtomMask(stdout, scalar->mask1, scalar->state);
	else
	  printAtomCompact2(stdout, scalar->atom1, scalar->state);

	fprintf(stdout, " & ");

	if (scalar->mask2)
	  printAtomMask(stdout, scalar->mask2, scalar->state);
	else
	  printAtomCompact2(stdout, scalar->atom2, scalar->state);

	fprintf(stdout, "\n");
      }
    }

    fprintf(stdout, "   AVERAGE: %8.4f (%.4f stddev)\n",
	    average, stddev);
    fprintf(stdout, "   INITIAL: %8.4f\n   FINAL:   %8.4f\n", 
	    scalar->value[0], scalar->value[scalar->totalFrames-1]);


    for (i=0; i < scalar->totalFrames; i++) {

      switch (scalar->mode) {

      case SCALAR_PUCKER:

	if (i == 0) {
	  for (j=0;j<10;j++) {
	    pucker_visits[j] = 0;
	    pucker_avg[j] = 0.0;
	    pucker_sd[j] = 0.0;
	    for (k=0;k<10;k++) {
	      pucker_transitions[j][k] = 0;
	    }
	  }
	}

	value = scalar->value[i];
	if (value < 0) value += 360.0;
	curbin = value / 36;

	pucker_visits[curbin]++;
	pucker_avg[curbin] += scalar->value[i];
	pucker_sd[curbin]  += scalar->value[i]*scalar->value[i];

	if (i > 0) {
	  if (curbin != prevbin) {
	    pucker_transitions[prevbin][curbin]++;
	  }
	}

	prevbin = curbin;

	break;


      case SCALAR_TORSION:

	if (i == 0) {
	  for (j=0;j<6;j++) {
	    torsion_visits[j] = 0;
	    torsion_avg[j] = 0.0;
	    torsion_sd[j] = 0.0;
	    for (k=0;k<6;k++) {
	      torsion_transitions[j][k] = 0;
	    }
	  }
	}

	value = scalar->value[i];
	if (value < 0) value += 360.0;
	curbin = (int) (value - 30.0) / 60;

	torsion_visits[curbin]++;
	value = scalar->value[i] + torsion_offset[curbin];
	/*
	 *  fix for trans averaging
	 */
	if (value < -150.0) value += 360.0;
	torsion_avg[curbin] += value;
	torsion_sd[curbin]  += value*value;

	if (i > 0) {
	  if (curbin != prevbin) {
	    torsion_transitions[prevbin][curbin]++;
	  }
	}

	prevbin = curbin;

	break;


      case SCALAR_DISTANCE:

	if (i == 0) {
	  for (j=0;j<6;j++) {
	    distance_visits[j] = 0;
	    distance_avg[j] = 0.0;
	    distance_sd[j] = 0.0;
	    for (k=0;k<6;k++) {
	      distance_transitions[j][k] = 0;
	    }
	  }
	}

	value = scalar->value[i];
	curbin = value - 1.5;
	if (curbin < 0)
	  curbin = 0;
	else if (curbin > 5)
	  curbin = 5;

	distance_visits[curbin]++;
	distance_avg[curbin] += value;
	distance_sd[curbin]  += value*value;

	if (i > 0) {
	  if (curbin != prevbin) {
	    distance_transitions[prevbin][curbin]++;
	  }
	}

	prevbin = curbin;

	if (scalar->type == SCALAR_TYPE_NOE) {

	  if (i == 0) {
	    fprintf(stdout, "   NOE SERIES: S < 2.9, M < 3.5, w < 5.0, blank otherwise.\n    |");
	    average = 0.0;
	    k = 0;
	    l = 0;
	  }
	  j = scalar->totalFrames / 50.0;
	  if (j < 1) j = 1;

	  if (scalar->value[i] < scalar->bound) k++;
	  if (scalar->value[i] < scalar->boundh) l++;
	  average += scalar->value[i];
	  if (j == 1 || i % j == 1) {
	    average /= j;
	    if (average < 2.9) {
	      printf("S");
	    } else if (average < 3.5) {
	      printf("M");
	    } else if (average < 5.0) {
	      printf("W");
	    } else {
	      printf(" ");
	    }
	    average = 0.0;
	  }
	  if (i+1 == scalar->totalFrames) {
	    printf("|\n");

	    if (scalar->bound > 0.0) {
	      fprintf(stdout, "   NOE < %.2f for %.2f%% of the time\n", 
		      scalar->bound, (double) k / scalar->totalFrames * 100.0);
	    }
	    if (scalar->boundh > 0.0) {
	      fprintf(stdout, "   NOE < %.2f for %.2f%% of the time\n", 
		      scalar->boundh, (double) l / scalar->totalFrames * 100.0);
	    }
	  }
	}


	break;

      } // END switch scalar->mode
    } // END loop over frames



    /*
     *  Try to provide a bit of analysis here based on type
     */
    switch( scalar->type) {

    case SCALAR_TYPE_UNDEFINED:
      break;

    case SCALAR_TYPE_PUCKER:

      fprintf(stdout, "\n   This is marked as a nucleic acid sugar pucker phase\n");
      break;

    case SCALAR_TYPE_HBOND:

      fprintf(stdout, "\n   This is marked as a hydrogen bond\n");
      if ( (torsion_visits[2] + torsion_visits[3] + torsion_visits[4] + torsion_visits[5] ) 
	   > (scalar->totalFrames * 0.1) )
	fprintf(stdout, "  *** hydrogen bond may have broken\n");
      
      break;
    }



    switch( scalar->mode ) {

    case SCALAR_PUCKER:


      fprintf(stdout, "\n            %s %s %s %s %s %s %s %s %s %s\n",
	      pucker_ss[0], pucker_ss[1], pucker_ss[2], pucker_ss[3], pucker_ss[4],
	      pucker_ss[5], pucker_ss[6], pucker_ss[7], pucker_ss[8], pucker_ss[9]);
      fprintf(stdout, "           -------------------------------------");
      fprintf(stdout, "------------------------------------------------------\n");

      for (j=0; j < 10; j++) {
	if (pucker_visits[j] > 0) {
	  pucker_avg[j] /= pucker_visits[j];
	  pucker_sd[j]  /= pucker_visits[j];
	  pucker_sd[j] = sqrt(pucker_sd[j] - pucker_avg[j]*pucker_avg[j]);
	}
      }

      fprintf(stdout, " %%occupied |");
      for (j=0; j < 10; j++) {
	if (pucker_visits[j] > 0) {
	  value = pucker_visits[j]*100.0/scalar->totalFrames;
	  fprintf(stdout," %6.1f |", value);
	} else
	  fprintf(stdout,"        |");
      }
      fprintf(stdout, "\n");
	  
      fprintf(stdout, " average   |");
      for (j=0; j < 10; j++) {
	if (pucker_visits[j] > 0) {
	  fprintf(stdout," %6.1f |", pucker_avg[j]);
	} else
	  fprintf(stdout,"        |");
      }
      fprintf(stdout, "\n");

      fprintf(stdout, " stddev    |");
      for (j=0; j < 10; j++) {
	if (pucker_visits[j] > 1) {
	  fprintf(stdout," %6.1f |", pucker_sd[j]);
	} else
	  fprintf(stdout,"        |");
      }
      fprintf(stdout, "\n           ----------------------------------------------------------");
      fprintf(stdout, "---------------------------------\n");

      if (prnlev > 0) {
      fprintf(stdout, "\nTRANSITIONS TABLE: (from/vertical to/horizontal)\n\n");
      fprintf(stdout, "           %s %s %s %s %s %s %s %s %s %s\n",
	      pucker_ss[0], pucker_ss[1], pucker_ss[2], pucker_ss[3], pucker_ss[4],
	      pucker_ss[5], pucker_ss[6], pucker_ss[7], pucker_ss[8], pucker_ss[9]);
      fprintf(stdout, "           ------------------------------------------");
      fprintf(stdout, "-------------------------------------------------\n");
      for (j=0; j<10; j++) {
	fprintf(stdout, "  %s |", pucker_ss[j]);
	for (k=0; k<10; k++) {
	  if (pucker_transitions[j][k] > 0)
	    fprintf(stdout, " %6i |", pucker_transitions[j][k]);
	  else
	    fprintf(stdout, "        |");
	}
	fprintf(stdout, "\n");
      }
      fprintf(stdout, "           ----------------------------------------------------------");
      fprintf(stdout, "---------------------------------\n\n");
      }
      break;


    case SCALAR_TORSION:

      fprintf(stdout, "\n               %s  %s  %s  %s  %s  %s\n",
	      torsion_ss[0], torsion_ss[1], torsion_ss[2],
	      torsion_ss[3], torsion_ss[4], torsion_ss[5]);
      fprintf(stdout, "           ---------------");
      fprintf(stdout, "----------------------------------------\n");

      for (j=0; j < 6; j++) {
	if (torsion_visits[j] > 0) {
	  torsion_avg[j] /= torsion_visits[j];
	  torsion_sd[j]  /= torsion_visits[j];
	  torsion_sd[j] = sqrt(torsion_sd[j] - torsion_avg[j]*torsion_avg[j]);
	  torsion_avg[j] -= torsion_offset[j];
	}
      }

      fprintf(stdout, " %%occupied |");
      for (j=0; j < 6; j++) {
	if (torsion_visits[j] > 0) {
	  value = torsion_visits[j]*100.0/scalar->totalFrames;
	  fprintf(stdout," %6.1f |", value);
	} else
	  fprintf(stdout,"        |");
      }
      fprintf(stdout, "\n");
	  
      fprintf(stdout, " average   |");
      for (j=0; j < 6; j++) {
	if (torsion_visits[j] > 0) {
	  fprintf(stdout," %6.1f |", torsion_avg[j]);
	} else
	  fprintf(stdout,"        |");
      }
      fprintf(stdout, "\n");

      fprintf(stdout, " stddev    |");
      for (j=0; j < 6; j++) {
	if (torsion_visits[j] > 1) {
	  fprintf(stdout," %6.1f |", torsion_sd[j]);
	} else
	  fprintf(stdout,"        |");
      }
      fprintf(stdout, "\n           --------------------------");
      fprintf(stdout, "-----------------------------\n");

      switch( scalar->type) {
	
      case SCALAR_TYPE_UNDEFINED:
	break;

      case SCALAR_TYPE_ALPHA:
	/*
	 *              "               g+       a+       t        a-       g-       c 
	 */
	fprintf(stdout, " ALPHA       minor             minor            canonical\n");
	fprintf(stdout, "\n   O3'-P-O5'-C5', SNB range is 270-300 deg (g-)\n");
	if ( (torsion_visits[0] + torsion_visits[1] + torsion_visits[2] + torsion_visits[5] ) 
	     > (scalar->totalFrames * 0.1) ) 
	  fprintf(stdout,"   *** > 10%% out of range population detected\n");

	break;

      case SCALAR_TYPE_BETA:

	/*
	 *              "               g+       a+       t        a-       g-       c 
	 */
	fprintf(stdout, " BETA                <-- canonical -->\n");

	fprintf(stdout, "\n   P-O5'-C5'-C4', SNB range is 130-200 deg (a+,t)\n");
	if ( (torsion_visits[0] + torsion_visits[3] + torsion_visits[4] + torsion_visits[5] ) 
	     > (scalar->totalFrames * 0.05) ) 
	  fprintf(stdout,"   *** > 5%% out of range population detected\n");

	break;
      
      case SCALAR_TYPE_GAMMA:

	/*
	 *              "               g+       a+       t        a-       g-       c 
	 */
	fprintf(stdout, " GAMMA     canonical           minor             minor\n");
	fprintf(stdout, "\n   O5'-C5'-C4'-C3', SNB range is 20-80 (g+)\n");
	if (torsion_visits[2] > (scalar->totalFrames* 0.1))
	  fprintf(stdout, "   *** GAMMA trans > 10%% detected!!!\n");

	break;
      

      case SCALAR_TYPE_DELTA:

	/*
	 *              "               g+       a+       t        a-       g-       c 
	 */
	fprintf(stdout, " DELTA      <------ canonical ------>\n");
	fprintf(stdout, "\n   C5'-C4'-C3'-O3', SNB range is 70-180\n");
	fprintf(stdout, "   DNA: ~128 with BI (a+), ~144 with BII (a+)\n");
	if ( (torsion_visits[0] + torsion_visits[3] + torsion_visits[4] + torsion_visits[5] ) 
	     > (scalar->totalFrames * 0.05) ) 
	  fprintf(stdout,"   *** > 5%% out of range population detected\n");

	break;
      
      case SCALAR_TYPE_EPSILON:

	/*
	 *              "               g+       a+       t        a-       g-       c 
	 */
	fprintf(stdout, " EPSILON                         BI       BII\n");
	fprintf(stdout, "\n   C4'-C3'-O3'-P, SNB range is 160-270\n");
	fprintf(stdout, "   BI = %6.2f%% (~184), BII = %6.2f%% (~246)\n",
		(torsion_visits[2]*100.0)/scalar->totalFrames,
		(torsion_visits[3]*100.0)/scalar->totalFrames);
	if ( (torsion_visits[0] + torsion_visits[1] + torsion_visits[4] + torsion_visits[5] ) 
	     > (scalar->totalFrames * 0.05) ) 
	  fprintf(stdout,"   *** > 5%% out of range population detected\n");

	break;
      
      case SCALAR_TYPE_ZETA:

	/*
	 *              "               g+       a+       t        a-       g-       c 
	 */
	fprintf(stdout, " ZETA                <----- BII ------------- BI ----->\n");
	fprintf(stdout, "\n   C3'-O3'-P-O5', SNB range is 130-300\n");
	fprintf(stdout, "   BI = %6.2f%% (~265, a-/g-), BII = %6.2f%% (~174, a+/t)\n",
		(torsion_visits[3]+torsion_visits[4])*100.0/scalar->totalFrames,
		(torsion_visits[1]+torsion_visits[2])*100.0/scalar->totalFrames);
	if ( (torsion_visits[0] + torsion_visits[5] ) 
	     > (scalar->totalFrames * 0.05) ) 
	  fprintf(stdout,"   *** > 5%% out of range population detected\n");

	break;

      case SCALAR_TYPE_CHI:

	/*
	 *              "               g+       a+       t        a-       g-       c 
	 */
	fprintf(stdout, " CHI                         <-------- anti ------->  <--syn---\n");
	fprintf(stdout, "\n   O4'-C1'-NX-CX, SNB range is 200-300\n");

	if ( (torsion_visits[0] + torsion_visits[5] ) 
	     > (scalar->totalFrames * 0.05) ) 
	  fprintf(stdout,"   *** CHI flips; > 5%% out of range populations detected (see table below)\n");
	if ( torsion_visits[1] > (scalar->totalFrames * 0.05) ) 
	  fprintf(stdout,"   *** Unexpected CHI population in a+ region, > 5%%\n");

	break;

      case SCALAR_TYPE_C2P:

	/*
	 *              "               g+       a+       t        a-       g-       c 
	 */
	fprintf(stdout, " C2' to base      in\n");
	fprintf(stdout, "\n   C2'-C1'-NX-CX\n\n");

	break;

      case SCALAR_TYPE_H1P:

	/*
	 *              "               g+       a+       t        a-       g-       c 
	 */
	fprintf(stdout, " H1'       below-plane                           above      in\n");
	fprintf(stdout, "\n   H1'-C1'-NX-CX, > 0 H1' below plane (check if sugar in plane)\n\n");

	break;
      }



      if (prnlev > 0) {
      fprintf(stdout, "\nTRANSITIONS TABLE: (from/vertical to/horizontal)\n\n");
      fprintf(stdout, "              %s  %s  %s  %s  %s  %s\n",
	      torsion_ss[0], torsion_ss[1], torsion_ss[2],
	      torsion_ss[3], torsion_ss[4], torsion_ss[5]);
      fprintf(stdout, "           -----------------------");
      fprintf(stdout, "--------------------------------\n");
      for (j=0; j<6; j++) {
	fprintf(stdout, "   %s |", torsion_ss[j]);
	for (k=0; k<6; k++) {
	  if (torsion_transitions[j][k] > 0)
	    fprintf(stdout, " %6i |", torsion_transitions[j][k]);
	  else
	    fprintf(stdout, "        |");
	}
	fprintf(stdout, "\n");
      }
      fprintf(stdout, "           ------------------");
      fprintf(stdout, "-------------------------------------\n\n");
      }
      break;


    case SCALAR_DISTANCE:

      fprintf(stdout, "\n              %s  %s  %s  %s  %s  %s\n",
	      distance_ss[0], distance_ss[1], distance_ss[2],
	      distance_ss[3], distance_ss[4], distance_ss[5]);
      fprintf(stdout, "           ---------------");
      fprintf(stdout, "----------------------------------------\n");

      for (j=0; j < 6; j++) {
	if (distance_visits[j] > 0) {
	  distance_avg[j] /= distance_visits[j];
	  distance_sd[j]  /= distance_visits[j];
	  distance_sd[j] = sqrt(distance_sd[j] - distance_avg[j]*distance_avg[j]);
	}
      }

      fprintf(stdout, " %%occupied |");
      for (j=0; j < 6; j++) {
	if (distance_visits[j] > 0) {
	  value = distance_visits[j]*100.0/scalar->totalFrames;
	  fprintf(stdout," %6.1f |", value);
	} else
	  fprintf(stdout,"        |");
      }
      fprintf(stdout, "\n");
	  
      fprintf(stdout, " average   |");
      for (j=0; j < 6; j++) {
	if (distance_visits[j] > 0) {
	  fprintf(stdout," %6.3f |", distance_avg[j]);
	} else
	  fprintf(stdout,"        |");
      }
      fprintf(stdout, "\n");

      fprintf(stdout, " stddev    |");
      for (j=0; j < 6; j++) {
	if (distance_visits[j] > 1) {
	  fprintf(stdout," %6.3f |", distance_sd[j]);
	} else
	  fprintf(stdout,"        |");
      }
      fprintf(stdout, "\n           --------------------------");
      fprintf(stdout, "-----------------------------\n");

      if (prnlev > 0) {
      fprintf(stdout, "\nTRANSITIONS TABLE: (from/vertical to/horizontal)\n\n");
      fprintf(stdout, "            %s  %s  %s  %s  %s  %s\n",
	      distance_ss[0], distance_ss[1], distance_ss[2],
	      distance_ss[3], distance_ss[4], distance_ss[5]);
      fprintf(stdout, "           -----------------------");
      fprintf(stdout, "--------------------------------\n");
      for (j=0; j<6; j++) {
	fprintf(stdout, "   %s |", distance_ss[j]);
	for (k=0; k<6; k++) {
	  if (distance_transitions[j][k] > 0)
	    fprintf(stdout, " %6i |", distance_transitions[j][k]);
	  else
	    fprintf(stdout, "        |");
	}
	fprintf(stdout, "\n");
      }
      fprintf(stdout, "           ------------------");
      fprintf(stdout, "-------------------------------------\n\n");
      }
      break;


    }

    if (analyze->iarg1 != 1) return 1;
  }

  return 1;

}

