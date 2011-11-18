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

#ifndef NO_PTRAJ_ANALYZE
/** ANALYZE ROUTINE *************************************************************
 *
 *  analyzeMatrix() --- diagonalizes matrices
 *
 *  Supplementary routines:
 *    dotprod (below)
 *    freeAnalyzeMatrixMemory (below)
 *
 ******************************************************************************/

   void
dotprod(int nelem, double *mat, double *vec, double *target){
  
  int i, j, ind;
  for(i = 0; i < nelem; i++){
    target[i] = 0.0;
  }

  for(i = 0; i < nelem; i++){
    for(j = i; j < nelem; j++){
      ind = nelem * i + j - (i * (i + 1)) / 2;
      target[i] += mat[ind] * vec[j];
      if(i != j)
        target[j] += mat[ind] * vec[i];
    }
  }
}

void freeAnalyzeMatrixMemory(analyzeInformation *analyze){

  modesInfo *modinfo;

  /* 
   * No freeing of analyze->carg1 because 
   * this is taken care of in action matrix 
   */
  if(analyze->carg2 != NULL)
    safe_free(analyze->carg2);

  if (analyze->carg4 != NULL)
    safe_free(analyze->carg4);

  modinfo = (modesInfo *) analyze->carg3;
  if(modinfo != NULL){
    /* 
     * modinfo->avg is cleaned by action matrix
     */
    if(modinfo->name != NULL)
      safe_free(modinfo->name);
    if(modinfo->freq != NULL)
      safe_free(modinfo->freq);
    if(modinfo->evec != NULL)
      safe_free(modinfo->evec);

    INITIALIZE_modesInfo(modinfo);
    safe_free(modinfo);
  }
}

   int
analyzeMatrix(analyzeInformation *analyze, stackType *sp, int mode)
{
  argStackType **argumentStackPointer;
  char * buffer;

  /*
   *  usage:
   *
   *    analyze matrix <matrixname> [out <filename>] [name <name>] [thermo | order] [vecs <vecs>] [reduce]
   *
   *  - This routine assumes that the matrix is represented as "the right upper triangle"
   *      including the main diagonal.
   *
   *  argument usage:
   *
   *    iarg1:
   *      0 -- no output of thermodynamic data
   *      1 -- output of thermodynamic data
   *    iarg2:
   *      x -- number of eigenvectors to calculate (must be >= 1, by default set to 1)
   *    iarg3:
   *      0 -- no reduction of vectors
   *      1 -- reduction of vectors
   *    carg1: the matrixInfo entry
   *    carg2: the outfile name
   *    carg3: the modesInfo entry
   *    carg4: order parameter filename
   */

  FILE *outfile;
  transformMatrixInfo *minfo;
  transformMatrixType type;
  modesInfo *modinfo;
  ptrajState *state;
  double tol, sigma, mass, temp, pressure, sum, tmp;
  double *vect, *mat, *workl, *workd, *eigval, *eigvali, *oparams, 
         *vout, *voutput, *resid, *vibn, *masses;
  int ithermo, reduce, mask1tot, nelem,
      lworkl, nevec, neval, nvec, ncv, ido, info, rvec, ldz,
      iparam[11], ipntr[11], *select, *mask1,
      i, j, k, n, crow, ind;
  char bmat[1], which[2], howmny[1], jobz[1], uplo[1];
  char *outfilename,*orderparamfilename;

  if (mode == PTRAJ_SETUP) {

    /*
     *  PTRAJ_SETUP:
     */

    argumentStackPointer = (argStackType **) analyze->carg1;
    analyze->carg1 = NULL;

    /*
     *  process matrix name
     */
    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL) {
      fprintf(stderr, "ptraj(), analyzeMatrix: missing the first name on\n");
      fprintf(stderr, "matrixStack to process...  Returning.\n");
      return -1;
    }
    minfo = matrixInfoStackGetName(&matrixStack, buffer); 
    if (minfo == NULL) {
      /*
       *  try searching for a match without a leading $ just in case the name
       *  was padded
       */
    minfo = matrixInfoStackGetName(&matrixStack, buffer+1); 
    }
    if (minfo == NULL) {
      fprintf(stderr, "ptraj(), analyzeMatrix: cannot find a match in the\n");
      fprintf(stderr, "matrixInfoStack for name (%s), returning.\n", buffer);
      safe_free(buffer);
      return -1;
    }
    analyze->carg1 = (void *) minfo;
    safe_free(buffer); 

    /*
     * process outfile name
     */
    analyze->carg2 = argumentStackKeyToString( argumentStackPointer, "out", NULL );

    // Order param file name
    analyze->carg4 = (char*) argumentStackKeyToString( argumentStackPointer,"orderparamfile",NULL);

    /*
     * process thermo and order flag
     */
    if( argumentStackContains( argumentStackPointer, "thermo" ) )
      analyze->iarg1 = 1;
    else if( argumentStackContains( argumentStackPointer, "order" ) )
      analyze->iarg1 = 2;  
    else
      analyze->iarg1 = 0;
    if (analyze->iarg1 == 1 && minfo->type != MATRIX_MWCOVAR) {
      fprintf(stderr, "ptraj(), analyzeMatrix: parameter \"thermo\" only\n");
      fprintf(stderr, "works for MATRIX_MWCOVAR, returning.\n");
      return -1;
    }
    else if (analyze->iarg1 == 2 && minfo->type != MATRIX_IRED) {
      fprintf(stderr, "ptraj(), analyzeMatrix: parameter \"order\" only\n");
      fprintf(stderr, "works for MATRIX_IRED, returning.\n");
      return -1;
    }


    /*
     * process nof eigenvectors - allow "0" only in the case of "thermo" flag given
     */
    nevec = argumentStackKeyToInteger(argumentStackPointer, "vecs", 0);
    if(analyze->iarg1 == 1){
      if(nevec < 0)
	nevec = 0;
    }
    else if(nevec <= 0){
        nevec = 1;
    }
    analyze->iarg2 = nevec;

    /*
     * process reduce flag
     */
    if( argumentStackContains( argumentStackPointer, "reduce" ) )
      analyze->iarg3 = 1;
    else
      analyze->iarg3 = 0;
    if (analyze->iarg3 > 0 && minfo->type != MATRIX_MWCOVAR &&
                              minfo->type != MATRIX_COVAR   &&
                              minfo->type != MATRIX_DISTCOVAR) {
      fprintf(stderr, "ptraj(), analyzeMatrix: parameter \"reduce\" only\n");
      fprintf(stderr, "works for MATRIX_MWCOVAR, ..._COVAR, ..._DISTCOVAR, returning.\n");
      return -1;
    }

    /*
     * Generate modesInfo struct
     */
    modinfo = (modesInfo *) safe_malloc(sizeof(modesInfo));
    INITIALIZE_modesInfo(modinfo);
    analyze->carg3 = (modesInfo *) modinfo;

    type = minfo->type;
    modinfo->type = (type == MATRIX_DIST ? MT_DIST :
                     (type == MATRIX_COVAR ? MT_COVAR :
                      (type == MATRIX_MWCOVAR ? MT_MWCOVAR :
                       (type == MATRIX_DISTCOVAR ? MT_DISTCOVAR :
                        (type == MATRIX_CORREL ? MT_CORREL :
                         (type == MATRIX_IDEA ? MT_IDEA :
                          (type == MATRIX_IRED ? MT_IRED : MT_UNKNOWN
                          )
                         )
                        )
                       )
                      )
                     )
                    );

    mask1tot = minfo->mask1tot;
    if(type == MATRIX_DIST || type == MATRIX_IDEA || type == MATRIX_IRED)
      nelem = mask1tot;
    else if(type == MATRIX_DISTCOVAR)
      nelem = mask1tot * (mask1tot - 1) / 2;
    else
      nelem = 3 * mask1tot;
    modinfo->navgelem = nelem;

    /*
     * process name for modesStack -> store modesInfo on modesStack
     */
    buffer = argumentStackKeyToString( argumentStackPointer, "name", NULL );
    if(buffer != NULL){
      modinfo->name = buffer;
      modinfo->source = MS_STACK;
      pushBottomStack(&modesStack, (void *) modinfo);
    }

    return 0;
  }
  else if (mode == PTRAJ_STATUS) {
    /*
     *  ACTION: PTRAJ_STATUS
     */
    fprintf(stdout,"  ANALYZE MATRIX: Analyzing matrix %s and dumping results to %s\n", 
                                      ((transformMatrixInfo *) analyze->carg1)->name, 
                                      (analyze->carg2 != NULL ? (char *) analyze->carg2 : "stdout"));
    fprintf(stdout,"      Calculating %i eigenvectors and %sthermodynamic data\n",
                                    analyze->iarg2, (analyze->iarg1 == 1 ? "" : "no "));
    fprintf(stdout,"      Calculating %i eigenvectors and %sorder parameters\n",
                                    analyze->iarg2, (analyze->iarg1 == 2 ? "" : "no "));
    if(analyze->iarg2 && analyze->iarg3)
      fprintf(stdout,"      Eigenvectors will be reduced\n");
    if(analyze->carg3 != NULL && ((modesInfo *) analyze->carg3)->name != NULL)
      fprintf(stdout,"      Storing modes on internal stack with name: %s\n", 
                            ((modesInfo *) analyze->carg3)->name);
     if(analyze->carg4!=NULL)
       fprintf(stdout,"      Order parameters will be written to %s\n",(char*) analyze->carg4);
  }
  else if (mode == PTRAJ_CLEANUP) {
    /*
     *  ACTION: PTRAJ_CLEANUP
     */
    freeAnalyzeMatrixMemory(analyze);
  }

  /*
   *  Other possible modes:
   *
   *  PTRAJ_PRINT   -- dump information to output files
   *  PTRAJ_NOOP    -- NULL mode
   */

  if (mode != PTRAJ_ACTION) return 0;

  /*
   *  Get data
   */
  minfo       = (transformMatrixInfo *) analyze->carg1;
  outfilename = (char *) analyze->carg2;
  orderparamfilename = (char*) analyze->carg4;
  modinfo     = (modesInfo *) analyze->carg3;
  ithermo     = analyze->iarg1;
  nevec       = analyze->iarg2;
  reduce      = analyze->iarg3;

  nelem       = modinfo->navgelem;
  mask1       = minfo->mask1;
  mask1tot    = minfo->mask1tot;
  type        = minfo->type;
  vect        = minfo->vect;
  mat         = minfo->mat;
  state       = minfo->state;

  /*
   *  Find eigenvalues and eigenvectors
   */

  if(nevec > nelem){
    nevec = nelem;
    fprintf(stderr,
	    "FYI: NEVEC > NELEM: Number of calculated evecs were reduced to %i\n", 
	    nevec);
  }
  if(nevec >= minfo->snap){
    nevec = minfo->snap;
    fprintf(stderr,
	    "FYI: NEVEC > SNAP: Number of calculated evecs were reduced to %i\n", 
	    nevec);
  }

  if(nevec == 0 || nelem == nevec){
    neval = nelem;
    eigval = (double *) safe_malloc(nelem * sizeof(double));
    /* scratch for the function dspev */
    workd  = (double *) safe_malloc(3 * nelem * sizeof(double));
    /* lower triangle of matrix is expected as input for dspev */
    uplo[0] = 'L';

    if (nevec == nelem){
    /* get all eigenvectors */
      /* for the eigenvectors */
      vout   = (double *) safe_malloc(nelem * nelem * sizeof(double));
      /* calculate eigenvalues and eigenvectors */
      jobz[0] = 'V';
      /* dimension of "vout" */
      ldz = nelem;
    }
    else if(nevec == 0){
    /* get only eigenvalues */
    /* only possible if thermo flag is set, */
    /* otherwise nevec is set to 1 per default. */
      /* for the eigenvectors */
      vout   = (double *) safe_malloc(nelem * sizeof(double));
      /* only calculate eigenvalues */
      jobz[0] = 'N';
      /* dimension of "vout" */
      ldz = 1;
    }

    dspev_(jobz, uplo, &nelem, mat, eigval, vout, &ldz, workd, &info);
    if(info != 0){
      fprintf(stderr,"  Warning in analyzeMatrix: dspev returned info = %i\n", info);
      return 0;
    }

    safe_free(workd);
  }
  else{
  /* get up to n-1 eigenvectors */
    neval = nevec;

    if(2*nevec <= nelem){
      ncv = 2*nevec;
      vout   = (double *) safe_malloc(nelem * ncv * sizeof(double));
    }
    else{
      ncv = nelem;
      if(reduce)
        vout   = (double *) safe_malloc(nelem * 2 * ncv * sizeof(double));
      else
        vout   = (double *) safe_malloc(nelem * ncv * sizeof(double));
    } 
    lworkl = ncv * (ncv + 8);
    eigval = (double *) safe_malloc(nelem * sizeof(double));
    workl  = (double *) safe_malloc(lworkl * sizeof(double));
    workd  = (double *) safe_malloc(3 * nelem * sizeof(double));
    resid  = (double *) safe_malloc(nelem * sizeof(double));
    select = (int *)    safe_malloc(ncv * sizeof(int));

    ido = 0;
    info = 0;
    iparam[0] = 1;
    iparam[2] = 300;
    iparam[3] = 1;
    iparam[6] = 1;
    tol = 0.0;
    bmat[0] = 'I';
    which[0] = 'L';
    which[1] = 'A';

    do{
      if(ido == -1 || ido == 1){
        dotprod(nelem, mat, workd + (ipntr[0] - 1), workd + (ipntr[1] - 1));
      }

      dsaupd_(&ido, bmat, &nelem, which, &nevec, &tol, resid, 
              &ncv, vout, &nelem, iparam, ipntr, workd, workl,
              &lworkl, &info);
    } while (ido == -1 || ido == 1);

    if(info != 0){
      fprintf(stderr,"  Warning in analyzeMatrix: dsaupd returned info = %i\n", info);
      fprintf(stderr,"  IPARAM(5) = %d\n", iparam[4]);
      return 0;
    }

    rvec = 1;
    howmny[0] = 'A';

    dseupd_(&rvec, howmny, select, eigval, vout, &nelem, &sigma,
            bmat, &nelem, which, &nevec, &tol, resid,
            &ncv, vout, &nelem, iparam, ipntr, workd, workl,
            &lworkl, &info);

    safe_free(workl);
    safe_free(workd);
    safe_free(resid);
    safe_free(select);
  }

  if(type == MATRIX_MWCOVAR){
    /*
     *  Convert eigenvalues to cm^-1
     */
    for(i = 0; i < neval; i++){
      if(eigval[i] > 0.0)
        eigval[i] = 108.587 * sqrt(0.6/eigval[i]); /* "0.6" is conversion of kT for 300K into kcal/mol(?) */
      else if (eigval[i] < 0.0)
        eigval[i] = -108.587 * sqrt(-0.6/eigval[i]);
      else{
        fprintf(stderr,"  Warning in analyzeMatrix: bad eigenvalue %i %6.2f\n", i, eigval[i]);
        return 0;
      }
    } 

    /*
     *  Mass-weight eigenvectors (and calc thermo chemistry)
     */
    crow = 0;
    for(i = 0; i < state->atoms; i++){
      if(mask1[i]){
        mass = 1.0 / sqrt(state->masses[i]);
        for(j = 0; j < nevec; j++){
          for(k = 0; k < 3; k++){
            ind = j * nelem + crow * 3 + k;
            vout[ind] = vout[ind] * mass;
          }
        }
        crow++;
      }
    }

    if(ithermo == 1){
      eigvali = (double *) safe_malloc(neval * sizeof(double));
      vibn = (double *) safe_malloc(4 * neval * sizeof(double));
      masses = (double *) safe_malloc(mask1tot * sizeof(double));

      for(i = 0; i < neval; i++)
        eigvali[i] = eigval[neval - 1 - i];

      crow = 0;
      for(i = 0; i < state->atoms; i++)
        if(mask1[i])
          masses[crow++] = state->masses[i];

      j = 1;
      temp = 298.15;
      pressure = 1.0;
      thermo_(&mask1tot, &neval, &j, vect, masses, eigvali, 
              vibn, vibn + 1*neval, vibn + 2*neval, vibn + 3*neval,
              &temp, &pressure);

      safe_free(eigvali);
      safe_free(vibn);
      safe_free(masses);
    }
  }

  /*
   * Calculation of S2 order parameters according to Prompers & Brüschweiler, JACS  124, 4522, 2002; added by A.N. Koller & H. Gohlke
   */
  if(type == MATRIX_IRED){
    if(ithermo == 2){
      if (orderparamfilename==NULL) {
        outfile=stdout;
      } else {
         //if ( openFile(&outfile,orderparamfilename,"w")==0 ) {
         if ( (outfile = safe_fopen(orderparamfilename,"w"))==NULL) {
           fprintf(stdout,"Error: Could not open %s.\n",orderparamfilename);
           fprintf(stdout,"Defaulting to STDOUT.\n");
           outfile=stdout;
         }
      }
      fprintf(outfile, "\n\t************************************\n\t- Calculated iRed order parameters -\n\t************************************\n\n");
      fprintf(outfile, "vector    S2\n----------------------\n");  

      for(i = 0; i < nelem; i++){ 			/* loop over all vector elements*/
          sum  = 0.0;                                   /* sum according to Eq. A22 in Prompers & Brüschweiler, JACS 124, 4522, 2002 */
          for(j = nevec - 6 ; j >= 0; j--){ 	        /* loop over all eigenvectors except the first five ones */
             sum += eigval[j] * vout[j * nelem + i] * vout[j * nelem + i];
          }
          fprintf(outfile, " %4i  %10.5f\n", i, 1.0 - sum);
      }
      if (outfile!=stdout) safe_fclose(outfile);
    }
  }

  /*
   * Reduction of eigenvectors (s. Abseher & Nilges, JMB 1998, 279, 911-920.)
   */
  if(reduce){
    voutput = vout + nevec * nelem;
    n = mask1tot;
    if(type == MATRIX_COVAR || type == MATRIX_MWCOVAR){
      for(i = 0; i < nevec; i++){
        for(j = 0; j < n; j++){
          voutput[i * n + j] = vout[i * nelem + j * 3    ] * vout[i * nelem + j * 3    ] + 
                               vout[i * nelem + j * 3 + 1] * vout[i * nelem + j * 3 + 1] +
                               vout[i * nelem + j * 3 + 2] * vout[i * nelem + j * 3 + 2];
        }
      }
    }
    else if(type == MATRIX_DISTCOVAR){
      for(i = 0; i < nevec; i++){
        for(j = 0; j < n; j++){
          voutput[i * n + j] = 0.0;
          for(k = 0; k < n; k++){
            if(k != j){
              ind = distindex(mask1tot, (j < k ? j : k), (j < k ? k : j));
              voutput[i * n + j] += vout[i * nelem + ind] * vout[i * nelem + ind];
            }
          }
        }
      }
    }
  }
  else{
    voutput = vout;
    n = nelem;
  }

  /*
   * Output average coordinates, eigenvectors, and eigenvalues
   */
  if(outfilename != NULL){
    outfile = safe_fopen(outfilename, "w");
    if(outfilename == NULL){
      fprintf(stderr, "WARNING in ptraj(), analyzeMatrix: error on opening %s for output\n",
	      outfilename);
      return 0;
    }

    /*
     * Average coordinates in the case of COVAR, MWCOVAR, CORREL;
     *   average distances in the case of DISTCOVAR;
     *   average of r*r/3 with r = "center of mass to atom" vector in the case of IDEA;
     *   average of Legendre polynomial P(cos(angle(ri, rj))) in the case of IRED;
     *   nothing in the case of DIST
     */
    fprintf(outfile, " %sEigenvector file: ",(reduce ? "Reduced " : ""));
    fprintf(outfile, "%s\n", (type == MATRIX_DIST ? "DIST" :
                              (type == MATRIX_COVAR ? "COVAR" :
                               (type == MATRIX_MWCOVAR ? "MWCOVAR" :
                                (type == MATRIX_DISTCOVAR ? "DISTCOVAR" :
                                 (type == MATRIX_CORREL ? "CORREL" :
                                  (type == MATRIX_IDEA ? "IDEA" :
                                   (type == MATRIX_IRED ? "IRED" : "UNKNOWN"
                                   )
                                  )
                                 )
                                )
                               )
                              )
                             ));
    fprintf(outfile, " %4i %4i\n", nelem, n);
    if(type != MATRIX_DIST){
      for(i = 0; i < nelem; i ++){
        fprintf(outfile, " %10.5f", vect[i]);
        if((i+1) % 7 == 0)
          fprintf(outfile, "\n");
      }
      if(nelem%7 != 0)
        fprintf(outfile, "\n");
    }

    /*
     * Eigenvectors and eigenvalues
     */
    for(i = neval - 1; i >= 0; i--){
      fprintf(outfile, " ****\n");
      fprintf(outfile, " %4i  %10.5f\n", neval - i, eigval[i]);
      if(nevec > 0){
        for(j = 0; j < n; j++){
          fprintf(outfile, " %10.5f", voutput[i * n + j]);
          if((j+1) % 7 == 0)
            fprintf(outfile, "\n");
        }
      }
      if(n % 7 != 0)
        fprintf(outfile, "\n");
    }

    safe_fclose(outfile);
  }

  /*
   * Get same ordering of arrays as in output if storing on stack requested
   */
  if(modinfo->source == MS_STACK){

    k = (int) (neval / 2);
    for(i = 0; i < k; i++){
      ind = neval - 1 - i;
      tmp = eigval[i];
      eigval[i] = eigval[ind];
      eigval[ind] = tmp;

      if(reduce){
        /*
         * Shift values from voutput(=vout+nevec*nelem) to vout and invert on the fly
         */
        if(nevec > 0)
          for(j = 0; j < n; j++)
            vout[i * n + j] = voutput[ind * n + j];
      }
      else{
        /*
         * Invert within vout(=voutput)
         */
        if(nevec > 0){
          for(j = 0; j < n; j ++){
            tmp = vout[i * n + j];
            vout[i * n + j] = vout[ind * n + j];
            vout[ind * n + j] = tmp;
          }
        }
      }
    }
  }

  /*
   * Store average coordinates, eigenvectors, and eigenvalues in modesInfo
   */
  modinfo->avg = vect;
  if(nevec > 0)
    modinfo->nvect = neval;
  else
    modinfo->nvect = 0;
  modinfo->nvectelem = n;
  modinfo->freq = eigval;
  modinfo->evec = vout;

  return 1;
}
#endif // ifndef NO_PTRAJ_ANALYZE

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

#ifndef NO_PTRAJ_ANALYZE
/** ANALYZE ROUTINE *************************************************************
 *
 *  analyzeTimecorr() --- auto and cross correlation functions 
 *
 *  Supplementary routines:
 *    corfdir (below)
 *    corffft (below)
 *    freeAnalyzeTimecorrMemory (below)
 *
 ******************************************************************************/
typedef enum _analyzeTimecorrMode {
  ATCM_UNKNOWN,
  ATCM_AUTO,
  ATCM_CROSS
} analyzeTimecorrMode;

typedef enum _analyzeTimecorrType {
  ATCT_UNKNOWN,
  ATCT_IRED,
  ATCT_NORMAL
} analyzeTimecorrType;

typedef struct _timecorrResults {
  analyzeTimecorrMode mode;
  analyzeTimecorrType type;
  int dplr;
  int norm;
  int drct;
  int relax;
  double tstep;
  double tcorr;
  double distnh;
  double freq;
  char *filename;
  char *noeFilename;
  int ndata;
  double *table;
  double *data1;
  double *data2;
  double *cf;
  double *cfinf;
  double *p2cf;
  double *rcf;
} timecorrResults;

#define INITIALIZE_timecorrResults(_p_) \
  _p_->mode     = ATCM_UNKNOWN;         \
  _p_->type     = ATCT_UNKNOWN;         \
  _p_->dplr     = 0;                    \
  _p_->norm     = 0;                    \
  _p_->drct     = 0;                    \
  _p_->relax    = 0;                    \
  _p_->tstep    = 1.0;                  \
  _p_->tcorr    = 10000.0;              \
  _p_->distnh   = 1.02;                 \
  _p_->freq     = -1.0;                 \
  _p_->filename = NULL;                 \
  _p_->noeFilename = NULL;              \
  _p_->ndata    = 0;                    \
  _p_->table    = NULL;                 \
  _p_->data1    = NULL;                 \
  _p_->data2    = NULL;                 \
  _p_->cf       = NULL;                 \
  _p_->cfinf    = NULL;                 \
  _p_->p2cf     = NULL;                 \
  _p_->rcf      = NULL;

   void
corfdir(int ndata, double *data1, double *data2, int nsteps, double *dtmp){

  /*
   * Calculates correlation functions
   * using the "direct" approach
   * (s. Comp. Sim. of Liquids, p.185)
   * - the result is not yet normalized by (no_of_discrete_data - t)**-1 (!)
   */

  int i, j, ndata2;
  int ind1, ind2;
  double dsum, dsumi;

  ndata2 = ndata / 2;

  if(data2 == NULL){
    for(i = 0; i < ndata2; i++){
      dsum = 0.0;
      for(j = i; j < ndata2; j++){
        ind1 = 2 * j;
        ind2 = 2 * (j-i);
        dsum += data1[ind1] * data1[ind2] + data1[ind1+1] * data1[ind2+1];
      }
      if(i < nsteps){
        ind1 = 2 * i;
        dtmp[ind1  ] = dsum;
        dtmp[ind1+1] = 0.0;
      }
      else{
        break;
      }
    }
  }
  else{    
    for(i = 0; i < ndata2; i++){
      dsum = 0.0;
      dsumi = 0.0;
      for(j = i; j < ndata2; j++){
        ind1 = 2 * j;
        ind2 = 2 * (j-i);
        dsum  += data2[ind1] * data1[ind2  ] + data2[ind1+1] * data1[ind2+1];
        dsumi += data2[ind1] * data1[ind2+1] - data2[ind1+1] * data1[ind2  ];
      }
      if(i < nsteps){
        ind1 = 2 * i;
        dtmp[ind1  ] = dsum;
        dtmp[ind1+1] = dsumi;
      }
      else{
        break;
      }
    }
  }

  for(i = 0; i < nsteps; i++){
    ind1 = 2 * i;
    data1[ind1  ] = dtmp[ind1  ];
    data1[ind1+1] = dtmp[ind1+1];
  }
}

   void
corffft(int ndata, double *data1, double *data2, double *table){

  /*
   * Calculates correlation functions
   * using the Wiener-Khinchin-Theorem
   * (s. Comp. Sim. of Liquids, p. 188)
   *
   * - ndata is the length of the arrays data1/2
   * - data1/2 is a real array of complex numbers: 
   *   data1/2(1)=real1/2(1), data1/2(2)=img1/2(1), ...
   * - it is recommended that the discrete data is appended
   *     with as much zeros to avoid spurious correlations
   * - in addition, data MUST have the dimension of power of 2
   *     (pad the "real" data (plus zeros) with additional 
   *     zeros up to the next power of 2)
   * - the result is not yet normalized by (no_of_discrete_data - t)**-1 (!)
   */

  int i, ndata2;
  double dtmp;
  
  ndata2 = ndata / 2;

  /*
   * FFT data
   */
  cfftf_(&ndata2, data1, table);
  if(data2 != NULL)
    cfftf_(&ndata2, data2, table);

  /*
   * Calc square modulus (in case of cross-correlation: calc [F(data1)]' * F(data2),
   *   where [F(data1)]' is the complex conjugate of F(data1)).
   */
  if(data2 == NULL)
    for(i = 0; i < ndata; i+=2){
      data1[i  ] = data1[i  ] * data1[i  ] + data1[i+1] * data1[i+1];
      data1[i+1] = 0.0;
    }
  else
    for(i = 0; i < ndata; i+=2){
      dtmp       = data1[i  ] * data2[i  ] + data1[i+1] * data2[i+1];
      data1[i+1] = data1[i  ] * data2[i+1] - data2[i  ] * data1[i+1];
      data1[i  ] = dtmp;
    }

  /*
   * Inverse FFT
   */
  cfftb_(&ndata2, data1, table);

  /*
   * Normalize with ndata/2 (since not done in inverse FFT routine)
   */
  dtmp = 1.0 / ((double) (ndata2));
  for(i = 0; i < ndata; i++)
    data1[i] *= dtmp;
}

   void
freeAnalyzeTimecorrMemory(analyzeInformation *analyze){

  timecorrResults *tcr = (timecorrResults *) analyze->carg1;

  if(tcr != NULL){
    if(tcr->filename != NULL) safe_free(tcr->filename);
    if(tcr->noeFilename!=NULL) safe_free(tcr->noeFilename);
    if(tcr->table != NULL) safe_free(tcr->table);
    if(tcr->data1 != NULL) safe_free(tcr->data1);
    if(tcr->data2 != NULL) safe_free(tcr->data2);
    if(tcr->cf    != NULL) safe_free(tcr->cf);
    if(tcr->cfinf != NULL) safe_free(tcr->cfinf);
    if(tcr->p2cf  != NULL) safe_free(tcr->p2cf);
    if(tcr->rcf   != NULL) safe_free(tcr->rcf);
    INITIALIZE_timecorrResults(tcr);
    safe_free(tcr);
  }
}

  double
calc_spectral_density(int nevec, int nelem, double *eigval, double *vout, double *taum, int i, double omega)
{
/*Calc spectral density (JACS 2002, 124, 4522; eq. A24)*/
  
  int j;
  double J;
                               
  J = 0.0;
  for(j = 0 ; j < nevec; j++){
    J += (eigval[j] * (vout[j * nelem + i] * vout[j * nelem + i])) * 2.0 * taum[j] /
         ( 1.0 + omega*omega * taum[j]*taum[j] );/*check order XXXX*/
  }                                                                                                                                  
  return J;
}
  

/*only for one vector*/
   int
analyzeTimecorr(analyzeInformation *analyze, stackType *sp, int mode)
{

  /* Variable declarations */
  argStackType **argumentStackPointer;
  char *buffer;

  //stackType *vectorStackTmp = NULL;
  //transformVectorInfo *vinfo;
  transformVectorInfo *vinfo1;
  transformVectorInfo *vinfo2;
  timecorrResults *timecorr;
  analyzeTimecorrMode tcmode;
  analyzeTimecorrType type;
  int i, j, k, ind1, ind2, ind3;
  int ndata, ndata2, nsteps;
  int nvect, nvectelem, mtot, frame; /*XXX*/
  int dplr, norm, drct, relax; /*XXX*/
  double tstep, tcorr, distnh, freq; /*XXX*/
  double rave1, r3iave1, r6iave1, rave2, r3iave2, r6iave2;
  double dnorm;
  double cfinfavgreal, cfinfavgimg;
  double *eigval, *vout; /*XXX*/
  double *table, *data1, *data2, *cf, *cf_cjt, *cfinf, *p2cf, *rcf, *taum; /*XXX*/
  double *avgcrd1, *cftmp1, *p2cftmp1, *rcftmp1,
         *avgcrd2, *cftmp2, *p2cftmp2, *rcftmp2;
  double *dpt1, *dpt2;
  char *filename;
  FILE *fp, *fp_cmt, *fp_cjt;
  
  const double factor = 0.8 * 3.141592654;  /* = 4/5*PI due to spherical harmonics addition theorem */
  /*
   *  USAGE:
   *
   *  analyze timecorr 
   *                   vec1 <vecname1> [vec2 <vecname2>]
   *		       [relax] [freq <hz>] [NHdist <distnh>]
   *                   tstep <tstep> tcorr <tcorr> out <filename>
   *
   *  analyze argument usage:
   *    carg1:
   *      pointer to timecorrResults
   *    carg2:
   *      pointer to vectorInfo1
   *    carg3:
   *      pointer to vectorInfo2
   */

  if (mode == PTRAJ_SETUP) {
    /*
     *  -------- ANALYZE: PTRAJ_SETUP
     */
    argumentStackPointer = (argStackType **) analyze->carg1;
    analyze->carg1 = NULL;

    /*
     * Generate timecorrResults object
     */
    timecorr = (timecorrResults *) safe_malloc(sizeof(timecorrResults));
    INITIALIZE_timecorrResults(timecorr);
    analyze->carg1 = (void *) timecorr;

    /*
     *  Get vectors
     */ 
    if((buffer = argumentStackKeyToString(argumentStackPointer, "vec1", NULL)) == NULL){
      fprintf(stderr,
              "WARNING in ptraj(), analyze timecorr: no vec1 given, ignoring command\n");
      freeAnalyzeTimecorrMemory(analyze);
      return -1;
    }
    else{
      vinfo1 = vectorInfoStackGetName(&vectorStack, buffer);
      analyze->carg2 = (void*) vinfo1;
/*      vinfo1 = NULL;
      for(vectorStackTmp = vectorStack;
          vectorStackTmp != NULL;
          vectorStackTmp = vectorStackTmp->next){
        vinfo = (transformVectorInfo *) vectorStackTmp->entry;
        if(strcmp(vinfo->name, buffer) == 0){
          vinfo1 = vinfo;
          analyze->carg2 = (void *) vinfo1;
          break;
        }
      }*/
      if(vinfo1 == NULL){
        fprintf(stderr,
                "WARNING in ptraj(), analyze timecorr: no vector with name %s found, ignoring command\n", buffer);
        freeAnalyzeTimecorrMemory(analyze);
        return -1;
      }
    }
    safe_free(buffer);
    vinfo2 = NULL;
    if((buffer = argumentStackKeyToString(argumentStackPointer, "vec2", NULL)) != NULL){
      vinfo2 = vectorInfoStackGetName(&vectorStack, buffer);
      analyze->carg3 = (void *) vinfo2;
/*      vinfo2 = NULL;
      for(vectorStackTmp = vectorStack;
          vectorStackTmp != NULL;
          vectorStackTmp = vectorStackTmp->next){
        vinfo = (transformVectorInfo *) vectorStackTmp->entry;
        if(strcmp(vinfo->name, buffer) == 0){
          vinfo2 = vinfo;
          analyze->carg3 = (void *) vinfo2;
          break;
        }
      }*/
      if(vinfo2 == NULL){
        fprintf(stderr,
                "WARNING in ptraj(), analyze timecorr: no vector with name %s found, ignoring command\n", buffer);
        freeAnalyzeTimecorrMemory(analyze);
        return -1;
      }
    }
    safe_free(buffer);
    /*
     *  Determine mode and type
     */
    if(vinfo1->mode == VECTOR_CORRIRED){
      timecorr->mode = tcmode = ATCM_AUTO;
      timecorr->type = type = ATCT_IRED;
      if(vinfo2 != NULL){
        fprintf(stderr,
                "WARNING in ptraj(), analyze timecorr: only calculating IRED corr for vec1, ignoring vec2\n");
        vinfo2 = NULL;
      }
    }
    else if(vinfo1->mode == VECTOR_CORR ||
            vinfo1->mode == VECTOR_CORRPLANE){
      timecorr->type = type = ATCT_NORMAL;
      if(vinfo2 == NULL){
        timecorr->mode = tcmode = ATCM_AUTO;
      }
      else{
        if((vinfo2->mode != VECTOR_CORR && 
            vinfo2->mode != VECTOR_CORRPLANE) ||
           vinfo1->frame != vinfo2->frame ||
           vinfo1->order != vinfo2->order){ 
          fprintf(stderr,
                  "WARNING in ptraj(), analyze timecorr: different attributes for vec1 and vec2, ignoring command\n");
          freeAnalyzeTimecorrMemory(analyze);
          return -1;
        }
        else{
          timecorr->mode = tcmode = ATCM_CROSS;
        }
      }
    }
    else{
      fprintf(stderr,
              "WARNING in ptraj(), analyze timecorr: wrong vector type given, ignoring command\n");
      freeAnalyzeTimecorrMemory(analyze);
      return -1;
    }
    /*
     * Get tstep, tcorr, filename
     */
    timecorr->tstep = (double) argumentStackKeyToFloat(argumentStackPointer, "tstep", 1.0);
    timecorr->tcorr = (double) argumentStackKeyToFloat(argumentStackPointer, "tcorr", 10000.0);
    buffer = argumentStackKeyToString(argumentStackPointer, "noefile",NULL);
    if (buffer!=NULL) {
      timecorr->noeFilename=(char*) safe_malloc( (strlen(buffer)+1) * sizeof(char));
      strcpy(timecorr->noeFilename,buffer);  
    }
    if((buffer = argumentStackKeyToString(argumentStackPointer, "out", NULL)) == NULL){
      fprintf(stderr,
              "WARNING in ptraj(), analyze timecorr: no outfile given, ignoring command\n");
      freeAnalyzeTimecorrMemory(analyze);
      return -1;
    }
    else{
      timecorr->filename = buffer;
    }

    /*
     * Get dplr, norm, drct, relax
     */
    if(argumentStackContains(argumentStackPointer, "dplr"))
      timecorr->dplr = 1;
    if(argumentStackContains(argumentStackPointer, "norm"))
      timecorr->norm = 1;
    if(argumentStackContains(argumentStackPointer, "drct"))
      timecorr->drct = 1;
    if(argumentStackContains(argumentStackPointer, "relax"))
      timecorr->relax = 1;
    
    if(timecorr->relax == 1){
    /*
     * Get freq, NH-distance
     */
      timecorr->freq = (double) argumentStackKeyToFloat(argumentStackPointer, "freq", -1.0);
      if(timecorr->freq == -1.0){
        fprintf(stderr,
        "WARNING in ptraj(), analyze timecorr: no frequency for the calculation of the relaxation params given, ignoring relax command\n");
        timecorr->relax = 0;
        timecorr->freq = 0.0;
      }
      timecorr->distnh = (double) argumentStackKeyToFloat(argumentStackPointer, "NHdist", 1.02); /*1.02 * 10**(-10) in Angstroem */
    }
    return 0;

  } else if (mode == PTRAJ_STATUS) {

    /*
     *  -------- ANALYZE: PTRAJ_STATUS
     */
    timecorr = (timecorrResults *) analyze->carg1;
    fprintf(stdout, "  ANALYZE TIMECORR: Calculating %s-correlation function of vector(s) %s %s,\n",
	    (timecorr->mode == ATCM_AUTO ? "auto" : 
             (timecorr->mode == ATCM_CROSS ? "cross" : "????"
             )
            ),
            ((transformVectorInfo *) analyze->carg2)->name,
            analyze->carg3 != NULL ? ((transformVectorInfo *) analyze->carg3)->name : "");
    fprintf(stdout, "      for a correlation time of %f, using a time step of %f.\n",
            timecorr->tcorr,
            timecorr->tstep);
    fprintf(stdout, "      Corr. func. are %s%snormalized\n",
            timecorr->dplr ? "for dipolar interactions and " : "",
            timecorr->norm ? "" : "not ");
    fprintf(stdout, "      Corr. func. are calculated using the %s approach\n",
            timecorr->drct ? "direct" : "FFT");
    fprintf(stdout, "      TauM, relaxation rates and NOEs are %s calculated using the iRED approach\n",
            timecorr->relax ? "" : "not");
    fprintf(stdout, "      %s using an NH distance of %f and a frequency of %f.\n",
            timecorr->relax ? "" : "not",
            timecorr->distnh,
            timecorr->freq
            );
    if (timecorr->noeFilename!=NULL)
      fprintf(stdout,"      NOEs and relaxation rates are written to %s\n",timecorr->noeFilename); 
    fprintf(stdout, "      Results are written to %s\n", timecorr->filename);

  } else if (mode == PTRAJ_CLEANUP) {
    /*
     *  -------- ANALYZE: PTRAJ_CLEANUP
     */
    freeAnalyzeTimecorrMemory(analyze);
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
  timecorr = (timecorrResults *) analyze->carg1;  
  filename = timecorr->filename;
  tcmode   = timecorr->mode;
  type     = timecorr->type;
  tstep    = timecorr->tstep;
  tcorr    = timecorr->tcorr;
  dplr     = timecorr->dplr;
  norm     = timecorr->norm;
  drct     = timecorr->drct;
  relax    = timecorr->relax;
  distnh   = timecorr->distnh;
  freq     = timecorr->freq;
  
  vinfo1   = (transformVectorInfo *) analyze->carg2;
  avgcrd1  = vinfo1->avgcrd;
  rave1    = vinfo1->rave;
  r3iave1  = vinfo1->r3iave;
  r6iave1  = vinfo1->r6iave;
  cftmp1   = vinfo1->cftmp;
  p2cftmp1 = vinfo1->p2cftmp;
  rcftmp1  = vinfo1->rcftmp;
  vinfo2   = (transformVectorInfo *) analyze->carg3;
  if(vinfo2 != NULL){
    avgcrd2  = vinfo2->avgcrd;
    rave2    = vinfo2->rave;
    r3iave2  = vinfo2->r3iave;
    r6iave2  = vinfo2->r6iave;
    cftmp2   = vinfo2->cftmp;
    p2cftmp2 = vinfo2->p2cftmp;
    rcftmp2  = vinfo2->rcftmp;
  }
 
  frame    = vinfo1->frame;
  if(((int) (tcorr / tstep) + 1) > frame){
    nsteps = frame;
  }
  else{
    nsteps = (int) (tcorr / tstep) + 1;
  }
  if(drct){
    ndata  = 2 * frame;
  }
  else{  
    ndata  = ldexp(1.0, (int) (log( (double) (4 * frame - 1)) / log(2.0)) + 1);
  }
  ndata2 = ndata / 2;
  mtot     = 2 * vinfo1->order + 1;

  /*
   *  Allocate memory and initialize arrays
   */
  if(type == ATCT_IRED){
    nvect   = vinfo1->modinfo->nvect;
    nvectelem = vinfo1->modinfo->nvectelem;
    eigval = vinfo1->modinfo->freq;
    vout = vinfo1->modinfo->evec;
    timecorr->cf = cf = (double *) safe_malloc(sizeof(double) * nvect * nsteps);
    cf_cjt = (double *) safe_malloc(sizeof(double) * nvect * nsteps);
    for(i = 0; i < nvect * nsteps; i++){
      cf[i] = 0.0;
      cf_cjt[i] = 0.0;
    }
    timecorr->cfinf = cfinf = (double *) safe_malloc(sizeof(double) * nvect);
    for(i = 0; i < nvect; i++)
      cfinf[i] = 0.0;
      
    if(relax==1){
      
      taum = (double *) safe_malloc(sizeof(double) * nvect);
      for(i = 0; i < nvect; i++)
        taum[i] = 0.0;
    }
  }
  else if(type == ATCT_NORMAL){
    timecorr->p2cf = p2cf = (double *) safe_malloc(sizeof(double) * nsteps);
    for(i = 0; i < nsteps; i++)
      p2cf[i] = 0.0;

    if(dplr){
      timecorr->cf = cf = (double *) safe_malloc(sizeof(double) * nsteps);
      timecorr->rcf = rcf  = (double *) safe_malloc(sizeof(double) * nsteps);
      for(i = 0; i < nsteps; i++){
        cf[i]   = 0.0;
        rcf[i]  = 0.0;
      }
    }
  }
  timecorr->data1 = data1 = (double *) safe_malloc(sizeof(double) * ndata);
  if(tcmode == ATCM_CROSS)
    timecorr->data2 = data2 = (double *) safe_malloc(sizeof(double) * ndata);
  if(drct){
    timecorr->table = table = (double *) safe_malloc(sizeof(double) * 2 * nsteps);
  }
  else{
    timecorr->table = table = (double *) safe_malloc(sizeof(double) * (2 * ndata + 15));
  }

  /*
   * Initialize FFT
   */
  if( ! drct)
    cffti_(&ndata2, table);
  if(type == ATCT_IRED){
    /*
     * Loop over all eigenvectors
     */
    for(i = 0; i < nvect; i++){
      ind1 = nsteps * i;
      /*
       * Loop over all m = -L, ..., L
       */
      for(j = 0; j < mtot; j++){
        ind2 = nvect * j;
        /*
         * Loop over all snapshots
         */
        cfinfavgreal = 0.0;
        cfinfavgimg = 0.0;
        for(k = 0; k < frame; k++){
          ind3 = 2 * (nvect * mtot * k + ind2 + i);
          data1[2*k  ] =  cftmp1[ind3  ];
          cfinfavgreal += cftmp1[ind3  ];
          data1[2*k+1] =  cftmp1[ind3+1];
          cfinfavgimg  += cftmp1[ind3+1];
        }
        cfinfavgreal /= frame;
        cfinfavgimg  /= frame;

	/*
	 * Calc plateau value of correlation function (= C(m,t->T) in Bruschweiler paper (A20))
	 */

        cfinf[i] += (cfinfavgreal * cfinfavgreal) + (cfinfavgimg * cfinfavgimg); 

        if(drct){
          /*
           * Calc correlation function (= C(m,l,t) in Bruschweiler paper) using direct approach
           */
          corfdir(ndata, data1, NULL, nsteps, table);
        }
        else{
          /*
           * Pad with zero's at the end 
           */
          for(k = 2 * frame; k < ndata; k++)
            data1[k] = 0.0;

          /*
           * Calc correlation function (= C(m,l,t) in Bruschweiler paper) using FFT
           */
          corffft(ndata, data1, NULL, table);
        }
 
        /*
         * Sum into cf (= C(m,t) in Bruschweiler paper)
         */
        for(k = 0; k < nsteps; k++){
          cf[ind1 + k] += data1[2 * k];
        }
      }
    }
   /*
    * Calculate correlation function for each vector:
    * Cj(t) according to eq. A23 in Prompers & Brüschweiler, JACS  124, 4522, 2002; added by A.N. Koller & H. Gohlke
    */
      for(i = 0; i < nvect; i++){ 
       for(k = 0; k < nsteps; k++){
       double sum = 0.0;
       for(j = 0; j < nvect; j++){ 
          sum+=(eigval[j] * (vout[j * nvectelem + i] * vout[j * nvectelem + i]))*(cf[nsteps * j + k] / (frame - k));
        }
        cf_cjt[nsteps * i + k]=sum;
      } 
    }
    
    if(relax == 1){
    /*
     * TauM calculation: eq. A19 from Prompers & Brüschweiler, JACS  124, 4522, 2002 is used.
     * Values are calculated from normalized correlation functions (Cm(t)s).
     * Added by Alrun N. Koller & H. Gohlke
     */
     double deltat = tstep * 1.0E-12; /* conversion to pico seconds (-> tstep) */
     if(nvectelem != nvect){
       fprintf(stderr,"  Warning in analyzeTimecorr: different number of eigen modes and eigenmode elements %i %i\n", nvect, nvectelem);
             return -1;
     }                                                      
     
     double new_nsteps = nsteps / 3.0; /* consider only first third of the correlation curve to avoid fitting errors */
           
     for(i = 0 ; i < nvect; i++){
        double integral = 0.0;
        double cfinfval = cfinf[i] / cf[nsteps * i] * frame;
	for ( j = 0 ; j < new_nsteps; j++ ) {
            double cfk=cf[nsteps * i + j] * frame / (cf[nsteps * i] * (frame - j));
            double cfk1=cf[nsteps * i + j + 1 ] * frame / (cf[nsteps * i] * (frame - (j+1)));
	    double fa = cfk - cfinfval;
	    double fb = cfk1 - cfinfval;
	    integral += tstep * 1.0E-12 * ( fa + fb ) * 0.5;
	}
        double taum_val = integral / ( 1.0 - cfinfval );
	taum[i]=taum_val;
     }	 
    
    /*
     * Relaxation calculation. Added by Alrun N. Koller & H. Gohlke
     */
      // DAN ROE - Using fp here - should be ok
      if (timecorr->noeFilename==NULL) {
        fp=stdout;
      } else {
        //if (openFile(&fp,timecorr->noeFilename,"w")==0) {
        if ( (fp = safe_fopen(timecorr->noeFilename,"w")) == NULL ) {
          fprintf(stdout,"Error: Could not open %s for writing.\n",timecorr->noeFilename);
          fprintf(stdout,"Redirecting to STDOUT\n");
          fp=stdout;
        }
      }
      fprintf(fp, "\n\t****************************************");
      fprintf(fp, "\n\t- Calculated relaxation rates and NOEs -");
      fprintf(fp, "\n\t****************************************\n\n");
      fprintf(fp, "vector   %10s   %10s   %10s\n","R1","R2","NOE");  
      
      double rnh = distnh * 1.0E-10; /* conversion from Angstrom to meter */
      
      /*
       * constants
       */                                                                                                             
      const double pi=3.141592654;
      const double mu_zero = 4.0 * pi * 1.0E-7;			/*in H m^-1; permeability */
      const double ha = 6.626176 * 1.0E-34;			/* in m^2 kg s1-1 ; Js ; Planck's constant */
      const double gamma_h = 2.6751987 * 1.0E8;   		/* in rad s^-1 T^-1 ; gyromagnetic ratio; T = absolut temperature in K */
      const double gamma_n = -2.7126 * 1.0E7;   		/* in rad s^-1 T^-1 ; gyromagnetic ratio */  
                                 
      const double csa = -170.0 * 1.0E-6;
                                                                            
      double spec_freq = freq * 2.0 * pi * 1.0E6;         	/* conversion from MHz to rad s^-1 */
      double b_zero = 2.0 * pi * freq * 1.0E6 / gamma_h;  	/* in T (Tesla) */
                                                                                                 
      const double lamfreqh = -1 * gamma_h * b_zero;		/* rad s^-1 */    
      const double lamfreqn = -1 * gamma_n * b_zero;		/* rad s^-1 */
      
      double c2 = lamfreqn*lamfreqn * csa*csa;
      double d2 = (mu_zero * ha * gamma_n * gamma_h)/( 8.0 * (pi*pi) * (rnh*rnh*rnh));
      d2 = d2*d2; /* fix from Junchao */
      
      for(i = 0; i < nvectelem; i++){ 			/* loop over all vector elements --> have only one element/vector here;  nvectelem = nelem*/
          
              /* Eq. A1 in Prompers & Brüschweiler, JACS  124, 4522, 2002 */
              double R1 = d2 / 20.0 * (
                                calc_spectral_density( nvect, nvectelem, eigval, vout, taum, i, lamfreqh - lamfreqn ) + 3.0 *
                                calc_spectral_density( nvect, nvectelem, eigval, vout, taum, i, lamfreqn ) + 6.0 *
                                calc_spectral_density( nvect, nvectelem, eigval, vout, taum, i, lamfreqh + lamfreqn)
                               ) + c2 / 15.0 * calc_spectral_density( nvect, nvectelem, eigval, vout, taum, i, lamfreqn );          
              
              /* Eq. A2 in Prompers & Brüschweiler, JACS  124, 4522, 2002 */                 
              double R2 = d2 / 40.0 * ( 4.0 *
                                 calc_spectral_density( nvect, nvectelem, eigval, vout, taum, i, 0.0 ) +
                                 calc_spectral_density( nvect, nvectelem, eigval, vout, taum, i, lamfreqh - lamfreqn ) + 3.0 *
                                 calc_spectral_density( nvect, nvectelem, eigval, vout, taum, i, lamfreqn ) + 6.0 *
                                 calc_spectral_density( nvect, nvectelem, eigval, vout, taum, i, lamfreqh ) + 6.0 *
                                 calc_spectral_density( nvect, nvectelem, eigval, vout, taum, i, lamfreqh + lamfreqn)
                               ) + c2 / 90.0 * ( 4.0 *
                                 calc_spectral_density( nvect, nvectelem, eigval, vout, taum, i, 0.0 ) + 3.0 *
                                 calc_spectral_density( nvect, nvectelem, eigval, vout, taum, i, lamfreqn));

              /* Eq. A3 in Prompers & Brüschweiler, JACS  124, 4522, 2002 */                   
              double Tj = d2 / 20.0 * ( 6.0 *
                                 calc_spectral_density( nvect, nvectelem, eigval, vout, taum, i, lamfreqh + lamfreqn ) -
                                 calc_spectral_density( nvect, nvectelem, eigval, vout, taum, i, lamfreqh + lamfreqn ));
                                             
              double Noe = 1.0 + ( gamma_h * 1.0/gamma_n ) * ( 1.0 / R1 ) * Tj;
                                                                              
          
           fprintf(fp, "%6i   %10.5f   %10.5f   %10.5f\n", i, R1, R2, Noe);
      }
      fprintf(fp, "\n\n");
      if (fp!=stdout) safe_fclose(fp);
    }
  }
  else if(type == ATCT_NORMAL){
    /*
     * Loop over all three correlation functions
     */
    for(i = 0; i < 3; i++){
      if(i == 0 && dplr){
        dpt1 = cftmp1;
        dpt2 = cftmp2;
      }
      else if(i == 1){
        dpt1 = p2cftmp1;
        dpt2 = p2cftmp2;
      }
      else if(i == 2 && dplr){
        dpt1 = rcftmp1;
        dpt2 = rcftmp2;
      }
      else{
        continue;
      }

      /*
       * Loop over all m=-L, ..., L
       */
      for(j = 0; j < mtot; j++){
        /*
         * Loop over all snapshots
         */
        for(k = 0; k < frame; k++){
          if(i < 2){
            ind1 = 2 * (mtot * k + j);
            ind2 = 2 * k;
            data1[ind2  ] = dpt1[ind1  ];
            data1[ind2+1] = dpt1[ind1+1];
            if(vinfo2 != NULL){
              data2[ind2  ] = dpt2[ind1  ];
              data2[ind2+1] = dpt2[ind1+1];
            }
          }
          else if(i == 2 && j == 0){
            ind1 = 2 * k;
            data1[ind1  ] = dpt1[ind1  ];
            data1[ind1+1] = dpt1[ind1+1];
            if(vinfo2 != NULL){
              data2[ind1  ] = dpt2[ind1  ];
              data2[ind1+1] = dpt2[ind1+1];
            }
          }
        }

        if(drct){
          /*
           * Calc correlation function using direct approach
           */
          if(vinfo2 == NULL)
            corfdir(ndata, data1, NULL, nsteps, table);
          else
            corfdir(ndata, data1, data2, nsteps, table);
        }
        else{
          /*
           * Pad with zero's at the end 
           */
          for(k = 2 * frame; k < ndata; k++){
            data1[k] = 0.0;
            if(vinfo2 != NULL)
              data2[k] = 0.0;
          }

          /*
           * Calc correlation function using FFT
           */
          if(vinfo2 == NULL)
            corffft(ndata, data1, NULL, table);
          else
            corffft(ndata, data1, data2, table);
        }
      
        /*
         * Sum into cf
         */
        for(k = 0; k < nsteps; k++){
          if(i == 0 && dplr)
            cf[k] += data1[2 * k];
          else if(i == 1)
            p2cf[k] += data1[2 * k];
          else if(i == 2 && j == 0 && dplr)
            rcf[k] += data1[2 * k];
          else
            break;
        }
      }
    }
  }
  /*
   * Normalize and output results
   */

  if(filename == NULL){
    fprintf(stderr, "WARNING in ptraj(), analyzeTimecorr: error on opening %s for output\n",
    filename);
    return 0;
  }
  
  if(type == ATCT_IRED){
    // Temporary filename for .cjt extension
    char *cjtfilename = (char *) malloc( (strlen(filename)+5) * sizeof(char));
    strcpy(cjtfilename, filename);
    strcat(cjtfilename, ".cjt");
    // Temporary filename for .cmt extension
    char *cmtfilename = (char*) malloc( (strlen(filename)+5) * sizeof(char));
    strcpy(cmtfilename, filename);
    strcat(cmtfilename, ".cmt");
    // Open output files  
    fp_cmt = safe_fopen(cmtfilename, "w");
    fp_cjt = safe_fopen(cjtfilename, "w");
    // Print header
    fprintf(fp_cmt,"%s-correlation functions Cm(t) for each eigenmode m, IRED type according to eq. A18 Prompers & Brüschweiler, JACS  124, 4522, 2002\n", 
          (tcmode == ATCM_AUTO? "Auto" :
           (tcmode == ATCM_CROSS? "Cross" : "????"
           )
          ));
    
    fprintf(fp_cjt,"%s-correlation functions Cj(t) for each ired vector j, IRED type according to eq. A23 Prompers & Brüschweiler, JACS  124, 4522, 2002\n", 
          (tcmode == ATCM_AUTO? "Auto" :
           (tcmode == ATCM_CROSS? "Cross" : "????"
           )
          ));

    
    fprintf(fp_cmt, "%12s", "XXX");
    fprintf(fp_cjt, "%12s", "XXX");
    for(i = 1; i <= nvect; i++){
      if(i < 10){
        fprintf(fp_cmt, "       Mode%i", i);
        fprintf(fp_cjt, "     Vector%i", i);
      }
      else if(i < 100){
        fprintf(fp_cmt, "      Mode%i", i);
        fprintf(fp_cjt, "    Vector%i", i);
      }
      else if(i < 1000){
        fprintf(fp_cmt, "     Mode%i", i);
        fprintf(fp_cjt, "   Vector%i", i);
      }
      else if(i < 10000){
        fprintf(fp_cmt, "    Mode%i", i);
        fprintf(fp_cjt, "  Vector%i", i);
      }
      else if(i < 100000){
        fprintf(fp_cmt, "   Mode%i", i);
        fprintf(fp_cjt, " Vector%i", i);
      }
    }
    fprintf(fp_cmt, "\n");
    fprintf(fp_cjt, "\n");

    /* Print cfinf */
    fprintf(fp_cmt, "%12s", "C(m,t->T)");
    for(i = 0; i < nvect; i++){
      if(norm)
	fprintf(fp_cmt, "%12.8f", cfinf[i] / cf[nsteps * i] * frame);
      else   
	fprintf(fp_cmt, "%12.8f", cfinf[i]);
    }
    fprintf(fp_cmt, "\n");

    if(relax ==1){
      /* Print Taum in ps*/
      fprintf(fp_cmt, "%12s", "Tau_m [ps]");
      for(i = 0; i < nvect; i++){
	fprintf(fp_cmt, "%12.6f", taum[i]* 1.0E12);
      }
      fprintf(fp_cmt, "\n");
    }

    /* Print cf */
    for(i = 0; i < nsteps; i++){
      fprintf(fp_cmt, "%12.8f", (double) i * tstep);
      fprintf(fp_cjt, "%12.8f", (double) i * tstep);
      for(j = 0; j < nvect; j++){
        if(norm){
          fprintf(fp_cmt, "%12.8f", cf[nsteps * j + i] * frame / (cf[nsteps * j] * (frame - i)));
          fprintf(fp_cjt, "%12.8f", cf_cjt[nsteps * j + i]/cf_cjt[nsteps * j]);
        }
        else{
          fprintf(fp_cmt, "%12.8f", factor * cf[nsteps * j + i] / (frame - i));
          fprintf(fp_cjt, "%12.8f", factor * cf_cjt[nsteps * j + i]);
        }
      }
      fprintf(fp_cmt, "\n");
      fprintf(fp_cjt, "\n");
    }
    safe_fclose(fp_cmt);
    safe_fclose(fp_cjt);
    safe_free(cjtfilename);
    safe_free(cmtfilename);
  } // END ATCT_IRED
  else if(type == ATCT_NORMAL){
    fp = safe_fopen(filename, "w");
    fprintf(fp,"%s-correlation functions, %s type\n", 
          (tcmode == ATCM_AUTO? "Auto" :
           (tcmode == ATCM_CROSS? "Cross" : "????"
           )
          ),
          (type == ATCT_IRED? "IRED" :
           (type == ATCT_NORMAL? "normal" : "????"
           )
          ));

    if(dplr){
      dnorm = 1.0 / ((double) frame);
      fprintf(fp, "***** Vector length *****\n");
      fprintf(fp, "%10s %10s %10s %10s\n", "<r>", "<rrig>", "<1/r^3>", "<1/r^6>");
      fprintf(fp, "%10.4f %10.4f %10.4f %10.4f\n",
              rave1 * dnorm,
              sqrt(avgcrd1[0]*avgcrd1[0] + avgcrd1[1]* avgcrd1[1] + avgcrd1[2]*avgcrd1[2]) *dnorm,
              r3iave1 * dnorm,
              r6iave1 * dnorm);
      if(vinfo2 != NULL){
        fprintf(fp, "%10.4f %10.4f %10.4f %10.4f\n",
                rave2 * dnorm,
                sqrt(avgcrd2[0]*avgcrd2[0] + avgcrd2[1]* avgcrd2[1] + avgcrd2[2]*avgcrd2[2]) *dnorm,
                r3iave2 * dnorm,
                r6iave2 * dnorm);
      }
    }

    fprintf(fp, "\n");
    fprintf(fp, "***** Correlation functions *****\n");
    if(dplr)
      fprintf(fp, "%10s %10s %10s %10s\n", "Time", "<C>", "<P2>", "<1/(r^3*r^3)>");
    else
      fprintf(fp, "%10s %10s\n", "Time", "<P2>");
    for(i = 0; i < nsteps; i++){
      if(dplr){
        if(norm){
          fprintf(fp, "%10.3f %10.4f %10.4f %10.4f\n",
                  (double) i * tstep,
                  cf[i]   * frame / (cf[0]   * (frame - i)),
                  p2cf[i] * frame / (p2cf[0] * (frame - i)),
                  rcf[i]  * frame / (rcf[0]  * (frame - i)));
        }
        else{
          fprintf(fp, "%10.3f %10.4f %10.4f %10.4f\n",
                  (double) i * tstep,
                  factor * cf[i]   / (frame - i),
                  factor * p2cf[i] / (frame - i),
                  rcf[i]  / (frame - i));
        }
      }
      else{
        if(norm){
          fprintf(fp, "%10.3f %10.4f\n",
                  (double) i * tstep,
                  p2cf[i] * frame / (p2cf[0] * (frame - i)));
        }
        else{
          fprintf(fp, "%10.3f %10.4f\n",
                  (double) i * tstep,
                  factor * p2cf[i] / (frame - i));
        }
      }
    }
    safe_fclose(fp);
  } // END ATCT_NORMAL

  return 1;
}
#endif // ifndef NO_PTRAJ_ANALYZE

