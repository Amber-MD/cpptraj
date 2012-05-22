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

// analyzeSet()
int analyzeSet(analyzeInformation *analyze, stackType *sp, int mode)
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

