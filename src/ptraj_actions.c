// ---------- CSTDLIB includes -------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

// ---------- Defines ----------------------------------------------------------
#define ACTION_MODULE
// ---------- PTRAJ includes ---------------------------------------------------
#include "ptraj_actions.h"
#include "ptraj_common.h" // scalarInfo
#include "ptraj_stack.h"
#include "ptraj_arg.h"
#include "ptraj_scalar.h"

// ---------- CPPTRAJ includes -------------------------------------------------
// Constants
#include "Constants.h"
// Distance routines
#include "DistRoutines.h"
// Torsion Routines
#include "TorsionRoutines.h"
// MPI worldrank and size
#include "MpiRoutines.h"

// ========== COMMON internal functions ========================================
#ifdef MPI
static void printError(char *actionName, char *fmt, ...) {
  va_list argp;
  va_start(argp, fmt);
//#ifdef MPI
  if (worldrank == 0) {
//#endif    
    printf("WARNING in ptraj(), %s: ", actionName);
    vprintf(fmt, argp);
//#ifdef MPI
  }
//#endif
  va_end(argp);
}
static void printParallelError(char *actionName) {
  printError(actionName, "Parallel implementation of action not supported.\nIgnoring command...\n");
}
#endif

#ifdef MPI
// printMPIerr()
/// Wrapper for MPI_Error string.
static void printMPIerr(int err, char *actionName) {
  int len,eclass,i;
  char buffer[BUFFER_SIZE];
  
  MPI_Error_string(err,buffer,&len);
  MPI_Error_class(err,&eclass);
  // Remove newlines from MPI error string
  for (i=0; i<len; i++) 
    if (buffer[i]=='\n') buffer[i]=':';
  fprintf(stdout,"[%i] MPI ERROR %d: %s: [%s]\n",worldrank,eclass,actionName,buffer);

  return;
}
#endif

/** ACTION ROUTINE *************************************************************
 *
 *  transformRandomizeIons() --- swap positions of ions and solvent randomly
 *
 ******************************************************************************/

   int
transformRandomizeIons(actionInformation *action, 
		       double *x, double *y, double *z, 
		       double *box, int mode)
{
  //char *name = "randomizeions";
  char *buffer;
  double distance, sx, sy, sz;
  int ion, i, j, w;
  int *around;
  int *solvent;
  argStackType **argumentStackPointer;

  double ucell[9], recip[9];

  /*
   *  USAGE:
   *
   *    randomizeions <mask> [around <mask> by <distance>] [overlap <value>] [noimage] [seed <value>]
   *
   *  action argument usage:
   *
   *  mask: the list of ions to be moved.  Each is assumed to be a single atom residue.
   *  iarg1: if 1, disable imaging
   *  iarg2: seed
   *  darg1: the minimum distance between ions (overlap)
   *  darg2: the minimum distance to the around mask
   *  carg1: the around mask (region of space to avoid)
   */

  if (mode == PTRAJ_SETUP) {

    /*
     *  ACTION: PTRAJ_SETUP
     */

#ifdef MPI

#endif

    argumentStackPointer = (argStackType **) action->carg1;
    action->carg1 = NULL;

    if (action->state->solventMolecules == 0) {
      fprintf(stdout, 
	      "WARNING in ptraj(), randomizeions: This command only works if solvent\n");
      fprintf(stdout, "information has been specified.  See the \"solvent\" command.\n");
      fprintf(stdout, "Ignoring this command.\n");
      return -1;
    }

    buffer = getArgumentString(argumentStackPointer, NULL);
    action->mask = processAtomMask(buffer, action->state);
    safe_free(buffer);

    if (action->mask == NULL) {
      fprintf(stdout, "WARNING in ptraj(), randomizeions: NULL mask for the ion specification\n");
      return -1;
    }

    /*
     *  check to see that each ion selected is only a single atom residue!
     */
    for (i=0; i < action->state->atoms; i++) {
      if (action->mask[i]) {

	j = atomToResidue(i,action->state->residues,action->state->ipres);

	if (prnlev > 6) {
	  printf("Atom %i is in residue %i which spans atoms %i to %i\n",
		 i+1, j+1, action->state->ipres[j], action->state->ipres[j+1]);
	}
		 
	if (action->state->ipres[j+1] - action->state->ipres[j] > 1) {
	  fprintf(stdout, 
		  "WARNING IN randomize ions: residue %i appears to contain more than 1 atom!\n",
		  j);
	}
      }
    }



    /*
     *  check the solvent information to make sure that each solvent listed has the
     *  same number of atoms in each molecule; otherwise a uniform trajectory is not
     *  possible and therefore this command will be ignored...
     */

    j = action->state->solventMoleculeStop[0] - action->state->solventMoleculeStart[0];
    for (i=1; i < action->state->solventMolecules; i++) {
      if (j != (action->state->solventMoleculeStop[i] - 
		action->state->solventMoleculeStart[i])) {
	fprintf(stdout, 
		"WARNING in ptraj(), randomizeions: the solvent molecules are not of uniform\n");
	fprintf(stdout, "size hence this command will be ignored.  [Try resetting the solvent\n");
	fprintf(stdout, "information with the \"solvent\" command...\n");
	return -1;
      }
    }

    action->iarg1 = argumentStackContains(argumentStackPointer, "noimage");
    action->iarg2 = argumentStackKeyToInteger(argumentStackPointer, "seed", -1);
    action->darg1 = argumentStackKeyToDouble(argumentStackPointer, "overlap",  3.5);
    action->darg2 = argumentStackKeyToDouble(argumentStackPointer, "by", 3.5);
    action->darg1 = action->darg1 * action->darg1;
    action->darg2 = action->darg2 * action->darg2;

    buffer = argumentStackKeyToString(argumentStackPointer, "around", NULL);
    if (buffer) {
      action->carg1 = processAtomMask(buffer, action->state);
      safe_free(buffer);
    } else
      action->carg1 = NULL;

  } else if (mode == PTRAJ_STATUS) {

    /*
     *  ACTION: PTRAJ_STATUS
     */

    fprintf(stdout, "  RANDOMIZEIONS: swapping the postions of the ions: ");
    printAtomMask(stdout, action->mask, action->state);
    fprintf(stdout, "\n");
    fprintf(stdout, "      with the solvent.  No ions can get closer than %5.2f angstroms to another ion\n",
	    sqrt(action->darg1));
    around = (int *) action->carg1;
    if (around != NULL) {
      fprintf(stdout, "      No ion can get closer than %5.2f angstroms to: ",
	    sqrt(action->darg2));
      printAtomMask(stdout, around, action->state);
      fprintf(stdout, "\n");
    }

    if (action->iarg1) {
      fprintf(stdout, "      Imaging of the coordinates will not be performed\n");
    }
    if (action->iarg2 > 0) {
      fprintf(stdout, "      Random number generator seed is %i\n", action->iarg2);
      srandom((unsigned) action->iarg2);
    }


  } else if (mode == PTRAJ_CLEANUP) {

    action->carg1 = NULL;

  }


  if (mode != PTRAJ_ACTION) return 0;

  /*
   *  ACTION: PTRAJ_ACTION
   */

  if (action->mask == NULL) return 0;


  if (action->iarg1 == 0 && box[3] == 0.0) {
    action->iarg1 = 1;
    fprintf(stdout, "  RANDOMIZEIONS: box angles are zero, disabling imaging!\n");
  }
  if (action->iarg1 == 0 && (box[3] != 90.0 || box[4] != 90.0 || box[5] != 90.0))
    boxToRecip(box, ucell, recip);

  around = (int *) action->carg1;

  solvent = (int *) safe_malloc(sizeof(int) * action->state->solventMolecules);
  memset(solvent, 0, sizeof(int) * action->state->solventMolecules);

  /*
   *  loop over all solvent molecules and mark those that are too close to the solute
   */
  for (j=0; j < action->state->solventMolecules; j++) {

    solvent[j] = 1; 

    /*
     *  is solvent molecule to near any atom in the around mask?
     */
    if (around != NULL) {

      for (i=0; i < action->state->atoms; i++) {
	if ( around[i] ) {
	  if (action->state->solventMoleculeStart[j] != i) {
	    distance = calculateDistance2(action->state->solventMoleculeStart[j], i, x, y, z, 
					  box, (double *) ucell, (double *) recip, 0.0, action->iarg1);
	    if (distance < action->darg2) {
	      solvent[j] = 0;
	      if (prnlev > 6) {
		fprintf(stdout, "  RANDOMIZEIONS: water %i is only %5.2f angstroms from atom %i\n",
			j+1, sqrt(distance), i+1);
		
	      }
	      i = action->state->atoms;
	    }
	  }
	}
      }
    }
  }

  if (prnlev > 4) {
    i = 0;
    if (prnlev > 6)
      fprintf(stdout, "RANDOMIZEIONS: The following waters are ACTIVE so far:\n");
    for (j=0; j < action->state->solventMolecules; j++) {
      if (solvent[j]) {
	i++;
	if (prnlev > 6) {
	  fprintf(stdout, " %5i ", j+1);
	  if (i%10 == 0) printf("\n");
	}
      }
    }
    fprintf(stdout, "  RANDOMIZEIONS: A total of %i waters (out of %i) are active\n", i, action->state->solventMolecules);
  }

  /*
   *  loop over all ions
   */
  for (ion=0; ion < action->state->atoms; ion++) {

    if (action->mask[ion]) {

      if (prnlev > 2) fprintf(stdout, "  RANDOMIZEIONS: Processing ion atom %i\n", ion+1);

	/* 
	 *  is a potential solvent molecule close to any of the ions (except this one)?
	 */
      for (j=0; j < action->state->solventMolecules; j++) {
	if ( solvent[j] ) {

	  /*
	   *  if this solvent is active, check distance to all other ions
	   */
	  for (i=0; i < action->state->atoms; i++) {
	  
	    if (action->mask[i] && ion != i) {
	      distance = calculateDistance2(action->state->solventMoleculeStart[j], i, x, y, z, 
					    box, (double *) ucell, (double *) recip, 0.0, action->iarg1);
	      if (distance < action->darg1) {
		i = action->state->atoms;
		solvent[j] = 0;
		if (prnlev > 6) {
		  fprintf(stdout, "  RANDOMIZEIONS: water %i is only %5.2f angstroms from (ion) atom %i\n",
			  j+1, sqrt(distance), i+1);
		}
	      }
	    }
	  }
	}
      }

      i = 1;
      while (i > 0 && i < 10000) {
	/*
	 *  Run the random number generator so that the same number is not produced
	 *  when the seed was set manually
	 */
#ifdef MPI	
	for (j = 0; j < worldsize; j++) {
	  if (j == worldrank)
	    w = random() % action->state->solventMolecules;
	  else
	    random();
	}
#else
	w = random() % action->state->solventMolecules;
#endif
	if (solvent[w] == 1) {
	  i = -1;
	} else {
	  i++;
	}
      }

      if (i > 0) {
	fprintf(stdout, "  RANDOMIZEIONS: warning tried 10000 random waters and couldn't meet criteria!  Skipping\n");
      }

      if (i < 0) {
	if (prnlev > 2) {
	  fprintf(stdout, "  RANDOMIZEIONS: Swaping solvent %i for ion %i\n",
		  w+1, ion+1);
	}


	i = action->state->solventMoleculeStart[w];
	sx = x[ion] - x[i];
	sy = y[ion] - y[i];
	sz = z[ion] - z[i];
	
	for (i = action->state->solventMoleculeStart[w]; i < action->state->solventMoleculeStop[w]; i++) {

	  x[i] += sx;
	  y[i] += sy;
	  z[i] += sz;

	}
	x[ion] -= sx;
	y[ion] -= sy;
	z[ion] -= sz;

      }
    }
  }
  safe_free(solvent);
  return 1;
}

/** ACTION ROUTINE *************************************************************
 *
 *  transformScale() --- Scale the coordinates by a specified amount
 *
 ******************************************************************************/


   int
transformScale(actionInformation *action, 
		   double *x, double *y, double *z,
		   double *box, int mode)
{
  //char *name = "scale";
  argStackType **argumentStackPointer;
  char *buffer;
  int i;

  /*
   *  USAGE:
   *
   *    scale [x <scalex>] [y <scaley>] [z <scalez>] [mask]
   *
   *  action argument usage:
   *
   *  mask : atom selection representing atoms to shift
   *  darg1: scalex
   *  darg2: scaley
   *  darg3: scalez 
   */

  if (mode == PTRAJ_SETUP) {
    /*
     *  ACTION: PTRAJ_SETUP
     */

#ifdef MPI
    printParallelError(name);
    return -1;
#endif

    argumentStackPointer = (argStackType **) action->carg1;
    action->carg1 = NULL;

    action->darg1 = argumentStackKeyToDouble(argumentStackPointer, "x", 0.0);
    action->darg2 = argumentStackKeyToDouble(argumentStackPointer, "y", 0.0);
    action->darg3 = argumentStackKeyToDouble(argumentStackPointer, "z", 0.0);

    buffer = safe_malloc(sizeof(char) * BUFFER_SIZE);
    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL) {
      action->mask = NULL;
    } else
      action->mask = processAtomMask(buffer, action->state);
    safe_free(buffer);

  } else if (mode == PTRAJ_STATUS) {

    /*
     *  ACTION: PTRAJ_STATUS
     */

    fprintf(stdout, "  SCALE coordinates: ");
    if (action->darg1 != 0.0) fprintf(stdout, "X by %.3f ", action->darg1);
    if (action->darg2 != 0.0) fprintf(stdout, "Y by %.3f ", action->darg2);
    if (action->darg3 != 0.0) fprintf(stdout, "Z by %.3f ", action->darg3);
    if (action->mask != NULL) {
      fprintf(stdout, " mask is ");
      printAtomMask(stdout, action->mask, action->state);
      fprintf(stdout, "\n");
    }
    if (action->mask == NULL) fprintf(stdout, "\n");

  }

  if (mode != PTRAJ_ACTION) return 0;

  /*
   *  ACTION: PTRAJ_ACTION
   */

  for (i=0; i < action->state->atoms; i++) 
    if (action->mask == NULL || action->mask[i]) {
      x[i] *= action->darg1;
      y[i] *= action->darg2;
      z[i] *= action->darg3;
    }

  return 1;

}

/** ACTION ROUTINE *************************************************************
 *
 *  transformTruncOct() --- trim/orient a box to make it a truncated octahedron
 *
 ******************************************************************************/
/* DISABLED FOR CPPTRAJ - relies on too much ingrained parm data in ptraj

   int
transformTruncOct(actionInformation *action, 
		  double *x, double *y, double *z,
		  double *box, int mode)
{
  char *name = "truncoct";
  argStackType **argumentStackPointer;
  char *buffer;
  ptrajState *state;
  int ii, i, j;
  int *mask;
  double cx, cy, cz;
  double total_mass, max_dist, dist;
  double sideDist,sideDist0,diagDist,diagDist2,diagCoord;
  double toDist, pdist,dnormCoord;
  int *zapMask, zap, zapTotal;
  double phi,cos1,sin1,cos2,sin2,tetra_angl;
  double t11,t12,t13,t21,t22,t23,t31,t32,t33,xx,yy;
  double ucell1[3];
  double ucell2[3];
  double ucell3[3];
  double gamma;
  int new_waters;
  //Parm *newparm, *tmpparm;
  FILE *fpout;


  //  USAGE:
  //
  //    truncoct <mask> <distance> prmtop <filename>
  //
  //  action argument usage:
  //
  //  mask: atom selection for solute
  //  iarg1: the index of the first solvent molecule
  //  darg1: the size of the truncated octahedron(?)

  if (mode == PTRAJ_SETUP) {
    //  ACTION: PTRAJ_SETUP

#ifdef MPI
    printParallelError(name);
    return -1;
#endif

    argumentStackPointer = (argStackType **) action->carg1;
    action->carg1 = NULL;

    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL) {
      fprintf(stdout, "WARNING in ptraj(), truncoct: No atom mask for the solute was\n");
      fprintf(stdout, "specified...  Ignoring command.\n");
      return -1;
    }
    action->mask = processAtomMask(buffer, action->state);
    safe_free(buffer);

    action->darg1 = getArgumentDouble(argumentStackPointer, -1.0);
    if (action->darg1 < 0) {
      fprintf(stdout, "WARNING in ptraj(), truncoct: The buffer distance specified is\n");
      fprintf(stdout, "out of range or was not specified.  Ignoring command\n");
      return -1;
    } 

    if (parm == NULL) {
      fprintf(stdout, "WARNING in ptraj(), truncoct: No AMBER prmtop file is present\n");
      fprintf(stdout, "This command only works with AMBER prmtop files, hence ignoring\n");
      fprintf(stdout, "command...\n");
      return -1;
    }

    if (action->state->solventMolecules == 0) {
      fprintf(stdout, "WARNING in ptraj(), truncoct: No solvent information has been\n");
      fprintf(stdout, "specified.  See the \"solvent\" command.  Ignoring...\n");
      return -1;
    }

    buffer = argumentStackKeyToString(argumentStackPointer, "prmtop", NULL);
    action->carg1 = (void *) buffer;


  } else if (mode == PTRAJ_STATUS) {

    //  ACTION: PTRAJ_STATUS

    fprintf(stdout, 
	    "  TRUNCATED OCTAHEDRON: will be created with minimum distance from solute\n");
    fprintf(stdout, "      to the sides of the truncated octahedron of %.3f angstroms\n",
	    action->darg1);
    buffer = (char *) action->carg1;
    if (buffer != NULL) {
      fprintf(stdout, "      Creating a prmtop named: %s\n", buffer);
    }

    fprintf(stdout, "      The solute mask is ");
    printAtomMask(stdout, action->mask, action->state);
    fprintf(stdout, "\n");

  }

  if (mode != PTRAJ_ACTION) return 0;

  //  ACTION: PTRAJ_ACTION

  state = (ptrajState *) action->state;

  //  update local state information
  for (i=0; i<6; i++)
    state->box[i] = box[i];

  //  FIRST CENTER of geometry of the solute at origin
  //  
  //  accumulate center of geometry...

  mask = action->mask;
  cx = 0.0;
  cy = 0.0;
  cz = 0.0;
  
  total_mass=0.;
  printf("\n***********************************************************\n");
  printf(  "*********  Truncated Octahedral Data                *******\n");
  printf(  "*********                                           *******\n");
  printf(  "***********************************************************\n");
  for (i=0; i < state->atoms; i++) {
      if (mask[i]) {
	  cx += x[i];
	  cy += y[i];
	  cz += z[i];
	  total_mass += 1.0;
      }
  }

  cx /= total_mass;
  cy /= total_mass;
  cz /= total_mass;
  max_dist=0.;
  printf("Center of geometry Offset     %lf %lf %lf\n",cx,cy,cz);
  for (i=0; i < state->atoms; i++) {
      x[i] -= cx;
      y[i] -= cy;
      z[i] -= cz;
      if (mask[i]) {
	  dist=x[i]*x[i]+y[i]*y[i]+z[i]*z[i];
	  if(dist>max_dist) max_dist=dist;	  
      }
  }
  max_dist=sqrt(max_dist);
  printf("max radius 0f solute is %lf\n",max_dist);

//     calculate the face distances
  toDist=action->darg1;
// printf("\n\nInside TruncOct toDist is %lf\n\n",toDist);
  diagDist=toDist+max_dist-0.5;
  diagCoord=diagDist/sqrt(3);
  dnormCoord=1./sqrt(3.);
  sideDist=diagCoord* 2.;
  sideDist0=diagCoord* 2.-0.5;
  if(sideDist > state->box[0]*0.5 || 
     sideDist > state->box[1]*0.5 ||
     sideDist > state->box[2]*0.5){
      printf("\nWARNING WARNING WARNING in truncoct: ");
      printf("Original box MAY not be big enough\n");
      printf("           ...... Continuing anyway ......\n\n");
  }
  printf("   TO cubic faces have dist %f while \n    orig. box sizes are %f %f %f\n\n",
	 2.*sideDist,state->box[0],state->box[1],state->box[2]);
// printf("Inside TruncOct side and diag are %lf %lf\n\n",sideDist,diagDist);

//    start removing solvent
//
//    NOTE: now we only keep track of solvent molecules we want to remove!

  zapTotal = 0;
  zapMask = (int *) safe_malloc(sizeof(int) * state->atoms);
  for (i=0; i < state->atoms; i++)
    zapMask[i] = 0;

  for(i=0; i < state->solventMolecules; i++) {

    zap = 0;
    for(ii=state->solventMoleculeStart[i]; ii < state->solventMoleculeStop[i]; ii++) {
      if(ABS(x[ii])>sideDist0 || ABS(y[ii])>sideDist0 || ABS(z[ii])>sideDist0) {
	zap=1;
      } else {
	pdist=(ABS(x[ii])+ABS(y[ii])+ABS(z[ii])-3.*diagCoord)*dnormCoord;
	if(pdist>0){
	  zap=1;
	}
      }
    }
    if (zap == 1) {
      for(ii=state->solventMoleculeStart[i]; ii < state->solventMoleculeStop[i]; ii++) {
	zapMask[ii] = 1;
      }
      zapTotal++;
    }
  }


  //  modify coordinates

  j = 0;
  for (i=0; i < state->atoms; i++) {
    if (zapMask[i] == 0) {
      x[j] = x[i];
      y[j] = y[i];
      z[j] = z[i];
      j++;
    }
  }

  //  modify the current state
  action->state = NULL;

  if (prnlev > 2) {
    printf("ZAP MASK IS: ");
    printAtomMask(stdout, zapMask, state);
    fprintf(stdout, "\n");
  }
  modifyStateByMask(&action->state, &state, zapMask, 1);
    
  safe_free(zapMask);
  zapMask = NULL;

  ptrajClearState(&state);
  state = action->state;

  new_waters=state->solventMolecules;
  printf("Number of waters in TO %d\n",new_waters);

  // *******************************************************
  //       Now Rotate the whole thing to line up the axes *
  // *******************************************************
  tetra_angl=2*acos(1./sqrt(3.));
  phi=PI/4.;
  cos1=cos(phi);
  sin1=sin(phi);
  phi=PI/2.-tetra_angl/2.;
  cos2=sqrt(2.)/sqrt(3.);
  sin2=1./sqrt(3.);
  
  // *****************************************************   
  //       45 around z axis, (90-tetra/2) around y axis,   
  //       90 around x axis                                 
  //                                                        
  //   (1  0  0)    (cos2  0 -sin2)    (cos1 -sin1  0)      
  //   (0  0 -1)    (   0  1     0)    (sin1  cos1  0)      
  //   (0  1  0)    (sin2  0  cos2)    (   0     0  1)      
  //                                                        
  //   cntr-clk       clock              clock              
  //   Looking down + axis of rotation toward origin        
  // *****************************************************   

  t11= cos2*cos1;
  t12=-cos2*sin1;
  t13=-sin2;
  t21=-sin2*cos1;
  t22= sin2*sin1;
  t23=-cos2;
  t31= sin1;
  t32= cos1;
  t33=0;

  for (i=0; i < state->atoms; i++) {
      xx = t11*x[i]+t12*y[i]+t13*z[i];
      yy = t21*x[i]+t22*y[i]+t23*z[i];
      z[i] = t31*x[i]+t32*y[i]+t33*z[i];
      x[i]=xx;
      y[i]=yy;
  }
  diagDist2=2.*diagDist;
  gamma=tetra_angl;
  ucell1[0] = diagDist2;
  ucell1[1] = 0.;
  ucell1[2] = 0.;
  ucell2[0] = diagDist2*cos(gamma);
  ucell2[1] = diagDist2*sin(gamma);
  ucell2[2] = 0.;
  ucell3[0] = diagDist2*cos(gamma);
  ucell3[1] = (diagDist2*diagDist2*cos(gamma)-
		 ucell3[0]*ucell2[0])/ucell2[1];
  ucell3[2] = sqrt( diagDist2*diagDist2 - ucell3[0]*ucell3[0] - 
		      ucell3[1]*ucell3[1] );

  if (prnlev > 2) {
    printf("TRUNCATED OCTAHEDRON GENERATION:\n");
    printf("UCELL %f %f %f \n",ucell1[0],ucell1[1],ucell1[2]);
    printf("UCELL %f %f %f \n",ucell2[0],ucell2[1],ucell2[2]);
    printf("UCELL %f %f %f \n",ucell3[0],ucell3[1],ucell3[2]);
  }
  printf("UCELL length for mdin file is %f, padded by 1.0 angstrom\n",ucell1[0]);

  state->box[0]=ucell1[0] + 1.0;
  state->box[1]=ucell1[0] + 1.0;
  state->box[2]=ucell1[0] + 1.0;
  state->box[3]=tetra_angl*RADDEG;
  state->box[4]=state->box[3];
  state->box[5]=state->box[3];
  for (i=0; i < 6; i++)
    box[i] = state->box[i];

  //  Dump out a new prmtop file if requested.
  buffer = (char *) action->carg1;
  action->carg1 = NULL;

  if (buffer) {
    fprintf(stdout,"Warning: truncoct Parm write disabled for Cpptraj.\n");
//
//  if ( (fpout=safe_fopen(buffer,"w")) == NULL ) {
//    fprintf(stdout, "WARNING in ptraj(), truncoct: Couldn't open prmtop file %s\n",
//            buffer);
//    return 1;
//  }

//  tmpparm = parm;
//  if (new_waters == 0) {
//    fprintf(stdout, "WARNING in ptraj(), truncoct: No waters were removed...\n");
//  } else {
//    newparm = modifyTIP3P(new_waters);
//    parm = newparm;
//  }
//  parm->IFBOX=2;
//  parm->box->beta  = state->box[4];
//  parm->box->box[0]= state->box[0];
//  parm->box->box[1]= state->box[0];
//  parm->box->box[2]= state->box[0];
//  writeParm( fpout, 1 ); 
//  parm = tmpparm;
//  safe_fclose(fpout);
//  safe_free(buffer);
  }  
  return 1;
}
*/

