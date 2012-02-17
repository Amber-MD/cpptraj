// ---------- CSTDLIB includes -------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
//#include <time.h> // for cluster
//#include <float.h>
//#include <ctype.h>

// ---------- Defines ----------------------------------------------------------
#define ACTION_MODULE
       /* Multiply the transpose of the 3x3 matrix times the 
        * coordinates specified in x, y and z.  xx, yy and zz 
        * are temporary variables. 
        * The x, y and z arrays are modified!
        */
#define VOP_3x3_TRANSPOSE_TIMES_COORDS(matrix, x, y, z, xx, yy, zz) \
  xx = matrix[0][0] * x + matrix[1][0] * y + matrix[2][0] * z;  \
  yy = matrix[0][1] * x + matrix[1][1] * y + matrix[2][1] * z;  \
  zz = matrix[0][2] * x + matrix[1][2] * y + matrix[2][2] * z;  \
  x = xx; y = yy; z = zz


#define VOP_3x3_TIMES_3x3(t, u, v) \
  t[0][0] = u[0][0] * v[0][0] + u[0][1] * v[1][0] + u[0][2] * v[2][0]; \
  t[0][1] = u[0][0] * v[0][1] + u[0][1] * v[1][1] + u[0][2] * v[2][1]; \
  t[0][2] = u[0][0] * v[0][2] + u[0][1] * v[1][2] + u[0][2] * v[2][2]; \
  t[1][0] = u[1][0] * v[0][0] + u[1][1] * v[1][0] + u[1][2] * v[2][0]; \
  t[1][1] = u[1][0] * v[0][1] + u[1][1] * v[1][1] + u[1][2] * v[2][1]; \
  t[1][2] = u[1][0] * v[0][2] + u[1][1] * v[1][2] + u[1][2] * v[2][2]; \
  t[2][0] = u[2][0] * v[0][0] + u[2][1] * v[1][0] + u[2][2] * v[2][0]; \
  t[2][1] = u[2][0] * v[0][1] + u[2][1] * v[1][1] + u[2][2] * v[2][1]; \
  t[2][2] = u[2][0] * v[0][2] + u[2][1] * v[1][2] + u[2][2] * v[2][2]

// ---------- PTRAJ includes ---------------------------------------------------
#include "ptraj_actions.h"
#include "ptraj_common.h" // scalarInfo
#include "ptraj_stack.h"
#include "ptraj_arg.h"
#include "ptraj_scalar.h"
//typedef struct _arrayType {
//  int length;
//  void *entry;
//} arrayType;
//# inc lude "../../ptraj/clusterLib.h"
//# inc lude "../../ptraj/cluster.h"

// ---------- CPPTRAJ includes -------------------------------------------------
// Constants
#include "Constants.h"
// Distance routines
#include "DistRoutines.h"
// Torsion Routines
#include "TorsionRoutines.h"
// MPI worldrank and size
#include "MpiRoutines.h"

// ========== TYPE Definitions =================================================
// coodType - Ptraj coordinate types
typedef enum _coordType {
  COORD_UNKNOWN,
  COORD_AMBER_TRAJECTORY,
  COORD_AMBER_RESTART,
  COORD_AMBER_NETCDF,
  COORD_AMBER_REMD,
  COORD_PDB,
  COORD_BINPOS,
  COORD_CHARMM_TRAJECTORY
} coordType;
// CoordinateInfo
typedef struct _coordinateInfo {
  FILE *file;      // File pointer
#ifdef MPI
  MPI_File *mfp;
#endif
  char *filename;  // File name
  int start;       // Frame to start processing
  int stop;        // Frame to end processing
  int Nframes;     // Total number of frames in the file.
  int offset;      // # of frames to skip
  int append;      // File will be appended to
  int isBox;       // File has box information
  int isVelocity;  
  int option1;
  int option2;
  int *mask;
  double *x;
  double *y;
  double *z;
  double *time;
  double *vx;
  double *vy;
  double *vz;
  char *title;
  char *application;
  char *program;
  char *version;
  void *info;          // Holds NETCDF or CHARMM trajectory info
//  netcdfTrajectoryInfo *NCInfo;  // Holds NETCDF trajectory info
  coordType type;      // Identify coordinate type
  int accessMode;      // 0 for read, 1 for write, 2 for append
  int compressType;    // 1 gzip 2 bzip 3 zip 
     //  LES information
/*  LesAction les_action;
  int nlescopy;
  LesStatus les_status;*/
     //  REMD Trajectory info
  struct _coordinateInfo **REMDtraj; // Hold information of other replica trajectories
  int isREMDTRAJ;       // 0 = normal trajectory, 1 = replica trajectory
  char *compressEXT;    /* If not NULL, contains compressed traj ext.*/
  int numREMDTRAJ;      /* How many replica trajectories are present */
  int firstREMDTRAJ;    /* Index of first replica                    */
  int EXTwidth;         /* Length of replica extension               */
  char* baseFilename;   /* To hold replica traj base filename        */
  int linesperset;      /* for fast scanthrough of other REMD files  */
  double remdtrajtemp;  /* Target temperature                        */
  /* Write file with MPI IO? */
  int isMPI;
  int isNetcdf;
   //  AMBER trajectory file information
  int seekable;
  int titleSize;
  int frameSize;
  int numBox;
  char *buffer;
} coordinateInfo;

#define INITIALIZE_coordinateInfo(_p_) \
  _p_->file = NULL; \
  _p_->filename = NULL; \
  _p_->start = 1; \
  _p_->stop = -1; \
  _p_->Nframes = 0; \
  _p_->offset = 1; \
  _p_->append = 0; \
  _p_->isBox = 0; \
  _p_->isVelocity = 0; \
  _p_->option1 = 0; \
  _p_->option2 = 0; \
  _p_->mask = NULL; \
  _p_->frameSize = 0; \
  _p_->x = NULL; \
  _p_->y = NULL; \
  _p_->z = NULL; \
  _p_->time = NULL; \
  _p_->vx = NULL; \
  _p_->vy = NULL; \
  _p_->vz = NULL; \
  _p_->title = NULL; \
  _p_->application = NULL; \
  _p_->program = NULL; \
  _p_->version = NULL; \
  _p_->info = NULL; \
  _p_->type = COORD_UNKNOWN; \
  _p_->accessMode = 0; \
  _p_->compressType = 0; \
  _p_->REMDtraj = NULL; \
  _p_->isREMDTRAJ = 0; \
  _p_->compressEXT = NULL; \
  _p_->numREMDTRAJ = 0; \
  _p_->firstREMDTRAJ = 0; \
  _p_->EXTwidth = 0; \
  _p_->baseFilename = NULL; \
  _p_->linesperset = 0; \
  _p_->remdtrajtemp = 0.0; \
  _p_->isMPI = 0; \
  _p_->isNetcdf = 0; \
  _p_->seekable = 0; \
  _p_->titleSize = 0; \
  _p_->frameSize = 0; \
  _p_->numBox = 0; \
  _p_->buffer = NULL; \

// ACTION: TRANSFORM_GRID, TRANSFORM_DIPOLE
typedef struct _transformGridInfo {
  double dx;
  double dy;
  double dz;
  int nx;
  int ny;
  int nz;
  int frames;
  float *grid;
  float *dipolex;
  float *dipoley;
  float *dipolez;
  char *filename;
} transformGridInfo;

#define INITIALIZE_transformGridInfo(_p_) \
  _p_->dx        = 0.0;  \
  _p_->dy        = 0.0;  \
  _p_->dz        = 0.0;  \
  _p_->nx        = 0;    \
  _p_->ny        = 0;    \
  _p_->nz        = 0;    \
  _p_->frames    = 0;    \
  _p_->grid      = NULL; \
  _p_->dipolex   = NULL; \
  _p_->dipoley   = NULL; \
  _p_->dipolez   = NULL; \
  _p_->filename  = NULL

// ========== GLOBAL variables =================================================
// prnlev - used to set overall debug level
//int prnlev;
// worldrank and worldsize - for MPI code
//const int worldrank = 0;
//const int worldsize = 1;
// from ptraj.h
//arrayType *attributeArray = NULL;
//arrayType *attributeArrayTorsion = NULL;
//stackType *vectorStack = NULL;
//stackType *matrixStack = NULL;
// Ref coords
coordinateInfo *referenceInfo = NULL;

// ========== Functions that access GLOBALs ====================================
// SetReferenceInfo()
// It seems that the functions that require a reference only need 
// x y and z set up.
void SetReferenceInfo(double *X, int natom) {
  int i3 = 0;
  int atom;

  referenceInfo = (coordinateInfo*) malloc(sizeof(coordinateInfo));
  INITIALIZE_coordinateInfo(referenceInfo);
  referenceInfo->x = (double*) malloc(natom * sizeof(double));
  referenceInfo->y = (double*) malloc(natom * sizeof(double));
  referenceInfo->z = (double*) malloc(natom * sizeof(double));
  for (atom = 0; atom < natom; atom++) {
    referenceInfo->x[atom] = X[i3++];
    referenceInfo->y[atom] = X[i3++];
    referenceInfo->z[atom] = X[i3++];
  }
}
// FreeReferenceInfo()
void FreeReferenceInfo() {
  if (referenceInfo!=NULL) {
    safe_free(referenceInfo->x);
    safe_free(referenceInfo->y);
    safe_free(referenceInfo->z);
    safe_free(referenceInfo);
    referenceInfo=NULL;
  }
}

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
// ptrajfprintf()
/// printf wrapper: If MPI use MPI_Write_ordered, otherwise use normal system calls.
static void ptrajfprintf(void *fp, char *fmt, ...) {
  va_list argp;
  va_start(argp, fmt);
#ifdef MPI
  char buffer[BUFFER_SIZE];
  int err;

  if (fp == stdout) 
    vprintf(fmt, argp);
  else if (fp == stderr) 
    vfprintf(stderr, fmt, argp);
  else {
    vsprintf(buffer, fmt, argp);
    err=MPI_File_write_ordered(*((MPI_File *) fp), buffer, strlen(buffer), MPI_CHAR, MPI_STATUS_IGNORE);
    if (err!=MPI_SUCCESS) printMPIerr(err,"ptrajfprintf"); 
  }
#else
  vfprintf((FILE *) fp, fmt, argp);
#endif
  va_end(argp);
  return;
}
// ptrajfprintfone()
/// Similar to ptrajfprintf but only print on master
static void ptrajfprintfone(void *fp, char *fmt, ...) {
  va_list argp;
  va_start(argp, fmt);
#ifdef MPI
  if (worldrank > 0) return;
  if (fp == stdout) {
    vprintf(fmt, argp);
    return;
  } else if (fp == stderr) {
    vfprintf(stderr, fmt, argp);
    return;
  }

  char buffer[BUFFER_SIZE];
  vsprintf(buffer, fmt, argp);
  MPI_File_write_shared(* ((MPI_File *) fp), buffer, strlen(buffer), MPI_CHAR, MPI_STATUS_IGNORE);
#else
  vfprintf((FILE *) fp, fmt, argp);
#endif
  va_end(argp);
}
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

// ========== STRING functions =================================================
/*
static char *toLowerCase( char *string_in ) {
  int i, length;

  length = strlen( string_in );
  for (i=0; i < length; i++ ) {
    string_in[i] = tolower( string_in[i] );
  }
  return string_in;
}
static int stringMatch(char *string, char *match) {
  char *lowerString, *lowerMatch;

  if (string == NULL || match == NULL) return 0;

  lowerString = safe_malloc( (size_t) sizeof(char) * (strlen(string)+1));
  strcpy(lowerString,string);
  lowerMatch  = safe_malloc( (size_t) sizeof(char) * (strlen(match)+1));
  strcpy(lowerMatch,match);
  lowerString = toLowerCase(lowerString);
  lowerMatch  = toLowerCase(lowerMatch);

   // this should be an exact, not partial match
   //    if (strncmp(lowerMatch, lowerString, strlen(lowerMatch)) == 0) {
  if (strcmp(lowerMatch, lowerString) == 0) {
    safe_free(lowerString);
    safe_free(lowerMatch);
    return 1;
  } else {
    safe_free(lowerString);
    safe_free(lowerMatch);
    return 0;
  }
}
*/
// ========== FILE IO functions ================================================
// ptrajOpenW()
/// Open a file using MPI calls if MPI defined, or regular system call if not.
static void *ptrajOpenW( char *filename ) {
#ifdef MPI
  MPI_File *fp;
  int err,errtotal;

  fp = (MPI_File*) safe_malloc(sizeof(MPI_File));
  err=MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, fp);
  if (err!=MPI_SUCCESS) printMPIerr(err,"ptrajOpenW");
  /* Check that all threads were able to open the file. If not, all will exit.
   * err is local error, errtotal is global error.
   */
  errtotal=0;
  MPI_Allreduce(&err,&errtotal,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if (errtotal>0) {
    safe_free(fp);
    return NULL;
  }
#else
  FILE *fp;
  fp = safe_fopen(filename, "w");
#endif  
  return (void *) fp;
}
// ptrajCloseFile()
/// Close an open file handle via appropriate system call.
// DAN ROE: The get position and set size routines end up causing an error,
// probably because they are supposed to be collective and are only called
// by rank 0. Removed for now.
static void ptrajCloseFile(void *fp) {
#ifdef MPI
/*  MPI_Offset offset;
  if (worldrank == 0) {
    MPI_File_get_position_shared(* ((MPI_File *) fp), &offset);
    MPI_File_set_size(* ((MPI_File *) fp), offset);
  }*/
  MPI_File_close((MPI_File *) fp);
  // DAN ROE - Free the memory used by the file pointer
  safe_free( fp );
#else
  safe_fclose((FILE *) fp);
#endif
}

// =============================================================================
/*
 *  The code in this file implements the various "actions" of ptraj to
 *  perform various coordinate manipulations or analyses.
 *
 *  The following routines (along with supplemental routines as necessary)
 *  are defined:
 *
 *    actionTest               --- the prototypical action routine (see comments below)
 *    transformAtomicFluct     --- compute atomic positional fluctuations or B-factors
 *    transformCheckOverlap    --- check for close contacts
 *    transformContacts        --- calculate number of nearest neighbors of an atom 
 *                                 (Holger Gohlke, JWGU)
 *    transformCorr            --- perform correlation analysis (Vickie Tsui, Scripps)
 *    transformDihedralCluster --- calculate torsion/dihedral angles (Dan Roe)
 *    transformDiffusion       --- calculate mean squared displacements vs. time
 *    transformDipole          --- bin dipoles (Jed Pitera, UCSF)
 *    transformDNAiontracker   --- track ions a' la Hambelberg/Wilson/Williams
 *    transformGrid            --- grid atomic densities
 *    transformMatrix          --- calculate matrices of covariances, cross-correlations, 
 *                                 distances (Gohlke)
 *    transformPrincipal       --- align coordinates along the principal axis
 *    transformScale           --- scale the coordinates by a specified amount
 *    transformStrip           --- strip coordinates
 *    transformTruncOct        --- trim/orient a box to make it a truncated octahedron
 *    transformUnwrap          --- unwrap atoms (InSuk Joung, Rutgers)
 *    transformVector          --- compute/store various vector quantities (IRED, 
 *                                 CORR, CORRIRED vectors: Gohlke)
 *    transformWatershell      --- calculate the number of waters in a given shell
 *
 *  Each of these routines performs its function by being called in a series of
 *  modes defined as follows:
 *
 *    PTRAJ_SETUP   --- perform initialization,
 *                      parse arguments off the argumentStack,
 *                      complete setup of the actionInformation structure.
 *                      NOTE: this works in conjunction with ptrajSetup()
 *                      in ptraj.c
 *
 *    PTRAJ_STATUS  --- print useful information about this action, such as 
 *                      a summary of the arguments.  This is called at
 *                      startup in ptraj() to print a summary.
 *
 *    PTRAJ_CLEANUP --- clean up this action, freeing any associated memory.
 *                      This mode is applied upon error detection or at the
 *                      and of the current round of trajectory processing.
 *
 *    PTRAJ_ACTION  --- perform the actual action on this set of coordinates.
 *                      This is called repeatedly for each set of coordinates.
 *
 *    PTRAJ_PRINT   --- Print out any data as necessary; this is currently called
 *                      after all of the coordinates have been processed.
 *
 *  The most complicated mode, and the only mode not directly called by
 *  ptraj() is the PTRAJ_SETUP mode, which involves a detailed obfuscation
 *  from user input file processing to the actual setup of the "transformActionStack"
 *  which contains the list of actions to be performed on each coordinate set.
 *
 *  This obfuscation is as follows:
 *
 *  dispatchToken() searches the ptrajTokenlist (both defined in dispatch.c
 *  for a match to a trigger string typed by the user.  If this is found, the 
 *  remaining text typed on the line is placed into a stack of white space
 *  separated strings called the argumentStack.  The subroutine ptrajSetup() 
 *  defined in ptraj.c is then called.  This routine allocates and initializes
 *  an actionInformation structure, sets this "action" structure to point to the
 *  appropriate function (defined in this file), places the argumentStack into
 *  the complex argument 1 slot of the action structure (action->carg1), and calls
 *  the function in the PTRAJ_SETUP mode.  In this mode, it is the functions 
 *  responsibility to parse the arguments off of the argument stack and to finish
 *  setup of the actionInformation structure.  Upon successful return (i.e. -1 is not
 *  returned), ptrajSetup() places this action onto the stack of actions that will
 *  be processed for every coordinate set, the globally accessible 
 *  "transformActionStack".
 *
 *  The other modes are called during coordinate processing by traversing 
 *  though this stack and executing each function.  This is performed in
 *  ptraj.c as follows:
 *
 *    for (actionStackTemp = transformActionStack;
 *         actionStackTemp != NULL;
 *	   actionStackTemp = actionStackTemp->next) {
 *
 *        action = (actionInformation *) actionStackTemp->entry;
 *        if (action->type != TRANSFORM_NOOP  &&
 *            action->type != TRANSFORM_TRANSFORM) {
 *
 *           action->fxn(action, X, Y, Z, box, PTRAJ_STATUS);
 *        }
 *    }
 *
 *  By reading the prototypical "action" routine actionTest below, or any of the
 *  other transformXXX routines defined in this file, this should all become
 *  clearer...
 *
 *  So, in summary if you would like to add a new "action" that involves
 *  processing of coordinate sets:
 *
 *   (1) add an entry to the ptrajTokenlist structure in dispatch.c using
 *       the same format as the current entries.  Set the expected number
 *       of arguments negative (to create the argumentStack) and set the
 *       function to be called on dispatch, i.e. ptrajSetup() in most cases
 *       (not the function you are defining).  You will also define an integer 
 *       enumerated type to correspond to this function (i.e. TRANSFORM_ACTION); 
 *       this will require modifying the actionType enumerated type defined in 
 *       actions.h 
 *
 *   (2) Modify ptrajSetup() in ptraj.c to recognize the new actionType and
 *       to set the action->fxn to point to the new routine you are defining.
 *
 *   (3) Write the actual routine with the same argument structure as the
 *       other action routines (i.e. action, X, Y, Z, box, mode) and 
 *       minimally setup code to handle the PTRAJ_SETUP mode 
 *       (for argumentStack processing) and the PTRAJ_ACTION mode for
 *       working on the coordinates.  See actionTest() below for detailed
 *       comments.
 *
 *  NOTES
 *
 *   o  transformDiffusion is currently only setup to image distances in
 *      orthorhombic lattices!!!
 */

/** ACTION ROUTINE *************************************************************
 *
 *  actionTest()   --- the prototypical action routine
 *
 *  This routine demonstrates how an action routine is used and information
 *  can be obtained from the actionInformation structure.  It is currently
 *  activated by the keyword "test" in the ptraj command list and therefore
 *  is functional.  Information in the actionInformation structure is setup
 *  and initialized in the routine ptrajSetup() in ptraj.c and filled out via
 *  a call to this function in the PTRAJ_SETUP mode.
 *
 ******************************************************************************/

   int
actionTest(actionInformation *action, 
	   double *x, double *y, double *z, 
	   double *box, int mode)
{
  //char *name = "test";
  argStackType **argumentStackPointer;
  ptrajState *state;
  char *buffer;

  /*  This is a prototype routine for processing the series of 
   *  coordinates.  Passed into this routine are:
   *
   *    action   -- a structure that contains global control information
   *                (such as the "state") and local information for this
   *                routine.  This gets used and setup in this routine.
   *
   *    x,y,z    -- the coordinates
   *
   *    box      -- the box information
   *
   *    mode     -- the current mode (see below).
   *
   *  The mode controls what is done within the routine and currently this
   *  routine will be called multiple times, each time with a different
   *  mode.  The current order is...
   *
   *      initially, call once at setup, each of:
   *
   *    PTRAJ_SETUP    -- process arguments and setup the action structure
   *    PTRAJ_STATUS   -- print a summary of what is going to happen
   *
   *      call repeatedly for each coordinate frame:
   *
   *    PTRAJ_ACTION   -- do the actual work on each coordinate
   *
   *      finally, at the end call once, each of:
   *
   *    PTRAJ_PRINT    -- print out any results
   *    PTRAJ_CLEANUP  -- clean up memory allocated, etc.
   *
   *
   *  The action (actionInformation *) structure is a complex beast.
   *
   *  It contains or will contain all the information necessary for 
   *  processing the trajectory files.  Initially on the first call
   *  (at PTRAJ_SETUP) it comes in partially setup via the following
   *  code in ptrajSetup():
   *
   *     statep = ptrajCurrentState();
   *     state = *statep;
   *
   *     action = (actionInformation *)
   *        safe_malloc(sizeof(actionInformation));
   *     INITIALIZE_actionInformation(action);
   *
   *     action->state = ptrajCopyState(ptrajCurrentState());
   *
   *     action->type = TRANSFORM_TEST;
   *     action->fxn  = (actionFunction) actionTest;
   *
   *  In summary, this:
   *
   *  -- sets the current state
   *  -- sets the typedef name for this transform function
   *  -- sets a function pointer to this routine.
   *
   *  Based on this and the processing of command line information, we can
   *  set up the rest of the action structure by setting the values of
   *  various placeholders:
   *
   *    iarg1   -- integer  argument 1
   *    iarg2   -- integer  argument 2
   *    iarg3   -- integer  argument 3
   *    iarg4   -- integer  argument 4
   *    darg1   -- double   argument 1
   *    darg2   -- double   argument 2
   *    darg3   -- double   argument 3
   *    darg4   -- double   argument 4
   *    carg1   -- (void *) argument 1 (i.e. a pointer to a structure)
   *    carg2   -- (void *) argument 2
   *    carg3   -- (void *) argument 3
   *    carg4   -- (void *) argument 4
   *
   *  Initially, command line arguments come in on action->carg1; these
   *  should be processed at the PTRAJ_SETUP stage as shown below.
   */

  state = (ptrajState *) action->state;
  if (prnlev > 2) {
    fprintf(stdout, "In actionTest: currently set to process %i frames\n",
	    state->maxFrames);
  }

  if (mode == PTRAJ_SETUP) {

    /*
     *  ACTION: PTRAJ_SETUP
     *
     *  Parse arguments off the stack and fill in any necessary
     *  information in the actionInformation "action" structure.
     *
     *  This mode is invoked by ptrajSetup().  The current
     *  argumentStack is passed in as action->carg1.
     */

#ifdef MPI
    printParallelError(name);
    return -1;
#endif

    argumentStackPointer = (argStackType **) action->carg1;
    action->carg1 = NULL;

    /*
     *  To process user input to this "action" command (which is in the form 
     *  of text typed in by the user that was placed onto the "argumentStack"
     *  as a series of white space separated strings) it is necessary to
     *  parse the information placed on the "argumentStack" that came into
     *  this routine.  There are a variety of routines that can aid in the 
     *  processing of these "strings" from the argumentStack (these are
     *  defined in dispatch.c):
     *
     *  (int )   argumentStackContains(&argumentStack, "foo");
     *
     *     Check to see if the string "foo" is a substring of any of the strings
     *     currently on the stack.  If it is, remove it from the stack and return
     *     1 (true), otherwise return 0 (false).
     *
     *  (char *) argumentStackKeyToString( &argumentStack, "foo", default);
     *  (int)    argumentStackKeyToInteger(&argumentStack, "foo", default);
     *  (float)  argumentStackKeyToFloat(  &argumentStack, "foo", default);
     *  (double) argumentStackKeyToDouble( &argumentStack, "foo", default);
     *
     *     If the string "foo" is a substring of any of the strings currently
     *     on the stack, remove this entry and grab the next entry converted
     *     to a string, integer, floating point or double precision value,
     *     as appropriate.  If "foo" is not found, return the default value.
     *
     *  (char *) getArgumentString( &argumentStack, default);
     *  (int)    getArgumentInteger(&argumentStack, default);
     *  (float)  getArgumentFloat(  &argumentStack, default);
     *  (double) getArgumentDouble( &argumentStack, default);
     *
     *     Grab the next string off of the argumentStack and return its
     *     value converted to a string, integer, float or double 
     *     depending on the function called.  If the string on the stack
     *     is null, return the default value.
     *
     *  As an example, lets say the possible commands to test are as
     *  follows:
     *
     *   test mask [alpha | beta | gamma] [sets 100]"
     *
     *  This can be parsed with the following, setting up the action "mask",
     *  iarg1 and iarg2.
     */

    buffer = getArgumentString(argumentStackPointer, NULL);
    action->mask = processAtomMask(buffer, action->state);
    safe_free(buffer);

    if (argumentStringContains(argumentStackPointer, "alpha"))
      action->iarg1 = 1;
    else if (argumentStringContains(argumentStackPointer, "beta"))
      action->iarg1 = 2;
    else if (argumentStringContains(argumentStackPointer, "gamma"))
      action->iarg1 = 3;

    action->iarg2 = argumentStackKeyToInteger(argumentStackPointer, "set", 0);
    
    /*
     *  For more examples, see the various transformSetup routines below.
     *
     *
     *  After argument processing, assuming successful return 
     *  (i.e. return 0 is performed by this routine), the action is
     *  placed on the transformActionStack by ptrajSetup() with
     *  the following code:
     *
     *     pushBottomStack( &transformActionStack, (void *) action );
     *
     *  If an error was encountered during argument processing, return -1
     *  and this action will be freed and not placed on the 
     *  transformActionStack
     */

    return 0;
  }


  if (mode == PTRAJ_STATUS) {

    /*
     *  Print out a summary of information about this action.
     *  Basically you would output a summary of what specific
     *  options were turned on, etc.
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
   *  Perform action on coordinates, etc.
   *
   *  To store information for later processing, some actions use stacks.
   *  For example transformAngle makes use of the scalarStack, which holds
   *  scalarInfo structures. The user passes in a name for the angle command,
   *  and this name is given to the scalarInfo structure, which is then placed
   *  on the scalarStack. The scalarInfo structure can be retrieved from the
   *  stack, which is done when the angle info is added, and when the angles
   *  are printed.
   */

  return 1;
}

/** ACTION ROUTINE *************************************************************
 *
 *  transformAtomicFluct() --- compute atomic positional fluctuations
 *
 ******************************************************************************/
   int
transformAtomicFluct(actionInformation *action, 
		     double *x, double *y, double *z,
		     double *box, int mode)
{
  //char *name = "atomicfluct";
  argStackType **argumentStackPointer;
  char *buffer;

  char *filename;
  coordinateInfo *info;
  double *xx, *yy, *zz, *results;
  double xi, yi, zi, fluct, bfactor;
  int i, j;

  /*
   *  USAGE:
   *
   *    atomicfluct [out filename] [<mask>] [start <start>] [stop <stop>] [offset <offset>]
   *                [byres | byatom | bymask] [bfactor]
   *
   *  action argument usage:
   *
   *    iarg1:
   *      0 -- by atom
   *      1 -- by residue
   *      2 -- by mask
   *    iarg2:
   *      the number of sets in the average
   *    iarg3:
   *      the number of visits
   *    iarg4:
   *      0 -- atomic positional fluctuations
   *      1 -- B factors
   *    carg1:
   *      a coordinate info structure
   *    carg2, carg3, carg4:
   *      the x, y and z coordinates (accumulated)
   *
   */

  if (mode == PTRAJ_SETUP) {

    /*
     *  ACTION: PTRAJ_SETUP
     */

#ifdef MPI

#endif

    argumentStackPointer = (argStackType **) action->carg1;
    action->carg1 = NULL;

    filename = argumentStackKeyToString( argumentStackPointer, "out", NULL );

    info = (coordinateInfo *) safe_malloc(sizeof(coordinateInfo));
    INITIALIZE_coordinateInfo(info);
    info->file = NULL;
    info->filename = filename;
    info->option1 = 0;
    info->option2 = 0;
    info->isVelocity = 0;
    info->info = NULL;
    info->mask = NULL;

    info->start = argumentStackKeyToInteger(argumentStackPointer, "start",  1);
    info->stop  = argumentStackKeyToInteger(argumentStackPointer, "stop",  -1);
    if (info->stop == -1) {
      info->stop  = argumentStackKeyToInteger(argumentStackPointer, "end",  -1);
    }
    info->offset = argumentStackKeyToInteger(argumentStackPointer, "offset", 1);

    action->iarg1 = 0;
    if ( argumentStackContains( argumentStackPointer, "byres" ) )
      action->iarg1 = 1;
    else if ( argumentStackContains( argumentStackPointer, "bymask" ) )
      action->iarg1 = 2;
    else if ( argumentStackContains( argumentStackPointer, "byatom" ) )
      action->iarg1 = 0;
    else if ( argumentStackContains( argumentStackPointer, "byatm" ) )
      action->iarg1 = 0;

    action->iarg4 = argumentStackContains( argumentStackPointer, "bfactor" );

    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL) {
      action->mask = processAtomMask( (char *) "*", action->state);
    } else {
      action->mask = processAtomMask(buffer, action->state);
      safe_free(buffer);
    }
    
    action->carg1 = (void *) info;

    xx = (double *) safe_malloc(sizeof(double) * action->state->atoms * 2);
    yy = (double *) safe_malloc(sizeof(double) * action->state->atoms * 2);
    zz = (double *) safe_malloc(sizeof(double) * action->state->atoms * 2);

    for (i=0; i < action->state->atoms*2; i++) {
      xx[i] = 0.0;
      yy[i] = 0.0;
      zz[i] = 0.0;
    }

    action->carg2 = (void *) xx;
    action->carg3 = (void *) yy;
    action->carg4 = (void *) zz;

    action->iarg2 = 0;
    action->iarg3 = worldrank + 1;

    return 0;

  } else if (mode == PTRAJ_STATUS) {

    /*
     *  ACTION: PTRAJ_STATUS
     */


    info = (coordinateInfo *) action->carg1;

    fprintf(stdout, "  ATOMICFLUCT: dumping %s %s %s",
	    (action->iarg4 ? "B factors" : "atomic positional fluctuations"),
	    (action->iarg1 == 2 ? "by mask" : (action->iarg1 == 1 ? "by residue" : "by atom")),
	    (info->filename == NULL ? "to standard output" : "to file "));
    if (info->filename != NULL)
      fprintf(stdout, "%s\n", info->filename);
    else
      fprintf(stdout, "\n");
    if (info->start != 1 && info->stop != -1 && info->offset != 1) {
      fprintf(stdout, "      start: %i", info->start);
      if (info->stop > 0)
	fprintf(stdout, "  stop: %i", info->stop);
      else
	fprintf(stdout, "  stop [at final frame]");
      fprintf(stdout, "  offset: %i\n", info->offset);
    }
    fprintf(stdout, "      Atom selection ");
    printAtomMask(stdout, action->mask, action->state);
    fprintf(stdout, "\n");

  } else if (mode == PTRAJ_PRINT) {

    /*
     *  ACTION: PTRAJ_PRINT
     */

    /*
     *  PF - multiptraj
     *  We need to combine the data using 3 MPI_Reduce's,
     *  passing them to rank 0, then let that rank print.
     *  First need to allocate memory for rank 0 recvbuf
     */

    int sets;

#ifdef MPI
    double *xx_local, *yy_local, *zz_local;

    xx_local = (double *) action->carg2;
    yy_local = (double *) action->carg3;
    zz_local = (double *) action->carg4;

    if (worldrank == 0) {
      xx = (double *) safe_malloc(sizeof(double) * action->state->atoms * 2);
      yy = (double *) safe_malloc(sizeof(double) * action->state->atoms * 2);
      zz = (double *) safe_malloc(sizeof(double) * action->state->atoms * 2);
    }

    MPI_Reduce(xx_local, xx, action->state->atoms * 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(yy_local, yy, action->state->atoms * 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(zz_local, zz, action->state->atoms * 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&action->iarg2, &sets, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    xx = (double *) action->carg2;
    yy = (double *) action->carg3;
    zz = (double *) action->carg4;
    sets = action->iarg2;
#endif    

    if (worldrank == 0) {
      results = (double *) safe_malloc(sizeof(double) * action->state->atoms);
      for (i=0; i < action->state->atoms; i++) {
	results[i] = 0.0;
      }

      info = (coordinateInfo *) action->carg1;

      for (i=0; i < action->state->atoms; i++) {

	xx[i] /= sets;
	yy[i] /= sets;
	zz[i] /= sets;
      }

      for (i=0; i < action->state->atoms; i++) {

	xi = xx[i+action->state->atoms]/sets - xx[i]*xx[i];
	yi = yy[i+action->state->atoms]/sets - yy[i]*yy[i];
	zi = zz[i+action->state->atoms]/sets - zz[i]*zz[i];
	xx[i] = xi;
	yy[i] = yi;
	zz[i] = zi;

      }

      if (info->filename)
	info->file = safe_fopen(info->filename, "w");
      else
	info->file = stdout;

      if (info->file == NULL) {
	fprintf(stdout, "WARNING in ptraj(), atomicfluct: error on opening %s for output\n",
		(info->filename == NULL ? "(stdout)" : info->filename));
	return 0;
      }

      fprintf(stdout, "PTRAJ ATOMICFLUCT: Dumping atomic positional fluctuations\n");
      if (info->file == stdout) {
	if (action->iarg1 == 0) 
	  fprintf(stdout, "  ATOM   FLUCTUATION\n");
	else if (action->iarg1 == 1) 
	  fprintf(stdout, "   RES   FLUCTUATION\n");
	else if (action->iarg1 == 2) 
	  fprintf(stdout, "  MASK   FLUCTUATION\n");
      }

      if (action->iarg4 > 0)
	bfactor = (8.0/3.0)*PI*PI;
      else
	bfactor = 1.0;

      for (i=0; i < action->state->atoms; i++) {
	fluct = xx[i] + yy[i] + zz[i];
	if (fluct > 0) {
	  if (action->iarg4 > 0)
	    /*
	     *  B-factors are (8/3)*PI*PI * <r>**2 hence we do not sqrt the fluctuations!!!
	     */
	    results[i] = bfactor * fluct;
	  else
	    results[i] = sqrt(fluct);
	} else
	  results[i] = 0.0;
      }

      if (action->iarg1 == 0) {
	/*
	 *  byatom print out
	 */
	j = 0;
	for (i=0; i < action->state->atoms; i++) {
	  if (action->mask[i]) {
	    if (action->iarg1 == 0) {  
	      fprintf(info->file, " %6i  %f\n", j+1, results[i]);
	      j++;
	    }
	  }
	}
      } else if (action->iarg1 == 1) {
	/*
	 *  byres print out
	 */
	for (i=0; i < action->state->residues; i++) {

	  xi = 0.0;
	  fluct = 0.0;
	  for (j=action->state->ipres[i]-1; j<action->state->ipres[i+1]-1; j++) {
	    if (action->mask[j]) {
	      xi += action->state->masses[j];
	      fluct += results[j] * action->state->masses[j];
	    }
	  }
	  if (xi > SMALL)
	    fprintf(info->file, " %6i  %f\n", i+1, fluct/xi);
	}
      } else if (action->iarg1 == 2) {
	/*
	 *  bymask print out
	 */
	xi = 0.0;
	fluct = 0.0;
	for (i=0; i < action->state->atoms; i++) {
	  if (action->mask[i]) {
	    xi += action->state->masses[i];
	    fluct += results[i] * action->state->masses[i];
	  }
	}
	if (xi > SMALL)
	  fprintf(info->file, " %6i  %f\n", 1, fluct/xi);
      }

      if (info->file != stdout) {
	safe_fclose(info->file);
	info->file = NULL;
      }
      safe_free(results);
#ifdef MPI
      safe_free(xx);
      safe_free(yy);
      safe_free(zz);
#endif      
    } // END if worldrank == 0

  } else if (mode == PTRAJ_CLEANUP) {

    /*
     *  ACTION: PTRAJ_CLEANUP
     */

    info = (coordinateInfo *) action->carg1;
    xx = (double *) action->carg2;
    yy = (double *) action->carg3;
    zz = (double *) action->carg4;

    safe_free(action->mask);
    safe_free(info->filename);
    safe_free(info);
    safe_free(xx);
    safe_free(yy);
    safe_free(zz);

  }

  if (mode != PTRAJ_ACTION) return 0;

  /*
   *  ACTION: PTRAJ_ACTION
   */

  info = (coordinateInfo *) action->carg1;

  xx = (double *) action->carg2;
  yy = (double *) action->carg3;
  zz = (double *) action->carg4;

  if (action->iarg3 >= info->start &&
      (info->stop < 0 || action->iarg3 <= info->stop) &&
      (action->iarg3 - info->start)%info->offset == 0) {

    action->iarg2++;

    for (i=0; i < action->state->atoms; i++) {
      xx[i] += x[i];
      yy[i] += y[i];
      zz[i] += z[i];
      xx[i+action->state->atoms] += x[i]*x[i];
      yy[i+action->state->atoms] += y[i]*y[i];
      zz[i+action->state->atoms] += z[i]*z[i];
    }
  }
  action->iarg3 += worldsize;

  return 1;
}

/** ACTION ROUTINE *************************************************************
 *
 *  transformAtomicFluct3D() --- compute atomic positional fluctuations in 3D
 *
 ******************************************************************************/
// windowInfoType
typedef struct _windowInfoType {
  int windowWidth;
  int windowOffset;
  int windows;
  int *windowStart;
  int *windowStop;
  int *windowVisits;
  double *x;
  double *y;
  double *z;
} windowInfoType;
// setupWindowInfo()
windowInfoType * setupWindowInfo(windowInfoType *windowInfo, 
		int windowWidth,   /* width between the windows */
		int windowOffset,  /* offset between windows */
		int windows,       /* number of windows to create or increase */
		int size)          /* amount of space to allocate per window */
{
  int i, sz;


  /*
   *  creates/reallocates a windowInfoType structure to running average 
   *  information over multiple windows of width "windowWidth" with an
   *  offset between windows of "windowOffset".
   *
   *  If windowInfo != NULL reallocate the structures to INCREASE
   *  the size by "windows"
   *
   *  If windowInfo == NULL do the initial allocation
   */

  if (windowInfo == NULL) {

    if (windowWidth < 1 || windowOffset < 1 || windows < 1) {
      warning("setupWindowInfo()", "Window offset (%i), width (%i), or #windows (%i) out of bounds!\n",
	      windowOffset, windowWidth, windows);
      return NULL;
    }

    windowInfo = safe_malloc(sizeof(windowInfoType));
    windowInfo->windowWidth = windowWidth;
    windowInfo->windowOffset = windowOffset;
    windowInfo->windows = windows;

    windowInfo->windowStart  = (int *) safe_malloc( sizeof(int) * windows );
    windowInfo->windowStop   = (int *) safe_malloc( sizeof(int) * windows );
    windowInfo->windowVisits = (int *) safe_malloc( sizeof(int) * windows );
    windowInfo->x = (double *) safe_malloc( sizeof(double) * windows * size );
    windowInfo->y = (double *) safe_malloc( sizeof(double) * windows * size );
    windowInfo->z = (double *) safe_malloc( sizeof(double) * windows * size );

    for (i=0; i < windows; i++) {
      windowInfo->windowStart[i] = i * windowOffset + 1;
      windowInfo->windowStop[i] = windowInfo->windowStart[i] + windowWidth;

      if (prnlev > 3) {
	fprintf(stdout, "Window %i from %i to %i\n", i+1, 
		windowInfo->windowStart[i], windowInfo->windowStop[i]);
      }
    }
    return windowInfo;
  }

  sz = sizeof(int);
  windowInfo->windowStart  = (int *)
    safe_realloc( (void *) windowInfo->windowStart, sz*windowInfo->windows, sz*windows );
  windowInfo->windowStop   = (int *)
    safe_realloc( (void *) windowInfo->windowStop, sz*windowInfo->windows, sz*windows );
  windowInfo->windowVisits = (int *)
    safe_realloc( (void *) windowInfo->windowVisits, sz*windowInfo->windows, sz*windows );

  sz = sizeof(double) * size;
  windowInfo->x = (double *)
    safe_realloc( (void *) windowInfo->x, sz*windowInfo->windows, sz*windows );
  windowInfo->y = (double *)
    safe_realloc( (void *) windowInfo->y, sz*windowInfo->windows, sz*windows );
  windowInfo->z = (double *)
    safe_realloc( (void *) windowInfo->z, sz*windowInfo->windows, sz*windows );
 
  for (i=0; i < windows; i++) {
    windowInfo->windowStart[i+windowInfo->windows] = 
      windowInfo->windowStart[windowInfo->windows-1] + (i+1)*windowInfo->windowOffset;
    windowInfo->windowStop[i+windowInfo->windows] = 
      windowInfo->windowStart[i+windowInfo->windows] + windowInfo->windowWidth;
    windowInfo->windowVisits[i+windowInfo->windows] = 0;

    if (prnlev > 3) {
      fprintf(stdout, "Window %i from %i to %i\n", i+windowInfo->windows+1, 
	      windowInfo->windowStart[i+windowInfo->windows], 
	      windowInfo->windowStop[i+windowInfo->windows]);
      }
  }
  windowInfo->windows = windowInfo->windows + windows;
  return(windowInfo);

}

// transformAtomicFluct3D()
int transformAtomicFluct3D(actionInformation *action, 
		       double *x, double *y, double *z,
		       double *box, int mode)
{
  //char *name = "atomicfluct3D";
  argStackType **argumentStackPointer;
  char *buffer;

  char *filename;
  coordinateInfo *info;
  //double *xx, *yy, *zz;
  double *results;
  double xi, yi, zi, fluct, bfactor;
  int i, j, w, ind, minwin, maxwin, windowWidth, windowOffset;
  windowInfoType *windowInfo;

  /*
   *  USAGE:
   *
   *    atomicfluct [out filename] [<mask>] [start <start>] [stop <stop>] [offset <offset>]
   *                [byres | byatom | bymask] [bfactor] [window <value> [by <value>]]
   *
   *  action argument usage:
   *
   *    iarg1:
   *      0 -- by atom
   *      1 -- by residue
   *      2 -- by mask
   *    iarg2:
   *      the number of sets in the average
   *    iarg3:
   *      the number of visits
   *    iarg4:
   *      0 -- atomic positional fluctuations
   *      1 -- B factors
   *    carg1:
   *      a coordinate info structure
   *    carg2, carg3, carg4:
   *      the x, y and z coordinates (accumulated)
   *
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

    filename = argumentStackKeyToString( argumentStackPointer, "out", NULL );

    info = (coordinateInfo *) safe_malloc(sizeof(coordinateInfo));
    INITIALIZE_coordinateInfo(info);
    info->file = NULL;
    info->filename = filename;
    info->option1 = 0;
    info->option2 = 0;
    info->isVelocity = 0;
    info->info = NULL;
    info->mask = NULL;

    info->start = argumentStackKeyToInteger(argumentStackPointer, "start",  1);
    info->stop  = argumentStackKeyToInteger(argumentStackPointer, "stop",  -1);
    if (info->stop == -1) {
      info->stop  = argumentStackKeyToInteger(argumentStackPointer, "end",  -1);
    }
    info->offset = argumentStackKeyToInteger(argumentStackPointer, "offset", 1);
    action->carg1 = (void *) info;

    action->iarg1 = 0;
    if ( argumentStackContains( argumentStackPointer, "byres" ) )
      action->iarg1 = 1;
    else if ( argumentStackContains( argumentStackPointer, "bymask" ) )
      action->iarg1 = 2;
    else if ( argumentStackContains( argumentStackPointer, "byatom" ) )
      action->iarg1 = 0;
    else if ( argumentStackContains( argumentStackPointer, "byatm" ) )
      action->iarg1 = 0;

    action->iarg4 = argumentStackContains( argumentStackPointer, "bfactor" );

    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL) {
      action->mask = processAtomMask( (char *) "*", action->state);
    } else {
      action->mask = processAtomMask(buffer, action->state);
      safe_free(buffer);
    }
    
    /*
     *  goal here is to create a 3D bfactor value that represents a series of
     *  moving windows in time.  The width of each window is specified by "window"
     *  and the offset for each window is specified by "by".  For example, if you
     *  specify "window 1000 by 100" this will create averages 0-1000, 100-1100, etc.
     */
    windowWidth  = argumentStackKeyToInteger(argumentStackPointer, "window", 1000);
    windowOffset = argumentStackKeyToInteger(argumentStackPointer, "by", 500);

    windowInfo = setupWindowInfo(NULL, windowWidth, windowOffset, 
				 100, 2*action->state->atoms);

    action->carg2 = (void *) windowInfo;
    action->iarg2 = 0;
    action->iarg3 = 0;
    
    return 0;

  } else if (mode == PTRAJ_STATUS) {

    /*
     *  ACTION: PTRAJ_STATUS
     */


    info = (coordinateInfo *) action->carg1;
    windowInfo = (windowInfoType *) action->carg2;

    fprintf(stdout, "  ATOMICFLUCT -3D- : dumping %s %s %s",
	    (action->iarg4 ? "B factors" : "atomic positional fluctuations"),
	    (action->iarg1 == 2 ? "by mask" : (action->iarg1 == 1 ? "by residue" : "by atom")),
	    (info->filename == NULL ? "to standard output" : "to file "));
    if (info->filename != NULL)
      fprintf(stdout, "%s\n", info->filename);
    else
      fprintf(stdout, "\n");
    fprintf(stdout, "    WINDOW width = %i, offset = %i\n", 
	    windowInfo->windowWidth, windowInfo->windowOffset);
    if (info->start != 1 && info->stop != -1 && info->offset != 1) {
      fprintf(stdout, "      start: %i", info->start);
      if (info->stop > 0)
	fprintf(stdout, "  stop: %i", info->stop);
      else
	fprintf(stdout, "  stop [at final frame]");
      fprintf(stdout, "  offset: %i\n", info->offset);
    }
    fprintf(stdout, "      Atom selection ");
    printAtomMask(stdout, action->mask, action->state);
    fprintf(stdout, "\n");

  } else if (mode == PTRAJ_PRINT) {

    /*
     *  ACTION: PTRAJ_PRINT
     */


    results = (double *) safe_malloc(sizeof(double) * action->state->atoms);
    info = (coordinateInfo *) action->carg1;
    
    windowInfo = (windowInfoType *) action->carg2;


    if (info->filename)
      info->file = safe_fopen(info->filename, "w");
    else
      info->file = stdout;

    if (info->file == NULL) {
      fprintf(stdout, "WARNING in ptraj(), atomicfluct: error on opening %s for output\n",
	      (info->filename == NULL ? "(stdout)" : info->filename));
      return 0;
    }


    fprintf(stdout, "PTRAJ ATOMICFLUCT: Dumping atomic positional fluctuations\n");
    if (info->file == stdout) {
      if (action->iarg1 == 0) 
	fprintf(stdout, "  ATOM   FLUCTUATION\n");
      else if (action->iarg1 == 1) 
	fprintf(stdout, "   RES   FLUCTUATION\n");
      else if (action->iarg1 == 2) 
	fprintf(stdout, "  MASK   FLUCTUATION\n");
    }

    if (action->iarg4 > 0)
      bfactor = (8.0/3.0)*PI*PI;
    else
      bfactor = 1.0;

    for (w=0; w < windowInfo->windows; w++) {
      if (windowInfo->windowVisits[w] > 0) {

	printf("WINDOW VISITS = %i for window %i\n", windowInfo->windowVisits[w], w);
	for (i=0; i < action->state->atoms; i++) {

	  ind = i+w*2*action->state->atoms;
	  windowInfo->x[ind] /= windowInfo->windowVisits[w];
	  windowInfo->y[ind] /= windowInfo->windowVisits[w];
	  windowInfo->z[ind] /= windowInfo->windowVisits[w];
	}

	for (i=0; i < action->state->atoms; i++) {
	  
	  ind = i+w*2*action->state->atoms;
	  xi = windowInfo->x[ ind + action->state->atoms] / windowInfo->windowVisits[w]
	    - windowInfo->x[ind] * windowInfo->x[ind];
	  yi = windowInfo->y[ ind + action->state->atoms] / windowInfo->windowVisits[w]
	    - windowInfo->y[ind] * windowInfo->y[ind];
	  zi = windowInfo->z[ ind + action->state->atoms] / windowInfo->windowVisits[w]
	    - windowInfo->z[ind] * windowInfo->z[ind];

	  windowInfo->x[ind] = xi;
	  windowInfo->y[ind] = yi;
	  windowInfo->z[ind] = zi;

	}

	for (i=0; i < action->state->atoms; i++) {
	  ind = i+w*2*action->state->atoms;
	  fluct = windowInfo->x[ind] + windowInfo->y[ind] + windowInfo->z[ind];

	  if (fluct > 0) {
	    if (action->iarg4 > 0)
	      /*
	       *  B-factors are (8/3)*PI*PI * <r>**2 hence we do not sqrt the fluctuations!!!
	       */
	      results[i] = bfactor * fluct;
	    else
	      results[i] = sqrt(fluct);
	  } else
	    results[i] = 0.0;
	}

	if (action->iarg1 == 0) {
	  /*
	   *  byatom print out
	   */
	  j = 0;
	  for (i=0; i < action->state->atoms; i++) {
	    if (action->mask[i]) {
	      if (action->iarg1 == 0) {  
		fprintf(info->file, " %6i %.2f %f\n", j+1, 
			(windowInfo->windowStart[w] + windowInfo->windowWidth/2.0), results[i]);
		j++;
	      }
	    }
	  }
	} else if (action->iarg1 == 1) {
	  /*
	   *  byres print out
	   */
	  for (i=0; i < action->state->residues; i++) {


	    xi = 0.0;
	    fluct = 0.0;
	    for (j=action->state->ipres[i]-1; j<action->state->ipres[i+1]-1; j++) {
	      if (action->mask[j]) {
		xi += action->state->masses[j];
		fluct += results[j] * action->state->masses[j];
	      }
	    }
	    if (xi > SMALL)
	      fprintf(info->file, " %6i %.2f %f\n", i+1, 
		      (windowInfo->windowStart[w] + windowInfo->windowWidth/2.0), 
		      fluct/xi);
	  }
	} else if (action->iarg1 == 2) {
	  /*
	   *  bymask print out
	   */
	  xi = 0.0;
	  fluct = 0.0;
	  for (i=0; i < action->state->atoms; i++) {
	    if (action->mask[i]) {
	      xi += action->state->masses[i];
	      fluct += results[i] * action->state->masses[i];
	    }
	  }
	  if (xi > SMALL)
	    fprintf(info->file, " %6i  %.2f %f\n", 1, 
		    (windowInfo->windowStart[w] + windowInfo->windowWidth/2.0), 
		    fluct/xi);
	}
      }
    }
    if (info->file != stdout) {
      safe_fclose(info->file);
      info->file = NULL;
    }

  } else if (mode == PTRAJ_CLEANUP) {

    /*
     *  ACTION: PTRAJ_CLEANUP
     */

    info = (coordinateInfo *) action->carg1;
    windowInfo = (windowInfoType *) action->carg2;

    safe_free(info);
    safe_free(windowInfo->windowStart);
    safe_free(windowInfo->windowStop);
    safe_free(windowInfo->windowVisits);
    safe_free(windowInfo->x);
    safe_free(windowInfo->y);
    safe_free(windowInfo->z);
    safe_free(windowInfo);

  }

  if (mode != PTRAJ_ACTION) return 0;

  /*
   *  ACTION: PTRAJ_ACTION
   */

  info = (coordinateInfo *) action->carg1;

  windowInfo = (windowInfoType *) action->carg2;

  action->iarg3++;
  if (action->iarg3 >= info->start &&
      (info->stop < 0 || action->iarg3 <= info->stop) &&
      (action->iarg3 - info->start)%info->offset == 0) {

    action->iarg2++;

    minwin = 0;
    maxwin = windowInfo->windows;
    if (maxwin > windowInfo->windows) {
      setupWindowInfo(windowInfo, -1, -1, 100, 2*action->state->atoms);
      maxwin = windowInfo->windows;
    }

    if (action->iarg2 > windowInfo->windowStop[minwin]) {
      minwin++;
      maxwin++;

      if (maxwin > windowInfo->windows)
	maxwin = windowInfo->windows;
    }
    
    for (j=minwin; j < maxwin; j++) {

      if ( action->iarg2 >= windowInfo->windowStart[j] &&
	   action->iarg2 < windowInfo->windowStop[j] ) {
	for (i=0; i < action->state->atoms; i++) {
	  ind=i+j*2*action->state->atoms;
	  windowInfo->x[ind] += x[i];
	  windowInfo->y[ind] += y[i];
	  windowInfo->z[ind] += z[i];
	  windowInfo->x[ind+action->state->atoms] += x[i]*x[i];
	  windowInfo->y[ind+action->state->atoms] += y[i]*y[i];
	  windowInfo->z[ind+action->state->atoms] += z[i]*z[i];
	}
	windowInfo->windowVisits[j] += 1;
      }
    }
  }

  return 1;
}

// -----------------------------------------------------------------------------
// ptrajCopyAction()
actionInformation *ptrajCopyAction(actionInformation **actioninp) {
  actionInformation *action, *actionin;
  /*
   *  Make a copy of the input action pointer.
   */
  actionin = *actioninp;
  action = (actionInformation*)safe_malloc(sizeof(actionInformation));
  INITIALIZE_actionInformation(action);
  action->fxn	= actionin->fxn;
  action->type 	= actionin->type;
  action->iarg1	= actionin->iarg1; 
  action->iarg2	= actionin->iarg2; 
  action->iarg3	= actionin->iarg3; 
  action->iarg4	= actionin->iarg4; 
  action->iarg5	= actionin->iarg5; 
  action->iarg6	= actionin->iarg6; 
  action->iarg7	= actionin->iarg7; 
  action->darg1	= actionin->darg1; 
  action->darg2	= actionin->darg2; 
  action->darg3	= actionin->darg3; 
  action->darg4	= actionin->darg4; 
  action->suppressProcessing = actionin->suppressProcessing;  
  action->performSecondPass	 = actionin->performSecondPass;
  /* The following members are pointers. Just copy from the old action. 
     So you need to allocate new space for the members. */
  action->state	= actionin->state;
  action->mask	= actionin->mask;
  action->carg1	= actionin->carg1; 
  action->carg2	= actionin->carg2; 
  action->carg3	= actionin->carg3; 
  action->carg4	= actionin->carg4; 
  action->carg5	= actionin->carg5; 
  action->carg6	= actionin->carg6; 
  action->carg7	= actionin->carg7; 
  return action;
}
// -----------------------------------------------------------------------------

/** ACTION ROUTINE ************************************************************
 *
 *  transformContacts()   --- perform referenced contact calculation
 *
 ******************************************************************************/
typedef enum _transformContactsType {
  CONTACTS_NULL = 0,
  CONTACTS_FIRST,
  CONTACTS_REFERENCE
} transformContactsType;

typedef struct _transformContactsInfo {
  double *refx, *refy, *refz;
  char *filename;
  void *outFile;
  void *outFile2;
  int currentFrame;
  int byResidue;
} transformContactsInfo;

#define INITIALIZE_transformContactsInfo(_p_) \
  _p_->refx      = NULL; \
  _p_->refy      = NULL; \
  _p_->refz      = NULL; \
  _p_->filename  = NULL; \
  _p_->outFile   = NULL; \
  _p_->outFile2  = NULL; \
  _p_->currentFrame = 0; \
  _p_->byResidue    = 0;

typedef struct _contactList {
  int index;
  char *name;
  struct _contactList *next;
} contactList;

#define INITIALIZE_contactList(_p_) \
  _p_->index = -1; \
  _p_->name  = NULL; \
  _p_->next  = NULL;

// transformContacts()
int transformContacts(actionInformation *action,
		  double *x, double *y, double *z,
		  double *box, int mode)
{
  //char *name = "contacts";
  argStackType **argumentStackPointer;
  char *buffer, *buffer2, *buffer3;
  transformContactsInfo *contactsInfo;
  char buffer2ptr[BUFFER_SIZE];
  char buffer3ptr[BUFFER_SIZE];

  double x1,y1,z1,x2,y2,z2,dist;
  int i, j;
  int *activeResidues;
  contactList *list;
  contactList *current;
  int *contactNumberList;
  int *nativeNumberList;
  int n, nativeNumber, contactNumber,resNum;

  nativeNumber = 0; contactNumber = 0;

  /*
   *  USAGE
   *
   *  contacts [first | reference] [byresidue]
   *      [out <filename>] [time <interval>] [distance <cutoff>] [<mask>]
   *
   *  action argument usage:
   *
   *  byresidue: calculate number of contacts for every specified atom and save result per residue
   *  iarg1: transformContactsType
   *    CONTACTS_FIRST     -- take first structure for reference contacts
   *    CONTACTS_REFERENCE -- take reference structure for reference contacts
   *  iarg2: frame counter
   *  iarg3: stores number of atoms
   *  darg1: time interval
   *  darg2: cutoff distance
   *  carg1: the contactsInfo structure
   *  carg2: list of all reference contacts
   *  carg3: list of active residues
   */

  /* Set up buffers for printing */

  //buffer2 = safe_malloc( BUFFER_SIZE * sizeof *buffer2);
  //buffer3 = safe_malloc( BUFFER_SIZE * sizeof *buffer3);
  //buffer2ptr = buffer2;
  //buffer3ptr = buffer3;
  buffer2 = buffer2ptr;
  buffer3 = buffer3ptr;

  if (mode == PTRAJ_SETUP) {

    /*
     *  ACTION PTRAJ_SETUP
     */

#ifdef MPI

#endif

    nativeNumberList = NULL;
    contactNumberList = NULL;

    argumentStackPointer = (argStackType **) action->carg1;
    action->carg1 = NULL;

    contactsInfo = (transformContactsInfo *) safe_malloc(sizeof(transformContactsInfo));
    INITIALIZE_transformContactsInfo(contactsInfo);

    action->iarg1 = (int) CONTACTS_FIRST;
    if (argumentStackContains(argumentStackPointer, "reference"))
      action->iarg1 = (int) CONTACTS_REFERENCE;
    else if (argumentStackContains(argumentStackPointer, "first"))
      action->iarg1 = (int) CONTACTS_FIRST;

    contactsInfo->filename = argumentStackKeyToString(argumentStackPointer, "out", NULL);
    if (contactsInfo->filename != NULL){
      contactsInfo->outFile = ptrajOpenW(contactsInfo->filename);
    }
    else{
      contactsInfo->outFile = stdout;
    }

    if (argumentStackContains(argumentStackPointer, "byresidue"))
      contactsInfo->byResidue = 1;

    action->iarg2 = 0;
    action->iarg3 = action->state->atoms;
    action->carg1 = contactsInfo;
    action->darg1 = argumentStackKeyToDouble(argumentStackPointer, "time", 1.0);
    action->darg2 = argumentStackKeyToDouble(argumentStackPointer, "distance", 7.0);

    /*
     * Get mask (here, everything else should have been processed from the argumentStack)
     */

    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL) {
      if (contactsInfo->byResidue) 
	action->mask = processAtomMask("@CA", action->state);
      else
	action->mask = processAtomMask("*", action->state);
    } else {
      action->mask = processAtomMask(buffer, action->state);
      safe_free(buffer);
    }

    /*
     * Do final setup for which mask information is needed
     */

    nativeNumber = 0;
    if (action->iarg1 == (int) CONTACTS_REFERENCE) {
      // Check that a reference structure has been defined
      if (referenceInfo == NULL) {
        fprintf(stdout,"Error: No reference coordinates defined.\n");
        return -1;
      }
      n = action->state->atoms;
      list = (contactList *) safe_malloc(sizeof(contactList) * n);
      for(i = 0; i < n; i++){
	list[i].index = -1;
	list[i].name = NULL;
	list[i].next = NULL;
      }

      for (i = 0; i < n; i++) {
        if (action->mask[i]) {
	  list[i].index = i;
	  list[i].name = action->state->atomName[i];
	  list[i].next = NULL;

          for (j = 0; j < n; j++) {
	    if (action->mask[j] && i != j) {
	      x1 = referenceInfo->x[i];
	      y1 = referenceInfo->y[i];
	      z1 = referenceInfo->z[i];
	      x2 = referenceInfo->x[j];
	      y2 = referenceInfo->y[j];
	      z2 = referenceInfo->z[j];
	      dist = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
	      if (dist < action->darg2) {
		current = (contactList *) safe_malloc(sizeof(contactList));
		INITIALIZE_contactList(current);
		current->name = action->state->atomName[i];
		current->index = j;
		current->next = list[i].next;
		list[i].next = current;
		nativeNumber++;
	      }
	    }
          }
        }
      }
      action->carg2 = (contactList *) list;
    }

    if(contactsInfo->byResidue){
      contactsInfo->outFile2 = ptrajOpenW(strncat(contactsInfo->filename, ".native", 7));
      sprintf(buffer2, "#time");
      buffer2 = buffer2ptr + strlen(buffer2ptr);

      activeResidues = safe_malloc(sizeof(int)*action->state->residues);
      for (i = 0; i < action->state->residues; i++)
	activeResidues[i] = 0;

      for (i = 0; i < action->state->atoms; i++) {
        if (action->mask[i]) {
	  resNum = atomToResidue(i+1, action->state->residues, action->state->ipres) - 1;
	  activeResidues[resNum] = 1;
        }
      }
      
      for (i = 0; i < action->state->residues; i++)
        if (activeResidues[i]) {
	  sprintf(buffer2, "\tresidue %d", i);
	  buffer2 = buffer2ptr + strlen(buffer2ptr);
	}
      
      sprintf(buffer2, "\n");
      buffer2 = buffer2ptr;
      ptrajfprintfone(contactsInfo->outFile, buffer2);
      ptrajfprintfone(contactsInfo->outFile2, buffer2);
      action->carg3 = (void *) activeResidues;
    }
    else {
      sprintf(buffer2, "#time\tContacts\tnative Contacts ");
      buffer2 = buffer2ptr + strlen(buffer2ptr);
      if (action->iarg1 == (int) CONTACTS_REFERENCE) {
	sprintf(buffer2, "(number of natives: %d)", nativeNumber);
	buffer2 = buffer2ptr + strlen(buffer2ptr);
      }
      sprintf(buffer2, "\n");
      buffer2 = buffer2ptr;
      ptrajfprintfone(contactsInfo->outFile, buffer2);
    }
  }
  else if (mode == PTRAJ_STATUS) {

    /*
     *  ACTION PTRAJ_STATUS
     */

    contactsInfo = (transformContactsInfo *) action->carg1;
    fprintf(stdout, "  CONTACTS: Calculating current contacts and comparing results to ");
    if (action->iarg1 == (int) CONTACTS_FIRST)
      fprintf(stdout, "first frame\n");
    else if (action->iarg1 == (int) CONTACTS_REFERENCE )
      fprintf(stdout, "reference structure\n");
    fprintf(stdout, "                Dumping results to %s,\n", contactsInfo->outFile == NULL ? "STDOUT" : contactsInfo->filename);
    fprintf(stdout, "                using a time interval of %f\n", action->darg1);
    fprintf(stdout, "                and a cutoff of %f A\n", action->darg2);
    if(contactsInfo->byResidue)
      fprintf(stdout, "                Results are output on a per-residue basis\n");
    fprintf(stdout, "                Atom selection follows ");
    printAtomMask(stdout, action->mask, action->state);
    fprintf(stdout, "\n");
  }
  else if (mode == PTRAJ_ACTION) {

    /*
     *  ACTION PTRAJ_ACTION
     */

    contactsInfo = (transformContactsInfo *) action->carg1;
    activeResidues = (int *) action->carg3;
    n = action->state->atoms;

    if (action->iarg1 == (int) CONTACTS_FIRST && action->iarg2 == 0) {
      list = (contactList * )safe_malloc(sizeof(contactList) * n);
      for(i = 0; i < n; i++){
	list[i].index = -1;
	list[i].name = NULL;
	list[i].next = NULL;
      }

      for (i = 0; i < n; i++) {
	if (action->mask[i]) {
	  list[i].index = i;
	  list[i].name = action->state->atomName[i];
	  list[i].next = NULL;

          for (j = 0; j < n; j++) {
	    if (action->mask[j] && i != j) {
	      dist = sqrt((x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]));
	      if (dist < action->darg2) {
		current = (contactList *) safe_malloc(sizeof(contactList));
		INITIALIZE_contactList(current);
		current->name = action->state->atomName[i];
		current->index = j;
		current->next = list[i].next;
		list[i].next = current;
	      }
	    }
          }
        }
      }
      action->carg2 = (contactList *) list;
    }

    list = (contactList *) action->carg2;

    contactNumberList = safe_malloc(sizeof(int) * action->state->residues);
    nativeNumberList  = safe_malloc(sizeof(int) * action->state->residues);
    for (i = 0; i < action->state->residues; i++) {
      contactNumberList[i] = 0;
      nativeNumberList[i]  = 0;
    }

    if (contactsInfo->byResidue) {
      for (i = 0; i < n; i++) {
        if (action->mask[i]) {
	  contactNumber = 0;
	  nativeNumber = 0;
	  resNum = atomToResidue(i+1, action->state->residues, action->state->ipres) - 1;
	  for (j = 0; j < n; j++) {
	    if(action->mask[j] && i != j) {
              dist = sqrt((x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]));
	      if (dist < action->darg2) {
		contactNumber++;
		contactNumberList[resNum]++;

		current = &list[i];
		current = current->next;
                while (current != NULL) {
		  if (current->index == j) {
                    nativeNumber++;
		    nativeNumberList[resNum]++;
                    break;
		  }
		  current = current->next;
		}
	      }
	    }
	  }
	}
      }

      sprintf(buffer2, "%10.2f", (double) (action->iarg2*worldsize+worldrank+1) * action->darg1);
      sprintf(buffer3, "%10.2f", (double) (action->iarg2*worldsize+worldrank+1) * action->darg1);
      buffer2 = buffer2ptr + strlen(buffer2ptr);
      buffer3 = buffer3ptr + strlen(buffer3ptr);
      for (i = 0; i < action->state->residues; i++)
	if (activeResidues[i]) {
          sprintf(buffer2, "\t%d", contactNumberList[i]);
          sprintf(buffer3, "\t%d", nativeNumberList[i]);
	  buffer2 = buffer2ptr + strlen(buffer2ptr);
	  buffer3 = buffer3ptr + strlen(buffer3ptr);
	}
      sprintf(buffer2, "%d\n", contactNumber);
      sprintf(buffer3, "%d\n", nativeNumber);
      buffer2 = buffer2ptr;
      buffer3 = buffer3ptr;
      ptrajfprintf(contactsInfo->outFile, buffer2);
      ptrajfprintf(contactsInfo->outFile2, buffer3);
    }
    else {
      contactNumber = nativeNumber = 0;
      for (i = 0; i < n; i++) {
        if (action->mask[i]) {
	  for (j = 0; j < n; j++) {
	    if(action->mask[j] && i != j) {
              dist = sqrt((x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]));
	      if (dist < action->darg2) {
                contactNumber++;

		current = &list[i];
		current = current->next;
		while (current != NULL) {
		  if (current->index == j) {
                    nativeNumber++;
		    break;
		  }
		  current = current->next;
		}
	      }
	    }
	  }
	}
      }
      sprintf(buffer2, "%10.2f\t%d\t%d\n", (action->iarg2*worldsize+worldrank+1) * action->darg1, contactNumber, nativeNumber);
      ptrajfprintf(contactsInfo->outFile, buffer2);
    }
    action->iarg2++;

    /* Clean up contactNumberList and nativeNumberList */
    safe_free(contactNumberList);
    safe_free(nativeNumberList);

  }
  else if (mode == PTRAJ_PRINT) {

    /*
     *  ACTION PTRAJ_PRINT
     */

    /* Nothing to do here */
  }
  else if (mode == PTRAJ_CLEANUP) {

    /*
     *  ACTION PTRAJ_CLEANUP
     */
    safe_free(action->mask);
    contactsInfo = (transformContactsInfo *) action->carg1;
    if(contactsInfo != NULL){
      if (contactsInfo->byResidue){
	if(contactsInfo->filename != NULL) {
	  ptrajCloseFile(contactsInfo->outFile2);
	}
	safe_free((int *) action->carg3);
      }
      if(contactsInfo->filename != NULL){
	ptrajCloseFile(contactsInfo->outFile);
	safe_free(contactsInfo->filename);
      }
      INITIALIZE_transformContactsInfo(contactsInfo);
      safe_free(contactsInfo);
    }

    /*
     *  This should be done each time in PTRAJ_ACTION
     *
     *  safe_free(contactNumberList);
     *  safe_free(nativeNumberList);
     */

    list = (contactList *) action->carg2;
    if(list != NULL){
      for(i = 0; i < action->iarg3; i++){
	current = list + i;
	current = current->next;
	while (current != NULL) {
	  contactList *next = current->next;
	  INITIALIZE_contactList(current);
	  safe_free(current);
	  current = next;
	}
      }
      safe_free(list);
    }
  }

  return 1;
}


/** ACTION ROUTINE *************************************************************
 *
 *  transformCorr()   --- perform correlation analysis (Vickie Tsui, Scripps)
 *
 *  Supplementary routines:
 *    compute_corr()
 *
 ******************************************************************************/
// compute_corr()
void compute_corr(char *outfile, double *x, double *y, double *z, 
                  int ts, int tc, int tm, int tf)
{
  typedef struct _complex {
    double real;
    double imaginary;
  } complex;

  int i,j, current_frame;
  double *dipcrd1, *dipcrd2, *dipcrd3;
  double rmag0, rmagt;
  double x0, y0, z0, xt, yt, zt;
  double *p2, *corr, *rcorr;
  double r6ave, r3ave, rave, avecrd[4], rrig;
  double dot, y2asy, y20;
  double th0, phi0;
  complex y21, y21c, y22, y22c;
  int *cfind, npts, ncorr, ind0, indt, ntot;
  int doit, jmax;
  FILE *ifp;


  /* allocate space */

  jmax = tm/ts + 1;
  if (tf+1 > jmax) jmax = tf+1;

  dipcrd1 = (double *) safe_malloc( sizeof(double) * jmax );
  dipcrd2 = (double *) safe_malloc( sizeof(double) * jmax );
  dipcrd3 = (double *) safe_malloc( sizeof(double) * jmax );

  p2     = (double *) safe_malloc( sizeof(double) * jmax );
  corr   = (double *) safe_malloc( sizeof(double) * jmax );
  rcorr  = (double *) safe_malloc( sizeof(double) * jmax );
  cfind  = (int *) safe_malloc( sizeof(int) * jmax );


/* initialize */
  for (i=tf; i>=1; --i)  {
    x[i]=x[i-1]; y[i]=y[i-1]; z[i]=z[i-1];
  }
  for (i=1; i<=tf; ++i)  {
    corr[i]=0.0;
    p2[i]=0.0;
    rcorr[i]=0.0;
    cfind[i]=i;
  }
  ntot=tm/ts;
  npts=ntot;  ncorr=0;  current_frame=ntot;
  for (i=1; i<=ntot; ++i)  {
    dipcrd1[i]=x[i];  dipcrd2[i]=y[i]; dipcrd3[i]=z[i];
  }
  jmax=ntot;
  r6ave=0.0; r3ave=0.0;
  avecrd[1]=0.0;  avecrd[2]=0.0;  avecrd[3]=0.0;
  rave=0.0; y2asy=0.0; y20=0.0;
  y21.real=0.0; y21.imaginary=0.0;
  y21c.real=0.0; y21c.imaginary=0.0;
  y22.real=0.0; y22.imaginary=0.0;
  y22c.real=0.0; y22c.imaginary=0.0;

/* main loop for calculating correlation functions */

  doit=1;
  while (doit > 0)  {
   ncorr=ncorr+1;
   ind0=cfind[1];

   rmag0=pow(dipcrd1[ind0],2)+pow(dipcrd2[ind0],2)+pow(dipcrd3[ind0],2);

   rmag0=sqrt(rmag0);
   x0=dipcrd1[ind0]/rmag0;
   y0=dipcrd2[ind0]/rmag0;
   z0=dipcrd3[ind0]/rmag0;

   r6ave=r6ave+1/pow(rmag0,6);
   r3ave=r3ave+1/pow(rmag0,3);
   rave=rave+rmag0;
   avecrd[1]=avecrd[1]+dipcrd1[ind0];
   avecrd[2]=avecrd[2]+dipcrd2[ind0];
   avecrd[3]=avecrd[3]+dipcrd3[ind0];

   th0=acos(z0);
   phi0=atan2(y0,x0);

   y22.real=y22.real+sqrt(3.0/4.0)*pow((sin(th0)),2)*(cos(2*phi0))/pow(rmag0,3);
   y22.imaginary=y22.imaginary+sqrt(3.0/4.0)*pow((sin(th0)),2)*(sin(2*phi0))/pow(rmag0,3);
   y22c.real=y22c.real+sqrt(3.0/4.0)*pow((sin(th0)),2)*(cos(2*phi0))/pow(rmag0,3);
   y22c.imaginary=y22c.imaginary+sqrt(3.0/4.0)*pow((sin(th0)),2)*(-sin(2*phi0))/pow(rmag0,3);
   y21.real=y21.real+sqrt(3.0)*cos(th0)*sin(th0)*cos(phi0)/pow(rmag0,3);
   y21.imaginary=y21.imaginary+sqrt(3.0)*cos(th0)*sin(th0)*sin(phi0)/pow(rmag0,3);
   y21c.real=y21c.real+sqrt(3.0)*cos(th0)*sin(th0)*cos(phi0)/pow(rmag0,3);
   y21c.imaginary=y21c.imaginary+sqrt(3.0)*cos(th0)*sin(th0)*(-sin(phi0))/pow(rmag0,3);
   y20=y20+(pow((3*(cos(th0))),2)-1)/(2.0*pow(rmag0,3));

   for (j=1; j<=jmax; ++j)  {
     indt=cfind[j];
     rmagt=pow(dipcrd1[indt],2)+pow(dipcrd2[indt],2)+pow(dipcrd3[indt],2);
     rmagt=sqrt(rmagt);
     xt=dipcrd1[indt]/rmagt;
     yt=dipcrd2[indt]/rmagt;
     zt=dipcrd3[indt]/rmagt;
     dot=(3*pow((x0*xt+y0*yt+z0*zt),2)-1)/2.0;
     corr[j]=corr[j]+dot/pow((rmag0*rmagt),3);
     p2[j]=p2[j]+dot;
     rcorr[j]=rcorr[j]+1/pow((rmag0*rmagt),3);
   }

   if (ncorr != npts)  {
     for (j=1; j<=jmax-1; ++j)  {
       cfind[j]=cfind[j+1];
     }
     cfind[jmax]=ind0;
     current_frame=current_frame+1;
     if (current_frame < tf)  {
       dipcrd1[current_frame]=x[current_frame];
       dipcrd2[current_frame]=y[current_frame];
       dipcrd3[current_frame]=z[current_frame];
       npts=npts+1;
     }
     else  {
       jmax=jmax-1;
     }
   }
   else  {
     doit=0;
   }
  }

/* normalize correlation functions */

  r6ave=r6ave/npts;
  r3ave=r3ave/npts;
  rave=rave/npts;
  avecrd[1]=avecrd[1]/npts;
  avecrd[2]=avecrd[2]/npts;
  avecrd[3]=avecrd[3]/npts;
  rrig=pow(avecrd[1],2)+pow(avecrd[2],2)+pow(avecrd[3],2);
  rrig=sqrt(rrig);

  y2asy=(y22.real*y22c.real+y21.real*y21c.real)+pow(y20,2);
  y2asy=y2asy/(npts*npts*r6ave);

  for (i=1; i<=ntot; ++i)  {
    corr[i]=corr[i]/((npts-i+1)*r6ave);
    rcorr[i]=rcorr[i]/((npts-i+1)*r6ave);
    p2[i]=p2[i]/(npts-i+1);
  }

/* output correlation functions */
  ifp=safe_fopen(outfile, "w");
  if (ifp == NULL) {
    warning("ptraj(), correlation: cannot open output file %s\n",
	    outfile);
  } else {
    fprintf(ifp, "# Rrigid= %lf  Rave= %lf \n", rrig, rave);
    fprintf(ifp, "# <1/r^3>= %lf  <1/r^6>= %lf\n", r3ave, r6ave);
    /*
    rfac = r6ave*pow(rave,6);
    qfac = y2asy*rfac;
    */
    fprintf(ifp, "#   time     C(t)      P2(t)      R(t)\n");
    i=tc/ts;
    for (j=1; j<=i; ++j)  {
      fprintf(ifp, "%d   %lf   %lf   %lf\n",
	      (j-1)*ts, corr[j], p2[j], rcorr[j]);
    }
    safe_fclose(ifp);
  }

/* deallocate space */

  safe_free(dipcrd1);
  safe_free(dipcrd2);
  safe_free(dipcrd3);
  safe_free(p2);
  safe_free(corr);
  safe_free(rcorr);
  safe_free(cfind);

}

// transformCorr()
int transformCorr(actionInformation *action, 
	      double *x, double *y, double *z, 
	      double *box, int mode)
{
  //char *name = "correlation";
  argStackType **argumentStackPointer;
  char *buffer;
  transformCorrInfo *corrInfo;
  ptrajState *state;
  int i;
  double cx, cy, cz, total_mass;
  double vx, vy, vz;


  /*
   *  USAGE:
   *
   *     correlation name mask1 mask2 tmin tcorr tmax [out filename]
   *
   *  action argument usage:
   *
   *  carg1:
   *     a transformCorrInfo structure
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

    /*
     *  set up complex argument
     */
    corrInfo = (transformCorrInfo *) safe_malloc(sizeof(transformCorrInfo));
    INITIALIZE_transformCorrInfo(corrInfo);
    corrInfo->totalFrames = -1;

    corrInfo->name = getArgumentString(argumentStackPointer, NULL);

    buffer = getArgumentString(argumentStackPointer, NULL);
    corrInfo->mask = processAtomMask(buffer, action->state);
    safe_free(buffer);

    buffer = getArgumentString(argumentStackPointer, NULL);
    corrInfo->mask2 = processAtomMask(buffer, action->state);
    safe_free(buffer);
    corrInfo->mode = VECTOR_MASK;

    corrInfo->tmin  = getArgumentInteger(argumentStackPointer, 1.0);
    corrInfo->tcorr = getArgumentInteger(argumentStackPointer, 1.0);
    corrInfo->tmax  = getArgumentInteger(argumentStackPointer, 1.0);

    /*
     *  assume "out" may be missing
     */

    corrInfo->filename = argumentStackKeyToString(argumentStackPointer, "out", NULL);
    if (corrInfo->filename == NULL) {
      corrInfo->filename = getArgumentString(argumentStackPointer, NULL);
      if (corrInfo->filename == NULL) {
	error("ptraj()", "correlation, no out file specified\n");
      }
    }
    if (corrInfo->name == NULL || corrInfo->mask == NULL ||
	corrInfo->mask2 == NULL) {
      error("ptraj()", "correlation arguments\n");
    }

    action->carg1 = (void *) corrInfo;

    return 0;
  }


  corrInfo = (transformCorrInfo *) action->carg1;

  if (mode == PTRAJ_STATUS) {

    /*
     *  ACTION: PTRAJ_STATUS
     */
    fprintf(stdout, "  CORRELATION: storage to array named %s",
            corrInfo->name);
    fprintf(stdout, " -- tmin: %i tcorr: %i tmax: %i\n",
	    corrInfo->tmin, corrInfo->tcorr, corrInfo->tmax);
    fprintf(stdout, "      Atom selection 1 is ");
    printAtomMask(stdout, corrInfo->mask, action->state);
    fprintf(stdout, "\n");
    fprintf(stdout, "      Atom selection 2 is ");
    printAtomMask(stdout, corrInfo->mask2, action->state);
    fprintf(stdout, "\n");

  } else if (mode == PTRAJ_PRINT) {

    /*
     *  ACTION: PTRAJ_PRINT
     */

    fprintf(stdout, "PTRAJ CORRELATION: calculating correlation %s\n",
	    corrInfo->name);
    if (corrInfo != NULL) {
      compute_corr(corrInfo->filename, corrInfo->vx, corrInfo->vy, corrInfo->vz,
		   corrInfo->tmin, corrInfo->tcorr, corrInfo->tmax, corrInfo->totalFrames);
    }
    return 0;

  } else if (mode == PTRAJ_CLEANUP) {

    /*
     *  ACTION: PTRAJ_CLEANUP
     */

    safe_free(corrInfo->cx);
    safe_free(corrInfo->cy);
    safe_free(corrInfo->cz);
    safe_free(corrInfo->vx);
    safe_free(corrInfo->vy);
    safe_free(corrInfo->vz);
    safe_free(corrInfo->mask);
    safe_free(corrInfo->mask2);
    safe_free(corrInfo->name);
    INITIALIZE_transformCorrInfo(corrInfo);

    safe_free(corrInfo);
  }


  if (mode != PTRAJ_ACTION) return 0;

  /*
   *  ACTION: PTRAJ_ACTION
   */

  state = (ptrajState *) action->state;

  if (corrInfo->totalFrames < 0) {
    corrInfo->totalFrames = state->maxFrames;
    corrInfo->cx = (double *)
      safe_malloc(sizeof(double) * corrInfo->totalFrames);
    corrInfo->cy = (double *)
      safe_malloc(sizeof(double) * corrInfo->totalFrames);
    corrInfo->cz = (double *)
      safe_malloc(sizeof(double) * corrInfo->totalFrames);
    corrInfo->vx = (double *)
      safe_malloc(sizeof(double) * (corrInfo->totalFrames+1));
    corrInfo->vy = (double *)
      safe_malloc(sizeof(double) * (corrInfo->totalFrames+1));
    corrInfo->vz = (double *)
      safe_malloc(sizeof(double) * (corrInfo->totalFrames+1));
  }


  if (corrInfo->frame > corrInfo->totalFrames) {
    warning("transformCorrr()", "Blowing array; too many frames!!\n");
    return 0;
  }

  corrInfo->mode = CORR_MASK;
  total_mass = 0.0;
  cx = 0.0;
  cy = 0.0;
  cz = 0.0;
  for (i=0; i < state->atoms; i++) {
    if (corrInfo->mask[i]) {
        cx += state->masses[i] * x[i];
        cy += state->masses[i] * y[i];
        cz += state->masses[i] * z[i];
        total_mass += state->masses[i];
    }
  }
  cx = cx / total_mass;
  cy = cy / total_mass;
  cz = cz / total_mass;

  total_mass = 0.0;
  vx = 0.0;
  vy = 0.0;
  vz = 0.0;
  for (i=0; i < state->atoms; i++) {
    if (corrInfo->mask2[i]) {
        vx += state->masses[i] * x[i];
        vy += state->masses[i] * y[i];
        vz += state->masses[i] * z[i];
        total_mass += state->masses[i];
    }
  }
  vx = vx / total_mass;
  vy = vy / total_mass;
  vz = vz / total_mass;

  corrInfo->vx[corrInfo->frame] = vx - cx;
  corrInfo->vy[corrInfo->frame] = vy - cy;
  corrInfo->vz[corrInfo->frame] = vz - cz;
  corrInfo->cx[corrInfo->frame] = cx;
  corrInfo->cy[corrInfo->frame] = cy;
  corrInfo->cz[corrInfo->frame] = cz;

  corrInfo->frame++;
  return 1;
}


/** ACTION ROUTINE ************************************************************
 *
 *  transformDihedralCluster()  --- Cluster trajectory by dihedral angles
 *
 *    Written by Dan Roe, Stonybrook U.
 *
 ******************************************************************************/
typedef struct _DCnodetype {
  struct _DCnodetype **branch;
  int *bin;
  long int numBranch;
  long int *count;
  double **frames;
} DCnodetype;

#define INITIALIZE_transformDihedralCluster(_p_) \
  _p_->branch = NULL; \
  _p_->bin = NULL; \
  _p_->numBranch = 0; \
  _p_->count = NULL; \
  _p_->frames = NULL ;

typedef struct _DCarray {
  int* Bins;
  long int count;
  double* frames;
} DCarray;
/*NOTE: Place checks after the realloc statements*/

/* DCbin: recursive function to place diedral bin combinations in a Tree
 * structure - this conserves memory and has better search performance than
 * a linear array.
 */
static int DCbin(int* Bins,int level, DCnodetype* T, int max, double FRAME,long int* numCluster)
{
  long int i,j;
  /*DEBUG Output*/
  if (prnlev>1) {
    if (level==0) printf("DIHEDRAL CLUSTER, FRAME %lf\n",FRAME);
    for (i=0; i<(level*2); i++)
      printf(" ");
    printf("Level %i Bin %i ",level,Bins[level]); 
  }
  /*Is this the first time bin has been called?*/
  if (T->numBranch==0) {
    if (prnlev>1) printf("Creating. ");
    T->numBranch=1;
    T->bin=(int*) safe_malloc(sizeof(int));
    T->bin[0]=Bins[level];
    /*If this is the lowest level allocate count and frames*/
    if (level==max) {
      T->branch=NULL;
      T->count=(long int*) safe_malloc(sizeof(long int));
      T->count[0]=1;
      (*numCluster)++;
      T->frames=(double**) safe_malloc(sizeof(double*));
      T->frames[0]=(double*) safe_malloc(sizeof(double));
      T->frames[0][0]=FRAME;
      if (prnlev>1) printf("Count= %li, Frame= %.0lf\n",T->count[0],T->frames[0][0]);
      return 0;
    /*Otherwise create a branch node for the next bin*/
    } else {
      T->count=NULL;
      T->frames=NULL;
      T->branch=(DCnodetype**) safe_malloc(sizeof(DCnodetype*));
      T->branch[0]=(DCnodetype*) safe_malloc(sizeof(DCnodetype));
      INITIALIZE_transformDihedralCluster(T->branch[0]);
      if (prnlev>1) printf("Next.\n");
      DCbin(Bins,level+1,T->branch[0],max,FRAME,numCluster);
      return 0;
    }
  } else {
  /*If not first call, does this Bin already exist?*/  
    if (prnlev>1) printf("Searching. ");
    for (i=0; i<T->numBranch; i++) {
      if (T->bin[i]==Bins[level]) {
        /*If it does and we're at lowest level, increment count, record frame*/ 
        if (level==max) {
          T->count[i]++;
          j=T->count[i];
          T->frames[i]=(double*) realloc(T->frames[i],j*sizeof(double));
          T->frames[i][j-1]=FRAME;
          if (prnlev>1) printf("Count= %li, Frame= %.0lf\n",j,T->frames[i][j-1]);
          return 0;
        } else {
        /*If it does and we're not at lowest level, continue search*/
          if (prnlev>1) printf("Next.\n");
          DCbin(Bins,level+1,T->branch[i],max,FRAME,numCluster);
          return 0;
        }
      }
    }
    /*Bin doesnt exist, create a new branch*/
    if (prnlev>1) printf("Newbranch. ");
    T->numBranch++;
    i=T->numBranch;
    T->bin=(int*) realloc(T->bin,i*sizeof(int));
    T->bin[i-1]=Bins[level];
    if (level==max) {
      /*If lowest level, set count and frame number*/
      T->branch=NULL;
      T->count=(long int*) realloc(T->count,i*sizeof(long int));
      T->count[i-1]=1;
      (*numCluster)++;
      T->frames=(double**) realloc(T->frames,i*sizeof(double*));
      T->frames[i-1]=(double*) safe_malloc(sizeof(double));
      T->frames[i-1][0]=FRAME;
      if (prnlev>1) printf("Count= %li, Frame= %.0lf\n",T->count[i-1],T->frames[i-1][0]);
      return 0;
    } else {
      /*Otherwise, continue down the branch*/
      T->count=NULL;
      T->frames=NULL;
      T->branch=(DCnodetype**) realloc(T->branch,i*sizeof(DCnodetype*));
      T->branch[i-1]=(DCnodetype*) safe_malloc(sizeof(DCnodetype));
      INITIALIZE_transformDihedralCluster(T->branch[i-1]);
      if (prnlev>1) printf("Next.\n");
      DCbin(Bins,level+1,T->branch[i-1],max,FRAME,numCluster);
      return 0;
    }
  }
  return 1;    
}

/* freeDCbin: recursive function to free memory allocated to Tree
 */
static void freeDCbin(DCnodetype* T) {
  long int i;

  if (T->branch!=NULL)
    for (i=0; i<T->numBranch; i++) { 
      freeDCbin(T->branch[i]);
      safe_free(T->branch[i]);
    }
  if (T->count!=NULL)
    safe_free(T->count);
  if (T->frames!=NULL) {
    for (i=0; i<T->numBranch; i++)
      safe_free(T->frames[i]);
    safe_free(T->frames);
  }
  safe_free(T->bin);
  safe_free(T->branch);
  return;
}

/* DCtree2array: recursive function to put contents of tree into array for sorting.
 */
static void DCtree2array(int* Bins, int level, DCnodetype* T,long int* numCluster,DCarray** A) {
  long int i,j,k;

  for (i=0; i<T->numBranch; i++) {
    Bins[level]=T->bin[i];
    if (T->branch!=NULL)
      DCtree2array(Bins,level+1,T->branch[i],numCluster,A);
    else {
      (*numCluster)++;
      k=(*numCluster);
      /*Store cluster info in array for sorting later*/
      if (prnlev>1) printf("DEBUG: Writing cluster %li\n",k);
      /*Assign Bins*/
      if (prnlev>1) printf("DEBUG: Bins= ");
      for (j=0; j<=level; j++) {
        A[k-1]->Bins[j]=Bins[j];
        if (prnlev>1) printf("%i ",A[k-1]->Bins[j]);
      }
      if (prnlev>1) printf("\n");
      /*Assign count*/
      A[k-1]->count=T->count[i];
      if (prnlev>1) printf("DEBUG: Count= %li\n",A[k-1]->count);
      /*Assign frames*/
      A[k-1]->frames=(double*) safe_malloc(T->count[i]*sizeof(double));
      if (prnlev>1) printf("DEBUG: frames= ");
      for (j=0; j<T->count[i]; j++){
        A[k-1]->frames[j]=T->frames[i][j];
        if (prnlev>1) printf("%.0lf ",A[k-1]->frames[j]);
      }
      if (prnlev>1) printf("\n");
      /*printf("DEBUG: Current Counts\n");
      for (j=0; j<k; j++)
        printf("  Count %i = %i\n",j,A[j]->count);*/
    }
  }
  return;
}

/*compare: used by the qsort function in sorting cluster array
*/
static int compare(const void *a, const void *b) {
  DCarray *A;
  DCarray *B;

  A=*(DCarray**) a;
  B=*(DCarray**) b;
  /*printf("    QSORT DEBUG: Comparing %i to %i\n",A->count,B->count);*/
  return ( B->count - A->count );
}  

// transformDihedralCluster
/** Usage: clusterdihedral
 *        [phibins N] number of bins in phi direction
 *        [psibins M] number of bins in psi direction
 *        Note: phibins and psibins only used if dihedralfile not specified
 *        [out <FILE>] file to print cluster information to
 *        [cut CUT] only print clusters with pop > CUT
 *        [framefile <FILE>] file to print frame-cluster info
 *        [clusterinfo <FILE>] print cluster info in format sander can read
 *                             for NB weighted RREMD
 *        [dihedralfile] read dihedral definitions from FILE. Format is 
 *                       ATOM#1 ATOM#2 ATOM#3 ATOM#4 BINS
 *                       Note: This functionality treats atom numbers as 
 *                       starting from 1 to be consistent with sander.
 *        [<MASK>] if not reading dihedrals from file Backbone dihedrals will
 *                 be searched for within MASK.

 * Action argument usage:
 * iarg1 = number of dihedral angles to keep track of
 * iarg2 = if 1, dihedral angles were read from a file, otherwise backbone
 *         dihedrals were searched for in MASK.
 * iarg4 = store CUT
 * darg3 = keep track of number of frames, used for memory allocation at end
 * darg4 = number of clusters
 * carg1 = int array; atom masks of each dihedral and the number of bins
 * carg2 = hold the DCnodetype tree
 * carg3 = int array; during run, store which bin each dihedral falls into
 *         Note: Even though no information from carg3 needs to be carried over
 *         from frame to frame, it is stored here to avoid reallocating the 
 *         memory each time.      
 * carg4 = char* array to hold output filenames
 */
int transformDihedralCluster(actionInformation *action,
                             double *x, double *y, double *z, double *box,
                             int mode)
{
  //char *name = "clusterdihedral";
  argStackType **argumentStackPointer;
  ptrajState *state;
  char *buffer = NULL;
  DCnodetype* T;
  DCarray** A;
  int** DCmasks = NULL;
  int* Bins;
  long int* framecluster;
  long int i,j,k;
  int C1,N2,CA,C2;
  int phibins,psibins;
  int numDihedral,CUT;
  long int numCluster;
  double FRAME;
  double phistep, PHI, temp;
  double cx1[3];
  double cx2[3];
  double cx3[3];
  double cx4[3];
  char** filenames;
  FILE* outfile;

  state = (ptrajState *) action->state;
  /*
   * ---=== PTRAJ_SETUP ===--- 
   */
  if (mode == PTRAJ_SETUP) {
#ifdef MPI
    printParallelError(name);
    return -1;
#endif
    argumentStackPointer = (argStackType **) action->carg1;
    action->carg1 = NULL;
    /*Parse Command Line*/
    phibins=argumentStackKeyToInteger(argumentStackPointer, "phibins", 10);
    if ((phibins>360)||(phibins<=1)) {
      fprintf(stdout,"clusterdihedral Error: phibins out of range 360 < x <= 1 (%i)\n",phibins);
      return -1;
    }
    psibins=argumentStackKeyToInteger(argumentStackPointer, "psibins", 10);
    if ((psibins>360)||(psibins<=1)) {
      fprintf(stdout,"clusterdihedral Error: psibins out of range 360 < x <= 1 (%i)\n",psibins);
      return -1;
    }
    /*Cluster cutoff*/
    CUT=argumentStackKeyToInteger(argumentStackPointer, "cut",0);
    /*Output Files*/
    filenames=(char**) safe_malloc(4*sizeof(char*));
    buffer=argumentStackKeyToString(argumentStackPointer, "out",NULL);
    if (buffer!=NULL) {
      filenames[0]=(char*) safe_malloc((strlen(buffer)+1)*sizeof(char));
      strcpy(filenames[0],buffer);
      safe_free(buffer);
    } else filenames[0]=NULL;
    buffer=argumentStackKeyToString(argumentStackPointer, "framefile",NULL);
    if (buffer!=NULL) {
      filenames[1]=(char*) safe_malloc((strlen(buffer)+1)*sizeof(char));
      strcpy(filenames[1],buffer);
      safe_free(buffer);
    } else filenames[1]=NULL;
    buffer=argumentStackKeyToString(argumentStackPointer, "clusterinfo",NULL);
    if (buffer!=NULL) {
      filenames[2]=(char*) safe_malloc((strlen(buffer)+1)*sizeof(char));
      strcpy(filenames[2],buffer);
      safe_free(buffer);
    } else filenames[2]=NULL;
    buffer=argumentStackKeyToString(argumentStackPointer, "clustervtime",NULL);
    if (buffer!=NULL) {
      filenames[3]=(char*) safe_malloc((strlen(buffer)+1)*sizeof(char));
      strcpy(filenames[3],buffer);
      safe_free(buffer);
    } else filenames[3]=NULL;
    action->carg4=(void*) filenames;
    /*Input Dihedral file*/
    action->iarg2=0;
    numDihedral=0;
    buffer=argumentStackKeyToString(argumentStackPointer, "dihedralfile",NULL);
    if (buffer!=NULL) {
      if ( (outfile=safe_fopen(buffer,"r"))==NULL ) {
        fprintf(stdout,"WARNING: Could not open dihedralfile %s",buffer);
      } else {
        /*Read Dihedrals and bins from a file*/
        printf("Reading Dihedrals from %s\n",buffer);
        DCmasks=(int**) safe_malloc(sizeof(int*));
        while ( fscanf(outfile,"%i %i %i %i %i",&C1,&N2,&CA,&C2,&phibins) != EOF) {
          numDihedral++;
          DCmasks=(int**) realloc(DCmasks,numDihedral*sizeof(int*));
          if (DCmasks==NULL) {
            fprintf(stdout,"clusterdihedral Error: Memory reallocation for masks failed.\n");
            return -1;
          }
          DCmasks[numDihedral-1]=(int*) safe_malloc(5*sizeof(int));
          /*amber atom numbers start from 1*/
          DCmasks[numDihedral-1][0]=C1-1;
          DCmasks[numDihedral-1][1]=N2-1;
          DCmasks[numDihedral-1][2]=CA-1;
          DCmasks[numDihedral-1][3]=C2-1;
          DCmasks[numDihedral-1][4]=phibins;
          printf("(%i)-(%i)-(%i)-(%i) Bins=%i\n",C1,N2,CA,C2,phibins);
        }
        printf("Read %i dihedrals\n",numDihedral);
        action->iarg2=1;
        safe_fclose(outfile);
      }
      safe_free(buffer);
    }
    /*Allocate Memory to hold cluster info*/
    action->darg4=0;
    T=(DCnodetype*) safe_malloc(sizeof(DCnodetype));
    INITIALIZE_transformDihedralCluster(T);
    FRAME=0;
    /* Process Mask */
    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL) {
      action->mask = processAtomMask("*", action->state);
    } else {
      action->mask = processAtomMask(buffer, action->state);
      safe_free(buffer);
    }
    /*Set up backbone dihedral angles if none were read*/
    if (numDihedral==0) {
      printf("Scanning for backbone dihedrals.\n");
      C1=-1; N2=-1; CA=-1; C2=-1;
      DCmasks=(int**) safe_malloc(sizeof(int*));
      for (i=0; i < state->atoms; i++) {
        if (action->mask[i]==1) {
          /*printf("  DEBUG: Atom %i: %s\n",i,state->atomName[i]);*/
          if (C2>-1) {
            /* If we have already found the last C in phi dihedral, this N is 
             * the last atom in psi dihedral - store both.
             */
            if ( strcmp(state->atomName[i],"N   ")==0 ) {
              numDihedral+=2;
              /* dynamically grow array */
              DCmasks=(int**) realloc(DCmasks,numDihedral*sizeof(int*));
              if (DCmasks==NULL) {
                fprintf(stdout,"clusterdihedral Error: Memory reallocation for masks failed.\n");
                return -1;
              }
              DCmasks[numDihedral-2]=(int*) safe_malloc(5*sizeof(int));
              DCmasks[numDihedral-2][0]=C1;
              DCmasks[numDihedral-2][1]=N2;
              DCmasks[numDihedral-2][2]=CA;
              DCmasks[numDihedral-2][3]=C2;
              DCmasks[numDihedral-2][4]=phibins;
              DCmasks[numDihedral-1]=(int*) safe_malloc(5*sizeof(int));
              DCmasks[numDihedral-1][0]=N2;
              DCmasks[numDihedral-1][1]=CA;
              DCmasks[numDihedral-1][2]=C2;
              DCmasks[numDihedral-1][3]=i;
              DCmasks[numDihedral-1][4]=psibins;
              if (prnlev>0) printf("DIHEDRAL PAIR FOUND: C1= %i, N2= %i, CA= %i, C2= %i, N3= %li\n",
                                   C1,N2,CA,C2,i);
              /* Since the carbonyl C/amide N probably starts a new dihedral,
               * reset to those
               */
              C1=C2;
              N2=i;
              C2=-1; CA=-1;
            }
          } else if (C1>-1) {
            /* If we've already found the first carbonyl, look for other atoms
             * in the dihedral pair.
             */
            if ( strcmp(state->atomName[i],"N   ")==0 )  N2=i;
            if ( strcmp(state->atomName[i],"CA  ")==0 ) CA=i;
            if ( strcmp(state->atomName[i],"C   ")==0 ) C2=i;
          } else if ( strcmp(state->atomName[i],"C   ")==0 ) C1=i; /*1st carbon*/
        }
      }
    }
    if ( numDihedral == 0 ) {
      fprintf(stdout,"clusterdihedral Error: No Backbone Dihedral Angles Found!\n");
      return -1;
    }
    if (prnlev>0) printf("FOUND %i DIHEDRAL ANGLES.\n",numDihedral);
    /* Allocate memory to store dihedral bins so we don't continuously
     * reallocate during the ptraj run. 
     */
    Bins=(int*) safe_malloc(numDihedral*sizeof(int));
    /* Assign action pointers */
    action->iarg1 = numDihedral;
    action->iarg4 = CUT;
    action->carg1 = (void*) DCmasks;
    action->carg2 = (void*) T;
    action->carg3 = (void*) Bins;
    action->darg3 = FRAME;
    return 0;
  }     
  /*
   * ---=== PTRAJ_CLEANUP ===--- 
   */
  if (mode == PTRAJ_CLEANUP) {
    numDihedral=action->iarg1;
    DCmasks=(int**) action->carg1;
    T=(DCnodetype*) action->carg2;
    Bins=(int*) action->carg3;
    filenames=(char**) action->carg4;
    if (filenames[0]!=NULL) safe_free(filenames[0]);
    if (filenames[1]!=NULL) safe_free(filenames[1]);
    if (filenames[2]!=NULL) safe_free(filenames[2]);
    if (filenames[3]!=NULL) safe_free(filenames[3]);
    safe_free(filenames);
    if (prnlev>1) printf("DIHEDRAL: Freeing %i dihedral masks.\n",numDihedral);
    for (i=0; i<numDihedral; i++)
      safe_free(DCmasks[i]);
    safe_free(DCmasks);
    freeDCbin(T);
    safe_free(T);
    safe_free(Bins);
    safe_free(action->mask);
    return 0;
  }
  /*
   * ---=== PTRAJ_STATUS ===--- 
   */
  if (mode == PTRAJ_STATUS) {
    numDihedral=action->iarg1;
    DCmasks=(int**) action->carg1;
    CUT=action->iarg4;
    filenames=(char**) action->carg4;
    fprintf(stdout,"\n  DIHEDRAL CLUSTERING on %i dihedral angles\n",numDihedral);
    if (filenames[0]!=NULL)
      fprintf(stdout,"  Cluster data will be output to %s\n",filenames[0]);
    if (CUT>0)
      fprintf(stdout,"  Only clusters with pop > %i will be printed.\n",CUT); 
    if (filenames[1]!=NULL)
      fprintf(stdout,"  Frame-Cluster data will be output to %s\n",filenames[1]);
    if (filenames[2]!=NULL)
      fprintf(stdout,"  Cluster information (pop. & ID) will be output to %s\n",filenames[2]);
    if (filenames[3]!=NULL)
      fprintf(stdout,"  Number of clusters v time will be output to %s\n",filenames[3]);
    if (action->iarg2==0) {
      fprintf(stdout, "  Looked for dihedrals within atom selection= ");
      printAtomMask(stdout, action->mask, action->state);
      fprintf(stdout, "\n");
    } else fprintf(stdout,"  Read dihedrals from an input file.\n");
    for (i=0; i<numDihedral; i++) {
      fprintf(stdout,"    %6li ",i);
      for (j=0; j<4; j++) {
        C1=DCmasks[i][j];
        /*amber atom nums start at 1*/
        fprintf(stdout,"%i(%-s)",C1+1,state->atomName[C1]);
      }
      fprintf(stdout," [Bins=%i]",DCmasks[i][4]);
      phistep=360/DCmasks[i][4];
      fprintf(stdout," [Step=%lf]\n",phistep);
    }
    fprintf(stdout,"\n");
    return 0;
  }
  /*
   * ---=== PTRAJ_PRINT ===--- 
   */
  if (mode == PTRAJ_PRINT) {
    DCmasks=(int**) action->carg1;
    filenames=(char**) action->carg4;
    numCluster=(long int) action->darg4;
    numDihedral=action->iarg1;
    CUT=action->iarg4;
    T=(DCnodetype*) action->carg2;
    Bins=(int*) action->carg3;
    /*Set up output file*/
    if (filenames[0]!=NULL) {
      if ( (outfile=safe_fopen(filenames[0],"w"))==NULL ) {
        fprintf(stdout,"Could not open %s, reverting to stdout.",filenames[0]);
        outfile=stdout;
      }
    } else outfile=stdout;
    /*Print Bin information*/
    printf("Printing Dihedral Clustering Results.\n");
    fprintf(outfile,"DIHEDRAL CLUSTER RESULTS");
    if (action->iarg2==0) {
      fprintf(outfile," for ");
      printAtomMask(outfile, action->mask, action->state);
      fprintf(outfile,"\n");
    } else fprintf(outfile,"\n");
    if (outfile!=stdout) {
      for (i=0; i<numDihedral; i++) {
        fprintf(outfile,"    %6li ",i);
        for (j=0; j<4; j++) {
          C1=DCmasks[i][j];
          /*shift atom numbers by 1, amber standard*/
          fprintf(outfile,"%-s(%i)",state->atomName[C1],C1+1);
        }
        fprintf(outfile," [Bins=%i]\n",DCmasks[i][4]);
      }
    }
    fprintf(outfile,"%li clusters.\n",numCluster);
    if (CUT>0) 
      fprintf(outfile,"Only printing clusters with pop > %i\n",CUT);
/*    fprintf(outfile,"Phi Bin Ranges:\n");
    for (i=0; i<phibins; i++) 
      fprintf(outfile,"  %8.2lf <= %3i < %8.2lf\n",(i*phistep)-180,i,((i+1)*phistep)-180);
    fprintf(outfile,"Psi Bin Ranges:\n");
    for (i=0; i<psibins; i++)
      fprintf(outfile,"  %8.2lf <= %3i < %8.2lf\n",(i*psistep)-180,i,((i+1)*psistep)-180);*/
    /*Allocate memory for array used in sorting*/
    A=(DCarray**) safe_malloc(numCluster*sizeof(DCarray*));
    for (j=0; j<numCluster; j++) {
      A[j]=(DCarray*) safe_malloc(sizeof(DCarray));
      A[j]->Bins=(int*) safe_malloc(numDihedral*sizeof(int));
    }
    /*Place tree elements into an array for sorting, then sort*/
    i=0;
    DCtree2array(Bins,0,T,&i,A);
    if (prnlev>1) printf("%li clusters written.\n",i);
    qsort(A,numCluster,sizeof(DCarray*),compare);
    /*Allocate memory for printing frame/cluster */
    i=(long int)action->darg3;
    framecluster=(long int*) safe_malloc((i+1)*sizeof(long int));
    /*Print sorted cluster array*/
    for (j=0; j<numCluster; j++) {
      if (A[j]->count>CUT) {
        fprintf(outfile,"Cluster %10li %10li [ ",j,A[j]->count);
        for (i=0; i<numDihedral; i++)
          fprintf(outfile,"%3i ",A[j]->Bins[i]);
        fprintf(outfile," ]\n");
        for (i=0; i<A[j]->count; i++) {
          fprintf(outfile,"%.0lf ",A[j]->frames[i]);
          /* store which cluster each frame belongs to. Not neccesary if user
           * didn't specify this option, but avoids a second loop if they did.
           */
          k=(long int) A[j]->frames[i];
          framecluster[k]=j;
        }
        fprintf(outfile,"\n");
      }
    }
    if (outfile!=NULL && outfile!=stdout) safe_fclose(outfile);
    /*cluster for each frame*/
    if (filenames[1]!=NULL) {
      printf("Printing cluster number for each frame.\n");
      if ( (outfile=safe_fopen(filenames[1],"w"))==NULL ) {
        fprintf(stdout,"WARNING: Could not open framefile %s",filenames[1]);
      } else {
        i=(long int) action->darg3;
        for (j=1; j<=i; j++) {
          C1=framecluster[j];
          fprintf(outfile,"%10li %10i %10li ",j,C1,A[C1]->count);
          for (k=0; k<numDihedral; k++)
            fprintf(outfile,"%03i",A[C1]->Bins[k]);
          fprintf(outfile,"\n");
        }
        safe_fclose(outfile);
      }
    }
    /*cluster information file*/
    if (filenames[2]!=NULL) {
      printf("Printing cluster information.\n");
      if ( (outfile=safe_fopen(filenames[2],"w"))==NULL ) {
        fprintf(stdout,"WARNING: Could not open clusterinfo file %s",filenames[2]);
      } else {
        fprintf(outfile,"%i\n",numDihedral);
        for (i=0; i<numDihedral; i++) {
          for (j=0; j<4; j++) {
            C1=DCmasks[i][j];
            /*shift atom numbers by 1, amber standard*/
            fprintf(outfile,"%10i ",C1+1);
          }
          fprintf(outfile,"%10i\n",DCmasks[i][4]);
        }
        fprintf(outfile,"%li\n",numCluster);
        for (i=0; i<numCluster; i++) {
          fprintf(outfile,"%10li %10li ",i,A[i]->count);
          for (k=0; k<numDihedral; k++)
            fprintf(outfile,"%03i",A[i]->Bins[k]);
          fprintf(outfile,"\n");
        }
        safe_fclose(outfile);
      }
    }
    /*Cleanup arrays*/
    for (j=0; j<numCluster; j++) {
      safe_free(A[j]->Bins);
      safe_free(A[j]->frames);
      safe_free(A[j]);
    }
    safe_free(A);
    safe_free(framecluster);
    return 0;
  }
  /*
   * ---=== PTRAJ_ACTION ===--- 
   */
  if (mode == PTRAJ_ACTION) {
    numDihedral=action->iarg1;
    DCmasks=(int**) action->carg1;
    T=(DCnodetype*) action->carg2;
    Bins=(int*) action->carg3;
    FRAME=action->darg3;
    FRAME++;
    numCluster=(long int) action->darg4;
    /*For each dihedral, calculate which bin it should go into and store bin#*/
    j=0;
    for (i=0; i<numDihedral; i++) {
      C1=DCmasks[i][0];
      N2=DCmasks[i][1];
      CA=DCmasks[i][2];
      C2=DCmasks[i][3];
      cx1[0]=x[C1]; cx1[1]=y[C1]; cx1[2]=z[C1];
      cx2[0]=x[N2]; cx2[1]=y[N2]; cx2[2]=z[N2];
      cx3[0]=x[CA]; cx3[1]=y[CA]; cx3[2]=z[CA];
      cx4[0]=x[C2]; cx4[1]=y[C2]; cx4[2]=z[C2];
      PHI=Torsion(cx1,cx2,cx3,cx4) * RADDEG; // Torsion is in radians
      if (prnlev>0) printf("%9s%8.2lf    ","Dihedral=",PHI);
      PHI+=180;
      phistep=360/DCmasks[i][4];
      PHI/=phistep;
      modf(PHI,&temp);
      phibins=temp;
      if (prnlev>0) printf("%4s%3i\n","Bin=",phibins);
      Bins[j]=phibins;
      j++;
    }
    /* At this point j=numDihedral */
    if (prnlev>0) {
      printf("[");
      for (i=0; i<j; i++)
        printf("%3i,",Bins[i]);
      printf("]\n");
    }
    /* Now place this combo in the tree, or if it is already there increment
     * the counter.
     */
    DCbin(Bins,0,T,j-1,FRAME,&numCluster);
    action->darg4=(double) numCluster;
    action->carg2=(void*) T;
    action->darg3=FRAME;
    return 1;
  }
  return 0;
}
/*DAN ROE*/

/** ACTION ROUTINE *************************************************************
 *
 *  transformDiffusion()   --- calculate mean squared displacements vs. time
 *
 ******************************************************************************/
typedef struct _transformDiffusionInfo {
  double *dx;           /* reference coordinates: MUST BE SET */
  double *dy;           /*  TO NULL IN transformSetup         */
  double *dz;
  double timePerFrame;  /*  in picoseconds */
  double *prevx;        /*  the previous x for active atoms */
  double *prevy;        /*  the previous y                  */
  double *prevz;        /*  the previous z                  */
  double *distancex;    /*  the active atoms distances      */
  double *distancey;    /*  the active atoms distances      */
  double *distancez;    /*  the active atoms distances      */
  double *distance;     /*  the active atoms distances      */
  double *deltax;
  double *deltay;
  double *deltaz;
  int activeAtoms;      /*  the total number of active atoms in the mask */
  int elapsedFrames;    /*  number of frames so far... */
  char *outputFilenameRoot;
  FILE *outputx;
  FILE *outputy;
  FILE *outputz;
  FILE *outputr;
  FILE *outputa;
  FILE *outputxyz;
} transformDiffusionInfo;

#define INITIALIZE_transformDiffusionInfo(_p_) \
  _p_->dx             = NULL; \
  _p_->dy             = NULL; \
  _p_->dz             = NULL; \
  _p_->timePerFrame   = 0.0;  \
  _p_->prevx          = NULL; \
  _p_->prevy          = NULL; \
  _p_->prevz          = NULL; \
  _p_->distancex      = NULL; \
  _p_->distancey      = NULL; \
  _p_->distancez      = NULL; \
  _p_->distance       = NULL; \
  _p_->deltax         = NULL; \
  _p_->deltay         = NULL; \
  _p_->deltaz         = NULL; \
  _p_->activeAtoms    = 0; \
  _p_->elapsedFrames  = 0; \
  _p_->outputFilenameRoot = NULL; \
  _p_->outputx        = NULL; \
  _p_->outputy        = NULL; \
  _p_->outputz        = NULL; \
  _p_->outputr        = NULL; \
  _p_->outputxyz      = NULL

// transformDiffusion()
int transformDiffusion(actionInformation *action, 
		   double *x, double *y, double *z,
		   double *box, int mode)
{
  //char *name = "diffusion";
  argStackType **argumentStackPointer;
  char *buffer;
  ptrajState *state;
  transformDiffusionInfo *diffusionInfo;
  int i;
  int currentAtom;
  double delx, dely, delz;
  double xx, yy, zz;
  double avgx, avgy, avgz;
  double average;
  double time;

  /*
   *  USAGE:
   *
   *     diffusion mask [average] [time <time per frame>]
   *
   *  action argument usage:
   *
   *  mask:
   *    atoms for which the diffusion is calculated
   *  iarg1: 
   *    0 -- default, print out average diffusion and diffusion values for
   *         each of the active atoms
   *    1 -- only print out averages
   *  carg1:
   *    the transformDiffusionInfo structure 
   *
   */
  xx = 0; yy = 0; zz = 0;

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

    diffusionInfo = safe_malloc(sizeof(transformDiffusionInfo));
    INITIALIZE_transformDiffusionInfo(diffusionInfo);

    buffer = getArgumentString(argumentStackPointer, NULL);
    action->mask = processAtomMask(buffer, action->state);
    safe_free(buffer);

    diffusionInfo->timePerFrame = getArgumentDouble(argumentStackPointer, 1.0);
    if (diffusionInfo->timePerFrame < 0) {
      error("ptraj()", "diffusion time per frame incorrectly specified\n");
    }

    action->iarg1 = argumentStackContains(argumentStackPointer, "average");

    diffusionInfo->outputFilenameRoot = getArgumentString(argumentStackPointer, "diffusion");
    action->carg1 = (void *) diffusionInfo;
    return 0;

  }


  diffusionInfo = (transformDiffusionInfo *) action->carg1;


  if (mode == PTRAJ_STATUS) {

    /*
     *  ACTION: PTRAJ_STATUS
     */

    fprintf(stdout, "  DIFFUSION\n");
    if (action->iarg1 == 1) {
      fprintf(stdout, "      Only the average results will ");
    } else {
      fprintf(stdout, "      The average and individual results will ");
    }
    fprintf(stdout, "be dumped to %s_?.xmgr\n", 
	    diffusionInfo->outputFilenameRoot);
    fprintf(stdout, "      The time between frames in psec is %5.3f.\n",
	    diffusionInfo->timePerFrame);
    fprintf(stdout, 
	    "      To calculated diffusion constants, calculate the slope of the lines(s)\n");
    fprintf(stdout, 
	    "      and multiply by 10.0/6.0; this will give units of 1x10**-5 cm**2/s\n");
    if (action->mask) {	    
      fprintf(stdout, "      The atoms in the calculation follow: ");
      printAtomMask(stdout, action->mask, action->state);
      fprintf(stdout, "\n");
    }
    return 0;

  } else if (mode == PTRAJ_CLEANUP) {

    /*
     *  ACTION: PTRAJ_CLEANUP
     */
    safe_free(action->mask); 
    safe_free(diffusionInfo->dx);
    safe_free(diffusionInfo->dy);
    safe_free(diffusionInfo->dz);
    safe_free(diffusionInfo->prevx);
    safe_free(diffusionInfo->prevy);
    safe_free(diffusionInfo->prevz);
    safe_free(diffusionInfo->distancex);
    safe_free(diffusionInfo->distancey);
    safe_free(diffusionInfo->distancez);
    safe_free(diffusionInfo->distance);
    safe_free(diffusionInfo->deltax);
    safe_free(diffusionInfo->deltay);
    safe_free(diffusionInfo->deltaz);
    safe_free(diffusionInfo->outputFilenameRoot);
    safe_fclose(diffusionInfo->outputx);
    safe_fclose(diffusionInfo->outputy);
    safe_fclose(diffusionInfo->outputz);
    safe_fclose(diffusionInfo->outputr);
    safe_fclose(diffusionInfo->outputa);
    INITIALIZE_transformDiffusionInfo(diffusionInfo);
    safe_free(diffusionInfo);

  }


  if (mode != PTRAJ_ACTION) return 0;

  /*
   *  ACTION: PTRAJ_ACTION
   */

  state = (ptrajState *) action->state;

  /*
   *  update local state information
   */
  for (i=0; i<6; i++)
    state->box[i] = box[i];

  diffusionInfo = (transformDiffusionInfo *) action->carg1;
  diffusionInfo->elapsedFrames++;

  /*
   *  load up initial frame if necessary
   */
  if ( diffusionInfo->dx == NULL ) {

    diffusionInfo->dx = safe_malloc(sizeof(double) * state->atoms);
    diffusionInfo->dy = safe_malloc(sizeof(double) * state->atoms);
    diffusionInfo->dz = safe_malloc(sizeof(double) * state->atoms);
  
    for (i=0; i < state->atoms; i++) {
      diffusionInfo->dx[i] = x[i];
      diffusionInfo->dy[i] = y[i];
      diffusionInfo->dz[i] = z[i];
    }
    for (i=0; i < state->atoms; i++) {
      if (action->mask == NULL || action->mask[i])
	diffusionInfo->activeAtoms++;
    }

    diffusionInfo->prevx = 
      safe_malloc(sizeof(double) * diffusionInfo->activeAtoms);
    diffusionInfo->prevy = 
      safe_malloc(sizeof(double) * diffusionInfo->activeAtoms);
    diffusionInfo->prevz = 
      safe_malloc(sizeof(double) * diffusionInfo->activeAtoms);

    currentAtom = 0;
    for (i=0; i < state->atoms; i++) {
      if (action->mask == NULL || action->mask[i]) {
	diffusionInfo->prevx[currentAtom] = x[i];
	diffusionInfo->prevy[currentAtom] = y[i];
	diffusionInfo->prevz[currentAtom] = z[i];
	currentAtom++;
      }
    }

    diffusionInfo->distancex = 
      safe_malloc(sizeof(double) * diffusionInfo->activeAtoms);
    diffusionInfo->distancey = 
      safe_malloc(sizeof(double) * diffusionInfo->activeAtoms);
    diffusionInfo->distancez = 
      safe_malloc(sizeof(double) * diffusionInfo->activeAtoms);
    diffusionInfo->distance = 
      safe_malloc(sizeof(double) * diffusionInfo->activeAtoms);
    
    diffusionInfo->elapsedFrames = 0;

    diffusionInfo->deltax = 
      safe_malloc(sizeof(double) * diffusionInfo->activeAtoms);
    diffusionInfo->deltay = 
      safe_malloc(sizeof(double) * diffusionInfo->activeAtoms);
    diffusionInfo->deltaz = 
      safe_malloc(sizeof(double) * diffusionInfo->activeAtoms);
    for (i=0; i < diffusionInfo->activeAtoms; i++) {
      diffusionInfo->deltax[i] = 0.0;
      diffusionInfo->deltay[i] = 0.0;
      diffusionInfo->deltaz[i] = 0.0;
    }

    buffer = (char *) safe_malloc(sizeof(char) * 
				  (strlen(diffusionInfo->outputFilenameRoot)+10));
    strcpy(buffer, diffusionInfo->outputFilenameRoot);
    strcat(buffer, "_x.xmgr");
    if ( ( diffusionInfo->outputx = fopen(buffer,"w") ) == NULL) {
      fprintf(stdout, "WARNING in ptraj(), diffusion: Cannot open diffusion output file\n");
    }

    strcpy(buffer, diffusionInfo->outputFilenameRoot);
    strcat(buffer, "_y.xmgr");
    if ( ( diffusionInfo->outputy = fopen(buffer,"w") ) == NULL) {
      fprintf(stdout, "WARNING in ptraj(), diffusion: Cannot open diffusion output file\n");
    }

    strcpy(buffer, diffusionInfo->outputFilenameRoot);
    strcat(buffer, "_z.xmgr");
    if ( ( diffusionInfo->outputz = fopen(buffer,"w") ) == NULL) {
      fprintf(stdout, "WARNING in ptraj(), diffusion: Cannot open diffusion output file\n");
    }

    strcpy(buffer, diffusionInfo->outputFilenameRoot);
    strcat(buffer, "_r.xmgr");
    if ( ( diffusionInfo->outputr = fopen(buffer,"w") ) == NULL) {
      fprintf(stdout, "WARNING in ptraj(), diffusion: Cannot open diffusion output file\n");
    }

    strcpy(buffer, diffusionInfo->outputFilenameRoot);
    strcat(buffer, "_a.xmgr");
    if ( ( diffusionInfo->outputa = fopen(buffer,"w") ) == NULL) {
      fprintf(stdout, "WARNING in ptraj(), diffusion: Cannot open diffusion output file\n");
    }
    if (prnlev > 2) {
      if ( ( diffusionInfo->outputxyz = fopen("diffusion_xyz.xmgr","w") ) == NULL) {
	fprintf(stdout, "WARNING in ptraj(), diffusion: Cannot open diffusion output file\n");
      }
    }
    safe_free(buffer);
    return 1;
  }


  currentAtom = 0;
  for (i=0; i < state->atoms; i++) {
    if ( action->mask == NULL || action->mask[i] == 1 ) {

      if ( currentAtom > diffusionInfo->activeAtoms )
	error("junk", "currentAtom out of bounds!\n");

      /*
       *  calculate distance to previous frames coordinates
       */
      delx = x[i] - diffusionInfo->prevx[currentAtom];
      dely = y[i] - diffusionInfo->prevy[currentAtom];
      delz = z[i] - diffusionInfo->prevz[currentAtom];

      /*
       *  if the particle moved more than half the box, assume
       *  it was imaged and adjust the distance of the total
       *  movement with respect to the original frame...
       */
      if ( state->box[0] > 0.0 ) {
	if ( delx > state->box[0]/2.0 )
	  diffusionInfo->deltax[currentAtom] -= state->box[0];
	else if ( delx < -state->box[0]/2.0 )
	  diffusionInfo->deltax[currentAtom] += state->box[0];
	if ( dely > state->box[1]/2.0 )
	  diffusionInfo->deltay[currentAtom] -= state->box[1];
	else if ( dely < -state->box[1]/2.0 )
	  diffusionInfo->deltay[currentAtom] += state->box[1];
	if ( delz > state->box[2]/2.0 ) 
	  diffusionInfo->deltaz[currentAtom] -= state->box[2];
	else if ( delz < -state->box[2]/2.0 )
	  diffusionInfo->deltaz[currentAtom] += state->box[2];
      }

      if (prnlev > 2) {
	fprintf(stdout, "ATOM: %5i %10.3f %10.3f %10.3f",
		i, x[i], delx,
		diffusionInfo->deltax[currentAtom]);
      }

      /*
       *  set the current x with reference to the un-imaged
       *  trajectory
       */
      xx = x[i] + diffusionInfo->deltax[currentAtom];
      yy = y[i] + diffusionInfo->deltay[currentAtom];
      zz = z[i] + diffusionInfo->deltaz[currentAtom];

      /*
       *  calculate the distance between this "fixed" coordinate
       *  and the reference (initial) frame
       */
      delx = xx - diffusionInfo->dx[i];
      dely = yy - diffusionInfo->dy[i];
      delz = zz - diffusionInfo->dz[i];

      if (prnlev > 2) 
	fprintf(stdout, " %10.3f\n", delx);


      /*
       *  store the distance for this atom
       */
      diffusionInfo->distancex[currentAtom] = delx*delx;
      diffusionInfo->distancey[currentAtom] = dely*dely;
      diffusionInfo->distancez[currentAtom] = delz*delz;
      diffusionInfo->distance[currentAtom] = delx*delx + dely*dely + delz*delz;

      /*
       *  update the previous coordinate set to match the
       *  current coordinates
       */
      diffusionInfo->prevx[currentAtom] = x[i];
      diffusionInfo->prevy[currentAtom] = y[i];
      diffusionInfo->prevz[currentAtom] = z[i];

      /*
       *  increment the current atom pointer
       */
      currentAtom++;
    }
  }

  /*
   *  accumulate averages
   */
  average = 0.0;
  avgx = 0.0;
  avgy = 0.0;
  avgz = 0.0;
  for (i=0; i < diffusionInfo->activeAtoms; i++) {
    average += diffusionInfo->distance[i];
    avgx += diffusionInfo->distancex[i];
    avgy += diffusionInfo->distancey[i];
    avgz += diffusionInfo->distancez[i];
  }
    
  average /= (double) diffusionInfo->activeAtoms;
  avgx /= (double) diffusionInfo->activeAtoms;
  avgy /= (double) diffusionInfo->activeAtoms;
  avgz /= (double) diffusionInfo->activeAtoms;


  /*
   *  dump output
   */
  time = diffusionInfo->elapsedFrames * diffusionInfo->timePerFrame;
  fprintf(diffusionInfo->outputx, "%8.3f  %8.3f", 
	  time,
	  avgx );
  fprintf(diffusionInfo->outputy, "%8.3f  %8.3f", 
	  time,
	  avgy );
  fprintf(diffusionInfo->outputz, "%8.3f  %8.3f", 
	  time,
	  avgz );
  fprintf(diffusionInfo->outputr, "%8.3f  %8.3f", 
	  time,
	  average);
  fprintf(diffusionInfo->outputa, "%8.3f  %8.3f", 
	  time,
	  sqrt(average));
  if (prnlev > 2)
    fprintf(diffusionInfo->outputxyz, "%8.3f  %8.3f  %8.3f  %8.3f", 
	    time, xx, yy, zz);

  /*
   *  dump individual values if requested
   */
  if (action->iarg1 == 0) {
    for (i = 0; i < diffusionInfo->activeAtoms; i++) {
      fprintf(diffusionInfo->outputx, "  %8.3f", 
	      diffusionInfo->distancex[i] );
      fprintf(diffusionInfo->outputy, "  %8.3f", 
	      diffusionInfo->distancey[i] );
      fprintf(diffusionInfo->outputz, "  %8.3f", 
	      diffusionInfo->distancez[i] );
      fprintf(diffusionInfo->outputr, "  %8.3f", 
	      diffusionInfo->distance[i] );
      fprintf(diffusionInfo->outputa, "  %8.3f", 
	      sqrt(diffusionInfo->distance[i]) );
    }
  }

  /*
   *  dump newlines 
   */
  fprintf(diffusionInfo->outputx, "\n");
  fprintf(diffusionInfo->outputy, "\n");
  fprintf(diffusionInfo->outputz, "\n");
  fprintf(diffusionInfo->outputr, "\n");
  fprintf(diffusionInfo->outputa, "\n");
  if (prnlev > 2)
    fprintf(diffusionInfo->outputxyz, "\n");

  fflush(diffusionInfo->outputx);
  fflush(diffusionInfo->outputy);
  fflush(diffusionInfo->outputz);
  fflush(diffusionInfo->outputr);
  fflush(diffusionInfo->outputa);
  return 1;
}

/** ACTION ROUTINE *************************************************************
 *
 *  transformDipole()   --- bin dipoles (Jed Pitera, UCSF)
 *
 ******************************************************************************/


   int
transformDipole(actionInformation *action, 
		double *x, double *y, double *z,
		double *box, int mode)
{
  //char *name = "dipole";
  argStackType **argumentStackPointer;
  char *buffer;
  ptrajState *state;
  transformGridInfo *dipoleInfo;
  int i, j, k;
  int i_solvent; 
  int isbox;
  double *xsol, *ysol, *zsol;
  double xcm, ycm, zcm, mass;
  int *mask;
  double dipolar_vector[3]; 
  int nx, ny, nz, index;
  double sx, sy, sz, max_density, density;
  FILE *outFile;

  /*
   *  USAGE:
   *
   *    dipole filename nx dx ny dx nz dz mask [box|origin] [max <%>]
   *
   *  action argument usage:
   *
   *  mask
   *     the waters for which the dipole is calculated
   *  carg1:
   *     a transformGridInfo structure
   *  iarg1:
   *     1 -- grid is centered at the origin
   *     0 -- grid is centered at the box center (if there is box information)
   *  iarg2:
   *    >0 -- only dump density >= iarg2% of the maximum density
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

    dipoleInfo = (transformGridInfo *)
      safe_malloc(sizeof(transformGridInfo));
    INITIALIZE_transformGridInfo(dipoleInfo);
	    
    dipoleInfo->filename = getArgumentString(argumentStackPointer, NULL);
    if (dipoleInfo->filename == NULL) {
      error("ptraj()", "No file was specified for grid/dipole output\n");
    }
    dipoleInfo->nx = getArgumentInteger(argumentStackPointer, -1);
    dipoleInfo->dx = getArgumentDouble( argumentStackPointer, -1);
    dipoleInfo->ny = getArgumentInteger(argumentStackPointer, -1);
    dipoleInfo->dy = getArgumentDouble( argumentStackPointer, -1);
    dipoleInfo->nz = getArgumentInteger(argumentStackPointer, -1);
    dipoleInfo->dz = getArgumentDouble( argumentStackPointer, -1);
    if (dipoleInfo->nx < 0 || dipoleInfo->ny < 0 || dipoleInfo->nz < 0 ||
	dipoleInfo->dx < 0 || dipoleInfo->dy < 0 || dipoleInfo->dz < 0) {
      error("ptraj()", "Specification of grid/dipole size or spacing\n");
    }

    if (dipoleInfo->nx % 2 == 1) {
      fprintf(stdout, "WARNING in ptraj(): grid/dipole -- number of grid points must be even!\n");
      dipoleInfo->nx++;
      fprintf(stdout, "Incrementing NX by 1 to %i\n", dipoleInfo->nx);
    }

    if (dipoleInfo->ny % 2 == 1) {
      fprintf(stdout, "WARNING in ptraj(): grid/dipole -- number of grid points must be even!\n");
      dipoleInfo->ny++;
      fprintf(stdout, "Incrementing NY by 1 to %i\n", dipoleInfo->ny);
    }

    if (dipoleInfo->nz % 2 == 1) {
      fprintf(stdout, "WARNING in ptraj(): grid/dipole -- number of grid points must be even!\n");
      dipoleInfo->nz++;
      fprintf(stdout, "Incrementing NZ by 1 to %i\n", dipoleInfo->nz);
    }

    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL) {
      fprintf(stdout, "WARNING in ptraj() dipole: No mask was specified for grid/dipole\n");
      return 0;
    }
    action->mask = processAtomMask(buffer, action->state);
    safe_free(buffer);

    if (argumentStackContains(argumentStackPointer, "box"))
      action->iarg1 = 1; 
    else if (argumentStackContains(argumentStackPointer, "origin"))
      action->iarg1 = 0;
    if (action->iarg1 == 1 &&
	(action->state->box[3] != 90.0 || box[3] != 90.0 ||
	 action->state->box[4] != 90.0 || box[4] != 90.0 ||
	 action->state->box[5] != 90.0 || box[5] != 90.0)) {
      fprintf(stdout, "WARNING in ptraj() dipole: Code to shift to the box center\n");
      fprintf(stdout, "is not implemented yet in transformDipole for non-orthorhomibic\n");
      fprintf(stdout, "unit cells.  Shifting to the origin instead!!!\n");
      action->iarg1 = 0;
    }

    if (argumentStackContains(argumentStackPointer, "negative"))
      action->iarg2 = 1; 
    else {
      action->iarg2 = argumentStackKeyToInteger(argumentStackPointer, "max", 0);
    }

    dipoleInfo->frames = 0;
    dipoleInfo->grid = (float *)
      safe_malloc(sizeof(float) * dipoleInfo->nx*dipoleInfo->ny*dipoleInfo->nz);
    memset(dipoleInfo->grid, 0,
	  sizeof(float) * dipoleInfo->nx*dipoleInfo->ny*dipoleInfo->nz);

    dipoleInfo->dipolex = (float *)
      safe_malloc(sizeof(float) * dipoleInfo->nx*dipoleInfo->ny*dipoleInfo->nz);
    memset(dipoleInfo->dipolex, 0,
	  sizeof(float) * dipoleInfo->nx*dipoleInfo->ny*dipoleInfo->nz);
    dipoleInfo->dipoley = (float *)
      safe_malloc(sizeof(float) * dipoleInfo->nx*dipoleInfo->ny*dipoleInfo->nz);
    memset(dipoleInfo->dipoley, 0,
	  sizeof(float) * dipoleInfo->nx*dipoleInfo->ny*dipoleInfo->nz);
    dipoleInfo->dipolez = (float *)
      safe_malloc(sizeof(float) * dipoleInfo->nx*dipoleInfo->ny*dipoleInfo->nz);
    memset(dipoleInfo->dipolez, 0,
	  sizeof(float) * dipoleInfo->nx*dipoleInfo->ny*dipoleInfo->nz);

    action->carg1 = (void *) dipoleInfo;
    return 0;
  }


  dipoleInfo = (transformGridInfo *) action->carg1;


  if (mode == PTRAJ_STATUS) {

    /*
     *  ACTION: PTRAJ_STATUS
     */

    fprintf(stdout, "  DIPOLE: grid at %s dumped to filename %s\n",
            (action->iarg1 == 0 ? "origin" : "box center"),
            dipoleInfo->filename);
    if (action->mask) {
      fprintf(stdout, "      Water selection is ");
      printAtomMask(stdout, action->mask, action->state);
      fprintf(stdout, "\n");
    }

    if (action->iarg2) {
      fprintf(stdout, "      Only keeping density >= to %i%% of the maximum density\n",
	      action->iarg2);
    }
    fprintf(stdout, "      dipole points : %5i %5i %5i\n",
            dipoleInfo->nx, dipoleInfo->ny, dipoleInfo->nz);
    fprintf(stdout, "      dipole spacing: %5.3f %5.3f %5.3f\n",
            dipoleInfo->dx, dipoleInfo->dy, dipoleInfo->dz);

  } else if (mode == PTRAJ_PRINT) {

    /*
     *  ACTION: PTRAJ_PRINT
     */

    fprintf(stdout, "PTRAJ DIPOLE: dumping dipole density\n");
    if ( dipoleInfo->filename == NULL ) {
      fprintf(stdout, "WARNING in ptraj(): dipole filename is NULL, not dumping dipole grid!\n");
    }

    if ( ( outFile = fopen(dipoleInfo->filename, "w") ) == NULL ) {
      fprintf(stdout, "WARNING in ptraj(), dipole: couldn't open file %s\n",
	      dipoleInfo->filename);
      return -1;
    }

    /*
     *      write out header, data here
     */
 
    fprintf(stdout, "PTRAJ: dipole, dumping data to output file %s\n", 
	    dipoleInfo->filename);
    fprintf(outFile, "field 8\n");
    fprintf(outFile, "size 1\n");
    fprintf(outFile, "nside 3\n");
    fprintf(outFile, "nlayer 1\n");
    fprintf(outFile, "directional\n");
    fprintf(outFile, "vector\n");
    fprintf(outFile, "data\n");

    /*
     *	sx, etc are center of grid in real coords
     */
    sx = (double) dipoleInfo->nx * dipoleInfo->dx/2.0;
    sy = (double) dipoleInfo->ny * dipoleInfo->dy/2.0;
    sz = (double) dipoleInfo->nz * dipoleInfo->dz/2.0;

    /*
     *  determine the maximum density
     */
    max_density = 0.0;
    for (k = 0; k < dipoleInfo->nz; k++) {
      for (j = 0; j < dipoleInfo->ny; j++) {
	for (i = 0; i < dipoleInfo->nx; i++) {
	  index = k * dipoleInfo->ny * dipoleInfo->nz +
	    j * dipoleInfo->nz + i;
	  density = (int) dipoleInfo->grid[index];
	  if ( max_density < density )
	    max_density = density;
	}
      }
    }
    fprintf(stdout, "        maximum density is %f\n", max_density);

    if ( action->iarg2 ) {
      max_density = action->iarg2 * max_density / 100.0;
      fprintf(stdout, "dumping density if >= to %f\n", max_density);
    } else
      max_density = 1.0;

    /*
     *  dump the data
     */
    for (k = 0; k < dipoleInfo->nz; k++) {
      for (j = 0; j < dipoleInfo->ny; j++) {
	for (i = 0; i < dipoleInfo->nx; i++) {
	  index = k * dipoleInfo->ny * dipoleInfo->nz +
	    j * dipoleInfo->nz + i;

	  density = (int) dipoleInfo->grid[index];

	  if ( density >= max_density ) {

	    /*
	     *	re-center coords
	     */
	    fprintf(outFile, "%8.3f %8.3f %8.3f", 
		    i*dipoleInfo->dx - sx,
		    j*dipoleInfo->dy - sy,
		    k*dipoleInfo->dz - sz);
	    /*
	     *	normalize dipoles by density
	     */
	    dipoleInfo->dipolex[index] =
	      (float) (dipoleInfo->dipolex[index]/density);
	    dipoleInfo->dipoley[index] =
	      (float) (dipoleInfo->dipoley[index]/density);
	    dipoleInfo->dipolez[index] =
	      (float) (dipoleInfo->dipolez[index]/density);
	    /*
	     *      writeout dipole components, length
	     */
	    fprintf(outFile, " %8.3f %8.3f %8.3f", 
		    dipoleInfo->dipolex[index],
		    dipoleInfo->dipoley[index],
		    dipoleInfo->dipolez[index]);
	    fprintf(outFile, " %8.3f %8.3f\n",
		    sqrt((dipoleInfo->dipolex[index])*
			 (dipoleInfo->dipolex[index]) +
			 (dipoleInfo->dipoley[index])*
			 (dipoleInfo->dipoley[index]) +
			 (dipoleInfo->dipolez[index])*
			 (dipoleInfo->dipolez[index])),
		    density); 
	  } 
	}
      }
    }
    safe_fclose(outFile);

  } else if (mode == PTRAJ_CLEANUP) {

    /*
     *  ACTION: PTRAJ_CLEANUP
     */

    safe_free(dipoleInfo->filename);
    safe_free(dipoleInfo->grid);
    safe_free(dipoleInfo->dipolex);
    safe_free(dipoleInfo->dipoley);
    safe_free(dipoleInfo->dipolez);
    INITIALIZE_transformGridInfo(dipoleInfo);
    safe_free(dipoleInfo);

  }


  if (mode != PTRAJ_ACTION) return 0;

  /*
   *  ACTION: PTRAJ_ACTION
   */

  state = (ptrajState *) action->state;
  /*
   *  update local state information
   */
  for (i=0; i<6; i++)
    state->box[i] = box[i];

  nx = dipoleInfo->nx;
  ny = dipoleInfo->ny;
  nz = dipoleInfo->nz;

  sx = (double) nx * dipoleInfo->dx/2.0;
  sy = (double) ny * dipoleInfo->dy/2.0;
  sz = (double) nz * dipoleInfo->dz/2.0;

  isbox = action->iarg1;

  if (state->IFBOX == 0)
    isbox = 0;


  mask = action->mask;

  /*
   *  traverse over solvent molecules to find out the 
   *  "largest" solvent molecule; allocate space for this
   *  many coordinates
   */
  j=0;
  for (i_solvent=0; i_solvent < state->solventMolecules; i_solvent++) {
    i = state->solventMoleculeStop[i_solvent] - state->solventMoleculeStart[i_solvent];
    if (i_solvent == 0 || j < i) j = i;
  }

  xsol = (double *) safe_malloc(sizeof(double) * j);
  ysol = (double *) safe_malloc(sizeof(double) * j);
  zsol = (double *) safe_malloc(sizeof(double) * j);

  /*
   *      traverse over solvent molecules
   */
  for (i_solvent=1; i_solvent < state->solventMolecules; i_solvent++) {

    /*
     *  get coordinates and shift to origin and then to appropriate spacing
     */

    i = 0;
    for (j = state->solventMoleculeStart[i_solvent];
	 j < state->solventMoleculeStop[i_solvent];
	 j++) {
      if (mask == NULL || mask[j]) {
	xsol[i] = x[j] + sx - (isbox ? state->box[0]/2.0 : 0.0);
	ysol[i] = y[j] + sy - (isbox ? state->box[1]/2.0 : 0.0);
	zsol[i] = z[j] + sz - (isbox ? state->box[2]/2.0 : 0.0);
	i++;
      }
    }
    
    /*
     *  calculate dipolar vector.  NOTE: the total charge on the solvent
     *  should be neutral for this to have any meaning...  The center
     *  of mass of the solvent is used as the "origin" for the vector
     */
    dipolar_vector[0] = 0.0;
    dipolar_vector[1] = 0.0;
    dipolar_vector[2] = 0.0;
    xcm = 0.0;
    ycm = 0.0;
    zcm = 0.0;
    mass = 0.0;
    i = 0;
    for (j = state->solventMoleculeStart[i_solvent];
	 j < state->solventMoleculeStop[i_solvent];
	 j++) {
    
      dipolar_vector[0] += state->charges[j] * xsol[i];
      dipolar_vector[1] += state->charges[j] * ysol[i];
      dipolar_vector[2] += state->charges[j] * zsol[i];
      xcm += state->masses[j] * xsol[i];
      ycm += state->masses[j] * ysol[i];
      zcm += state->masses[j] * zsol[i];
      mass += state->masses[j];
      i++;
    }
    xcm /= mass;
    ycm /= mass;
    zcm /= mass;
    
    /*
     *  determine indices into arrays
     */
    i = (int) (xcm / dipoleInfo->dx);
    j = (int) (ycm / dipoleInfo->dy);
    k = (int) (zcm / dipoleInfo->dz);

    /*
     *  check bounds and increment grid, dipole if appropriate
     */

    if (prnlev > 2)
      printf("Dipole: index is %i %i %i, water %i\n",
	     i, j, k, i_solvent);

    if (i > 0 && i < nx && j > 0 && j < ny && k > 0 && k < nz) {
      dipoleInfo->grid[ k*ny*nz + j*nz + i ] += 1.0;
      dipoleInfo->dipolex[ k*ny*nz + j*nz + i] += dipolar_vector[0];
      dipoleInfo->dipoley[ k*ny*nz + j*nz + i] += dipolar_vector[1];
      dipoleInfo->dipolez[ k*ny*nz + j*nz + i] += dipolar_vector[2];
    }
  }

  if (prnlev > 2)
    printf("Dipole: finished traversing waters\n");

  return 1;
}

/** ACTION ROUTINE *************************************************************
 *
 *  transformDNAiontracker
 *
 ******************************************************************************/
// setupDistance()
static void setupDistance(double *X, double *Y, double *Z, double *x, double *y, double *z,
                          ptrajState *state, int *mask1, int *mask2, int atom1, int atom2)
{
  double total_mass1;
  double total_mass2;
  double atommass;
  int i;

  X[0] = 0.0;
  Y[0] = 0.0;
  Z[0] = 0.0;
  total_mass1 = 0.0;
  X[1] = 0.0;
  Y[1] = 0.0;
  Z[1] = 0.0;
  total_mass2 = 0.0;
  atommass=1.0;

  if (atom1 == -1) {

    for (i=0; i < state->atoms; i++) {
      if (mask1[i]) {
        if (state->masses[i]!=0.0) atommass=state->masses[i];
        X[0] += atommass * x[i];
        Y[0] += atommass * y[i];
        Z[0] += atommass * z[i];
        total_mass1 += atommass;
      }
    }
    X[0] /= total_mass1;
    Y[0] /= total_mass1;
    Z[0] /= total_mass1;

  } else {

    X[0] = x[atom1];
    Y[0] = y[atom1];
    Z[0] = z[atom1];

  }

  if (atom2 == -1) {

    for (i=0; i < state->atoms; i++) {
      if (mask2[i]) {
        if (state->masses[i]!=0.0) atommass=state->masses[i];
        X[1] += atommass * x[i];
        Y[1] += atommass * y[i];
        Z[1] += atommass * z[i];
        total_mass2 += atommass;
      }
    }
    X[1] /= total_mass2;
    Y[1] /= total_mass2;
    Z[1] /= total_mass2;

  } else {

    X[1] = x[atom2];
    Y[1] = y[atom2];
    Z[1] = z[atom2];
  }
}

// transformDNAiontracker()
int transformDNAiontracker(actionInformation *action, 
		       double *x, double *y, double *z, 
		       double *box, int mode)
{
  //char *name = "dnaiontracker";
  argStackType **argumentStackPointer;
  char *buffer;
  scalarInfo *distanceInfo;
  ptrajState *state;
  int i, mask1tot, mask2tot, mask3tot, mask4tot;
  double X[2], Y[2], Z[2], ucell[9], recip[9];
  double pp_centroidx, pp_centroidy, pp_centroidz;
  double d_p1ion, d_p2ion, d_baseion, d_cut;
  double d_min, d_tmp;
  double d_pp, poffset, d_pbase;
  int bound, boundLower, boundUpper;

  FILE *outFile;

  /*
   *  USAGE:
   *
   *    dnaiontracker name mask_p1 mask_p2 mask_base mask_ions \
   *      [poffset <value>] [out <filename>] [time <interval>] [noimage] [shortest | count]
   *
   *  action argument usage:
   *
   *  iarg1: 1 implies don't image
   *  iarg3: flag to determine if distance (shortest) or count is saved
   *  darg1: time interval in ps (for output)
   *  darg2: poffset (perpendicular offset)
   *  carg1:
   *     a scalarInfo structure
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

       /*
        *  set up the information necessary to place this on the scalarStack
        */
    distanceInfo = (scalarInfo *) safe_malloc(sizeof(scalarInfo));
    INITIALIZE_scalarInfo(distanceInfo);
    distanceInfo->mode = SCALAR_DISTANCE;
    distanceInfo->totalFrames = -1;

    distanceInfo->name = getArgumentString(argumentStackPointer, NULL);
    if (distanceInfo->name == NULL) {
      fprintf(stdout, "WARNING: ptraj(), dnaiontracker: It is necessary to specify a unique name\n");
      fprintf(stdout, "for each specified tracking.  Ignoring command...\n");
      safe_free(distanceInfo);
      return -1;
    } /*else if ( scalarStackGetName(&scalarStack, distanceInfo->name) != NULL ) {
      fprintf(stdout, "WARNING: ptraj(), dnaiontracker: The chosen name (%s) has already been used.\n",
	      distanceInfo->name);
      fprintf(stdout, "Ignoring command...\n");
      safe_free(distanceInfo);
      return -1;
    }*/
    distanceInfo->state = action->state;

       /*
        *  grab the filename
        */
    distanceInfo->filename = argumentStackKeyToString(argumentStackPointer, "out", NULL);

       /*
        *  grab the perpendicular offset (poffset)
        */
    action->darg2 = argumentStackKeyToDouble(argumentStackPointer, "poffset", 5.0);


       /*
        *  decide whether to bin the shortest distances seen (to phosphates or base
	*  centroid) or whether to simply bin count or counttopcone or countbottomcone
	*/
    action->iarg3 = 0;
    if (argumentStackContains(argumentStackPointer, "shortest"))
      action->iarg3 = 1;
    else if (argumentStackContains(argumentStackPointer, "counttopcone"))
      action->iarg3 = 2;
    else if (argumentStackContains(argumentStackPointer, "countbottomcone"))
      action->iarg3 = 3;
    else if (argumentStackContains(argumentStackPointer, "count"))
      action->iarg3 = 0;


       /*
        *  push the distance info on to the distance stack
        */
      fprintf(stdout,"Warning: scalarStack disabled for Cpptraj\n");
//    pushBottomStack(&scalarStack, (void *) distanceInfo);

       /*
        *  grab a time interval between frames in ps (for output)
        */
    action->darg1 = argumentStackKeyToDouble(argumentStackPointer, "time", 1.0);

       /*
        *  check to see if we want imaging disabled
        */
    action->iarg1 = argumentStackContains(argumentStackPointer, "noimage");

       /*
        *  process the atom masks
        */
    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL) {
      fprintf(stdout, 
	      "WARNING in ptraj(), dnaiontracker: Error in specification of the first phosphate mask\n");
      fprintf(stdout, "Ignoring command\n");
      safe_free(distanceInfo);
      return -1;

    } else {
      distanceInfo->mask1 = processAtomMask(buffer, action->state);
      safe_free(buffer);
    }

    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL) {
      fprintf(stdout, 
	      "WARNING in ptraj(), dnaiontracker: Error in specification of the second phosphate mask\n");
      fprintf(stdout, "Ignoring command\n");
      safe_free(distanceInfo);
      return -1;

    } else {
      distanceInfo->mask2 = processAtomMask(buffer, action->state);
      safe_free(buffer);
    }

    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL) {
      fprintf(stdout, 
	      "WARNING in ptraj(), dnaiontracker: Error in specification of the base centroid mask\n");
      fprintf(stdout, "Ignoring command\n");
      safe_free(distanceInfo);
      return -1;

    } else {
      distanceInfo->mask3 = processAtomMask(buffer, action->state);
      safe_free(buffer);
    }

    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL) {
      fprintf(stdout, 
	      "WARNING in ptraj(), dnaiontracker: Error in specification of the ion mask\n");
      fprintf(stdout, "Ignoring command\n");
      safe_free(distanceInfo);
      return -1;

    } else {
      distanceInfo->mask4 = processAtomMask(buffer, action->state);
      safe_free(buffer);
    }


    /*
     *  check to see if each mask only represents a single atom or not (to save on
     *  memory)
     */
    mask1tot = 0; distanceInfo->atom1 = -1;
    mask2tot = 0; distanceInfo->atom2 = -1;
    mask3tot = 0; distanceInfo->atom4 = -1;
    mask4tot = 0; distanceInfo->atom3 = -1;
    for (i=0; i < action->state->atoms; i++) {
      if (distanceInfo->mask1[i] == 1) {
	mask1tot++;
	distanceInfo->atom1 = i;
      }
      if (distanceInfo->mask2[i] == 1) {
	mask2tot++;
	distanceInfo->atom2 = i;
      }
      if (distanceInfo->mask3[i] == 1) {
	mask3tot++;
	distanceInfo->atom3 = i;
      }
      if (distanceInfo->mask4[i] == 1) {
	mask4tot++;
	distanceInfo->atom4 = i;
      }
    }

    if (mask1tot == 0) {
      fprintf(stdout, 
	      "WARNING in ptraj(), distance: No atoms selected in mask1, ignoring command\n");
      safe_free(distanceInfo->mask1);
      safe_free(distanceInfo);
      return -1;
    } else if (mask1tot == 1) {
      safe_free(distanceInfo->mask1);
      distanceInfo->mask1 = NULL;
    } else
      distanceInfo->atom1 = -1;

    if (mask2tot == 0) {
      fprintf(stdout, 
	      "WARNING in ptraj(), distance: No atoms selected in mask2, ignoring command\n");
      safe_free(distanceInfo->mask2);
      safe_free(distanceInfo);
      return -1;
    } else if (mask2tot == 1) {
      safe_free(distanceInfo->mask2);
      distanceInfo->mask2 = NULL;
    } else
      distanceInfo->atom2 = -1;

    if (mask3tot == 0) {
      fprintf(stdout, 
	      "WARNING in ptraj(), distance: No atoms selected in mask3, ignoring command\n");
      safe_free(distanceInfo->mask3);
      safe_free(distanceInfo);
      return -1;
    } else if (mask3tot == 1) {
      safe_free(distanceInfo->mask3);
      distanceInfo->mask3 = NULL;
    } else
      distanceInfo->atom3 = -1;

    if (mask4tot == 0) {
      fprintf(stdout, 
	      "WARNING in ptraj(), distance: No atoms selected in mask4, ignoring command\n");
      safe_free(distanceInfo->mask4);
      safe_free(distanceInfo);
      return -1;
    } else if (mask4tot == 1) {
      safe_free(distanceInfo->mask4);
      distanceInfo->mask4 = NULL;
    } else
      distanceInfo->atom4 = -1;

    action->carg1 = (void *) distanceInfo;

    return 0;
  }


  distanceInfo = (scalarInfo *) action->carg1;


  if (mode == PTRAJ_STATUS) {

    /*
     *  ACTION: PTRAJ_STATUS
     */

    fprintf(stdout, "  DNAIONTRACKER: Data representing the ");
    if (action->iarg3 == 0)
      fprintf(stdout, "count within the cone will be\n");
    else if (action->iarg3 == 1)
      fprintf(stdout, "shortest distance to a phosphate or base centroid will be\n");
    else if (action->iarg3 == 2)
      fprintf(stdout, "count in the top half of the cone (and sort-of bound) will be\n");
    else if (action->iarg3 == 3)
      fprintf(stdout, "count in the bottom half of the cone will be\n");
    fprintf(stdout, "      saved to array named %s\n", distanceInfo->name);
    fprintf(stdout, "      Perpendicular offset for cone is %5.2f angstroms\n", action->darg2);
    if (action->iarg1)
      fprintf(stdout, "      Imaging has been disabled\n");
    if (distanceInfo->atom1 == -1) {
      fprintf(stdout, "      Atom selection 1 is ");
      printAtomMask(stdout, distanceInfo->mask1, action->state);
      fprintf(stdout, "\n");
    } else {
      fprintf(stdout, "      Atom selection 1 is :%i@%s\n",
	      atomToResidue(distanceInfo->atom1+1, action->state->residues, action->state->ipres),
	      action->state->atomName[distanceInfo->atom1]);
    }
    if (distanceInfo->atom2 == -1) {
      fprintf(stdout, "      Atom selection 2 is ");
      printAtomMask(stdout, distanceInfo->mask2, action->state);
      fprintf(stdout, "\n");
    } else {
      fprintf(stdout, "      Atom selection 2 is :%i@%s\n",
	      atomToResidue(distanceInfo->atom2+1, action->state->residues, action->state->ipres),
	      action->state->atomName[distanceInfo->atom2]);
    }
    if (distanceInfo->atom3 == -1) {
      fprintf(stdout, "      Atom selection 3 is ");
      printAtomMask(stdout, distanceInfo->mask3, action->state);
      fprintf(stdout, "\n");
    } else {
      fprintf(stdout, "      Atom selection 3 is :%i@%s\n",
	      atomToResidue(distanceInfo->atom3+1, action->state->residues, action->state->ipres),
	      action->state->atomName[distanceInfo->atom3]);
    }
    if (distanceInfo->atom4 == -1) {
      fprintf(stdout, "      Atom selection 4 is ");
      printAtomMask(stdout, distanceInfo->mask4, action->state);
      fprintf(stdout, "\n");
    } else {
      fprintf(stdout, "      Atom selection 4 is :%i@%s\n",
	      atomToResidue(distanceInfo->atom4+1, action->state->residues, action->state->ipres),
	      action->state->atomName[distanceInfo->atom4]);
    }


    if (distanceInfo->filename != NULL) {
      fprintf(stdout, "      Data will be dumped to a file named %s\n",
	      distanceInfo->filename);
    }

  } else if (mode == PTRAJ_PRINT) {

    /*
     *  ACTION: PTRAJ_PRINT
     */

    if ( ( outFile = fopen(distanceInfo->filename, "w") ) == NULL ) {
      fprintf(stdout, "WARNING in ptraj(), dnaiontracker: couldn't open file %s\n",
	      distanceInfo->filename);
      return 0;
    }
    if (prnlev > 2)
      fprintf(stdout, "PTRAJ DNAIONTRACKER dumping distance %s\n",
	      distanceInfo->name);
    for (i=0; i < action->state->maxFrames; i++) {
      fprintf(outFile, "%10.2f %f\n", (i+1) * action->darg1, distanceInfo->value[i]);
    }
    safe_fclose(outFile);

  } else if (mode == PTRAJ_CLEANUP) {

    /*
     *  ACTION: PTRAJ_CLEANUP
     */

    safe_free(distanceInfo->name);
    safe_free(distanceInfo->filename);
    safe_free(distanceInfo->mask1);
    safe_free(distanceInfo->mask2);
    safe_free(distanceInfo->mask3);
    safe_free(distanceInfo->mask4);
    safe_free(distanceInfo->value);
    INITIALIZE_scalarInfo(distanceInfo);
    safe_free(distanceInfo);

  }



  if (mode != PTRAJ_ACTION) return 0;


  /*
   *  ACTION: PTRAJ_ACTION
   */


  state = (ptrajState *) action->state;

  /*
   *  update local state information
   */
  for (i=0; i<6; i++)
    state->box[i] = box[i];

  if (distanceInfo->totalFrames < 0) {
    distanceInfo->totalFrames = state->maxFrames;
    distanceInfo->value = (double *) 
      safe_malloc(sizeof(double) * distanceInfo->totalFrames);
  }

  if (distanceInfo->frame > distanceInfo->totalFrames) {
    warning("transformDNAiontracker()", "Blowing array; too many frames!!\n");
    return 0;
  }


  /*
   *  setup for imaging if necessary
   */

  if (box[3] <= 0.0 && action->iarg1 == 0) {
    action->iarg1 = 1;
    fprintf(stdout, "  DNAIONTRACKER: box angles are zero, disabling imaging!\n");
  }
  if (action->iarg1 == 0 && (box[3] != 90.0 || box[4] != 90.0 || box[5] != 90.0))
    boxToRecip(box, ucell, recip);

  /*
   *  calculate distances
   */

     /*
      *  P -- P distance (as specified in masks1 and masks2)
      */
  setupDistance(X, Y, Z, x, y, z, state, distanceInfo->mask1, distanceInfo->mask2,
		distanceInfo->atom1, distanceInfo->atom2);

  pp_centroidx = (X[0] + X[1]) / 2.0;
  pp_centroidy = (Y[0] + Y[1]) / 2.0;
  pp_centroidz = (Z[0] + Z[1]) / 2.0;

  d_pp = calculateDistance2(0, 1, X, Y, Z, 
			    box, ucell, recip, 0.0, action->iarg1);
  d_pp = sqrt(d_pp);


     /*
      *  perpendicular offset
      */
  poffset = action->darg2;

     /*
      *  P -- base centroid to median point
      */

  setupDistance(X, Y, Z, x, y, z, state, NULL, distanceInfo->mask3,
		1, distanceInfo->atom3);

  X[0] = pp_centroidx;
  Y[0] = pp_centroidy;
  Z[0] = pp_centroidz;

  d_pbase = calculateDistance2(0, 1, X, Y, Z, 
			       box, ucell, recip, 0.0, action->iarg1);
  d_pbase = sqrt(d_pbase);

  /*
   *  loop over ion positions
   */
  d_min = 9999999999.0;
  if (action->iarg3 == 1)
    distanceInfo->value[distanceInfo->frame] = d_min;

  for (i=0; i < state->atoms; i++) {

    if (distanceInfo->mask4[i] == 1) {

      setupDistance(X, Y, Z, x, y, z, state, NULL, distanceInfo->mask1,
		    i, distanceInfo->atom1);
      d_p1ion = calculateDistance2(0, 1, X, Y, Z, 
				   box, ucell, recip, 0.0, action->iarg1);
      d_p1ion = sqrt(d_p1ion);



      setupDistance(X, Y, Z, x, y, z, state, NULL, distanceInfo->mask2,
		    i, distanceInfo->atom2);
      d_p2ion = calculateDistance2(0, 1, X, Y, Z, 
				   box, ucell, recip, 0.0, action->iarg1);
      d_p2ion = sqrt(d_p2ion);

      setupDistance(X, Y, Z, x, y, z, state, NULL, distanceInfo->mask3,
		    i, distanceInfo->atom3);
      d_baseion = calculateDistance2(0, 1, X, Y, Z, 
				     box, ucell, recip, 0.0, action->iarg1);
      d_baseion = sqrt(d_baseion);

      printf("DEBUG: ion atom %i to P1 is %f\n", i+1, d_p1ion);
      printf("DEBUG: ion atom %i to P2 is %f\n", i+1, d_p2ion);
      printf("DEBUG: ion atom %i to base is %f\n", i+1, d_baseion);

      d_cut = sqrt( (d_pp*d_pp*0.25 + poffset*poffset) );

      printf("DEBUG: d_pp is %f, poffset is %f, d_cut is %f\n", d_pp, poffset, d_cut);

      bound = 0;
      boundLower = 0;
      boundUpper = 0;

      if (d_p1ion < d_cut && d_p2ion < d_cut)
	bound = 1;
      if (d_baseion < d_pbase)
	boundLower = 1;

      if (d_p1ion < d_p2ion)
	d_tmp = d_p1ion;
      else
	d_tmp = d_p2ion;
      if (d_tmp > d_baseion)
	d_tmp = d_baseion;

      if (d_tmp > d_min)
	d_min = d_tmp;

      if (bound && boundLower == 0)
	boundUpper = 1;


      if (action->iarg3 == 0)
	distanceInfo->value[distanceInfo->frame] += bound;
      else if (action->iarg3 == 2)
	distanceInfo->value[distanceInfo->frame] += boundUpper;
      else if (action->iarg3 == 3)
	distanceInfo->value[distanceInfo->frame] += boundLower;
      else if (action->iarg3 == 1)
	if (distanceInfo->value[distanceInfo->frame] > d_min)
	  distanceInfo->value[distanceInfo->frame] = d_min;

    }

  }


  distanceInfo->frame++;

  return 1;
}


/** ACTION ROUTINE *************************************************************
 *
 *  transformEcho
 *
 ******************************************************************************/

   int
transformEcho(actionInformation *action, 
	      double *x, double *y, double *z, 
	      double *box, int mode)
{
  //char *name = "echo";
  argStackType **argumentStackPointer;

  /*
   *  USAGE:
   *
   *    echo "string" ["string"] ...
   *
   */


  if (mode == PTRAJ_SETUP) {
    /*
     *  ACTION: PTRAJ_SETUP
     */

    argumentStackPointer = (argStackType **) action->carg1;
    action->carg1 = NULL;

       /*
        *  set up the information necessary to place this on the scalarStack
        */

    action->carg1 = getArgumentString(argumentStackPointer, NULL);
    action->carg2 = getArgumentString(argumentStackPointer, NULL);

    fprintf(stdout, "ECHO (PTRAJ_SETUP): %s %s\n", 
	    (action->carg1 != NULL ? (char *) action->carg1 : ""), 
	    (action->carg2 != NULL ? (char *) action->carg2 : ""));

    return 0;
  }


  if (mode == PTRAJ_STATUS) {

    /*
     *  ACTION: PTRAJ_STATUS
     */


    fprintf(stdout, "ECHO (PTRAJ_STATUS): %s %s\n", 
	    (action->carg1 != NULL ? (char *) action->carg1 : ""), 
	    (action->carg2 != NULL ? (char *) action->carg2 : ""));



  } else if (mode == PTRAJ_PRINT) {

    /*
     *  ACTION: PTRAJ_PRINT
     */

    fprintf(stdout, "ECHO (PTRAJ_PRINT): %s %s\n", 
	    (action->carg1 != NULL ? (char *) action->carg1 : ""), 
	    (action->carg2 != NULL ? (char *) action->carg2 : ""));




  } else if (mode == PTRAJ_CLEANUP) {

    /*
     *  ACTION: PTRAJ_CLEANUP
     */

    fprintf(stdout, "ECHO (PTRAJ_CLEANUP): %s %s\n", 
	    (action->carg1 != NULL ? (char *) action->carg1 : ""), 
	    (action->carg2 != NULL ? (char *) action->carg2 : ""));

  }


  if (mode != PTRAJ_ACTION) return 0;


  /*
   *  ACTION: PTRAJ_ACTION
   */

  fprintf(stdout, "ECHO (PTRAJ_ACTION): %s %s\n", 
	  (action->carg1 != NULL ? (char *) action->carg1 : ""), 
	  (action->carg2 != NULL ? (char *) action->carg2 : ""));

  return 1;
}


/** ACTION ROUTINE *************************************************************
 *
 *  transformEnergy         --- WORK IN PROGRESS!!!
 *
 ******************************************************************************/
/*
   int
transformEnergy(actionInformation *action, 
		double *x, double *y, double *z, 
		double *box, int mode)
{
  char *name = "energy";
  stackType **argumentStackPointer;
  ptrajState *state;
  char *filename;
  double bondE;

  rtfInfo *rtf;
  int i;

  //  USAGE:
  //
  //    energy [prmtop filename | psf filename] 
  //
  //  action argument usage:
  //
  //  iarg1:
  //  iarg2:
  //  darg1:
  //  carg1:
  //


  if (mode == PTRAJ_SETUP) {
     // ACTION: PTRAJ_SETUP

#ifdef MPI
    printParallelError(name);
    return -1;
#endif

    argumentStackPointer = (stackType **) action->carg1;
    action->carg1 = NULL;

       //  set up the information necessary to place this on the scalarStack

    filename = argumentStackKeyToString(argumentStackPointer, "prmtop", NULL);
    if (filename == NULL) {
      action->iarg1 = 1;
      filename = argumentStackKeyToString(argumentStackPointer, "psf", NULL);
    }

    if (action->iarg1 == 0) {
      rtf = loadEnergyInfoFromPrmtop(filename);
    }


    action->carg1 = (void *) rtf;

    return 0;
  }


  rtf = (rtfInfo *) action->carg1;


  if (mode == PTRAJ_STATUS) {

    //  ACTION: PTRAJ_STATUS

    fprintf(stdout, "  ENERGY: \n");

  } else if (mode == PTRAJ_PRINT) {

    //  ACTION: PTRAJ_PRINT


  } else if (mode == PTRAJ_CLEANUP) {

    //  ACTION: PTRAJ_CLEANUP

  }


  if (mode != PTRAJ_ACTION) return 0;


  //  ACTION: PTRAJ_ACTION

  bondE = calculateBondEnergy(rtf, x, y, z);
  printf("The BOND energy is %8.4f\n", bondE);

  state = (ptrajState *) action->state;

  //  update local state information
  for (i=0; i<6; i++)
    state->box[i] = box[i];

  return 1;
}
*/




/** ACTION ROUTINE *************************************************************
 *
 *  transformGrid()   --- grid atomic densities
 *
 ******************************************************************************/


   int
transformGrid(actionInformation *action, 
	      double *x, double *y, double *z,
	      double *box, int mode)
{
  //char *name = "grid";
  argStackType **argumentStackPointer;
  char *buffer;
  ptrajState *state;
  transformGridInfo *gridInfo;
  int n, i, j, k, c;
  int isbox, negative;
  int *mask;
  double xx, yy, zz;
  double sx, sy, sz, gridMax;
  int nx, ny, nz, index;
  FILE *outFile;

  /*
   *  USAGE:
   *
   *     grid <filename> nx dx ny dy nz dz mask [origin] [negative] [max fraction] [smoothdensity value] [invert]
   *
   *  action argument usage:
   *
   *  mask
   *     the atoms for which the positions are binned (atomic density)
   *  carg1:
   *     a transformGridInfo structure
   *  iarg1:
   *     1 -- grid is centered at the origin
   *     0 -- grid is centered at the box center (if there is box information)
   *  iarg2:
   *     1 -- dump as positive density
   *     0 -- dump as negative density
   *  iarg3:
   *     1 -- invert density
   *     0 -- normal
   *  darg1: percent of max to dump
   *  darg2: madura +/- option
   *  darg3: smoothdensity option
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

    gridInfo = (transformGridInfo *) 
      safe_malloc(sizeof(transformGridInfo));
    INITIALIZE_transformGridInfo(gridInfo);
    
    gridInfo->filename = getArgumentString(argumentStackPointer, NULL);
    if (gridInfo->filename == NULL) {
      error("ptraj()", "No file was specified for grid/dipole output\n");
    }
    gridInfo->nx = getArgumentInteger(argumentStackPointer, -1);
    gridInfo->dx = getArgumentDouble( argumentStackPointer, -1);
    gridInfo->ny = getArgumentInteger(argumentStackPointer, -1);
    gridInfo->dy = getArgumentDouble( argumentStackPointer, -1);
    gridInfo->nz = getArgumentInteger(argumentStackPointer, -1);
    gridInfo->dz = getArgumentDouble( argumentStackPointer, -1);
    if (gridInfo->nx < 0 || gridInfo->ny < 0 || gridInfo->nz < 0 ||
	gridInfo->dx < 0 || gridInfo->dy < 0 || gridInfo->dz < 0) {
      error("ptraj()", "Specification of grid/dipole size or spacing\n");
    }

    if (gridInfo->nx % 2 == 1) {
      fprintf(stdout, "WARNING in ptraj(): grid -- number of grid points must be even!\n");
      gridInfo->nx++;
      fprintf(stdout, "Incrementing NX by 1 to %i\n", gridInfo->nx);
    }

    if (gridInfo->ny % 2 == 1) {
      fprintf(stdout, "WARNING in ptraj(): grid -- number of grid points must be even!\n");
      gridInfo->ny++;
      fprintf(stdout, "Incrementing NY by 1 to %i\n", gridInfo->ny);
    }

    if (gridInfo->nz % 2 == 1) {
      fprintf(stdout, "WARNING in ptraj(): grid -- number of grid points must be even!\n");
      gridInfo->nz++;
      fprintf(stdout, "Incrementing NZ by 1 to %i\n", gridInfo->nz);
    }

    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer == NULL) {
      error("ptraj()", "No mask was specified for grid/dipole\n");
    }
    action->mask = processAtomMask(buffer, action->state);
    safe_free(buffer);

    if (argumentStackContains(argumentStackPointer, "box"))
      action->iarg1 = 1; 
    else if (argumentStackContains(argumentStackPointer, "origin"))
      action->iarg1 = 0;

    if (action->iarg1 == 1 &&
	(action->state->box[3] != 90.0 || box[3] != 90.0 ||
	 action->state->box[4] != 90.0 || box[4] != 90.0 ||
	 action->state->box[5] != 90.0 || box[5] != 90.0)) {
      fprintf(stdout, "WARNING in ptraj() dipole: Code to shift to the box center\n");
      fprintf(stdout, "is not implemented yet in transformDipole for non-orthorhomibic\n");
      fprintf(stdout, "unit cells.  Shifting to the origin instead!!!\n");
      action->iarg1 = 0;
    }

    action->iarg2 = 0;
    if (argumentStackContains(argumentStackPointer, "negative"))
      action->iarg2 = 1; 

    action->darg1 = argumentStackKeyToDouble(argumentStackPointer, "max", 0.80);
    action->darg2 = argumentStackKeyToDouble(argumentStackPointer, "madura", 0.0);
    action->darg3 = argumentStackKeyToDouble(argumentStackPointer, "smoothdensity", 0.0);

    if (argumentStackContains(argumentStackPointer, "invert"))
      action->iarg3 = 1; 

    gridInfo->frames = 0;
    gridInfo->grid = (float *)
      safe_malloc(sizeof(float) * gridInfo->nx*gridInfo->ny*gridInfo->nz);
    memset(gridInfo->grid, 0,
	  sizeof(float) * gridInfo->nx*gridInfo->ny*gridInfo->nz);

    action->carg1 = (void *) gridInfo;
    return 0;
  }



  gridInfo = (transformGridInfo *) action->carg1;


  if (mode == PTRAJ_STATUS) {

    /*
     *  ACTION: PTRAJ_STATUS
     */

    fprintf(stdout, "  GRID at %s will be dumped to filename %s as %s density\n",
	    (action->iarg1 == 0 ? "origin" : "box center"),
	    gridInfo->filename,
	    (action->iarg2 == 0 ? "positive" : "negative"));
    fprintf(stdout, "      Atom selection is ");
    printAtomMask(stdout, action->mask, action->state);
    fprintf(stdout, "\n");
    fprintf(stdout, "      grid points : %5i %5i %5i\n",
	    gridInfo->nx, gridInfo->ny, gridInfo->nz);
    fprintf(stdout, "      grid spacing: %5.3f %5.3f %5.3f\n",
	    gridInfo->dx, gridInfo->dy, gridInfo->dz);
    return 0;


  } else if (mode == PTRAJ_PRINT) {

    /*
     *  ACTION: PTRAJ_PRINT
     */

    if ( gridInfo->filename == NULL ) {
      fprintf(stdout, "WARNING in ptraj(), grid: filename is NULL, not dumping grid data\n");
      return 0;
    }

    fprintf(stdout, "PTRAJ GRID dumping atomic density\n");
    if ( ( outFile = fopen(gridInfo->filename, "w") ) == NULL ) {
      fprintf(stdout, "WARNING in ptraj(), grid: couldn't open file %s\n",
	      gridInfo->filename);
      return 0;
    }

    fprintf(outFile, "This line is ignored\n");
    fprintf(outFile, "%8i\n", 1);
    fprintf(outFile, "rdparm generated grid density\n");
    fprintf(outFile, "%8i%8i%8i",
	    gridInfo->nx, -gridInfo->nx/2+1, gridInfo->nx/2);
    fprintf(outFile, "%8i%8i%8i",
	    gridInfo->ny, -gridInfo->ny/2+1, gridInfo->ny/2);
    fprintf(outFile, "%8i%8i%8i\n",
	    gridInfo->nz, -gridInfo->nz/2+1, gridInfo->nz/2);
    fprintf(outFile, "%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n",
	    (double) gridInfo->nx * gridInfo->dx,
	    (double) gridInfo->ny * gridInfo->dy,
	    (double) gridInfo->nz * gridInfo->dz,
	    90.0, 90.0, 90.0);
    fprintf(outFile, "ZYX\n");

    gridMax = 0.0;
    for (k = 0; k < gridInfo->nz; k++) {
      fprintf(outFile, "%8i\n", k - gridInfo->nz/2+1);
      for (j = 0; j < gridInfo->ny; j++) {
	c = 1;
	for (i = 0; i < gridInfo->nx; i++) {
	  
	  index = i * gridInfo->ny * gridInfo->nz +
	    j * gridInfo->nz + k;


	  if (action->darg3 > 0.0 && action->iarg3) {

	    yy = gridInfo->grid[index] - action->darg3;
	    xx = yy*yy / (0.2 * action->darg3 * action->darg3);
	    xx = exp ( -xx );

	    if (gridInfo->grid[index] > action->darg3) { 
	      gridInfo->grid[index] = -5.0;
	    } else {
	      gridInfo->grid[index] = gridInfo->grid[index] - gridInfo->grid[index] * xx;
	    }

	    /*
	    if (gridInfo->grid[index] < action->darg3) {
	      gridInfo->grid[index] = 0.0;
	    }
	    */

	    if (gridInfo->grid[index] >= 0.0) {
	      gridInfo->grid[index] = action->darg3 - gridInfo->grid[index];
	    }

	  }


	  if (action->darg3 > 0.0 && action->iarg3 == 0) {

	    yy = gridInfo->grid[index] - action->darg3;
	    xx = yy*yy / (0.2 * action->darg3 * action->darg3);
	    xx = exp ( -xx );

	    if (gridInfo->grid[index] < action->darg3) { 
	      gridInfo->grid[index] = 0.0;
	    } else {
	      gridInfo->grid[index] = gridInfo->grid[index] - gridInfo->grid[index] * xx;
	    }

	    if (gridInfo->grid[index] < action->darg3) {
	      gridInfo->grid[index] = 0.0;
	    }

	  }


	  
	  if (action->darg2 > 0.0 &&
	      gridInfo->grid[index] > 0.0 &&
	      gridInfo->grid[index] < action->darg2 )

	    /*
	     *  do the madura negative option to expose low density
	     */
	    fprintf(outFile, "%12.5f", -gridInfo->grid[index]);

	  else

	    fprintf(outFile, "%12.5f", gridInfo->grid[index]);


	  if (c && c%6 == 0)
	    fprintf(outFile, "\n");
	  c++;

	  if ( gridInfo->grid[index] > gridMax )
	    gridMax = gridInfo->grid[index];
	}
	if ( (c-1) % 6 != 0 )   /* unless a newline was just written.. */
	  fprintf(outFile, "\n");
      }
    } 
	  
    safe_fclose(outFile);

    c = 1;
    printf("PTRAJ GRID: grid max is %5.3f\n", gridMax);
    printf("            dumping a pseudo-pdb representing all points > %5.3f\n",
	   0.80 * gridMax);
	  
    for (k = 0; k < gridInfo->nz; k++) {
      for (j = 0; j < gridInfo->ny; j++) {
	for (i = 0; i < gridInfo->nx; i++) {
	  index = i * gridInfo->ny * gridInfo->nz +
	    j * gridInfo->nz + k;

	  xx = (double) i*gridInfo->dx - gridInfo->nx*gridInfo->dx/2.0
	    + 0.5 * gridInfo->dx;
	  yy = (double) j*gridInfo->dy - gridInfo->ny*gridInfo->dy/2.0
	    + 0.5 * gridInfo->dy;
	  zz = (double) k*gridInfo->dz - gridInfo->nz*gridInfo->dz/2.0
	    + 0.5 * gridInfo->dz;

	  if ( gridInfo->grid[index] > (action->darg1 * gridMax) ) {
	    fprintf(stdout, "HETATM %4i  XX  XXX   %3i     %7.3f %7.3f %7.3f\n",
		    c, c, xx, yy, zz);
	    c++;
	  }
	}
      }
    }

  } else if (mode == PTRAJ_CLEANUP) {

    /*
     *  ACTION: PTRAJ_CLEANUP
     */
    safe_free(action->mask);
    safe_free(gridInfo->filename);
    safe_free(gridInfo->grid);
    INITIALIZE_transformGridInfo(gridInfo);
    safe_free(gridInfo);
  }



  if (mode != PTRAJ_ACTION) return 0;


  /*
   *  ACTION: PTRAJ_ACTION
   */

  state = (ptrajState *) action->state;

  /*
   *  update local state information
   */
  for (i=0; i<6; i++)
    state->box[i] = box[i];

  nx = gridInfo->nx;
  ny = gridInfo->ny;
  nz = gridInfo->nz;

  sx = (double) nx * gridInfo->dx/2.0;
  sy = (double) ny * gridInfo->dy/2.0;
  sz = (double) nz * gridInfo->dz/2.0;

  isbox = action->iarg1;
  negative = action->iarg2;
  if (state->IFBOX == 0)
    isbox = 0;

  mask = action->mask;

  if (prnlev > 2) {
    printf("GRID DEBUGGING\n");
    printf("NX, NY, NZ (%i %i %i)\n", nx, ny, nz);
    printf("DX, DY, DZ (%5.3f %5.3f %5.3f)\n", gridInfo->dx, gridInfo->dy, gridInfo->dz);
    printf("Half grid is %5.3f %5.3f %5.3f\n", sx, sy, sz);
  }

  for (n=0; n < state->atoms; n++) {
    /*
     *  get coordinates and shift to origin then to half the grid,
     *  later we need to shift back to an origin reference
     */
    if (mask[n]) {
      if ( isbox ) {
	xx = x[n] - state->box[0]/2.0 + sx;
	yy = y[n] - state->box[1]/2.0 + sy;
	zz = z[n] - state->box[2]/2.0 + sz;
      } else {
	xx = x[n] + sx;
	yy = y[n] + sy;
	zz = z[n] + sz;
      }
      
      /*
       *  determine indices into grid
       */

      i = (int) (xx / gridInfo->dx ) - 1;
      j = (int) (yy / gridInfo->dy ) - 1;
      k = (int) (zz / gridInfo->dz ) - 1;

      /*
       *  check bounds and increment grid if appropriate
       */
      if (prnlev > 2)
	printf("Coords are (%5.2f %5.2f %5.2f), indices are (%i %i %i), atom %i\n", 
	       xx, yy, zz, i, j, k, n+1);

      if (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz) {
	if (negative)
	  gridInfo->grid[ i*ny*nz + j*nz + k ] -= 1.0;
	else
	  gridInfo->grid[ i*ny*nz + j*nz + k ] += 1.0;
      }
    }
  }

  return 1;

}


/** ACTION ROUTINE *************************************************************
 *
 *  transformMatrix() --- calculate matrices of covariance, correlation, distance
 *
 *  Supplementary routines:
 *    halfmatindex (below)
 *    distindex (ptraj_scalar.c)
 *    free_matrix_memory (below)
 *    lenpl (below)
 ******************************************************************************/

   int
halfmatindex(int mask1tot, int i, int j){

  /* Assure in call that i <= j.
   * Returns index for "upper right half" matrix,
   *   including main diagonal.
   */

  return (i * mask1tot - (i * (i-1) / 2) + (j - i));
}

   void
freeTransformMatrixMemory(actionInformation *action){

  transformMatrixInfo *minfo;
  coordinateInfo *cinfo;

  cinfo = (coordinateInfo *) action->carg1;
  if(cinfo != NULL){
    if(cinfo->filename != NULL)
      safe_free(cinfo->filename);
    INITIALIZE_coordinateInfo(cinfo);
    safe_free(cinfo);
  }

  minfo = (transformMatrixInfo *) action->carg2;
  if(minfo != NULL){
    if(minfo->name != NULL)
      safe_free(minfo->name);
    if(minfo->vect != NULL)
      safe_free(minfo->vect);
    if(minfo->vect2 != NULL)
      safe_free(minfo->vect2);
    if(minfo->mat != NULL)
      safe_free(minfo->mat);
    if(minfo->mask1 != NULL)
      safe_free(minfo->mask1);
    if(minfo->mask2 != NULL)
      safe_free(minfo->mask2);
    INITIALIZE_transformMatrixInfo(minfo);
    safe_free(minfo);
  }

  if(action->carg3)
    safe_free(action->carg3);

}

   double
lenpl(int order, double val){
  /* Calc Legendre polynomials
   *   (see Num. Rec. in C, p. 178 and 680)
   */

  int i;
  double twox, f2, f1, d;
  double pNplus1, pN, pNminus1;

  if(order == 0)
    return 1.0;
  else if(order == 1)
    return val;
  else{
    pNminus1 = 1.0;
    pN = val;
    twox = 2.0 * val;
    f2 = val;
    d = 1.0;

    for(i=2; i<=order; i++){
      f1 = d++;
      f2 += twox;
      pNplus1 = (f2 * pN - f1 * pNminus1) / d;
      pNminus1 = pN;
      pN = pNplus1;
    }
    return pN;
  }
}

   int
transformMatrix(actionInformation *action, 
  	        double *x, double *y, double *z,
		double *box, int mode)
{
  //char *name = "matrix";
  argStackType **argumentStackPointer;
  char *buffer;

  char *filename;
  coordinateInfo *cinfo;
  transformMatrixInfo *minfo;
  int *itmp, *mask1, *mask2, *maskA, *maskB;
  int mask1tot, mask2tot, vectsize, matsize;
  int i, j, k, l, iend, jstart, lend, toprint, toprint2;
  int crow, crowold, ccol, ccolold, ind, ind2, ind3;
  int atcnt1, atcnt2, atcnt3, atcnt4;
  int order, snap;
  double *vect, *vect2, *mat;
  double val, valnorm, mass, totmass;
  double val1, val2, val3;
  double dist1, dist2;
  double ri[3], rj[3];
  double xi, yi, zi, xj, yj, zj;
  double xk, yk, zk, xl, yl, zl;
  stackType *vectorStackTmp  = NULL; // NOTE: STACK TYPE CHANGE
  stackType *vectorStackTmp2 = NULL; // NOTEL STACK TYPE CHANGE
  transformVectorInfo *vInfo1, *vInfo2;
 
  l = 0;
  /*
   *  USAGE:
   *
   *  matrix dist|covar|mwcovar|distcovar|correl|idea|ired
   *                                             [name <name>] [order <order>]
   *                                             [<mask1>] [<mask2>] [out <filename>] 
   *                                             [start <start>] [stop <stop>] [offset <offset>]
   *                                             [byatom|byres|bymask] [mass]
   * 
   *  - If MATRIX_IRED, mask1 and mask2 are ignored and the number of matrix elements
   *      is determined by the number of vector definitions given PRIOR to the
   *      matrix command. Here, only the "upper right half" of the matrix is allocated.
   *  - Otherwise:
   *    - Upon input, ||mask1|| >= ||mask2||; this is checked below.
   *    - If only mask1 (or none) is given, only the "upper right half" of the matrix
   *        is allocated, including the main diagonal.
   *        Non-squared elements ii are contained in "vect", squared are in "vect2".
   *        This is done to be consistent if mask1 and mask2 is given and mask1 != mask2.
   *        In the case of MATRIX_DISTCOVAR and MATRIX_IRED, only (mask1tot * (mask1tot - 1)/2) resp. mask1tot
   *          elements of "vect" are used; "vect2" acts as a temporary array to store distances resp.
   *          vector lengths for each snapshot.
   *    - If both mask1 and mask2 are given, the full matrix is allocated, assuming that atoms
   *        in both masks do not necessarily correspond. (To generate full, symmetric matrices, call the 
   *        function with mask1 == mask2 upon input.)
   *    - The matrix will be stored internally with the name "name" on the matrixStack for later
   *      processing (w/ the "analyze matrix" command) ONLY if mask1 (or none) is given.
   * 
   *  - For "covar, mwcovar, distcovar, idea, ired", only "byatom" output may be chosen.
   *  - Since "distcovar, idea, ired" is mainly intended for subsequent analysis with "analyze matrix",
   *      only input of mask1 (or none) is possible.
   *
   *  action argument usage:
   *
   *    iarg1:
   *      0 -- by atom
   *      1 -- by residue
   *      2 -- by mask
   *    iarg2:
   *      unused
   *    iarg3:
   *      the number of visits
   *    iarg4:
   *      order of Legendre polynomials for MATRIX_IRED
   *    carg1:
   *      a coordinate info structure
   *    carg2:
   *      a matrix info structure
   *    carg3:
   *      0 -- no mass weighting
   *      1 -- mass weighting
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

    /*
     *  Get filename, start, stop, offset
     */

    filename = argumentStackKeyToString( argumentStackPointer, "out", NULL );

    cinfo = (coordinateInfo *) safe_malloc(sizeof(coordinateInfo));
    INITIALIZE_coordinateInfo(cinfo);
    cinfo->file = NULL;
    cinfo->filename = filename;
    cinfo->option1 = 0;
    cinfo->option2 = 0;
    cinfo->isVelocity = 0;
    cinfo->info = NULL;
    cinfo->mask = NULL;
    cinfo->start = argumentStackKeyToInteger(argumentStackPointer, "start", 1);
    cinfo->stop  = argumentStackKeyToInteger(argumentStackPointer, "stop", -1);
    if (cinfo->stop == -1) {
      cinfo->stop  = argumentStackKeyToInteger(argumentStackPointer, "end", -1);
    }
    cinfo->offset= argumentStackKeyToInteger(argumentStackPointer, "offset", 1);
    action->carg1 = (void *) cinfo;

    /*
     *  Get order
     */
    action->iarg4 = argumentStackKeyToInteger(argumentStackPointer, "order", 1);
    if(action->iarg4 <= 0){
      fprintf(stdout,
              "WARNING in ptraj(), matrix: order parameter <= 0, ignoring command\n");
      freeTransformMatrixMemory(action);
      return -1;
    }

    /*
     *  Get byatom, byres, bymask
     */

    action->iarg1 = 0;
    if ( argumentStackContains( argumentStackPointer, "byres" ) )
      action->iarg1 = 1;
    else if ( argumentStackContains( argumentStackPointer, "bymask" ) )
      action->iarg1 = 2;
    else if ( argumentStackContains( argumentStackPointer, "byatom" ) )
      action->iarg1 = 0;
    else if ( argumentStackContains( argumentStackPointer, "byatm" ) )
      action->iarg1 = 0;

    /*
     *  Get mass
     */

    itmp = (int *) safe_malloc(sizeof(int));
    *itmp = 0;
    if ( argumentStackContains( argumentStackPointer, "mass" ) )
      *itmp = 1;
    action->carg3 = (void *) itmp;

    /*
     *  Get name
     */ 

    minfo = (transformMatrixInfo *) safe_malloc(sizeof(transformMatrixInfo));
    INITIALIZE_transformMatrixInfo(minfo);

    minfo->name = argumentStackKeyToString( argumentStackPointer, "name", NULL );
    
    /*
     *  Get covar, correl, dist ...
     */

    minfo->type = MATRIX_DIST;
    if ( argumentStackContains( argumentStackPointer, "distcovar" ) )
      minfo->type = MATRIX_DISTCOVAR;
    else if ( argumentStackContains( argumentStackPointer, "mwcovar" ) )
      minfo->type = MATRIX_MWCOVAR;
    else if ( argumentStackContains( argumentStackPointer, "dist" ) )
      minfo->type = MATRIX_DIST;
    else if ( argumentStackContains( argumentStackPointer, "covar" ) )
      minfo->type = MATRIX_COVAR;
    else if ( argumentStackContains( argumentStackPointer, "correl" ) )
      minfo->type = MATRIX_CORREL;
    else if ( argumentStackContains( argumentStackPointer, "idea" ) )
      minfo->type = MATRIX_IDEA;
    else if ( argumentStackContains( argumentStackPointer, "ired" ) )
      minfo->type = MATRIX_IRED;

    /*
     *  Get mask(s)
     */

    mask1tot = 0;
    mask2tot = 0;
    buffer = getArgumentString(argumentStackPointer, NULL);
    if(minfo->type == MATRIX_IRED){
      if(buffer != NULL){
        fprintf(stdout,
	        "WARNING in ptraj(), matrix: mask input does not work with ired, ignoring command\n");
        freeTransformMatrixMemory(action);
        return -1;
      }
      else{
        for(vectorStackTmp = vectorStack;
            vectorStackTmp != NULL;
            vectorStackTmp = vectorStackTmp->next){
          vInfo1 = (transformVectorInfo *) vectorStackTmp->entry;
          if(vInfo1->mode == VECTOR_IRED)
            mask1tot++;
        }
        if(mask1tot == 0){
          fprintf(stdout,
	          "WARNING in ptraj(), matrix: no vector defined for IRED, ignoring command\n");
          freeTransformMatrixMemory(action);
          return -1;
        }
      }
    }
    else{
      if (buffer == NULL) {
        minfo->mask1 = processAtomMask( (char *) "*", action->state);
      } else {
        minfo->mask1 = processAtomMask(buffer, action->state);
        safe_free(buffer);
      }
    
      buffer = getArgumentString(argumentStackPointer, NULL);
      if (buffer == NULL) {
        minfo->mask2 = NULL;
      } else {
        minfo->mask2 = processAtomMask(buffer, action->state);
        safe_free(buffer);
      }

      for(i=0; i < action->state->atoms; i++){
        if(minfo->mask1[i])
          mask1tot++;
        if(minfo->mask2 != NULL && minfo->mask2[i])
          mask2tot++;
      }
    }

    if(mask1tot < mask2tot){
      fprintf(stdout,
	      "WARNING in ptraj(), matrix: # of atoms in mask1 < # of atoms in mask2, ignoring command\n");
      freeTransformMatrixMemory(action);
      return -1;
    }
    else{
      minfo->mask1tot = mask1tot;
      minfo->mask2tot = mask2tot;
    }

    if(minfo->name != NULL && minfo->mask2 != NULL){
      fprintf(stdout,
	      "WARNING in ptraj(), matrix: matrix only stored if no mask2,\n");
      fprintf(stdout,
              "ignoring command\n");
      freeTransformMatrixMemory(action);
      return -1;
    }

    if((minfo->type == MATRIX_COVAR || 
	minfo->type == MATRIX_MWCOVAR || 
	minfo->type == MATRIX_IRED) && action->iarg1 != 0){
      fprintf(stdout,
	      "WARNING in ptraj(), matrix: for COVAR, MWCOVAR, or IRED matrix only byatom output possible,\n");
      fprintf(stdout,
              "ignoring command\n");
      freeTransformMatrixMemory(action);
      return -1;
    }

    if((minfo->type == MATRIX_DISTCOVAR || minfo->type == MATRIX_IDEA) && 
       (minfo->mask2 != NULL || action->iarg1 != 0)){
      fprintf(stdout,
	      "WARNING in ptraj(), matrix: DISTCOVAR or IDEA matrix only generated if no mask2 and byatom output,\n");
      fprintf(stdout,
              "ignoring command\n");
      freeTransformMatrixMemory(action);
      return -1;
    }

    /*
     *  Alloc matrix memory, initialize values
     */ 

    if(minfo->type != MATRIX_NULL && minfo->type != MATRIX_DIST){ 
      /* No vector necessary for distance matrix */
      if(minfo->type == MATRIX_DISTCOVAR)
        vectsize = mask1tot * (mask1tot - 1) / 2;
      else
        vectsize = mask1tot + mask2tot;
      vect  = (double *) safe_malloc(sizeof(double) * vectsize * 3);
      vect2 = (double *) safe_malloc(sizeof(double) * vectsize * 3);

      for(i = 0; i < vectsize; i++){
        vect[i*3    ]  = 0.0;
        vect[i*3 + 1]  = 0.0;
        vect[i*3 + 2]  = 0.0;
        vect2[i*3    ] = 0.0;
        vect2[i*3 + 1] = 0.0;
        vect2[i*3 + 2] = 0.0;
      }
      minfo->vect  = vect;
      minfo->vect2  = vect2;
      minfo->vectsize = vectsize;
    }

    if(mask2tot == 0){
      /* "upper right half" matrix, including main diagonal */
      if(minfo->type == MATRIX_DISTCOVAR)
        matsize = mask1tot * (mask1tot - 1) * (mask1tot * (mask1tot - 1) / 2 + 1) / 4;
      else if(minfo->type == MATRIX_COVAR || minfo->type == MATRIX_MWCOVAR)
        matsize = 9 * mask1tot * (mask1tot + 1) / 2;
      else /* MATRIX_DIST || MATRIX_CORREL || MATRIX_IDEA || MATRIX_IRED*/
        matsize = mask1tot * (mask1tot + 1) / 2; 
    }
    else{
      /* full matrix -> no MATRIX_DISTCOVAR, MATRIX_IDEA, or MATRIX_IRED possible */
      if(minfo->type == MATRIX_COVAR || minfo->type == MATRIX_MWCOVAR)
        matsize = 9 * mask1tot * mask2tot;
      else /* MATRIX_DIST || MATRIX_CORREL */
        matsize = mask1tot * mask2tot;           
    }
    mat = (double *) safe_malloc(sizeof(double) * matsize);
    for(i = 0; i < matsize; i++){
      mat[i] = 0.0;
    }
    minfo->mat = mat;
    minfo->matsize = matsize;

    minfo->state = action->state;
    if(minfo->name){
      pushBottomStack(&matrixStack, (void *) minfo);
    }

    action->carg2 = (void *) minfo;
    action->iarg3 = 0;
    
    return 0;

  } else if (mode == PTRAJ_STATUS) {

    /*
     *  ACTION: PTRAJ_STATUS
     */

    cinfo = (coordinateInfo *)      action->carg1;
    minfo = (transformMatrixInfo *) action->carg2;

    fprintf(stdout, "  MATRIX: Calculating %s %s%s",
	    (minfo->type == MATRIX_DIST ? "distance matrix" : 
              (minfo->type == MATRIX_COVAR ? "covar matrix" : 
                (minfo->type == MATRIX_MWCOVAR ? "mass weighted covar matrix" :
                  (minfo->type == MATRIX_CORREL ? "correlation matrix" :
                    (minfo->type == MATRIX_DISTCOVAR ? "distance covar matrix" :
                      (minfo->type == MATRIX_IDEA ? "idea matrix" :
                        (minfo->type == MATRIX_IRED ? "ired matrix" : "Should not occur")
                      )
                    )
                  )
                )
              )
            ),
	    (action->iarg1 == 2 ? "by mask" : (action->iarg1 == 1 ? "by residue" : "by atom")),
	    (cinfo->filename == NULL ? "" : ", dumping to file "));
    if (cinfo->filename != NULL)
      fprintf(stdout, "%s, using ", cinfo->filename);
    else
      fprintf(stdout, " using ");
    fprintf(stdout, "%s\n",
            (*((int *) action->carg3) == 0 ? "no mass weighting" : "mass weighting"));
    if(minfo->type == MATRIX_IRED)
      fprintf(stdout, "      Order of Legendre polynomials: %i\n", action->iarg4);
    if(minfo->name != NULL)
      fprintf(stdout, "      Storing matrix on internal stack with name: %s\n",minfo->name);
    if (cinfo->start != 1 || cinfo->stop != -1 || cinfo->offset != 1) {
      fprintf(stdout, "      start: %i", cinfo->start);
      if (cinfo->stop > 0)
	fprintf(stdout, "  stop: %i", cinfo->stop);
      else
	fprintf(stdout, "  stop [at final frame]");
      fprintf(stdout, "  offset: %i\n", cinfo->offset);
    }
    if(minfo->type != MATRIX_IRED){
      fprintf(stdout, "      Atom selection 1 follows ");
      printAtomMask(stdout, minfo->mask1, action->state);
      fprintf(stdout, "\n");
      if(minfo->mask2 != NULL){
        fprintf(stdout, "      Atom selection 2 follows ");
        printAtomMask(stdout, minfo->mask2, action->state);
	fprintf(stdout, "\n");
      }
    }

  } else if (mode == PTRAJ_PRINT) {

    /*
     *  ACTION: PTRAJ_PRINT
     */

    cinfo = (coordinateInfo *) action->carg1;
    minfo = (transformMatrixInfo *) action->carg2;
    vect  = minfo->vect;
    vect2  = minfo->vect2;
    vectsize = minfo->vectsize;
    mat   = minfo->mat;
    matsize = minfo->matsize;
    mask1 = minfo->mask1;
    mask1tot = minfo->mask1tot;
    if(minfo->mask2 != NULL){
      mask2 = minfo->mask2;
      mask2tot = minfo->mask2tot;
    }
    else{
      mask2 = mask1;
      mask2tot = minfo->mask1tot;
    }
    snap = minfo->snap;

    /*
     *  Calc average over number of sets
     */
    if(vect != NULL){
      for(i=0; i < vectsize; i++){
        vect[i*3  ]  /= (double) snap;
        vect[i*3+1]  /= (double) snap;
        vect[i*3+2]  /= (double) snap;
        vect2[i*3  ] /= (double) snap;
        vect2[i*3+1] /= (double) snap;
        vect2[i*3+2] /= (double) snap;
      }
    }
    for(i=0; i < matsize; i++){
      mat[i] /= (double) snap;
    }

    if(minfo->type == MATRIX_DIST || minfo->type == MATRIX_IRED){
      /*
       * Nothing to do in case of distance or ired matrix
       */
    }
    else if(minfo->type == MATRIX_IDEA){
      for(i=0; i < vectsize; i++){
        vect[i*3  ]  /= 3.0;
        vect[i*3+1]  /= 3.0;
        vect[i*3+2]  /= 3.0;
        vect2[i*3  ] /= 3.0;
        vect2[i*3+1] /= 3.0;
        vect2[i*3+2] /= 3.0;
      }
      for(i=0; i < matsize; i++){
        mat[i] /= 3.0;
      }
    }
    else if(minfo->type == MATRIX_DISTCOVAR){
      /*
       * Calc distance covariance matrix
       */      
      atcnt1 = 0;
      for(i=0; i < action->state->atoms; i++){
        if(mask2[i]){
          atcnt2 = atcnt1 + 1;
          for(j=i+1; j < action->state->atoms; j++){
            if(mask2[j]){
              ind = distindex(mask1tot, atcnt1, atcnt2);
              atcnt3 = atcnt1;
              for(k=i; k < action->state->atoms; k++){
                if(mask1[k]){
                  atcnt4 = (k>=j ? atcnt3+1 : atcnt2);
                  for(l=(k>=j ? k+1 : j); l < action->state->atoms; l++){
                    if(mask1[l]){
                      ind2 = distindex(mask1tot, atcnt3, atcnt4);
                      ind3 = halfmatindex(mask1tot * (mask1tot - 1) / 2, ind, ind2);
                      /*                      
                      printf("%i(%i) %i(%i) -> %i ||| %i(%i) %i(%i) -> %i ||| %i\n", 
                             atcnt1, i, atcnt2, j, ind,
                             atcnt3, k, atcnt4, l, ind2,
                             ind3);
                      */
                      mat[ind3] -= vect[ind] * vect[ind2];
                      atcnt4++;
                    }
                  } /* end for l */
                  atcnt3++;
                }
              } /* end for k */
              atcnt2++;
            }
          } /* end for j */
          atcnt1++;
        }
      } /* end for i */
    }
    else{ /* MATRIX_COVAR || MATRIX_MWCOVAR || MATRIX_CORREL */
      /*
       * Calc covariance or correlation matrix
       */

      /* Calc <riri> - <ri><ri> */
      for(i = 0; i < vectsize; i++){
        vect2[i*3  ] -= vect[i*3  ] * vect[i*3  ];
        vect2[i*3+1] -= vect[i*3+1] * vect[i*3+1];
        vect2[i*3+2] -= vect[i*3+2] * vect[i*3+2];
      }

      /* Calc <rirj> - <ri><rj> */
      ind  = 0;
      ind2 = (mask1 == mask2 ? 0 : mask1tot);
      lend = (minfo->type == MATRIX_COVAR || minfo->type == MATRIX_MWCOVAR ? 3 : 1);
      for(i = 0; i < mask2tot; i++){
        for(l = 0; l < lend; l++){
          for(j = 0; j < mask1tot; j++){
            if((mask1 == mask2 && j >= i) || mask1 != mask2){
              if(minfo->type == MATRIX_COVAR || minfo->type == MATRIX_MWCOVAR){
                if(mask1 == mask2 && i == j){                
                  if(l == 0){
                    mat[ind++] -= vect[(ind2 + i)*3  ] * vect[j*3  ];
                    mat[ind++] -= vect[(ind2 + i)*3  ] * vect[j*3+1];
                    mat[ind++] -= vect[(ind2 + i)*3  ] * vect[j*3+2];
                  }
                  else if(l == 1){
                    mat[ind++] -= vect[(ind2 + i)*3+1] * vect[j*3+1];
                    mat[ind++] -= vect[(ind2 + i)*3+1] * vect[j*3+2];
                  }
                  else if(l == 2){
                    mat[ind++] -= vect[(ind2 + i)*3+2] * vect[j*3+2];
                  }
                }
                //else if(mask1 == mask2 && i < j || mask1 != mask2){
                else if(((mask1 == mask2) && (i < j)) || (mask1 != mask2)){
                  if(l == 0){
                    mat[ind++] -= vect[(ind2 + i)*3  ] * vect[j*3  ];
                    mat[ind++] -= vect[(ind2 + i)*3  ] * vect[j*3+1];
                    mat[ind++] -= vect[(ind2 + i)*3  ] * vect[j*3+2];
                  }
                  else if(l == 1){
                    mat[ind++] -= vect[(ind2 + i)*3+1] * vect[j*3  ];
                    mat[ind++] -= vect[(ind2 + i)*3+1] * vect[j*3+1];
                    mat[ind++] -= vect[(ind2 + i)*3+1] * vect[j*3+2];
                  }
                  else if(l == 2){
                    mat[ind++] -= vect[(ind2 + i)*3+2] * vect[j*3  ];
                    mat[ind++] -= vect[(ind2 + i)*3+2] * vect[j*3+1];
                    mat[ind++] -= vect[(ind2 + i)*3+2] * vect[j*3+2];
                  }
                }
              }
              else if(minfo->type == MATRIX_CORREL){ 
                mat[ind] -= (vect[j*3  ] * vect[(ind2 + i)*3  ] +
                             vect[j*3+1] * vect[(ind2 + i)*3+1] +
                             vect[j*3+2] * vect[(ind2 + i)*3+2]);
                /* Normalize */
                mat[ind] /= sqrt((vect2[j*3] + vect2[j*3+1] + vect2[j*3+2]) *
                                 (vect2[(ind2 + i)*3] + vect2[(ind2 + i)*3+1] + vect2[(ind2 + i)*3+2]));
                ind++;
              }
            }
          }
        }
      }

      if(minfo->type == MATRIX_MWCOVAR){
        ind = 0;
        crow = 0;
        for (i=0; i < action->state->atoms; i++) {
          if(mask2[i]){
            for(k = 0; k < 3; k++){
              ccol = 0;          
              for (j=0; j < action->state->atoms; j++) {
                if(mask1[j]){
                  if(i == j){
                    mass = action->state->masses[i];
                    vect2[ccol*3]   *= mass;
                    vect2[ccol*3+1] *= mass;
                    vect2[ccol*3+2] *= mass;
                  }
                
                  if(mask1 == mask2 && j >= i){ 
                    /*** "upper right half" matrix ***/
                    mass = sqrt(action->state->masses[i] * action->state->masses[j]);
                    if(crow*3+k <= ccol*3){
                      ind = halfmatindex(mask1tot * 3, crow*3+k, ccol*3);
                      mat[ind] *= mass;
                    }
                    if(crow*3+k <= ccol*3+1){
                      ind = halfmatindex(mask1tot * 3, crow*3+k, ccol*3+1);
                      mat[ind] *= mass;
                    }
                    if(crow*3+k <= ccol*3+2){
                      ind = halfmatindex(mask1tot * 3, crow*3+k, ccol*3+2);
                      mat[ind] *= mass;
                    }
                  }
                  else if(mask1 != mask2){ 
                    /*** full matrix ***/
                    mass = sqrt(action->state->masses[i] * action->state->masses[j]);
                    mat[ind++] *= mass;
                    mat[ind++] *= mass;
                    mat[ind++] *= mass;
                  }
                  /* Next column in matrix */
                  ccol++;
                }
              }
            }
            /* Next row in matrix */
            crow++;
          }
        }
      }
    }
    
    if (cinfo->filename){
      cinfo->file = safe_fopen(cinfo->filename, "w");
      if (cinfo->file == NULL) {
        fprintf(stdout, "WARNING in ptraj(), matrix: error on opening %s for output\n",
	        cinfo->filename);
        return 0;
      }
      fprintf(stdout, "PTRAJ MATRIX: Dumping matrix values\n");
    }

    if (action->iarg1 == 0) {
      /*
       *  byatom print out
       */
      if(minfo->type == MATRIX_DISTCOVAR){
        atcnt1 = 0;
        for(i=0; i < action->state->atoms; i++){
          if(mask2[i]){
            atcnt2 = atcnt1 + 1;
            for(j=i+1; j < action->state->atoms; j++){
              if(mask2[j]){
                ind = distindex(mask1tot, (atcnt1 <= atcnt2 ? atcnt1 : atcnt2),
                                          (atcnt1 <= atcnt2 ? atcnt2 : atcnt1));
                atcnt3 = 0;
                for(k=0; k < action->state->atoms; k++){
                  if(mask1[k]){
                    atcnt4 = atcnt3 + 1;
                    for(l=k+1; l < action->state->atoms; l++){
                      if(mask1[l]){
                        ind2 = distindex(mask1tot, (atcnt3 <= atcnt4 ? atcnt3 : atcnt4),
                                                   (atcnt3 <= atcnt4 ? atcnt4 : atcnt3));
                        ind3 = halfmatindex(mask1tot * (mask1tot - 1) / 2, 
                                            (ind <= ind2 ? ind : ind2),
                                            (ind <= ind2 ? ind2 : ind));
                        /*
                        printf("%i(%i) %i(%i) -> %i ||| %i(%i) %i(%i) -> %i ||| %i\n", 
                               atcnt1, i, atcnt2, j, ind,
                               atcnt3, k, atcnt4, l, ind2,
                               ind3);
                        */
                        if(cinfo->file != NULL)
                          fprintf(cinfo->file,"%6.2f ",mat[ind3]);
                        atcnt4++;
                      }
                    } /* end for l */
                    atcnt3++;
                  }
                } /* end for k */
                if(cinfo->file != NULL)
                  fprintf(cinfo->file,"\n");
                atcnt2++;
              }
            } /* end for j */
            atcnt1++;
          }
        } /* end for i */
        if(cinfo->file != NULL)
          fprintf(cinfo->file,"\n");

      }
      else if(minfo->type == MATRIX_IRED){
        crow = 0;
        for(vectorStackTmp = vectorStack;
            vectorStackTmp != NULL;
            vectorStackTmp = vectorStackTmp->next){
          vInfo1 = (transformVectorInfo *) vectorStackTmp->entry;
          if(vInfo1->mode == VECTOR_IRED){
            ccol = 0;          
            for(vectorStackTmp2 = vectorStack;
                vectorStackTmp2 != NULL;
                vectorStackTmp2 = vectorStackTmp2->next){
              vInfo2 = (transformVectorInfo *) vectorStackTmp2->entry;
              if(vInfo2->mode == VECTOR_IRED){
                ind  = halfmatindex(mask1tot,
                                    (crow <= ccol ? crow : ccol),
                                    (crow <= ccol ? ccol : crow));
                val = mat[ind];
                if(cinfo->file != NULL)
                  fprintf(cinfo->file,"%6.3f ", val);

                /* Next column in matrix */
                ccol++;
              }
            }
            if(cinfo->file != NULL)
              fprintf(cinfo->file,"\n");
          }
          /* Next row in matrix */
          crow++;
        }
        if(cinfo->file != NULL)
          fprintf(cinfo->file,"\n");
      }
      else{
        ind = 0;
        crow = 0;
        lend = (minfo->type == MATRIX_COVAR || minfo->type == MATRIX_MWCOVAR ? 3 : 1);
        for (i=0; i < action->state->atoms; i++) {
          if(mask2[i]){
            for(l = 0; l < lend; l++){
              ccol = 0;          
              for (j=0; j < action->state->atoms; j++) {
                if(mask1[j]){
                  if(mask1 == mask2){ 
                    if(minfo->type == MATRIX_COVAR || minfo->type == MATRIX_MWCOVAR){
                      ind  = halfmatindex(mask1tot * 3,
                                          (crow*3+l <= ccol*3 ? crow*3+l : ccol*3),
                                          (crow*3+l <= ccol*3 ? ccol*3   : crow*3+l));
                      ind2 = halfmatindex(mask1tot * 3,
                                          (crow*3+l <= ccol*3+1 ? crow*3+l : ccol*3+1),
                                          (crow*3+l <= ccol*3+1 ? ccol*3+1 : crow*3+l));
                      ind3 = halfmatindex(mask1tot * 3,
                                          (crow*3+l <= ccol*3+2 ? crow*3+l : ccol*3+2),
                                          (crow*3+l <= ccol*3+2 ? ccol*3+2 : crow*3+l));
                      if(cinfo->file != NULL)
                        fprintf(cinfo->file,"%6.3f %6.3f %6.3f ", mat[ind], mat[ind2], mat[ind3]);
                    }
                    else{
                      ind  = halfmatindex(mask1tot,
                                          (crow <= ccol ? crow : ccol),
                                          (crow <= ccol ? ccol : crow));
                      val = mat[ind];
                      if(cinfo->file != NULL)
                        fprintf(cinfo->file,"%6.3f ", val);
                    }
                  }
                  else{ 
                    /*** full matrix ***/
                    if(minfo->type == MATRIX_COVAR || minfo->type == MATRIX_MWCOVAR){
                      val1 = mat[ind++];
                      val2 = mat[ind++];
                      val3 = mat[ind++];
                      if(cinfo->file != NULL)
                        fprintf(cinfo->file,"%6.3f %6.3f %6.3f ", val1, val2, val3);
                    }
                    else{
                      val = mat[ind++];
                      if(cinfo->file != NULL)
                        fprintf(cinfo->file,"%6.3f ", val);
                    }
                  }
                  /* Next column in matrix */
                  ccol++;
                }
              }
              if(cinfo->file != NULL)
                fprintf(cinfo->file,"\n");
            }
            /* Next row in matrix */
            crow++;
          }
        }
        if(cinfo->file != NULL)
          fprintf(cinfo->file,"\n");
      }
    } 
    else if (action->iarg1 == 1) {
      /*
       *  byres print out
       */
      crow = 0;
      for (i=0; i < action->state->residues; i++) {
        toprint = 0;
        /* Store actual row value */
        crowold = crow; 
        /* Init column value */
        ccol = 0;          
        for (k=0; k < action->state->residues; k++) {
          toprint2 = 0;
          val = 0.0;
          valnorm = 0.0;
          /* Restore row value */
          crow = crowold;
          /* Store actual column value */
          ccolold = ccol;
  	  for (j=action->state->ipres[i]-1; j<action->state->ipres[i+1]-1; j++) {
            if(mask2[j]){
              /* Restore column value */
              ccol = ccolold;
              
   	      for (l=action->state->ipres[k]-1; l<action->state->ipres[k+1]-1; l++) {
                if(mask1[l]){
                  mass = *((int *) action->carg3) == 0 ? 1.0 : 
                           action->state->masses[j] * action->state->masses[l];
                  valnorm += mass;
                  toprint = toprint2 = 1;
                  if(mask1 == mask2){ 
                    /*** "upper right half" matrix ***/
                    ind = halfmatindex(mask1tot,
                                       (crow <= ccol ? crow : ccol),
                                       (crow <= ccol ? ccol : crow));
                    val += mat[ind] * mass;
                  }
                  else{ 
                    /*** full matrix ***/
                    val += mat[crow * mask1tot + ccol] * mass;
                  }
                  /* Next column in matrix */
                  ccol++;
                }
              }
              /* Next row in matrix */
              crow++;
            }
          }
          if(toprint2 && cinfo->file != NULL)
            fprintf(cinfo->file,"%6.2f ",val / valnorm);
        }
        if(toprint && cinfo->file != NULL)
          fprintf(cinfo->file,"\n");
      }
      if(cinfo->file != NULL)
        fprintf(cinfo->file,"\n");
    } 
    else if (action->iarg1 == 2) {
      /*
       *  bymask print out
       *
       *  if(mask1 == mask2):
       *    internal average over mask1
       *  else:
       *    mask1/mask1 mask1/mask2 mask2/mask2
       */
      if(mask1 == mask2){
        iend = 1;
      }
      else{
        iend = 3;
      }            
      maskA = mask1;
      maskB = mask1;

      for(i = 0; i < iend; i++){
        if(i > 0){
          maskA = maskB;
          maskB = mask2;
        }
        val = 0.0;
        valnorm = 0.0;
        crow = 0;
        for (j=0; j < action->state->atoms; j++) {
          if(maskB[j]){
            ccol = 0;          
            for (k=0; k < action->state->atoms; k++) {
              if(maskA[k]){
                mass = *((int *) action->carg3) == 0 ? 1.0 : 
                         action->state->masses[j] * action->state->masses[l];
                valnorm += mass;
                if(mask1 == mask2){ 
                  /*** "upper right half" matrix ***/
                  ind = halfmatindex(mask1tot,
                                     (crow <= ccol ? crow : ccol),
                                     (crow <= ccol ? ccol : crow));
                  val += mat[ind] * mass;
                }
                else{ 
                  /*** full matrix ***/
                  val += mat[crow * mask1tot + ccol] * mass;
                }
                /* Next column in matrix */
                ccol++;
              }
            }
            /* Next row in matrix */
            crow++;
          }
        }
        if(cinfo->file != NULL)
          fprintf(cinfo->file,"%6.2f ",val / valnorm);
      }
      if(cinfo->file != NULL)
        fprintf(cinfo->file,"\n");
    }

    if (cinfo->file != NULL) {
      safe_fclose(cinfo->file);
      cinfo->file = NULL;
    }

  }
  else if (mode == PTRAJ_CLEANUP) {

    /*
     *  ACTION: PTRAJ_CLEANUP
     */

    freeTransformMatrixMemory(action);
  }

  if (mode != PTRAJ_ACTION) return 0;

  /*
   *  ACTION: PTRAJ_ACTION
   */

  cinfo  = (coordinateInfo *)      action->carg1;
  minfo  = (transformMatrixInfo *) action->carg2;
  order  = action->iarg4;
  vect   = minfo->vect;
  vect2  = minfo->vect2;
  vectsize = minfo->vectsize;
  mat    = minfo->mat;
  matsize = minfo->matsize; 
  mask1  = minfo->mask1;
  if(minfo->mask2 != NULL){
    mask2 = minfo->mask2;
  }
  else{
    mask2 = mask1;
  }
  mask1tot = minfo->mask1tot;
  mask2tot = minfo->mask2tot;

  action->iarg3++;
  if (action->iarg3 >= cinfo->start &&
      (cinfo->stop < 0 || action->iarg3 <= cinfo->stop) &&
      (action->iarg3 - cinfo->start)%cinfo->offset == 0) {

    minfo->snap++;

    if(minfo->type == MATRIX_DIST){
      /*
       * Calc distance matrix
       */
      ind = 0;
      for (i=0; i < action->state->atoms; i++) {
        if(mask2[i]){
	  xi = x[i]; 
	  yi = y[i]; 
          zi = z[i];
          jstart = (mask1 == mask2 ? i : 0);
          for (j=jstart; j < action->state->atoms; j++) {
            if(mask1[j]){
              xj = x[j];
              yj = y[j]; 
              zj = z[j];
              mat[ind] += sqrt( (xi - xj)*(xi - xj) + 
                                (yi - yj)*(yi - yj) +
                                (zi - zj)*(zi - zj) );
              ind++;
            }
          }
        }
      }
    }
    else if(minfo->type == MATRIX_COVAR || minfo->type == MATRIX_MWCOVAR || minfo->type == MATRIX_CORREL){
      /*
       * Calc covariance matrix or correlation matrix
       */
      ind  = 0;
      ind2 = 0;
      ind3 = 0;
      k    = 1;      
      lend = (minfo->type == MATRIX_COVAR || minfo->type == MATRIX_MWCOVAR ? 3 : 1);
      for (i=0; i < action->state->atoms; i++) {
        if(mask2[i]){
 	  xi = x[i]; 
	  yi = y[i]; 
          zi = z[i];
          if(mask1 != mask2){
            vect[(mask1tot + ind2)*3  ] += xi;
            vect[(mask1tot + ind2)*3+1] += yi;
            vect[(mask1tot + ind2)*3+2] += zi;
            vect2[(mask1tot + ind2)*3  ] += xi*xi;
            vect2[(mask1tot + ind2)*3+1] += yi*yi;
            vect2[(mask1tot + ind2)*3+2] += zi*zi;
            ind2++;
          }
          for(l = 0; l < lend; l++){
            for (j=0; j < action->state->atoms; j++) {
              if(mask1[j]){
                xj = x[j];
                yj = y[j]; 
                zj = z[j];
                if(k == 1){  
                  vect[ind3*3  ] += xj;
                  vect[ind3*3+1] += yj;
                  vect[ind3*3+2] += zj;
                  vect2[ind3*3  ] += xj*xj;
                  vect2[ind3*3+1] += yj*yj;
                  vect2[ind3*3+2] += zj*zj;
                  ind3++;
                }
                if((mask1 == mask2 && j >= i) || mask1 != mask2){
                  if(minfo->type == MATRIX_COVAR || minfo->type == MATRIX_MWCOVAR){
                    if(mask1 == mask2 && i == j){
                      if(l == 0){
                        mat[ind++] += xi*xj;
                        mat[ind++] += xi*yj;
                        mat[ind++] += xi*zj;
                      }
                      else if(l == 1){
                        mat[ind++] += yi*yj;
                        mat[ind++] += yi*zj;
                      }
                      else if(l == 2)
                        mat[ind++] += zi*zj;
                    }
                    else if((mask1 == mask2 && i < j) || mask1 != mask2){
                      if(l == 0){
                        mat[ind++] += xi*xj;
                        mat[ind++] += xi*yj;
                        mat[ind++] += xi*zj;
                      }
                      else if(l == 1){
                        mat[ind++] += yi*xj;
                        mat[ind++] += yi*yj;
                        mat[ind++] += yi*zj;
                      }
                      else if(l == 2){
                        mat[ind++] += zi*xj;
                        mat[ind++] += zi*yj;
                        mat[ind++] += zi*zj;
                      }
                    }
                  }
                  else{
                    mat[ind++] += xi*xj + yi*yj + zi*zj;
                  }
                }
              }
            }
            k = 0;
          }
        }
      }
    }
    else if(minfo->type == MATRIX_IDEA){
      /*
       * Calc isotropically distributed ensemble matrix
       *   (see Proteins 2002, 46, 177; eq. 7)
       */

      /*
       * Find center of mass coordinates
       */
      val1 = val2 = val3 = 0.0;
      totmass = 0.0;
      for (i=0; i < action->state->atoms; i++) {
        if(mask2[i]){
          mass = action->state->masses[i];
          totmass += mass;
          val1 += mass * x[i];
          val2 += mass * y[i];
          val3 += mass * z[i];
        }      
      }
      val1 /= totmass;
      val2 /= totmass;
      val3 /= totmass;

      /*
       * Get ri, rj and calc ri * rj
       */
      ind = ind2 = 0;
      for (i=0; i < action->state->atoms; i++) {
        if(mask2[i]){
          ri[0] = x[i] - val1;
          ri[1] = y[i] - val2;
          ri[2] = z[i] - val3;
          for (j=i; j < action->state->atoms; j++) {
            if(mask1[j]){
              rj[0] = x[j] - val1;
              rj[1] = y[j] - val2;
              rj[2] = z[j] - val3;
              val = ri[0] * rj[0] + ri[1] * rj[1] + ri[2] * rj[2];
              mat[ind++] += val;
              if(j == i){
                vect[ind2] += val;
                vect2[ind2] += (val * val);
                ind2++;
              }
            }
          }
        }
      }  
    }
    else if(minfo->type == MATRIX_IRED){
      /*
       * Calc isotropic reorientational eigenmode dynamics
       *   (see JACS 2002, 124, 4522, eq. A14;
       *    CAVEAT: omegaK-omegaL is not "just" the intra
       *            molecular angle there)
       */    
     
      /*
       * Store length of vectors in vect2
       */
      ind = 0;
      for(vectorStackTmp = vectorStack;
          vectorStackTmp != NULL;
          vectorStackTmp = vectorStackTmp->next){
        vInfo1 = (transformVectorInfo *) vectorStackTmp->entry;
        if(vInfo1->mode == VECTOR_IRED){
          if(ind >= mask1tot){
            /* 
             * This can happen if IRED vectors are defined
             *   after IRED matrix command
             */
            fprintf(stdout, "WARNING in ptraj(), matrix: IRED vectors defined after IRED matrix command\n");
            return 0;
          }
          vect2[ind++] = sqrt(vInfo1->vx[0] * vInfo1->vx[0] +
                              vInfo1->vy[0] * vInfo1->vy[0] +
                              vInfo1->vz[0] * vInfo1->vz[0]);
        }
      }

      ind = ind2 = 0;
      for(vectorStackTmp = vectorStack;
          vectorStackTmp != NULL;
          vectorStackTmp = vectorStackTmp->next){
        vInfo1 = (transformVectorInfo *) vectorStackTmp->entry;
        if(vInfo1->mode == VECTOR_IRED){
          val1 = vect2[ind2];
          ind3 = ind2;
          for(vectorStackTmp2 = vectorStackTmp;
              vectorStackTmp2 != NULL;
              vectorStackTmp2 = vectorStackTmp2->next){
            vInfo2 = (transformVectorInfo *) vectorStackTmp2->entry;
            if(vInfo2->mode == VECTOR_IRED){
              val2 = vect2[ind3++];
              val3 = lenpl(order, (vInfo1->vx[0]*vInfo2->vx[0] +
                                   vInfo1->vy[0]*vInfo2->vy[0] +
                                   vInfo1->vz[0]*vInfo2->vz[0]) /
                                  (val1 * val2));
              mat[ind++] += val3;
              if(vectorStackTmp2 == vectorStackTmp)
                vect[ind2++] += val3; 
            }
          }
        }
      }
    }
    else if(minfo->type == MATRIX_DISTCOVAR){
      /*
       * Calc distance covariance matrix
       */
      atcnt1 = 0;
      for(i=0; i < action->state->atoms; i++){
        if(mask2[i]){
 	  xi = x[i]; 
	  yi = y[i]; 
          zi = z[i];
          atcnt2 = atcnt1 + 1;
          for(j=i+1; j < action->state->atoms; j++){
            if(mask2[j]){
              xj = x[j];
              yj = y[j]; 
              zj = z[j];
              ind = distindex(mask1tot, atcnt1, atcnt2);
              if(atcnt2 > 1){
                dist1 = vect2[ind];
              }
              else{
                dist1 = sqrt((xi - xj) * (xi - xj) +
                             (yi - yj) * (yi - yj) +
                             (zi - zj) * (zi - zj));
              }
              
              atcnt3 = atcnt1;
              for(k=i; k < action->state->atoms; k++){
                if(mask1[k]){
    	          xk = x[k]; 
	          yk = y[k]; 
                  zk = z[k];
                  atcnt4 = (k>=j ? atcnt3+1 : atcnt2);
                  for(l=(k>=j ? k+1 : j); l < action->state->atoms; l++){
                    if(mask1[l]){
    	              xl = x[l]; 
	              yl = y[l]; 
                      zl = z[l];
                      ind2 = distindex(mask1tot, atcnt3, atcnt4);
                      if(atcnt2 > 1){
                        dist2 = vect2[ind2];
                      }
                      else{
                        dist2 = sqrt((xk - xl) * (xk - xl) +
                                     (yk - yl) * (yk - yl) +
                                     (zk - zl) * (zk - zl));
                        vect2[ind2] = dist2;
                      }

                      ind3 = halfmatindex(mask1tot * (mask1tot - 1) / 2, ind, ind2);

                      /*
                      printf("%i(%i) %i(%i) -> %i: %f ||| %i(%i) %i(%i) -> %i: %f ||| %i\n", 
                             atcnt1, i, atcnt2, j, ind,  dist1,
                             atcnt3, k, atcnt4, l, ind2, dist2,
                             ind3);
                      */
                      
                      mat[ind3] += dist1 * dist2;
                      if(ind == ind2)
                        vect[ind] += dist1;

                      atcnt4++;
                    }
                  } /* end for l */
                  atcnt3++;
                }
              } /* end for k */
              atcnt2++;
            }
          } /* end for j */
          atcnt1++;
        }
      } /* end for i */
    }

  }

  return 1;
}

/** ACTION ROUTINE *************************************************************
 *
 *  transformPrincipal()   --- align coordinates along the principal axis
 *
 *  Supplementary routines:
 *
 *    jacobi() -- diagonalization
 *    jabobiCheckChirality()
 *    calculatePrincipalAxis()
 *
 ******************************************************************************/

#define MAX_ITERATIONS 50
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
          a[k][l]=h+s*(g-h*tau);

void jacobi(double a[3][3],int n,double d[3], double v[3][3])
{
  int j,iq,ip,i,nrot;
  double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
  
  b=safe_malloc(sizeof(double) * n);
  z=safe_malloc(sizeof(double) * n);

  for (ip=1;ip<=n;ip++) {
    for (iq=1;iq<=n;iq++) v[ip-1][iq-1]=0.0;
    v[ip-1][ip-1]=1.0;
  }
  for (ip=1;ip<=n;ip++) {
    b[ip-1]=d[ip-1]=a[ip-1][ip-1];
    z[ip-1]=0.0;
  }
  nrot=0;
  for (i=1;i<=MAX_ITERATIONS;i++) {
    sm=0.0;
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++)
	sm += fabs(a[ip-1][iq-1]);
    }

    if (sm == 0.0) {
      safe_free(b);
      safe_free(z);
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++) {
	g=100.0*fabs(a[ip-1][iq-1]);
	if (i > 4 && fabs(d[ip-1])+g == fabs(d[ip-1])
	    && fabs(d[iq-1])+g == fabs(d[iq-1]))
	  a[ip-1][iq-1]=0.0;
	else if (fabs(a[ip-1][iq-1]) > tresh) {
	  h=d[iq-1]-d[ip-1];
	  if (fabs(h)+g == fabs(h))
	    t=(a[ip-1][iq-1])/h;
	  else {
	    theta=0.5*h/(a[ip-1][iq-1]);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*a[ip-1][iq-1];
	  z[ip-1] -= h;
	  z[iq-1] += h;
	  d[ip-1] -= h;
	  d[iq-1] += h;
	  a[ip-1][iq-1]=0.0;
	  for (j=1;j<=ip-1;j++) {
	    ROTATE(a,j-1,ip-1,j-1,iq-1)
	    }
	  for (j=ip+1;j<=iq-1;j++) {
	    ROTATE(a,ip-1,j-1,j-1,iq-1)
	    }
	  for (j=iq+1;j<=n;j++) {
	    ROTATE(a,ip-1,j-1,iq-1,j-1)
	    }
	  for (j=1;j<=n;j++) {
	    ROTATE(v,j-1,ip-1,j-1,iq-1)
	    }
	  ++nrot;
	}
      }
    }
    for (ip=1;ip<=n;ip++) {
      b[ip-1] += z[ip-1];
      d[ip-1]=b[ip-1];
      z[ip-1]=0.0;
    }
  }
  error("jacobi()", "Too many iterations in routine JACOBI");
}

#undef ROTATE
#undef MAX_ITERATIONS




/*  The jacobi diagonalization procedure can sometimes result
 *  in eigenvectors which when applied to transform the coordinates
 *  result in a a chiral inversion about the Y axis.  This code catches
 *  this case, reversing the offending eigenvectors.
 *  
 *  NOTE: the idea of rotating the coordinate basis vectors came from 
 *  some code posted to the computational chemistry mailing list 
 *  (chemistry@osc) in a summary of methods to perform principal axis 
 *  alignment...
 */

   int
jacobiCheckChirality(double evalue[3], double ev[3][3])
{
  double points[3][3], result[3][3];
  double transform[3][3];
  double xtemp, ytemp, ztemp;
  double r;

  /* transform the coordinate basis vectors (identity matrix) 
   * to check for chiral inversion...
   */
  points[0][0] = 1.0; points[0][1] = 0.0; points[0][2] = 0.0;
  points[1][0] = 0.0; points[1][1] = 1.0; points[1][2] = 0.0;
  points[2][0] = 0.0; points[2][1] = 0.0; points[2][2] = 1.0;

  VOP_3x3_TRANSPOSE_TIMES_COORDS(ev, 
				 points[0][0], points[1][0], points[2][0],
				 xtemp, ytemp, ztemp);
  VOP_3x3_TRANSPOSE_TIMES_COORDS(ev, 
				 points[0][1], points[1][1], points[2][1],
				 xtemp, ytemp, ztemp);
  VOP_3x3_TRANSPOSE_TIMES_COORDS(ev, 
				 points[0][2], points[1][2], points[2][2],
				 xtemp, ytemp, ztemp);

  /* rotate vector three into XZ plane */
  r = sqrt( points[0][2] * points[0][2] + points[1][2] * points[1][2] );
  transform[0][0] = points[0][2] / r;
  transform[1][1] = points[0][2] / r;
  transform[0][1] = points[1][2] / r;
  transform[1][0] = -points[1][2] / r;
  transform[2][2] = 1.0;
  transform[0][2] = 0.0;
  transform[1][2] = 0.0;
  transform[2][0] = 0.0;
  transform[2][1] = 0.0;
  VOP_3x3_TIMES_3x3(result, transform, points);

  /* rotate vector three into Z axis */
  r = sqrt( result[0][2] * result[0][2] + result[2][2] * result[2][2] );
  transform[0][0] = result[2][2] / r;
  transform[2][2] = result[2][2] / r;
  transform[0][2] = -result[0][2] / r;
  transform[2][0] = result[0][2] / r;
  transform[1][1] = 1.0;
  transform[0][1] = 0.0;
  transform[1][0] = 0.0;
  transform[1][2] = 0.0;
  transform[2][1] = 0.0;
  VOP_3x3_TIMES_3x3(points, transform, result);

  /* rotate vector one into XZ */
  r = sqrt( points[0][0] * points[0][0] + points[1][0] * points[1][0] );
  transform[0][0] = points[0][0] / r;
  transform[1][1] = points[0][0] / r;
  transform[0][1] = points[1][0] / r;
  transform[1][0] = -points[1][0] / r;
  transform[2][2] = 1.0;
  transform[0][2] = 0.0;
  transform[1][2] = 0.0;
  transform[2][0] = 0.0;
  transform[2][1] = 0.0;
  VOP_3x3_TIMES_3x3(result, transform, points);

  /* rotate vector one into X */
  r = sqrt( result[0][0] * result[0][0] + result[0][2] * result[0][2] );
  transform[0][0] = result[0][0] / r;
  transform[2][2] = result[0][0] / r;
  transform[2][0] = result[0][2] / r;
  transform[0][2] = -result[0][2] / r;
  transform[1][1] = 1.0;
  transform[0][1] = 0.0;
  transform[1][0] = 0.0;
  transform[1][2] = 0.0;
  transform[2][1] = 0.0;
  VOP_3x3_TIMES_3x3(points, transform, result);

  /* has Y changed sign? */
  if ( points[1][1] < 0 ) {
    ev[0][1] = -ev[0][1];
    ev[1][1] = -ev[1][1];
    ev[2][1] = -ev[2][1];
    return 1;
  }
  return 0;
}



   double *
calculatePrincipalAxis(ptrajState *state, int *mask,
		       double *x, double *y, double *z,
		       int com, int doRotation, int doReturn)
{
  double cx, cy, cz;
  double *amass, total_mass;
  double Ixx, Iyy, Izz;
  double Ixy, Iyz, Ixz;
  double xtemp, ytemp, ztemp;
  double inertia[3][3];
  double evalue[3], evector[3][3];
  double *returnValue;
  int i, j, activeAtoms;
  int i1, i2, i3;

  /*
   *  allocate return value
   */
  if (doReturn)
    returnValue = (double *) safe_malloc(sizeof(double) * 12);
  else
    returnValue = NULL;

  /*
   *  determine the total mass and the center of mass/geometry 
   *  for the selected atoms
   */
  amass = safe_malloc(sizeof(double) * state->atoms);
  /* amass is zeroed */

  activeAtoms = 0;
  total_mass = 0.0;
  cx = 0.0;
  cy = 0.0;
  cz = 0.0;
  for (i=0; i < state->atoms; i++) {
    if ( mask && mask[i] ) {
      activeAtoms++;
      amass[i] = (com ? state->masses[i] : 1.0);
      total_mass += amass[i];
      cx += amass[i] * x[i];
      cy += amass[i] * y[i];
      cz += amass[i] * z[i];
    }
  }
  cx = cx / total_mass;
  cy = cy / total_mass;
  cz = cz / total_mass;

  /*
   *  calculate the moments of inertia and products of
   *  inertia
   */
  Ixx = 0.0; Iyy = 0.0; Izz = 0.0;
  Ixy = 0.0; Ixz = 0.0; Iyz = 0.0;

  /*
  fprintf(stdout, "IN CALCULATE PRINCIPAL AXIS TOTAL NUMBER OF ACTIVE IS %i\n",activeAtoms);
  fprintf(stdout, "%8.3f %8.3f %8.3f      %8.3f %8.3f %8.3f\n", x[0], y[0], z[0], cx, cy, cz);
  fprintf(stdout, "%8.3f %8.3f %8.3f\n", x[1], y[1], z[1]);
  fprintf(stdout, "%8.3f %8.3f %8.3f\n", x[2], y[2], z[2]);
  fprintf(stdout, "%8.3f %8.3f %8.3f\n", x[3], y[3], z[3]);
  fprintf(stdout, "%8.3f %8.3f %8.3f\n", x[4], y[4], z[4]);
  fprintf(stdout, "%8.3f %8.3f %8.3f\n", x[5], y[5], z[5]);
  fprintf(stdout, "---------------\n");
  */

  for (j=0; j < state->atoms; j++) {
    /*
     *  The formulas for the moments and products of
     *  inertia are:
     *
     *  Ixx =  SUM (amass[j] * ((y[j] - cy) * (y[j] - cy) + 
     *                          (z[j] - cz) * (z[j] - cz)));
     *
     *  Ixy = -SUM (amass[j] * (x[j] - cx) * (y[j] - cy));
     */
    
    xtemp = x[j] - cx;
    ytemp = y[j] - cy;
    ztemp = z[j] - cz;

    Ixx += amass[j] * ( ytemp * ytemp + ztemp * ztemp );
    Iyy += amass[j] * ( xtemp * xtemp + ztemp * ztemp );
    Izz += amass[j] * ( xtemp * xtemp + ytemp * ytemp );
    Ixy -= amass[j] * xtemp * ytemp;
    Iyz -= amass[j] * ytemp * ztemp;
    Ixz -= amass[j] * xtemp * ztemp;
  }

  inertia[0][0] = Ixx;
  inertia[0][1] = Ixy;
  inertia[0][2] = Ixz;
  
  inertia[1][0] = Ixy;
  inertia[1][1] = Iyy;
  inertia[1][2] = Iyz;
  
  inertia[2][0] = Ixz;
  inertia[2][1] = Iyz;
  inertia[2][2] = Izz;
  
  jacobi(inertia, 3, evalue, evector);


  /*
   *  reorder according to the absolute value of the 
   *  eigenvalues; the maximal one comes first...
   */

  i1 = 0; i2 = 1; i3 = 2;

  if (evalue[0] < 0) evalue[0] = -evalue[0];
  if (evalue[1] < 0) evalue[1] = -evalue[1];
  if (evalue[2] < 0) evalue[2] = -evalue[2];

  if (evalue[0] > evalue[1] && 
      evalue[0] > evalue[2]) {
    if (evalue[1] > evalue[2]) {
      i1 = 0; i2 = 1; i3 = 2;
    } else {
      i1 = 0; i2 = 2; i3 = 1;
    }
  } else if (evalue[1] > evalue[0] && 
	     evalue[1] > evalue[2]) {
    if (evalue[0] > evalue[2]) {
      i1 = 1; i2 = 0; i3 = 2;
    } else {
      i1 = 1; i2 = 2; i3 = 0;
    }
  } else if (evalue[0] > evalue[1]) {
    i1 = 2; i2 = 0; i3 = 1;
  } else {
    i1 = 2; i2 = 1; i3 = 0;
  }

  /*
   *  swap around eigenvectors
   */


  if (i1 != 0 || i2 != 1) {
    
    if (prnlev > 2) {

      fprintf(stdout, "PRINCIPAL, EIGENVECTORS/VALUES ARE SWAPPED: %i %i %i\n", 
	      i1, i2, i3);

    }
  }

  inertia[0][0] = evector[0][i1];
  inertia[0][1] = evector[0][i2];
  inertia[0][2] = evector[0][i3];
  inertia[1][0] = evector[1][i1];
  inertia[1][1] = evector[1][i2];
  inertia[1][2] = evector[1][i3];
  inertia[2][0] = evector[2][i1];
  inertia[2][1] = evector[2][i2];
  inertia[2][2] = evector[2][i3];

  evector[0][0] = inertia[0][0];
  evector[0][1] = inertia[0][1];
  evector[0][2] = inertia[0][2];
  evector[1][0] = inertia[1][0];
  evector[1][1] = inertia[1][1];
  evector[1][2] = inertia[1][2];
  evector[2][0] = inertia[2][0];
  evector[2][1] = inertia[2][1];
  evector[2][2] = inertia[2][2];

  /*
   *  invert eigenvalue signs post swap to avoid chiral inversion
   */
  if (i1 == 0 && i2 == 2 && i3 == 1) {
    evector[0][1] = -evector[0][1];
    evector[1][1] = -evector[1][1];
    evector[2][1] = -evector[2][1];
  } else if (i1 == 2 && i2 == 0 && i3 == 1) {
    evector[0][0] = -evector[0][0];
    evector[1][0] = -evector[1][0];
    evector[2][0] = -evector[2][0];
    evector[0][1] = -evector[0][1];
    evector[1][1] = -evector[1][1];
    evector[2][1] = -evector[2][1];
    evector[0][2] = -evector[0][2];
    evector[1][2] = -evector[1][2];
    evector[2][2] = -evector[2][2];
  }


  /*
   *  swap eigenvalues
   */


  inertia[0][0] = evalue[i1];
  inertia[0][1] = evalue[i2];
  inertia[0][2] = evalue[i3];

  evalue[0] = inertia[0][0];
  evalue[1] = inertia[0][1];
  evalue[2] = inertia[0][2];



  if (prnlev > 2) {
    fprintf(stdout, "\nJACOBI\n\n");
    fprintf(stdout, "Ixx = %10.3f Iyy = %10.3f Izz = %10.3f\n", Ixx, Iyy, Izz);
    fprintf(stdout, "Ixy = %10.3f Ixz = %10.3f Iyz = %10.3f\n\n", Ixy, Ixz, Iyz);

    fprintf(stdout, "ATOM      1 XXX  XX      1    %8.3f%8.3f%8.3f\n",
	    cx, cy, cz);
    fprintf(stdout, "ATOM      2 XXX  XX      1    %8.3f%8.3f%8.3f\n",
	    cx+evector[0][0], cy+evector[0][1], cz+evector[0][2]);

    fprintf(stdout, "ATOM      1 YYY  YY      2    %8.3f%8.3f%8.3f\n",
	    cx, cy, cz);
    fprintf(stdout, "ATOM      2 YYY  YY      2    %8.3f%8.3f%8.3f\n",
	    cx+evector[1][0], cy+evector[1][1], cz+evector[1][2]);

    fprintf(stdout, "ATOM      1 ZZZ  ZZ      3    %8.3f%8.3f%8.3f\n",
	    cx, cy, cz);
    fprintf(stdout, "ATOM      2 ZZZ  ZZ      3    %8.3f%8.3f%8.3f\n",
	    cx+evector[2][0], cy+evector[2][1], cz+evector[2][2]);

    fprintf(stdout,"\n");

    fprintf(stdout, "          Eigenvalues %10.3f %10.3f %10.3f\n", 
	    evalue[0], evalue[1], evalue[2]);
  }
 

  /*
   *  check for chiral inversion!
   */

  if ( jacobiCheckChirality(evalue, evector) == 1 && prnlev > 0 ) {
    fprintf(stdout, "\nPRINCIPAL: WARNING!!! CHECK CHIRALITY: vectors swapped!\n");
  }



  /* 
   *  Perform rotation if requested, i.e. multiply the transpose of the
   *  evector matrix by the coords...
   */
  if ( doRotation ) {

    for (j=0; j < state->atoms; j++) {


     VOP_3x3_TRANSPOSE_TIMES_COORDS(evector, x[j], y[j], z[j], 
				     xtemp, ytemp, ztemp);

     /* alternate code to avoid macro
      xtemp = evector[0][0] * x[j] + evector[1][0] * y[j] + evector[2][0] * z[j];  
      ytemp = evector[0][1] * x[j] + evector[1][1] * y[j] + evector[2][1] * z[j];  
      ztemp = evector[0][2] * x[j] + evector[1][2] * y[j] + evector[2][2] * z[j];  
      x[j] = xtemp; 
      y[j] = ytemp;
      z[j] = ztemp;
     */

    }
  }

  safe_free(amass);

  if (doReturn) { 

    /*
     *  order according to the absolute value of the 
     *  eigenvalues; the maximal one comes first...
     */
    returnValue[0] = evector[i1][0];
    returnValue[1] = evector[i1][1];
    returnValue[2] = evector[i1][2];
    returnValue[3] = evector[i2][0];
    returnValue[4] = evector[i2][1];
    returnValue[5] = evector[i2][2];
    returnValue[6] = evector[i3][0];
    returnValue[7] = evector[i3][1];
    returnValue[8] = evector[i3][2];
    returnValue[9] = cx;
    returnValue[10] = cy;
    returnValue[11] = cz;
  } 
  return( returnValue );

}


   int
transformPrincipal(actionInformation *action, 
		   double *x, double *y, double *z, 
		   double *box, int mode)
{
  //char *name = "principal";
  argStackType **argumentStackPointer;
  char *buffer;

  /*
   *  USAGE:
   *
   *    principal mask [dorotation] [mass]
   *
   *  action argument usage:
   *
   *  mask -- atoms for which to calculate the principal axis
   *  iarg1 -- = 1 if coordinates will be modified, 0 by default
   */

  if (mode == PTRAJ_SETUP) {

    /*
     *  ACTION: PTRAJ_SETUP
     */

#ifdef MPI

#endif

    argumentStackPointer = (argStackType **) action->carg1;
    action->carg1 = NULL;

    buffer = getArgumentString(argumentStackPointer, NULL);
    action->mask = processAtomMask(buffer, action->state);
    safe_free(buffer);

    action->iarg1 = argumentStackContains(argumentStackPointer, "dorotation");
    action->iarg2 = argumentStackContains(argumentStackPointer, "mass");
    

  } else if (mode == PTRAJ_STATUS) {

    fprintf(stdout, "  PRINCIPAL %s rotation by center of %s, atom selection is: ", 
	    (action->iarg1 ? "with" : "without"),
	    (action->iarg2 ? "mass" : "geometry"));
    printAtomMask(stdout, action->mask, action->state);
    fprintf(stdout, "\n");

  }

  if (mode != PTRAJ_ACTION) return 0;


  calculatePrincipalAxis(action->state, action->mask,
			 x, y, z, action->iarg2, action->iarg1, 0);

  return 1;

}


/** ACTION ROUTINE *************************************************************
 *
 *  transformProjection() --- project snapshots on normal modes
 *
 *  Supplementary routines:
 *    freeTransformProjectionMemory (below)
 *    readEvecFile (in evec.h/.c)
 *
 ******************************************************************************/
typedef struct _transformProjectionInfo {
  char *outfile;
  FILE *fp;
  int ibeg;
  int iend;
  int start;
  int stop;
  int offset;
  double *sqrtmasses;
} transformProjectionInfo;

#define INITIALIZE_transformProjectionInfo(_p_) \
  _p_->outfile    = NULL; \
  _p_->fp         = NULL; \
  _p_->ibeg       = 1;    \
  _p_->iend       = 2;    \
  _p_->start      = 1;    \
  _p_->stop       = -1;   \
  _p_->offset     = 1;    \
  _p_->sqrtmasses = NULL;

// freeTransformProjectionMemory()
static void freeTransformProjectionMemory(actionInformation *action) { 
  transformProjectionInfo *pinfo;
  modesInfo *modinfo;

  pinfo = (transformProjectionInfo *) action->carg1;
  if(pinfo != NULL){
    if(pinfo->outfile != NULL)
      safe_free(pinfo->outfile);
    if(pinfo->sqrtmasses != NULL)
      safe_free(pinfo->sqrtmasses);
    INITIALIZE_transformProjectionInfo(pinfo);
    safe_free(pinfo);
  }

  modinfo = (modesInfo *) action->carg2;
  if(modinfo != NULL){
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
}

// transformProjection()
int transformProjection(actionInformation *action, 
   	            double *x, double *y, double *z,
		    double *box, int mode)
{
  //char *name = "projection";
  argStackType **argumentStackPointer;
  ptrajState *state;
  char *buffer;

  transformProjectionInfo *pinfo;
  modesInfo *modinfo;
  modesType type;
  FILE *fp;
  int i, j, nvect, natoms, cnt, ind1, ind2;
  int start, stop, offset, ibeg, iend, masktot;
  int *mask;
  double sqrtmass, proj, proj1, proj2, proj3;
  double *sqrtmasses, *avg, *evec; 

  /*
   *  USAGE:
   *
   *  projection
   *             modes <modesfile> out <outfile>
   *             [beg <beg>] [end <end>] [<mask>]
   *             [start <start>] [stop <stop>] [offset <offset>]
   *  
   *  action argument usage:
   *    carg1:
   *      pointer to transformProjectionInfo
   *    carg2:
   *      pointer to modesInfo
   *
   *    iarg1:
   *      counter of snapshots
   *
   */

  if (mode == PTRAJ_SETUP) {

    /*
     *  -------- ACTION: PTRAJ_SETUP
     */

#ifdef MPI
    printParallelError(name);
    return -1;
#endif

    argumentStackPointer = (argStackType **) action->carg1;
    action->carg1 = NULL;

    /*
     *  Alloc / init transformProjectionInfo
     */ 
    pinfo = (transformProjectionInfo *) safe_malloc(sizeof(transformProjectionInfo));
    INITIALIZE_transformProjectionInfo(pinfo);
    action->carg1 = (void *) pinfo;

    /*
     * Get ibeg, iend, start, stop, offset
     */
    pinfo->ibeg = ibeg = argumentStackKeyToInteger(argumentStackPointer, "beg", 1);
    pinfo->iend = iend = argumentStackKeyToInteger(argumentStackPointer, "end", 2);
    pinfo->start       = argumentStackKeyToInteger(argumentStackPointer, "start", 1);
    pinfo->stop        = argumentStackKeyToInteger(argumentStackPointer, "stop", -1);
    pinfo->offset      = argumentStackKeyToInteger(argumentStackPointer, "offset", 1);

    /*
     *  Get modes
     */
    buffer = argumentStackKeyToString(argumentStackPointer, "modes", NULL);
    if(buffer == NULL){
      fprintf(stdout,
              "WARNING in ptraj(), transformProjection: no modesfile given, ignoring command\n");
      freeTransformProjectionMemory(action);
      return -1;
    }
    else{
      /*
       *  Allocate modesInfo structure
       */
      modinfo = (modesInfo *) safe_malloc(sizeof(modesInfo));
      action->carg2 = (void *) modinfo;
      INITIALIZE_modesInfo(modinfo);
      modinfo->name = buffer;
      modinfo->type = MT_UNKNOWN;
      modinfo->source = MS_FILE;

      /*
       *  Read evec file
       */
      fp = safe_fopen(buffer, "r");
      if(fp == NULL){
        fprintf(stdout,
                "WARNING in ptraj(), transformProjection: file %s not opened, ignoring command\n", buffer);
        freeTransformProjectionMemory(action);
        return -1;
      }
      if(readEvecFile(fp, ibeg, iend, modinfo)){
        fprintf(stdout,
                "WARNING in ptraj(), transformProjection: error while reading %s, ignoring command\n", buffer);
        freeTransformProjectionMemory(action);
        return -1;
      }
      if(modinfo->nvect != (iend - ibeg + 1)){
        fprintf(stdout,
                "FYI: Number of read evecs is %i, number of requested evecs is %i\n", 
                modinfo->nvect, iend - ibeg + 1);
      }
      safe_fclose(fp);

      if(modinfo->type != MT_COVAR && 
         modinfo->type != MT_MWCOVAR &&
         modinfo->type != MT_IDEA){
        fprintf(stdout,
             "WARNING in ptraj(), transformProjection: evecs not of type COVAR, MWCOVAR, or IDEA, ignoring command\n");
        freeTransformProjectionMemory(action);
        return -1;
      }
    }

    /*
     *  Get outfile
     */
    pinfo->outfile = buffer = argumentStackKeyToString(argumentStackPointer, "out", NULL);
    if(buffer == NULL){
      fprintf(stdout,
              "WARNING in ptraj(), transformProjection: no outfile given, ignoring command\n");
      freeTransformProjectionMemory(action);
      return -1;
    }
    else{
      /*
       * Open outfile and store fp in pinfo
       */
      fp = safe_fopen(buffer, "w");
      if(fp == NULL){
        fprintf(stdout,
                "WARNING in ptraj(), transformProjection: file %s not opened, ignoring command\n", buffer);
        freeTransformProjectionMemory(action);
        return -1;
      }
      else{
        pinfo->fp = fp;

        /* 
         * Write header line
         */
        fprintf(fp, "Projection of snapshots onto modes\n");
        fprintf(fp, "%10s", "Snapshot");
        for(i = ibeg; i < modinfo->nvect+ibeg; i++){
          if(i < 10)
            fprintf(fp, "     Mode%i", i);
          else if(i < 100)
            fprintf(fp, "    Mode%i", i);
          else if(i < 1000)
            fprintf(fp, "   Mode%i", i);
          else if(i < 10000)
            fprintf(fp, "  Mode%i", i);
          else if(i < 100000)
            fprintf(fp, " Mode%i", i);

          if(modinfo->type == MT_IDEA)
            fprintf(fp, "                              ");
        }
        fprintf(fp, "\n");
      }
    }

    /*
     *  Get mask
     */
    buffer = getArgumentString(argumentStackPointer, NULL);
    if(buffer == NULL){
      action->mask = processAtomMask( (char *) "*", action->state);
    }else{
      action->mask = processAtomMask(buffer, action->state);
      safe_free(buffer);
    }
    
    masktot = 0;
    for(i=0; i < action->state->atoms; i++)
      if(action->mask[i])
        masktot++;

    if(modinfo->type == MT_COVAR ||
       modinfo->type == MT_MWCOVAR){
      /*
       *  Check if (3 * number of atoms in mask) and nvectelem agree
       */ 
      if(3 * masktot != modinfo->navgelem ||
         3 * masktot != modinfo->nvectelem){
         fprintf(stdout,
           "WARNING in ptraj(), transformProjection: no. of atom coords does not agree with no. of vect elementsignoring command\n");
         freeTransformProjectionMemory(action);
         return -1;
      }
    }
    else if(modinfo->type == MT_IDEA){
      /*
       *  Check if (number of atoms in mask) and nvectelem agree
       */ 
      if(masktot != modinfo->navgelem ||
         masktot != modinfo->nvectelem){
         fprintf(stdout,
           "WARNING in ptraj(), transformProjection: no. of atom coords does not agree with no. of vect elementsignoring command\n");
         freeTransformProjectionMemory(action);
         return -1;
      }
    }

    /*
     *  Precalc sqrt of mass for each coordinate
     */
    sqrtmasses = (double *) safe_malloc(sizeof(double) * masktot);
    pinfo->sqrtmasses = sqrtmasses;
    cnt = 0;
    for(i=0; i < action->state->atoms; i++){
      if(action->mask[i]){
        if(modinfo->type == MT_MWCOVAR)
          sqrtmasses[cnt++] = sqrt(action->state->masses[i]);
        else /* MT_COVAR - no mass-weighting necessary */
          sqrtmasses[cnt++] = 1.0;
      }
    }

    /*
     *  Init variables
     */ 
    action->iarg1 = 0;
    
    return 0;

  } else if (mode == PTRAJ_STATUS) {

    /*
     *  -------- ACTION: PTRAJ_STATUS
     */

    pinfo = (transformProjectionInfo *) action->carg1;
    modinfo = (modesInfo *) action->carg2;

    fprintf(stdout, "  PROJECTION: Calculating projection using modes %i to %i of file %s\n",
            pinfo->ibeg,
            pinfo->iend,
	    modinfo->name);
    fprintf(stdout, "                  Results are written to %s\n", pinfo->outfile);
    if (pinfo->start != 1 || pinfo->stop != -1 || pinfo->offset != 1) {
      fprintf(stdout, "                 Start: %i", pinfo->start);
      if (pinfo->stop > 0)
	fprintf(stdout, "   Stop: %i", pinfo->stop);
      else
	fprintf(stdout, "   Stop: at final frame");
      fprintf(stdout, "   Offset: %i\n", pinfo->offset);
    }
    fprintf(stdout, "                  Atom selection follows ");
    printAtomMask(stdout, action->mask, action->state);
    fprintf(stdout, "\n");

  } else if (mode == PTRAJ_PRINT) {

    /*
     *  -------- ACTION: PTRAJ_PRINT
     */

    /*
     *  Nothing to do here - output happens in PTRAJ_ACTION, to avoid allocating results vector
     */
  }
  else if (mode == PTRAJ_CLEANUP) {

    /*
     *  -------- ACTION: PTRAJ_CLEANUP
     */

    safe_fclose(((transformProjectionInfo *) action->carg1)->fp);
    freeTransformProjectionMemory(action);

  }

  if (mode != PTRAJ_ACTION) return 0;

  /*
   *  -------- ACTION: PTRAJ_ACTION
   */

  state      = action->state;
  natoms     = state->atoms;
  mask       = action->mask;

  pinfo      = (transformProjectionInfo *) action->carg1;
  start      = pinfo->start;
  stop       = pinfo->stop;
  offset     = pinfo->offset;
  sqrtmasses = pinfo->sqrtmasses;
  fp         = pinfo->fp;

  modinfo    = (modesInfo *) action->carg2;
  type       = modinfo->type;
  avg        = modinfo->avg;
  nvect      = modinfo->nvect;
  evec       = modinfo->evec;


  action->iarg1++;
  if(action->iarg1 >= start && 
     (stop < 0 || action->iarg1 <= stop) &&
     (offset == 1 || (action->iarg1 - start) % offset == 0)){

    fprintf(fp, "%10i", action->iarg1);

    /*
     * Project snapshots on modes
     */
    if(type == MT_COVAR ||
       type == MT_MWCOVAR){
      ind1 = 0;
      for(i = 0; i < nvect; i++){
        proj = 0.0;
        ind2 = 0;
        cnt  = 0;
        for(j = 0; j < natoms; j++){
          if(mask[j]){
            sqrtmass = sqrtmasses[cnt];
            proj += (x[j] - avg[ind2  ]) * sqrtmass * evec[ind1  ];
            proj += (y[j] - avg[ind2+1]) * sqrtmass * evec[ind1+1];
            proj += (z[j] - avg[ind2+2]) * sqrtmass * evec[ind1+2];
            ind1 += 3;
            ind2 += 3;
            cnt++;
          }
        }
        /* Output proj */
        fprintf(fp, " %9.3f", proj);
      }
      /* Newline */
      fprintf(fp, "\n");
    }
    else if(type == MT_IDEA){
      ind1 = 0;
      for(i = 0; i < nvect; i++){
        proj  = 0.0;
        proj1 = 0.0;
        proj2 = 0.0;
        proj3 = 0.0;
        for(j = 0; j < natoms; j++){
          if(mask[j]){
            proj1 += x[j] * evec[ind1];
            proj2 += y[j] * evec[ind1];
            proj3 += z[j] * evec[ind1];
            ind1++;
          }
        }
        /* Output proj */
        fprintf(fp, " %9.3f %9.3f %9.3f %9.3f", 
                    proj1, proj2, proj3, 
                    sqrt(proj1*proj1 + proj2*proj2 + proj3*proj3));
      }
      /* Newline */
      fprintf(fp, "\n");
    }
  }   

  return 1;
}



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

/** ACTION ROUTINE *************************************************************
 *
 *  transformUnwrap() --- opposite of image 
 *
 ******************************************************************************/
typedef struct _transformUnwrapInfo {
  double *refx;
  double *refy;
  double *refz;
} transformUnwrapInfo;
// transformUnwrap()
int transformUnwrap( actionInformation *action,
                 double *x, double *y, double *z,
                 double *box, int mode )
{
  //char *name = "unwrap";
  argStackType **argumentStackPointer;
  transformUnwrapInfo *unwrapInfo; 
  char *buffer;
  ptrajState *state;
  coordinateInfo *refInfo;
  int *mask;
  int orthog;
  int natoms;
  int i; 
  double *rx, *ry, *rz;
  double ucell[9], recip[9];
  double dx, dy, dz;
  double cx, cy, cz;
  double ccx, ccy, ccz;
  //double boxXtrans, boxYtrans, boxZtrans;
  int ix, iy, iz;
  double newX, newY, newZ;
  double minX, minY, minZ;
  double distanceSquare, minDistanceSquare;
  /*
   *  USAGE:
   *
   *    unwrap [reference] [mask]
   *
   *  action argument usage:
   *
   *  iarg1:
   *     0 -- first frame is not unwrapped
   *     1 -- first frame is unwrapped wrt the reference
   *  carg1:
   *     the transformUnwrapInfo structure
   *  carg2:
   *     the current reference structure
   *
   */

  if ( mode == PTRAJ_SETUP ) {
     /*
      *  ACTION: PTRAJ_SETUP
      */

#ifdef MPI
    printParallelError(name);
    return -1;
#endif

    argumentStackPointer = (argStackType **) action->carg1;
    action->carg1 = NULL;

    unwrapInfo = (transformUnwrapInfo *) safe_malloc( sizeof(transformUnwrapInfo) );
    // Initialize unwrapInfo
    unwrapInfo->refx = NULL;
    unwrapInfo->refy = NULL;
    unwrapInfo->refz = NULL;

    if ( argumentStackContains( argumentStackPointer, "reference" ) ) {
      action->iarg1 = 1;
      if ( referenceInfo == NULL ) {
        fprintf( stdout, "WARNING in ptraj(), unwrap to reference: missing reference structure.\n" );
        fprintf( stdout, "Set this prior to this unwrap with the command \"reference\"\n" );
        fprintf( stdout, "Ignoring reference...\n" );
        action->iarg1 = 0;
      }
    }
    else {
      action->iarg1 = 0;
    }

    if ( action->iarg1 ) {
      natoms = action->state->atoms;
      rx = (double *) safe_malloc( sizeof(double) * natoms );
      ry = (double *) safe_malloc( sizeof(double) * natoms );
      rz = (double *) safe_malloc( sizeof(double) * natoms );
      for ( i = 0; i < natoms; i++ ) {
        rx[i] = referenceInfo->x[i];
        ry[i] = referenceInfo->y[i];
        rz[i] = referenceInfo->z[i];
      }
      unwrapInfo->refx = rx;
      unwrapInfo->refy = ry;
      unwrapInfo->refz = rz;
    }

    buffer = getArgumentString(argumentStackPointer, NULL);
    if ( buffer == NULL ) {
      action->mask = processAtomMask( "*", action->state );
    }
    else {
      action->mask = processAtomMask( buffer, action->state );
      safe_free(buffer);
    }

    action->carg1 = (void *) unwrapInfo;
    action->carg2 = (void *) referenceInfo;
  }

  else if ( mode == PTRAJ_STATUS ) {
     /*
      *  ACTION: PTRAJ_STATUS
      */

    refInfo = (coordinateInfo *) action->carg2;

    fprintf( stdout, "  UNWRAP\n" );
    if ( action->iarg1 == 1 ) {
      fprintf( stdout, "      First frame is unwrapped to reference structure (%s).\n",
               refInfo->filename );
    }
    fprintf( stdout, "      The atoms in the calculation follow: ");
    printAtomMask( stdout, action->mask, action->state );
    fprintf( stdout, "\n" );
  }

  else if ( mode == PTRAJ_CLEANUP ) {
     /*
      *  ACTION: PTRAJ_CLEANUP
      */
    unwrapInfo = (transformUnwrapInfo *) action->carg1;
    safe_free(unwrapInfo->refx);
    safe_free(unwrapInfo->refy);
    safe_free(unwrapInfo->refz);
    unwrapInfo->refx = NULL;
    unwrapInfo->refy = NULL;
    unwrapInfo->refz = NULL;
    safe_free(unwrapInfo);
  }

  if ( mode != PTRAJ_ACTION ) return 0;

  /*
   *  ACTION: PTRAJ_ACTION
   */

  state = (ptrajState *) action->state;
  natoms = state->atoms;

  /*
   *  update local state information
   */
  for (i=0; i<6; i++)
    state->box[i] = box[i];

  if ( prnlev > 4 ) {
    fprintf( stdout, "  UNWRAP: box size is %8.3f %8.3f %8.3f  %f %f %f\n",
             box[0], box[1], box[2], box[3], box[4], box[5] );
  }

  /*
   *  process arguments
   */
  mask = action->mask;
  if ( mask == NULL ) return 0;

  unwrapInfo = (transformUnwrapInfo *) action->carg1;

  if ( box[3] != 90.0 || box[4] != 90.0 || box[5] != 90.0 )
    orthog = 0;
  else
    orthog = 1;

  if ( unwrapInfo->refx == NULL || unwrapInfo->refy == NULL || unwrapInfo->refz == NULL ) {
    // previous coordinates do not exist
    rx = (double *) safe_malloc( sizeof(double) * natoms );
    ry = (double *) safe_malloc( sizeof(double) * natoms );
    rz = (double *) safe_malloc( sizeof(double) * natoms );
    for ( i = 0; i < natoms; i++ ) {
      rx[i] = x[i];
      ry[i] = y[i];
      rz[i] = z[i];
    }
    unwrapInfo->refx = rx;
    unwrapInfo->refy = ry;
    unwrapInfo->refz = rz;

    return 0;
  }

  if ( box[0] <= 0.0 || box[1] <= 0.0 || box[2] <= 0.0 || box[3] <= 0.0 ) {
    if (prnlev > 1)
      fprintf(stdout, "  UNWRAP: box size is out of range (%8.2f %8.2f %8.2f, alpha = %8.2f), returning\n", box[0], box[1], box[2], box[3]);
    return 0;
  }

  rx = unwrapInfo->refx;
  ry = unwrapInfo->refy;
  rz = unwrapInfo->refz;

  if ( orthog ) {
    for ( i = 0; i < natoms; i++ ) {
      if ( ! mask[i] ) continue; 

      dx = x[i] - rx[i];
      dy = y[i] - ry[i];
      dz = z[i] - rz[i];

      rx[i] = x[i] = x[i] - floor( dx / box[0] + 0.5 ) * box[0];
      ry[i] = y[i] = y[i] - floor( dy / box[1] + 0.5 ) * box[1];
      rz[i] = z[i] = z[i] - floor( dz / box[2] + 0.5 ) * box[2];
    }
  }
  else {
    if (prnlev > 0 && box[3] <= 0.0) {
      fprintf( stdout, "  UNWRAP: Warning, box angles are <= 0.0 (%8.2f %8.2f %8.2f)\n",
               box[3], box[4], box[5] );
    }
    boxToRecip( box, ucell, recip );

    if (prnlev > 4) {
      printf("UNWRAP TRICLINIC\n");
      printf("BOX is   %f %f %f %f %f %f\n", box[0],box[1],box[2],box[3],box[4],box[5]);
      printf("UCELL is %f %f %f\n", ucell[0], ucell[1], ucell[2]);
      printf("         %f %f %f\n", ucell[3], ucell[4], ucell[5]);
      printf("         %f %f %f\n", ucell[6], ucell[7], ucell[8]);

      printf("RECIP is %f %f %f\n", recip[0], recip[1], recip[2]);
      printf("RECIP is %f %f %f\n", recip[3], recip[4], recip[5]);
      printf("RECIP is %f %f %f\n", recip[6], recip[7], recip[8]);
    }

    for ( i = 0; i < natoms; i++ ) {
      if ( ! mask[i] ) continue;

      dx = x[i] - rx[i];
      dy = y[i] - ry[i];
      dz = z[i] - rz[i];

      cx=floor( dx*recip[0] + dy*recip[1] + dz*recip[2] );
      cy=floor( dx*recip[3] + dy*recip[4] + dz*recip[5] );
      cz=floor( dx*recip[6] + dy*recip[7] + dz*recip[8] );

      minDistanceSquare = pow(dx,2) + pow(dy,2) + pow(dz,2);
      minX = x[i];
      minY = y[i];
      minZ = z[i];

      for ( ix = -1; ix <= 1; ix++ ) {
        for ( iy = -1; iy <= 1; iy++ ) {
          for ( iz = -1; iz <= 1; iz++ ) {
            ccx = cx + (double) ix;
            ccy = cy + (double) iy;
            ccz = cz + (double) iz;
     
            newX = x[i] - ( ccx * ucell[0] + ccy * ucell[3] + ccz * ucell[6] );
            newY = y[i] - ( ccx * ucell[1] + ccy * ucell[4] + ccz * ucell[7] );
            newZ = z[i] - ( ccx * ucell[2] + ccy * ucell[5] + ccz * ucell[8] );

            distanceSquare =   pow( rx[i] - newX, 2 )
                             + pow( ry[i] - newY, 2 )
                             + pow( rz[i] - newZ, 2 );

            if ( minDistanceSquare > distanceSquare ) {
                minDistanceSquare = distanceSquare;
                minX = newX;
                minY = newY;
                minZ = newZ;
            }
          }
        }
      }

      rx[i] = x[i] = minX;
      ry[i] = y[i] = minY;
      rz[i] = z[i] = minZ;
    }
  }

  return 0;
}

/*
int nint(double x){
    int i;
    i = (x >= 0) ? (int) (x+.5) : (int) (x-.5);
    return (i);
}
*/
    
void cross(double *x, double *y, double *z){
  z[0]=x[1]*y[2]-x[2]*y[1];
  z[1]=-x[0]*y[2]+x[2]*y[0];
  z[2]=x[0]*y[1]-x[1]*y[0];
  return;
}
double dot(double *x, double *y){
  return(x[0]*y[0]+x[1]*y[1]+x[2]*y[2]);
}



/** ACTION ROUTINE *************************************************************
 *
 *  transformVector()   --- compute/store various vector quantities
 *
 *  Supplementary routines:
 *    freeTransformVectorInfo (below)
 *    spherharm (below)
 *    lsqplane (below)
 *    solvecubiceq (below)
 *    cmpdouble (below)
 *
 ******************************************************************************/

   int
cmpdouble(const void *v1, const void *v2){

  double d1, d2;
  d1 = *((double *) v1);
  d2 = *((double *) v2);

  if(d1 < d2)
    return -1;
  else if(d1 > d2)
    return 1;
  else
    return 0;
}

   void
solve_cubic_eq(double a, double b, double c, double d, double *droot){

  /* 
   * Solves a cubic equation
   * ax^3 + bx^2 + cx + d = 0
   * using "Cardan's formula"
   * (see: Bronstein, S.131f)
   */

  //const double PI = 3.141592654;
  const double one3 = 1.0 / 3.0;
  const double one27 = 1.0 / 27.0;

  double r, s, t;
  double p, q, rho, phi;
  double D, u, v;
  double dtmp[3];

  /* Coeff. for normal form x^3 + rx^2 + sx + t = 0 */
  r = b / a;
  s = c / a;
  t = d / a;
  
  /* Coeff. for red. eq. y^3 + py + q = 0 with y = x + r/3 bzw. (x = y - r/3) */
  p = s - r * r * one3;
  q = 2.0 * r * r * r * one27 - r * s * one3 + t;

  /* Dummy variables */
  rho = sqrt(-p * p * p * one27);
  phi = acos(-q / (2.0 * rho));

  /* Discriminante(?) */
  D = pow((p * one3),3) + q * q * 0.25;
  
  if(D > 0){ /* x real -> one real solution */
    u = pow(-q * 0.5 + sqrt(D), one3);
    v = -p / u * one3;
    *droot = (u + v) - r * one3;
  }
  else if(D <= 0){ /* three real solutions (d < 0) | one real solution + one real double solution or 
                                                     one real triple solution (d = 0) */
    dtmp[0] = 2.0 * pow(rho, one3) * cos(phi * one3) - r * one3;
    dtmp[1] = 2.0 * pow(rho, one3) * cos((phi + 2.0 * PI) * one3) - r * one3;
    dtmp[2] = 2.0 * pow(rho, one3) * cos((phi + 4.0 * PI) * one3) - r * one3;

    qsort((void *) dtmp, (size_t) 3, sizeof(double), cmpdouble);
    *droot = dtmp[0];
  }
}

   void 
lsqplane(int n,
         double *cx, double *cy, double *cz,
         double *a, double *b, double *c){

  /*
   * Calcs (least-squares best) plane through a series of points
   * relative to their center of geom. (the latter has to be done outside this routine), 
   * returns (normalized) coeff. for plane eq. ax + by + cz = 0
   * following: Crystal Structure Analysis for Chem. and Biol.,
   * Glusker, Lewis, Rossi, S. 460ff
   */

  int i;
  double dSumXX, dSumYY, dSumZZ, dSumXY, dSumXZ, dSumYZ;
  double o, p, q, root;
  double dnorm;
  double x1, y1, z1, x2, y2, z2;

  root = 0;

  if(n == 3){
    x1 = cx[1] - cx[0];
    y1 = cy[1] - cy[0];
    z1 = cz[1] - cz[0];
    x2 = cx[2] - cx[1];
    y2 = cy[2] - cy[1];
    z2 = cz[2] - cz[1];

    *a = y1 * z2 - z1 * y2;
    *b = z1 * x2 - x1 * z2;
    *c = x1 * y2 - y1 * x2;
  }
  else{
    /* Calc Var. */
    dSumXX = 0.0;
    dSumYY = 0.0;
    dSumZZ = 0.0;
    dSumXY = 0.0;
    dSumXZ = 0.0;
    dSumYZ = 0.0;
  
    for(i = 0; i < n; i++){
      dSumXX += cx[i] * cx[i];
      dSumYY += cy[i] * cy[i];
      dSumZZ += cz[i] * cz[i];

      dSumXY += cx[i] * cy[i];
      dSumXZ += cx[i] * cz[i];
      dSumYZ += cy[i] * cz[i];    
    }

    /* Calc coeff. for -l^3 + o * l^2 + p * l + q = 0 */
    o = dSumXX + dSumYY + dSumZZ;
    p = pow(dSumXY,2) + pow(dSumXZ,2) + pow(dSumYZ,2) - 
        (dSumXX * dSumYY + dSumXX * dSumZZ + dSumYY * dSumZZ);
    q = dSumXX * dSumYY * dSumZZ + 2.0 * dSumXY * dSumXZ * dSumYZ -
      (dSumXX * dSumYZ * dSumYZ + dSumYY * dSumXZ * dSumXZ + dSumZZ * dSumXY * dSumXY);

    /* Solve cubic eq. */
    solve_cubic_eq(-1.0, o, p, q, &root);

    /* Calc determinantes */
    *a = (dSumYY - root) * dSumXZ - dSumXY * dSumYZ;
    *b = (dSumXX - root) * dSumYZ - dSumXY * dSumXZ;
    *c =  dSumXY         * dSumXY - (dSumYY - root) * (dSumXX - root);

  }

  /* Normalize */
  dnorm = 1.0 / sqrt((*a) * (*a) + (*b) * (*b) + (*c) * (*c));
  *a *= dnorm;
  *b *= dnorm;
  *c *= dnorm;
}

   void
spherharm(int l, int m, double x, double y, double z, double r, 
          double *dreal, double *dimg){

    
  /*
   * Calc spherical harmonics of order l=0,1,2
   * and -l<=m<=l with cartesian coordinates as input
   * (see e.g. Merzbacher, Quantum Mechanics, p. 186)
   */

  const double SH00=0.28209479;
  const double SH10=0.48860251;
  const double SH11=0.34549415;
  const double SH20=0.31539157;
  const double SH21=0.77254840;
  const double SH22=0.38627420;

  double ri;

  *dreal = 0.0;
  *dimg = 0.0;
  ri = 1.0 / r;

  if(l == 0 && m == 0){
    *dreal = SH00;
  }
  else if(l == 1){
    if(m == 0){
      *dreal = SH10 * z * ri;
    }
    else{
      *dreal = -m * SH11 * x * ri;
      *dimg  = -    SH11 * y * ri;
    }
  }
  else if(l == 2){
    if(m == 0){
      *dreal = SH20 * (2.0*z*z - x*x - y*y) * ri * ri;
    }
    else if(fabs(m) == 1){
      *dreal = -m * SH21 * x * z * ri * ri;
      *dimg  = -    SH21 * y * z * ri * ri;
    }
    else{
      *dreal = SH22 * (x*x - y*y) * ri * ri;
      *dimg  = m * SH22 * x * y * ri * ri;
    }
  }
}

   void
freeTransformVectorMemory(actionInformation *action){

  transformVectorInfo *vinfo;
  modesInfo *modinfo;

  vinfo = (transformVectorInfo *) action->carg1;
  if(vinfo != NULL){
    if(vinfo->mode == VECTOR_CORRIRED && vinfo->master == 1){
      modinfo = (modesInfo *) vinfo->modinfo;
      if(modinfo != NULL){
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
    }

    if(vinfo->name != NULL)
      safe_free(vinfo->name);
    if(vinfo->filename != NULL)
      safe_free(vinfo->filename);
    if(vinfo->mask != NULL)
      safe_free(vinfo->mask);
    if(vinfo->mask2 != NULL)
      safe_free(vinfo->mask2);
    if(vinfo->vx != NULL)
      safe_free(vinfo->vx);
    if(vinfo->vy != NULL)
      safe_free(vinfo->vy);
    if(vinfo->vz != NULL)
      safe_free(vinfo->vz);
    if(vinfo->cx != NULL)
      safe_free(vinfo->cx);
    if(vinfo->cy != NULL)
      safe_free(vinfo->cy);
    if(vinfo->cz != NULL)
      safe_free(vinfo->cz);

    if(vinfo->avgcrd != NULL)
      safe_free(vinfo->avgcrd);
    if( (((vinfo->mode == VECTOR_CORRIRED) && (vinfo->master == 1)) ||
         (vinfo->mode != VECTOR_CORRIRED)) &&
        (vinfo->cftmp != NULL) )
      safe_free(vinfo->cftmp);
    if(vinfo->p2cftmp != NULL)
      safe_free(vinfo->p2cftmp);
    if(vinfo->rcftmp != NULL)
      safe_free(vinfo->rcftmp);

    INITIALIZE_transformVectorInfo(vinfo);
    safe_free(vinfo);
  }
}

   int
transformVector(actionInformation *action, 
		double *x, double *y, double *z, 
		double *box, int mode)
{
  //char *name = "vector";
  argStackType **argumentStackPointer;
  stackType *vectorStackTmp = NULL; // NOTE: STACK TYPE CHANGE
  char *buffer;
  ptrajState *state;

  transformVectorInfo *vectorInfo, *vectorInfoTmp;
  vectorMode vmode;
  modesInfo *modinfo;
  int i, j, n, frame;
  int ibeg, iend, npair;
  int indtot, indsnap, indplus, indminus;
  int nvect, nvectelem, order;
  double cx, cy, cz, vx, vy, vz, total_mass;
  double r, r3, r3i;
  double dplusreal, dplusimg, dminusreal, dminusimg, q;
  double *avgcrd, *cftmp, *p2cftmp, *rcftmp; 
  double *principal, *evec;
  FILE *outFile, *fp;

  npair = 0;
  nvect = 0;
  nvectelem = 0;
  avgcrd = NULL;
  p2cftmp = NULL;
  rcftmp = NULL;
  evec = NULL;
  r3i = 0;
  indplus = 0;
  indminus = 0;
  dminusreal = 0;
  dminusimg = 0;

  /*
   *  USAGE:
   *
   *    vector 
   *           name mask 
   *           [principal | dipole | box | corrplane | ired mask2 | corr mask2 | corrired mask2 | mask2] 
   *           [order <order>] 
   *           [modes <modesfile>] [beg <beg>] [end <end>] [npair <npair>]
   *           [out <filename>]
   *
   *  action argument usage:
   *    carg1:
   *      transformVectorInfo
   *
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

    /*
     *  set up complex argument
     */
    vectorInfo = (transformVectorInfo *)
      safe_malloc(sizeof(transformVectorInfo));
    INITIALIZE_transformVectorInfo(vectorInfo);
    vectorInfo->totalFrames = -1;
    vectorInfo->frame = 0;

    /*
     *  get vector name
     */
    vectorInfo->name = getArgumentString(argumentStackPointer, NULL);

    /*
     *  get atom mask for this vector
     */
    buffer = getArgumentString(argumentStackPointer, NULL);
    vectorInfo->mask = processAtomMask(buffer, action->state);
    safe_free(buffer);

    /*
     *  check to see if there are any files to be output
     */
    vectorInfo->filename = argumentStackKeyToString(argumentStackPointer, "out", NULL);
    
    /*
     *  get order for Legendre polynomials
     */
    vectorInfo->order = argumentStackKeyToInteger(argumentStackPointer, "order", 2);
    if(vectorInfo->order < 0 || vectorInfo->order > 2){
      fprintf(stdout,
              "FYI: Order was given out of bounds (<0 or >2), resetting it to 2\n");
      vectorInfo->order = 2;
    }

    /*
     *  process the remaining arguments
     */
    if ( argumentStackContains(argumentStackPointer, "principal")  ) {
      vectorInfo->mode = VECTOR_PRINCIPAL_X;
      // DRR - Check the next argument for x, y, or z. 
      //buffer = getArgumentString(argumentStackPointer, NULL);
      //if (buffer!=NULL) {
        if (argumentStringContains(argumentStackPointer,"x")) {
          vectorInfo->mode = VECTOR_PRINCIPAL_X; 
        } else if (argumentStringContains(argumentStackPointer,"y")) {
	  vectorInfo->mode = VECTOR_PRINCIPAL_Y;
        } else if (argumentStringContains(argumentStackPointer,"z")) {
	  vectorInfo->mode = VECTOR_PRINCIPAL_Z;
        } //else
          //pushBottomStack(argumentStackPointer, buffer); 
      //}
    } else if (argumentStackContains(argumentStackPointer, "dipole")) 
      vectorInfo->mode = VECTOR_DIPOLE;
    else if (argumentStackContains(argumentStackPointer, "box")) 
      vectorInfo->mode = VECTOR_BOX;
    else if (argumentStackContains(argumentStackPointer, "corrplane"))
      vectorInfo->mode = VECTOR_CORRPLANE;
    else if (argumentStackContains(argumentStackPointer, "corrired"))
      vectorInfo->mode = VECTOR_CORRIRED;
    else if (argumentStackContains(argumentStackPointer, "corr"))
      vectorInfo->mode = VECTOR_CORR;
    else if (argumentStackContains(argumentStackPointer, "ired"))
      vectorInfo->mode = VECTOR_IRED;
    else
      vectorInfo->mode = VECTOR_MASK;

    /*
     *  for VECTOR_CORRIRED
     */
    if(vectorInfo->mode == VECTOR_CORRIRED){
      /*
       *  get pair number
       */
      vectorInfo->npair = argumentStackKeyToInteger(argumentStackPointer, "npair", 0);
      if(vectorInfo->npair == 0){
        fprintf(stdout, "WARNING in ptraj(), vector: no npair information given, ignoring command\n");
        freeTransformVectorMemory(action);
        return -1;
      }

      if((buffer = argumentStackKeyToString(argumentStackPointer, "modes", NULL)) == NULL){
        fprintf(stdout, "WARNING in ptraj(), vector: no modes file given, ignoring command\n");
        freeTransformVectorMemory(action);
        return -1;
      }
      else{
        vectorInfo->ibeg = ibeg = argumentStackKeyToInteger(argumentStackPointer, "beg", 1);
        vectorInfo->iend = iend = argumentStackKeyToInteger(argumentStackPointer, "end", 50);

        /*
         *  See if modes info is already available 
         */
        for(vectorStackTmp = vectorStack;
            vectorStackTmp != NULL;
            vectorStackTmp = vectorStackTmp->next){

          vectorInfoTmp = (transformVectorInfo *) vectorStackTmp->entry;

          if(vectorInfoTmp->mode == vectorInfo->mode &&
             vectorInfoTmp->ibeg == ibeg &&
             vectorInfoTmp->iend == iend &&
             vectorInfoTmp->modinfo != NULL &&
             strcmp(vectorInfoTmp->modinfo->name, buffer) == 0){

            /*
             * Link to already known modes; set vinfo->master=0 to indicate
             *   in freeTransformVectorMemory that this action doesn't need to
             *   clean up modes
             */
            vectorInfo->modinfo = vectorInfoTmp->modinfo;
            vectorInfo->master = 0;
            break;
          }
        }

        if(vectorInfo->modinfo == NULL){
          /*
           * Need to load modes from file; set vinfo->master=1 to indicate
           *   that this action has to clean up modes in freeTransformVectorMemory
           */
          vectorInfo->master = 1;

          /*
           *  Allocate modesInfo structure
           */
          vectorInfo->modinfo = modinfo = (modesInfo *) safe_malloc(sizeof(modesInfo));
          INITIALIZE_modesInfo(modinfo);
          modinfo->name = buffer;
          modinfo->type = MT_UNKNOWN;
          modinfo->source = MS_FILE;

          /*
           *  Read evec file
           */
          fp = safe_fopen(buffer, "r");
          if(fp == NULL){
            fprintf(stdout,
                    "WARNING in ptraj(), vector: file %s not opened, ignoring command\n", buffer);
            freeTransformVectorMemory(action);
            return -1;
          }
          if(readEvecFile(fp, ibeg, iend, modinfo)){
            fprintf(stdout,
                    "WARNING in ptraj(), vector: error while reading %s, ignoring command\n", buffer);
            freeTransformVectorMemory(action);
            return -1;
          }
          if(modinfo->nvect != (iend - ibeg + 1)){
            fprintf(stdout,
                    "FYI: Number of read evecs is %i, number of requested evecs is %i\n", 
                    modinfo->nvect, iend - ibeg + 1);
          }
          safe_fclose(fp);
        }
      }
    }

    /*
     * for VECTOR_CORRPLANE
     */

    if(vectorInfo->mode == VECTOR_CORRPLANE){
      n = 0;
      for(i = 0; i < action->state->atoms; i++)
        if(vectorInfo->mask[i])
          n++;
      if(n < 3){
        fprintf(stdout, "WARNING in ptraj(), vector: < 3 atoms given for vector corrplane, ignoring command\n");
        freeTransformVectorMemory(action);
        return -1;
      }
    }

    /*
     * Get second mask if necessary
     */
    if(vectorInfo->mode == VECTOR_IRED ||
       vectorInfo->mode == VECTOR_CORR ||
       vectorInfo->mode == VECTOR_CORRIRED ||
       vectorInfo->mode == VECTOR_MASK){
      buffer = getArgumentString(argumentStackPointer, NULL);
      if (buffer==NULL) {
        fprintf(stdout,"Error: vector: specified vector mode requires a second mask.\n");
        freeTransformVectorMemory(action);
        return -1;
      }
      vectorInfo->mask2 = processAtomMask(buffer, action->state);
      safe_free(buffer);
    }

    action->carg1 = (void *) vectorInfo;

    /*
     *  push the vector info onto the vector stack and store it in action object
     */

    if(vectorInfo->mode == VECTOR_IRED){
      /*
       *  Inverse vector storage necessary for IRED; 
       *  otherwise IRED matrix will be (N,N)->(1,1) instead of (1,1)->(N,N) 
       */
      pushBottomStack(&vectorStack, (void *) vectorInfo);
    }
    else{
      pushStack(&vectorStack, (void *) vectorInfo);
    }

    return 0;
  }

  vectorInfo = (transformVectorInfo *) action->carg1;

  if (mode == PTRAJ_STATUS) {

    /*
     *  ACTION: PTRAJ_STATUS
     */

    fprintf(stdout, "  VECTOR: storage to array named %s\n", vectorInfo->name);

    switch (vectorInfo->mode) {

    case VECTOR_DIPOLE:
      fprintf(stdout, "      The dipole moment vector with respect to the center of mass\n");
      fprintf(stdout, "      will be processed for mask atoms: ");
      printAtomMask(stdout, vectorInfo->mask, action->state);
      fprintf(stdout, "\n");
      break;

    case VECTOR_PRINCIPAL_X:
    case VECTOR_PRINCIPAL_Y:
    case VECTOR_PRINCIPAL_Z:
      fprintf(stdout, "      The principal axis vector ");
      if (vectorInfo->mode == VECTOR_PRINCIPAL_X) {
	fprintf(stdout, "(%c) with respect to the\n", 'X');
      } else if (vectorInfo->mode == VECTOR_PRINCIPAL_Y) {
	fprintf(stdout, "(%c) with respect to the\n", 'Y');
      } else {
	fprintf(stdout, "(%c) with respect to the\n", 'Z');
      }
      fprintf(stdout, "      center of mass of the mask atoms will be dumped: ");
      printAtomMask(stdout, vectorInfo->mask, action->state);
      fprintf(stdout, "\n");
      break;

    case VECTOR_MASK:
    case VECTOR_IRED:
    case VECTOR_CORR:
    case VECTOR_CORRIRED:
      fprintf(stdout, 
         "      Calculate the vector between the center of mass of the two atom selections\n");
      fprintf(stdout, 
	      "      which follow (with the origin at the center of mass of the first)\n");
      fprintf(stdout, "      Atom selection 1 is ");
      printAtomMask(stdout, vectorInfo->mask, action->state);
      fprintf(stdout, "\n");
      fprintf(stdout, "      Atom selection 2 is ");
      printAtomMask(stdout, vectorInfo->mask2, action->state);
      fprintf(stdout, "\n");
      break;

    case VECTOR_CORRPLANE:
      fprintf(stdout, 
         "      Calculate the vector perpendicular to the least squares best plane\n");
      fprintf(stdout, 
	      "      through the atom selection which follows\n");
      fprintf(stdout, "      Atom selection is ");
      printAtomMask(stdout, vectorInfo->mask, action->state);
      fprintf(stdout, "\n");
      break;

    case VECTOR_BOX:
      fprintf(stdout, "      The box lengths will be treated as a vector\n");
      break;
    case VECTOR_NOOP:
      fprintf(stdout, "Error: vector: No operation.\n");
      return -1;
    }

    if(vectorInfo->mode == VECTOR_CORRPLANE ||
       vectorInfo->mode == VECTOR_CORR ||
       vectorInfo->mode == VECTOR_CORRIRED)
      fprintf(stdout, "      The order of Legendre polynomials is %i\n", vectorInfo->order);

    if(vectorInfo->mode == VECTOR_CORRIRED){
      fprintf(stdout, "      IRED modes are read from %s,\n", vectorInfo->modinfo->name);
      fprintf(stdout, "      and the pair %i is considered\n", vectorInfo->npair);
    }     

    if(vectorInfo->mode == VECTOR_CORRPLANE ||
       vectorInfo->mode == VECTOR_CORR ||
       vectorInfo->mode == VECTOR_CORRIRED ||
       vectorInfo->mode == VECTOR_IRED) {
      
      if (vectorInfo->filename) {
	fprintf(stdout, "      Warning: Output of corr, ired, corrired or corrplane vectors is not yet supported!\n");
	safe_free(vectorInfo->filename);
	vectorInfo->filename = NULL;
      }
    }

    if (vectorInfo->filename != NULL) {
      fprintf(stdout, "      Output will be dumped to a file, %s\n",
	      vectorInfo->filename);
    }

  } else if (mode == PTRAJ_PRINT) {

    /*
     *  ACTION: PTRAJ_PRINT
     */

    if (vectorInfo->filename == NULL) return 0;

    if ( ( outFile = fopen(vectorInfo->filename, "w") ) == NULL ) {
      fprintf(stdout, "WARNING in ptraj(), vector: couldn't open file %s\n",
	      vectorInfo->filename);
      return 0;
    }
    fprintf(stdout, "PTRAJ VECTOR: dumping vector information %s\n",
	    vectorInfo->name);
    fprintf(outFile, 
	    "# FORMAT: frame vx vy vz cx cy cz cx+vx cy+vy cz+vz\n");
    fprintf(outFile,
	    "# FORMAT where v? is vector, c? is center of mass...\n");
    for (i=0; i < action->state->maxFrames; i++) {
      fprintf(outFile, "%i %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
	      i+1, 
	      vectorInfo->vx[i], vectorInfo->vy[i], vectorInfo->vz[i],
	      vectorInfo->cx[i], vectorInfo->cy[i], vectorInfo->cz[i],
	      vectorInfo->cx[i]+vectorInfo->vx[i], 
	      vectorInfo->cy[i]+vectorInfo->vy[i], 
	      vectorInfo->cz[i]+vectorInfo->vz[i]);
    }
    safe_fclose(outFile);

  } else if (mode == PTRAJ_CLEANUP) {

    /*
     *  ACTION: PTRAJ_CLEANUP
     */
    freeTransformVectorMemory(action);
  }

  if (mode != PTRAJ_ACTION) return 0;

  /*
   *  ACTION: PTRAJ_ACTION
   */

  state = (ptrajState *) action->state;

  /*
   *  update local state information
   */
  for (i=0; i<6; i++)
    state->box[i] = box[i];

  if (vectorInfo->totalFrames < 0) {
    if(vectorInfo->mode == VECTOR_IRED)
      vectorInfo->totalFrames = 1;
    else
      vectorInfo->totalFrames = state->maxFrames;

    if(vectorInfo->mode == VECTOR_CORRIRED){
      if(vectorInfo->master == 0){
        /*
         * Find master vector for this "slave"
         */
        for(vectorStackTmp = vectorStack;
            vectorStackTmp != NULL;
            vectorStackTmp = vectorStackTmp->next){

          vectorInfoTmp = (transformVectorInfo *) vectorStackTmp->entry;

          if(vectorInfoTmp->mode    == VECTOR_CORRIRED &&
             vectorInfoTmp->master  == 1 &&
             vectorInfoTmp->modinfo == vectorInfo->modinfo){
            /*
             * Link to already known vector
             */
            vectorInfo->cftmp = vectorInfoTmp->cftmp;
            break;
          }
        }
      }
      else{
        /*
         * I'm the master; allocate memory for cftmp
         */
        n = 2 * 
            vectorInfo->totalFrames * 
            (2 * vectorInfo->order + 1) *
            vectorInfo->modinfo->nvect;
        vectorInfo->cftmp = (double *) safe_malloc(sizeof(double) * n);
        for(i = 0; i < n; i++)
          vectorInfo->cftmp[i] = 0.0;
      }
    }
    else if(vectorInfo->mode == VECTOR_CORRPLANE || 
            vectorInfo->mode == VECTOR_CORR){
      n = 3;
      vectorInfo->avgcrd = (double *) safe_malloc(sizeof(double) * n);
      for(i = 0; i < n; i++)
        vectorInfo->avgcrd[i] = 0.0;
      n = 2 * 
          vectorInfo->totalFrames * 
          (2 * vectorInfo->order + 1);
      vectorInfo->cftmp = (double *) safe_malloc(sizeof(double) * n);
      vectorInfo->p2cftmp = (double *) safe_malloc(sizeof(double) * n);
      for(i = 0; i < n; i++){
        vectorInfo->cftmp[i] = 0.0;
        vectorInfo->p2cftmp[i] = 0.0;
      }
      n = 2 * 
          vectorInfo->totalFrames;
      vectorInfo->rcftmp = (double *) safe_malloc(sizeof(double) * n);
      for(i = 0; i < n; i++)
        vectorInfo->rcftmp[i] = 0.0;

      if(vectorInfo->mode == VECTOR_CORRPLANE){
        n = 0;
        for (i=0; i < state->atoms; i++)
          if (vectorInfo->mask[i])
            n++;
        vectorInfo->cx = (double *) 
          safe_malloc(sizeof(double) * n);
        vectorInfo->cy = (double *) 
          safe_malloc(sizeof(double) * n);
        vectorInfo->cz = (double *) 
          safe_malloc(sizeof(double) * n);
      }
    }
    else{
      vectorInfo->cx = (double *) 
        safe_malloc(sizeof(double) * vectorInfo->totalFrames);
      vectorInfo->cy = (double *) 
        safe_malloc(sizeof(double) * vectorInfo->totalFrames);
      vectorInfo->cz = (double *) 
        safe_malloc(sizeof(double) * vectorInfo->totalFrames);
      vectorInfo->vx = (double *) 
        safe_malloc(sizeof(double) * vectorInfo->totalFrames);
      vectorInfo->vy = (double *) 
        safe_malloc(sizeof(double) * vectorInfo->totalFrames);
      vectorInfo->vz = (double *) 
        safe_malloc(sizeof(double) * vectorInfo->totalFrames);
    }
  }

  if (vectorInfo->frame > vectorInfo->totalFrames) {
    warning("transformVector()", "Blowing array; too many frames!!\n");
    return 0;
  }

  switch( vectorInfo->mode ) {

  case VECTOR_CORRIRED:
  case VECTOR_CORR:
  case VECTOR_CORRPLANE:
    vmode = vectorInfo->mode;
    frame = vectorInfo->frame;
    order = vectorInfo->order;
    indsnap = (2*order + 1) * frame;
    cftmp = vectorInfo->cftmp;
    if(vmode == VECTOR_CORRIRED){
      nvect = vectorInfo->modinfo->nvect;
      nvectelem = vectorInfo->modinfo->nvectelem;
      evec  = vectorInfo->modinfo->evec;
      npair = vectorInfo->npair - 1;
      indsnap *= nvect;
    }
    else if(vmode == VECTOR_CORR || 
            vmode == VECTOR_CORRPLANE){
      avgcrd  = vectorInfo->avgcrd;
      p2cftmp = vectorInfo->p2cftmp;
      rcftmp  = vectorInfo->rcftmp;
    }
 
    /*
     * Calc center of mass of masks and vector v(x,y,z)
     */
    total_mass = 0.0;
    cx = 0.0;
    cy = 0.0;
    cz = 0.0;
    for (i=0; i < state->atoms; i++) {
      if (vectorInfo->mask[i]) {
	cx += state->masses[i] * x[i];
	cy += state->masses[i] * y[i];
	cz += state->masses[i] * z[i];
	total_mass += state->masses[i];
      }
    }
    cx /= total_mass;
    cy /= total_mass;
    cz /= total_mass;

    if(vectorInfo->mode == VECTOR_CORR ||
       vectorInfo->mode == VECTOR_CORRIRED){
      total_mass = 0.0;
      vx = 0.0;
      vy = 0.0;
      vz = 0.0;
      for (i=0; i < state->atoms; i++) {
        if (vectorInfo->mask2[i]) {
  	  vx += state->masses[i] * x[i];
	  vy += state->masses[i] * y[i];
	  vz += state->masses[i] * z[i];
	  total_mass += state->masses[i];
        }
      }
      vx /= total_mass;
      vy /= total_mass;
      vz /= total_mass;

      vx -= cx;
      vy -= cy;
      vz -= cz;
    }
    else if(vectorInfo->mode == VECTOR_CORRPLANE){
      n = 0;
      for (i=0; i < state->atoms; i++) {
        if (vectorInfo->mask[i]) {
          vectorInfo->cx[n] = x[i] - cx;      
          vectorInfo->cy[n] = y[i] - cy;      
          vectorInfo->cz[n] = z[i] - cz;
          n++;
        }
      }

      lsqplane(n,
               vectorInfo->cx, vectorInfo->cy, vectorInfo->cz,
               &vx, &vy, &vz);
    }
    
    /*
     * Calc vector length
     */
    r = sqrt(vx*vx + vy*vy + vz*vz);

    /*
     * Update avgcrd, rave, r3iave, r6iave for VECTOR_CORR, VECTOR_CORRPLANE
     */
    if(vmode == VECTOR_CORR ||
       vmode == VECTOR_CORRPLANE){
      avgcrd[0]          += vx;
      avgcrd[1]          += vy;
      avgcrd[2]          += vz;
      vectorInfo->rave   += r;
      r3                  = r*r*r;
      r3i                 = 1.0 / r3;
      vectorInfo->r3iave += r3i;
      vectorInfo->r6iave += r3i*r3i;
    }

    /*
     * Loop over m=0, ..., +L
     */
    for(i = 0; i <= order; i++){

      /*
       * Calc spherical harmonics
       */
      spherharm(order, i, vx, vy, vz, r, &dplusreal, &dplusimg);
      if(vmode == VECTOR_CORRIRED)
        indplus = nvect * (order + i);
      else if(vmode == VECTOR_CORR || 
              vmode == VECTOR_CORRPLANE)
        indplus = order + i;

      if(i > 0){                
        spherharm(order, -i, vx, vy, vz, r, &dminusreal, &dminusimg);
        if(vmode == VECTOR_CORRIRED)
          indminus = nvect * (order - i);
        else if(vmode == VECTOR_CORR || 
                vmode == VECTOR_CORRPLANE)
          indminus = order - i;
      }
     
      if(vmode == VECTOR_CORRIRED){
        /*
         * Loop over all eigenvectors
         */
        for(j = 0; j < nvect; j++){
          q = evec[j * nvectelem + npair];
          indtot = 2 * (indsnap + indplus + j);
          cftmp[indtot  ] += (q * dplusreal);
          cftmp[indtot+1] += (q * dplusimg);
          if(i > 0){
            indtot = 2 * (indsnap + indminus + j);
            cftmp[indtot  ] += q * dminusreal;
            cftmp[indtot+1] += q * dminusimg;
          }
        }
      }
      else if(vmode == VECTOR_CORR || 
              vmode == VECTOR_CORRPLANE){
        indtot = 2 * (indsnap + indplus);
        cftmp[indtot  ] += r3i * dplusreal;
        cftmp[indtot+1] += r3i * dplusimg;
        p2cftmp[indtot  ] += dplusreal;
        p2cftmp[indtot+1] += dplusimg;
        if(i > 0){
          indtot = 2 * (indsnap + indminus);
          cftmp[indtot  ] += r3i * dminusreal;
          cftmp[indtot+1] += r3i * dminusimg;
          p2cftmp[indtot  ] += dminusreal;
          p2cftmp[indtot+1] += dminusimg;
        }
        else if(i == 0){
          indtot = 2 * frame;
          rcftmp[indtot  ] += r3i;
          rcftmp[indtot+1]  = 0.0;
        }
      }
    }

    break;
    
  case VECTOR_DIPOLE:

    cx = 0.0;
    cy = 0.0;
    cz = 0.0;
    vx = 0.0;
    vy = 0.0;
    vz = 0.0;
    total_mass = 0.0;
    for (i=0; i < state->atoms; i++) {
      if (vectorInfo->mask[i]) {
	cx += state->masses[i] * x[i];
	cy += state->masses[i] * y[i];
	cz += state->masses[i] * z[i];
	total_mass += state->masses[i];

	vx += x[i] * state->charges[i];
	vy += y[i] * state->charges[i];
	vz += z[i] * state->charges[i];

      }
    }
    cx = cx / total_mass;
    cy = cy / total_mass;
    cz = cz / total_mass;

    vectorInfo->vx[vectorInfo->frame] = vx;
    vectorInfo->vy[vectorInfo->frame] = vy;
    vectorInfo->vz[vectorInfo->frame] = vz;
    vectorInfo->cx[vectorInfo->frame] = cx;
    vectorInfo->cy[vectorInfo->frame] = cy;
    vectorInfo->cz[vectorInfo->frame] = cz;

    if (prnlev > 2) {
      fprintf(stdout, "\nDipole...\n");
      fprintf(stdout, "ATOM      1 DD1  DD      1    %8.3f%8.3f%8.3f\n",
	      cx, cy, cz);
      fprintf(stdout, "ATOM      2 DD2  DD      1    %8.3f%8.3f%8.3f\n",
	      cx+vx, cy+vy, cz+vz);
    }

    break;

  case VECTOR_PRINCIPAL_X:
  case VECTOR_PRINCIPAL_Y:
  case VECTOR_PRINCIPAL_Z:

    principal = calculatePrincipalAxis(state, vectorInfo->mask, 
				       x, y, z, 1, 0, 1);

    if (vectorInfo->mode == VECTOR_PRINCIPAL_X) { 
      vectorInfo->vx[vectorInfo->frame] = principal[0];
      vectorInfo->vy[vectorInfo->frame] = principal[1];
      vectorInfo->vz[vectorInfo->frame] = principal[2];
    } else if (vectorInfo->mode == VECTOR_PRINCIPAL_Y) {
      vectorInfo->vx[vectorInfo->frame] = principal[3];
      vectorInfo->vy[vectorInfo->frame] = principal[4];
      vectorInfo->vz[vectorInfo->frame] = principal[5];
    } else {
      vectorInfo->vx[vectorInfo->frame] = principal[6];
      vectorInfo->vy[vectorInfo->frame] = principal[7];
      vectorInfo->vz[vectorInfo->frame] = principal[8];
    }
    vectorInfo->cx[vectorInfo->frame] = principal[9];
    vectorInfo->cy[vectorInfo->frame] = principal[10];
    vectorInfo->cz[vectorInfo->frame] = principal[11];

    if (prnlev > 2) {
      fprintf(stdout, "\ntransformVector PRINCIPAL AXIS:\n");
      fprintf(stdout, "ATOM      1 PP1  PP      1    %8.3f%8.3f%8.3f\n",
	      vectorInfo->cx[vectorInfo->frame], 
	      vectorInfo->cy[vectorInfo->frame], 
	      vectorInfo->cz[vectorInfo->frame]);

      fprintf(stdout, "ATOM      2 PP2  PP      1    %8.3f%8.3f%8.3f\n",
	      vectorInfo->cx[vectorInfo->frame]+
	      vectorInfo->vx[vectorInfo->frame],
	      vectorInfo->cy[vectorInfo->frame]+
	      vectorInfo->vy[vectorInfo->frame],
	      vectorInfo->cz[vectorInfo->frame]+
	      vectorInfo->vz[vectorInfo->frame]);
    }

    safe_free(principal);

    break;

  case VECTOR_MASK:
  case VECTOR_IRED:

    total_mass = 0.0;
    cx = 0.0;
    cy = 0.0;
    cz = 0.0;
    for (i=0; i < state->atoms; i++) {
      if (vectorInfo->mask[i]) {
	cx += state->masses[i] * x[i];
	cy += state->masses[i] * y[i];
	cz += state->masses[i] * z[i];
	total_mass += state->masses[i];
      }
    }
    cx = cx / total_mass;
    cy = cy / total_mass;
    cz = cz / total_mass;

    total_mass = 0.0;
    vx = 0.0;
    vy = 0.0;
    vz = 0.0;
    for (i=0; i < state->atoms; i++) {
      if (vectorInfo->mask2[i]) {
	vx += state->masses[i] * x[i];
	vy += state->masses[i] * y[i];
	vz += state->masses[i] * z[i];
	total_mass += state->masses[i];
      }
    }
    vx = vx / total_mass;
    vy = vy / total_mass;
    vz = vz / total_mass;

    vectorInfo->vx[vectorInfo->frame] = vx - cx;
    vectorInfo->vy[vectorInfo->frame] = vy - cy;
    vectorInfo->vz[vectorInfo->frame] = vz - cz;
    vectorInfo->cx[vectorInfo->frame] = cx;
    vectorInfo->cy[vectorInfo->frame] = cy;
    vectorInfo->cz[vectorInfo->frame] = cz;

    if (prnlev > 2) {
      fprintf(stdout, "\nMASK...\n");
      fprintf(stdout, "ATOM      1 MM1  MM      1    %8.3f%8.3f%8.3f\n",
	      cx, cy, cz);
      fprintf(stdout, "ATOM      1 MM2  MM      1    %8.3f%8.3f%8.3f\n",
	      vx, vy, vz);
    }

    break;

  case VECTOR_BOX:

    vectorInfo->vx[vectorInfo->frame] = state->box[0];
    vectorInfo->vy[vectorInfo->frame] = state->box[1];
    vectorInfo->vz[vectorInfo->frame] = state->box[2];
    vectorInfo->cx[vectorInfo->frame] = 0.0;
    vectorInfo->cy[vectorInfo->frame] = 0.0;
    vectorInfo->cz[vectorInfo->frame] = 0.0;

    if (prnlev > 2) {
      fprintf(stdout, "\nBOX %8.3f %8.3f %8.3f\n", 
	      vectorInfo->vx[vectorInfo->frame],
	      vectorInfo->vy[vectorInfo->frame], 
	      vectorInfo->vz[vectorInfo->frame]);
    }

    break;
  case VECTOR_NOOP:
    return 0;

  }

  if(vectorInfo->mode != VECTOR_IRED)
    vectorInfo->frame++;

  return 1;

}



/** ACTION ROUTINE *************************************************************
 *
 *  transformWatershell()  --- calculate the number of waters in a given shell
 *
 ******************************************************************************/
typedef struct _transformShellInfo {
  int *soluteMask; 
  int *solventMask;
  int *activeResidues;
  int *lower;
  int *upper;
  int visits;
  double lowerCutoff;
  double upperCutoff;
  char *filename;

} transformShellInfo;

#define INITIALIZE_transformShellInfo(_p_) \
  _p_->soluteMask     = NULL; \
  _p_->solventMask    = NULL; \
  _p_->activeResidues = NULL; \
  _p_->lower          = NULL; \
  _p_->upper          = NULL; \
  _p_->visits         = 0;    \
  _p_->lowerCutoff    = 0.0;  \
  _p_->upperCutoff    = 0.0;  \
  _p_->filename       = NULL

   int
transformWatershell(actionInformation *action, 
		    double *x, double *y, double *z,
		    double *box, int mode)
{
  //char *name = "watershell";
  argStackType **argumentStackPointer;
  char *buffer;//, buffer2[BUFFER_SIZE];
  ptrajState *state;
  transformShellInfo *info;
  int i, j, jj;
  double distance;
  double ucell[9], recip[9];
  void *outFile;

  /*
   *  USAGE:
   *
   *    watershell mask filename [lower <lower cut>] [upper <upper cut>] [noimage]
   *
   *  action argument usage:
   *
   *  carg1: a transformShellInfo structure
   *  iarg1: disable imaging?
   *
   */


  if (mode == PTRAJ_SETUP) {

    /*
     *  ACTION: PTRAJ_SETUP
     */

#ifdef MPI

#endif

    // DEBUG
    //fprintf(stdout,"DEBUG:\taction address (watershell): %x\n",action);
    //fprintf(stdout,"DEBUG:\taction->carg1 address (watershell): %x\n",action->carg1);

    argumentStackPointer = (argStackType **) action->carg1;
    action->carg1 = NULL;
    // DEBUG
    //fprintf(stdout,"DEBUG:\targumentStack address (watershell): %x\n",*argumentStackPointer);
    //printArgumentStack( argumentStackPointer );

    info = (transformShellInfo *)
      safe_malloc(sizeof(transformShellInfo));
    INITIALIZE_transformShellInfo(info);
    info->lowerCutoff = 3.4;
    info->upperCutoff = 5.0;

    action->iarg1 = argumentStackContains(argumentStackPointer, "noimage");

    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer==NULL) {
      fprintf(stdout,"ERROR: WATERSHELL: Solute mask must be specified.\n");
      return -1;
    }

    info->soluteMask = processAtomMask(buffer, action->state);
    if (info->soluteMask==NULL) {
      fprintf(stdout,"ERROR: WATERSHELL: Solute mask %s corresponds to 0 atoms.\n",buffer);
      safe_free(buffer);
      return -1;
    }
    safe_free(buffer);

    info->filename = getArgumentString(argumentStackPointer, NULL);
    if (info->filename==NULL) {
      fprintf(stdout,"ERROR: WATERSHELL: Output filename must be specified.\n");
      return -1;
    }

    info->lowerCutoff = 
      argumentStackKeyToDouble(argumentStackPointer, "lower", info->lowerCutoff);
    info->upperCutoff = 
      argumentStackKeyToDouble(argumentStackPointer, "upper", info->upperCutoff);

    buffer = getArgumentString(argumentStackPointer, NULL);
    if (buffer != NULL) { 
      //fprintf(stdout,"DEBUG:\tWATERSHELL: Using %s as solvent mask.\n",buffer);
      info->solventMask = processAtomMask(buffer, action->state);
    } else {
      //fprintf(stdout,"DEBUG:\tWATERSHELL: Using :WAT as solvent mask.\n");
      info->solventMask = processAtomMask(":WAT", action->state);
    }
    if (info->solventMask==NULL) {
      if (buffer!=NULL)
        fprintf(stdout,"ERROR: WATERSHELL: Solvent mask %s corresponds to 0 atoms.\n",buffer);
      else {
        fprintf(stdout,"ERROR: WATERSHELL: Default solvent mask :WAT corresponds to 0 atoms.\n");
        fprintf(stdout,
          "                   Solvent mask can be specified as the third argument.\n");
      }
      safe_free(buffer);
      return -1;
    }
    safe_free(buffer);

    action->carg1 = (void *) info;

    return 0;

  }


  info = (transformShellInfo *) action->carg1;


  if (mode == PTRAJ_STATUS) {

    /*
     *  ACTION: PTRAJ_STATUS
     */

    fprintf(stdout, "  WATER SHELL: Output to %s\n", info->filename);
    if (action->iarg1)
      fprintf(stdout, "      Imaging is disabled\n");
    fprintf(stdout, 
	    "      The first shell will contain water < %5.3f angstroms from\n",
	    info->lowerCutoff);
    fprintf(stdout,
	    "      the solute; the second shell < %5.3f angstroms...\n",
	    info->upperCutoff);
    fprintf(stdout, "      The solute atoms are ");
    printAtomMask(stdout, info->soluteMask, action->state);
    fprintf(stdout, "\n");
    fprintf(stdout, "      The solvent atoms are ");
    printAtomMask(stdout, info->solventMask, action->state);
    fprintf(stdout, "\n");

  } else if (mode == PTRAJ_PRINT) {

    /*
     *  ACTION: PTRAJ_PRINT
     */

    outFile = ptrajOpenW(info->filename);
    if ( outFile == NULL ) {
      fprintf(stdout, "WARNING in ptraj(), watershell: couldn't open output file %s\n",
	      info->filename);
      return 0;
    }

    fprintf(stdout, "PTRAJ WATERSHELL: dumping data to output file\n");
    for (i=0; i < action->state->maxFrames/worldsize; i++) {
      ptrajfprintf(outFile, "%i %i %i\n",
	      i*worldsize+worldrank+1, info->lower[i], info->upper[i]);
    }
    ptrajCloseFile(outFile);

  } else if (mode == PTRAJ_CLEANUP) {

    /*
     *  ACTION: PTRAJ_CLEANUP
     */

    safe_free(info->filename);
    safe_free(info->upper);
    safe_free(info->lower);
    safe_free(info->activeResidues);
    safe_free(info->solventMask);
    safe_free(info->soluteMask);
    INITIALIZE_transformShellInfo(info);
    safe_free(info);

  }



  if (mode != PTRAJ_ACTION) return 0;

  /*
   *  ACTION: PTRAJ_ACTION
   */

  state = (ptrajState *) action->state;
     /*
      *  update local state information
      */
  for (i=0; i<6; i++)
    state->box[i] = box[i];

     /*
      *  set up information for imaging non-orthorhombic if necessary
      */
  if (action->iarg1 == 0 &&
      (box[3] != 90.0 || box[4] != 90.0 || box[5] != 90.0))
    boxToRecip(box, ucell, recip);

     /*
      *  allocate space for saving results if this is the first visit
      */
  if ( info->activeResidues == NULL ) {
    info->activeResidues = (int *)
      safe_malloc(sizeof(int) * state->residues);
  }
  if ( info->lower == NULL ) {
    info->lower = (int *)
      safe_malloc(sizeof(int) * state->maxFrames);
    info->upper = (int *)
      safe_malloc(sizeof(int) * state->maxFrames);
    for (i=0; i < state->maxFrames; i++) {
      info->lower[i] = 0.0;
      info->upper[i] = 0.0;
    }
    info->visits = 0;
  }


  for (j=0; j < state->residues; j++) {
    info->activeResidues[ j ] = 0;
  }

     /*
      *  loop over all active solute atoms
      */
  for (i=0; i < state->atoms; i++) {

    if ( info->soluteMask[i] > 0 ) {

         /*
          *  loop over solvent atoms by residue
          */

      for (j=0; j < state->residues; j++) {
	
	for (jj=state->ipres[j]-1; jj < state->ipres[j+1]-1; jj++) {
	  if (jj > state->atoms)
	    printf("WARNING in ptraj(), watershell: Blew atom arrays\n");

	  if ( info->solventMask[jj] > 0 ) { 

	    distance = calculateDistance2(i, jj, x, y, z, box, ucell, recip, 0.0, action->iarg1);
	    distance = sqrt(distance);
            // DEBUG
            //fprintf(stderr,"%i %i %i %lf\n",info->visits, i, jj, distance);

	    if (distance < info->upperCutoff && 
		info->activeResidues[j] == 0) {
	      info->activeResidues[j] = 1;
	    }
	    if ( distance < info->lowerCutoff ) {
	      info->activeResidues[j] = 2;
	      break;
	    }
	  }
	}
      }
    }
  }


  info->lower[info->visits] = 0;
  info->upper[info->visits] = 0;
  for (i=0; i < state->residues; i++) {
    if ( info->activeResidues[i] == 2 ) {
      info->lower[info->visits] += 1;
      info->upper[info->visits] += 1;
    } else if ( info->activeResidues[i] == 1 ) {
      info->upper[info->visits] += 1;
    }
  }

  info->visits++;

  return 1;
}




