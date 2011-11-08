#ifndef INC_PTRAJ_ACTIONS_H
#define INC_PTRAJ_ACTIONS_H
#ifdef __cplusplus
extern "C" {
#endif
// This file was adapted from ptraj/actions.h for cpptraj
// Dan Roe, Rutgers, November 2011
// INCLUDES
#include "ptraj_evec.h" // modesInfo in transformVectorInfo
#include "Name.h" // NAME in ptrajState

// ---------- Argument Stack routines ------------------------------------------
// argStackType - just used to set up arglist
typedef struct _argStackType {
  int nargs;
  char **arglist;
  char *marked;
} argStackType;
void printArgumentStack(argStackType **);

// ---------- Routines that access global ptraj_actions vars -------------------
void SetReferenceInfo(double*,int);
void FreeReferenceInfo();
void SetPrnlev(int);

// ---------- ptrajMode --------------------------------------------------------
// originally from ptraj_local.h
typedef enum _ptrajMode {
  PTRAJ_NOOP,  PTRAJ_ACTION, PTRAJ_FIRSTPASS, PTRAJ_SECONDPASS,
  PTRAJ_SETUP, PTRAJ_STATUS, PTRAJ_PRINT,     PTRAJ_CLEANUP 
} ptrajMode;

// ---------- Ptraj State ------------------------------------------------------
// originally from ptraj_local.h
// Redefine Name to match definition in Ptraj
typedef NAME Name;
typedef struct _ptrajState {
  double box[6];             // box lengths and angles 
  double *masses;            // atom masses 
  double *charges;           // atom charges 
  int atoms;                 // number of atoms 
  int residues;              // number of residues 
  int *ipres;                // atoms in each residue (atom #s start from 1)
  int *ipres_mask;           // ipres for mask parser, atom #s start from 0 
  int IFBOX;                 // is there box information? 
  int boxfixed[6];           // equal to 1 if this box coordinate is fixed 
  int molecules;             // total number of molecules 
  int *moleculeInfo;         // number of atoms in each molecule 
  int *solventMask;          // atoms in the solvent 
  int solventMolecules;      // number of solvent molecules 
  int *solventMoleculeStart; // pointer into solventMask for first atom of each solvent 
  int *solventMoleculeStop;  // pointer into solventMask for last atom of each solvent 
  int solventAtoms;          // number of solvent atoms 
  Name *atomName;            // atom names 
  Name *residueName;         // residue names 
  int maxFrames;             // number of useful frames in 1 pass of trajin's 
  double temp0;              // DAN TEST: for writing out netcdf temp trajs 
} ptrajState;

#define INITIALIZE_ptrajState( _p_ ) \
  _p_->box[0] = 0.0; _p_->box[1] = 0.0; _p_->box[2] = 0.0; \
  _p_->box[3] = 90.0; _p_->box[4] = 90.0; _p_->box[5] = 90.0; \
  _p_->masses = NULL; \
  _p_->charges = NULL; \
  _p_->atoms = 0; \
  _p_->residues = 0; \
  _p_->ipres = NULL; \
  _p_->IFBOX = 0; \
  _p_->boxfixed[0] = 0; _p_->boxfixed[1] = 0; _p_->boxfixed[2] = 0; \
  _p_->boxfixed[3] = 0; _p_->boxfixed[4] = 0; _p_->boxfixed[5] = 0; \
  _p_->molecules = 0; \
  _p_->moleculeInfo = NULL; \
  _p_->solventMask = NULL; \
  _p_->solventMolecules = 0; \
  _p_->solventMoleculeStart = NULL; \
  _p_->solventMoleculeStop = NULL; \
  _p_->solventAtoms = 0; \
  _p_->atomName = NULL; \
  _p_->residueName = NULL; \
  _p_->maxFrames = 0; \
  _p_->temp0 = 0.0

// ---------- Ptraj actions ----------------------------------------------------
// possible ptraj actions
typedef enum _actionType {
  TRANSFORM_NOOP, TRANSFORM_ACCEPTOR, TRANSFORM_ANALYZE, TRANSFORM_ANGLE,
  TRANSFORM_ATOMICFLUCT, TRANSFORM_ATOMICFLUCT3D, TRANSFORM_AVERAGE,
  TRANSFORM_BENCHMARK, TRANSFORM_BOX, TRANSFORM_CENTER, TRANSFORM_CHECKDNA,
  TRANSFORM_CHECKOVERLAP, TRANSFORM_CLOSESTWATERS, TRANSFORM_CLUSTER,
  TRANSFORM_CLUSTERATTRIBUTE, TRANSFORM_CONTACTS, TRANSFORM_DIFFUSION,
  TRANSFORM_DIHEDRAL, TRANSFORM_DIHEDRALCLUSTER, TRANSFORM_DIPOLE,
  TRANSFORM_DISTANCE, TRANSFORM_DONOR, TRANSFORM_DNAIONTRACKER, TRANSFORM_ECHO,
  TRANSFORM_ENERGY, TRANSFORM_GRID, TRANSFORM_HBOND, TRANSFORM_IMAGE,
  TRANSFORM_MATRIX, TRANSFORM_PRINCIPAL, TRANSFORM_PRNLEV, TRANSFORM_PROJECTION,
  TRANSFORM_PUCKER, TRANSFORM_RADIAL, TRANSFORM_RADIUSOFGYRATION,
  TRANSFORM_RANDOMIZEIONS, TRANSFORM_REFERENCE, TRANSFORM_RMS,
  TRANSFORM_RUNNINGAVERAGE, TRANSFORM_SCALE, TRANSFORM_SECONDARYSTRUCT,
  TRANSFORM_STRIP, TRANSFORM_SMARTIMAGE, TRANSFORM_SOLVENT, TRANSFORM_TRAJIN,
  TRANSFORM_TRAJOUT, TRANSFORM_TRANSFORM, TRANSFORM_TRANSLATE,
  TRANSFORM_TRUNCOCT, TRANSFORM_TEST, TRANSFORM_UNWRAP, TRANSFORM_VECTOR,
  TRANSFORM_CORRELATION, TRANSFORM_WATERSHELL, TRANSFORM_2DRMS
} actionType;

// ACTION FUNCTION TYPE DEFINITION
// action->fxn(action, X, Y, Z, box, ptrajMode)
typedef int (*actionFunction)( void *,double *,double *,double *,double *,int);

// ACTION INFORMATION (passed in as the first argument 
// to the action function)
typedef struct _actionInformation {
  actionFunction fxn;
  actionType type;
  ptrajState *state;
  int suppressProcessing;
  int performSecondPass;
  int *mask;
  int iarg1;
  int iarg2;
  int iarg3;
  int iarg4;
  int iarg5;
  int iarg6;
  int iarg7;
  int iarg8;
  double darg1;
  double darg2;
  double darg3;
  double darg4;
  double darg5;
  double darg6;
  double darg7;
  double darg8;
  void *carg1;
  void *carg2;
  void *carg3;
  void *carg4;
  void *carg5;
  void *carg6;
  void *carg7;
  void *carg8;
} actionInformation;

#define INITIALIZE_actionInformation(_p_) \
  _p_->fxn = NULL;            \
  _p_->type = TRANSFORM_NOOP; \
  _p_->state = NULL;          \
  _p_->suppressProcessing = 0;\
  _p_->performSecondPass = 0; \
  _p_->mask = NULL;           \
  _p_->iarg1 = 0;             \
  _p_->iarg2 = 0;             \
  _p_->iarg3 = 0;             \
  _p_->iarg4 = 0;             \
  _p_->iarg5 = 0;             \
  _p_->iarg6 = 0;             \
  _p_->iarg7 = 0;             \
  _p_->iarg8 = 0;             \
  _p_->darg1 = 0.0;           \
  _p_->darg2 = 0.0;           \
  _p_->darg3 = 0.0;           \
  _p_->darg4 = 0.0;           \
  _p_->darg5 = 0.0;           \
  _p_->darg6 = 0.0;           \
  _p_->darg7 = 0.0;           \
  _p_->darg8 = 0.0;           \
  _p_->carg1 = NULL;          \
  _p_->carg2 = NULL;          \
  _p_->carg3 = NULL;          \
  _p_->carg4 = NULL;          \
  _p_->carg5 = NULL;          \
  _p_->carg6 = NULL;          \
  _p_->carg7 = NULL;          \
  _p_->carg8 = NULL

// -----------------------------------------------------------------------------
//  ACTION: TRANSFORM_MATRIX_INFO -- this structure needs to be global since it is
//                                   accessed by analyze.c functions
typedef enum _transformMatrixType {
  MATRIX_NULL,   MATRIX_DIST,      MATRIX_COVAR, MATRIX_MWCOVAR,
  MATRIX_CORREL, MATRIX_DISTCOVAR, MATRIX_IDEA,  MATRIX_IRED
} transformMatrixType;

typedef struct _transformMatrixInfo {
  ptrajState *state;
  char *name;
  transformMatrixType type;
  double *vect;
  double *vect2;
  double *mat;
  int vectsize;
  int matsize;
  int *mask1;
  int *mask2;
  int mask1tot;
  int mask2tot;
  int snap;
} transformMatrixInfo;

#define INITIALIZE_transformMatrixInfo(_p_) \
  _p_->state         = NULL;   \
  _p_->name          = NULL;   \
  _p_->type          = MATRIX_NULL; \
  _p_->vect          = NULL;   \
  _p_->vect2         = NULL;   \
  _p_->mat           = NULL;   \
  _p_->vectsize      = 0;      \
  _p_->matsize       = 0;      \
  _p_->mask1         = NULL;   \
  _p_->mask2         = NULL;   \
  _p_->mask1tot      = 0;      \
  _p_->mask2tot      = 0;      \
  _p_->snap          = 0;

// -----------------------------------------------------------------------------
//  ACTION: TRANSFORM_VECTOR
typedef enum _vectorMode { 
  VECTOR_NOOP,      VECTOR_PRINCIPAL_X, VECTOR_PRINCIPAL_Y, VECTOR_PRINCIPAL_Z, 
  VECTOR_DIPOLE,    VECTOR_BOX,         VECTOR_MASK,        VECTOR_IRED,
  VECTOR_CORRPLANE, VECTOR_CORR,        VECTOR_CORRIRED
} vectorMode;

typedef struct _transformVectorInfo {
  char *name;
  char *filename;
  int totalFrames;
  int frame;
  vectorMode mode;
  int *mask;
  int *mask2;
  double *cx, *cy, *cz;
  double *vx, *vy, *vz;

  int master;
  modesInfo *modinfo;
  int ibeg;
  int iend;
  int order;
  int npair;
  double *avgcrd;
  double rave;
  double r3iave;
  double r6iave;
  double *cftmp;
  double *p2cftmp;
  double *rcftmp;
} transformVectorInfo;

#define INITIALIZE_transformVectorInfo(_p_) \
  _p_->name        = NULL;                \
  _p_->filename    = NULL;                \
  _p_->totalFrames = 0;                   \
  _p_->frame       = 0;                   \
  _p_->mode        = VECTOR_NOOP;         \
  _p_->mask        = NULL;                \
  _p_->mask2       = NULL;                \
  _p_->cx          = NULL;                \
  _p_->cy          = NULL;                \
  _p_->cz          = NULL;                \
  _p_->vx          = NULL;                \
  _p_->vy          = NULL;                \
  _p_->vz          = NULL;                \
  _p_->master      = 0;                   \
  _p_->modinfo     = NULL;                \
  _p_->ibeg        = 1;                   \
  _p_->iend        = 50;                  \
  _p_->order       = 2;                   \
  _p_->npair       = -1;                  \
  _p_->avgcrd      = NULL;                \
  _p_->rave        = 0.0;                 \
  _p_->r3iave      = 0.0;                 \
  _p_->r6iave      = 0.0;                 \
  _p_->cftmp       = NULL;                \
  _p_->p2cftmp     = NULL;                \
  _p_->rcftmp      = NULL;

// -----------------------------------------------------------------------------
// ACTION: TRANSFORM_CORRELATION
typedef struct _transformCorrInfo {
  char *name;
  char *filename;
  int totalFrames;
  int frame;
  int mode;
  int *mask;
  int *mask2;
  double *cx, *cy, *cz;
  double *vx, *vy, *vz;
  int tmin;
  int tcorr;
  int tmax;
} transformCorrInfo;

enum _corrInfo { CORR_BOX, CORR_MASK };

#define INITIALIZE_transformCorrInfo(_p_) \
  _p_->name        = NULL;                \
  _p_->filename    = NULL;                \
  _p_->totalFrames = 0;                   \
  _p_->frame       = 0;                   \
  _p_->mode        = 0;                   \
  _p_->mask        = NULL;                \
  _p_->mask2       = NULL;                \
  _p_->cx          = NULL;                \
  _p_->cy          = NULL;                \
  _p_->cz          = NULL;                \
  _p_->vx          = NULL;                \
  _p_->vy          = NULL;                \
  _p_->vz          = NULL;                \
  _p_->tmin        = 0;                   \
  _p_->tcorr       = 0;                   \
  _p_->tmax        = 0

// ---------- Functions --------------------------------------------------------
int distindex(int, int, int);

int actionTest(actionInformation *, double *, double *, double *, double *, int);
int transformAtomicFluct(actionInformation *, double *, double *, double *, double *, int);
int transformAtomicFluct3D(actionInformation *, double *, double *, double *, double *, int);
int transformCheckOverlap(actionInformation *, double *, double *, double *, double *, int);
//int transformCluster(actionInformation *, double *, double *, double *, double *, int);
//int transformClusterAttribute(actionInformation*,double*,double*,double*, double *, int);
int transformContacts(actionInformation *, double *, double *, double *, double *, int);
int transformCorr(actionInformation *, double *, double *, double *, double *, int);
int transformDiffusion(actionInformation *, double *, double *, double *, double *, int);
int transformDihedralCluster(actionInformation*,double*,double *, double *, double *, int);
int transformDipole(actionInformation *, double *, double *, double *, double *, int);
int transformDNAiontracker(actionInformation *, double *, double *, double *, double *, int);
int transformEcho(actionInformation *, double *, double *, double *, double *, int);
//int transformEnergy(actionInformation *, double *, double *, double *, double *, int);
int transformGrid(actionInformation *, double *, double *, double *, double *, int);
int transformMatrix(actionInformation *, double *, double *, double *, double *, int);
int transformPrincipal(actionInformation *, double *, double *, double *, double *, int);
int transformProjection(actionInformation *, double *, double *, double *, double *, int);
int transformRandomizeIons(actionInformation *, double *, double *, double *, double *, int);
int transformRunningAverage(actionInformation*,double*,double *, double *, double *, int);
int transformScale(actionInformation *, double *, double *, double *, double *, int);
//int transformTruncOct(actionInformation *, double *, double *, double *, double *, int);
int transformUnwrap(actionInformation *, double *, double *, double *, double *, int);
int transformVector(actionInformation *, double *, double *, double *, double *, int);
int transformWatershell(actionInformation *, double *, double *, double *, double *, int);
actionInformation* ptrajCopyAction(actionInformation**);

#ifdef __cplusplus
}
#endif
#endif
