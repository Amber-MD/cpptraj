#ifndef INC_PTRAJ_ACTIONS_H
#define INC_PTRAJ_ACTIONS_H
#ifdef __cplusplus
extern "C" {
#endif
/*! \file ptraj_actions.h
    \brief Actions originally from Ptraj

 This file has ptraj actions originally from actions.c that have not yet
 been incorporated into the cpptraj framework. They have been modified
 to use Cpptraj mask, distance, and torsion routines and are called from
 Cpptraj via the PtrajAction class. 
 Ptraj originally written by Thomas E. Cheatham III et al.
 Adapted for cpptraj by Dan Roe, Rutgers, November 2001.
 See $AMBERHOME/AmberTools/src/ptraj/contributors.h for more author info.
 */

// INCLUDES
#include "ptraj_state.h" // ptrajState and mask parsing routines, modes

// ---------- Routines that access global ptraj_actions vars -------------------
void SetReferenceInfo(double*,int);
void FreeReferenceInfo();

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
int actionTest(actionInformation *, double *, double *, double *, double *, int);
int transformAtomicFluct(actionInformation *, double *, double *, double *, double *, int);
int transformAtomicFluct3D(actionInformation *, double *, double *, double *, double *, int);
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
