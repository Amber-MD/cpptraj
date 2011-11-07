#ifndef PTRAJ_ACTIONHEADER_H
#define PTRAJ_ACTIONHEADER_H
// Contains definitions for ptraj-like actions
#ifdef __cplusplus
extern "C" {
#endif

// Redefine Name to match definition in Ptraj
#include "Name.h"
typedef NAME Name;

// argStackType - just used to set up arglist
typedef struct _argStackType {
  int nargs;
  char **arglist;
  char *marked;
} argStackType;

// Ptraj State - originally from ptraj_local.h
typedef struct _ptrajState {
  double box[6];             // box lengths and angles 
  double *masses;            // atom masses 
  double *charges;           // atom charges 
  int atoms;                 // number of atoms 
  int residues;              // number of residues 
  int *ipres;                // atoms in each residue 
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

// possible ptraj actions
typedef enum _actionType {
  TRANSFORM_NOOP,
  TRANSFORM_ACCEPTOR,
  TRANSFORM_ANALYZE,
  TRANSFORM_ANGLE,
  TRANSFORM_ATOMICFLUCT,
  TRANSFORM_ATOMICFLUCT3D,
  TRANSFORM_AVERAGE,
  TRANSFORM_BENCHMARK,
  TRANSFORM_BOX,
  TRANSFORM_CENTER,
  TRANSFORM_CHECKDNA,
  TRANSFORM_CHECKOVERLAP,
  TRANSFORM_CLOSESTWATERS,
  TRANSFORM_CLUSTER,
  TRANSFORM_CLUSTERATTRIBUTE,
  TRANSFORM_CONTACTS,
  TRANSFORM_DIFFUSION,
  TRANSFORM_DIHEDRAL,
  TRANSFORM_DIHEDRALCLUSTER,
  TRANSFORM_DIPOLE,
  TRANSFORM_DISTANCE,
  TRANSFORM_DONOR,
  TRANSFORM_DNAIONTRACKER,
  TRANSFORM_ECHO,
  TRANSFORM_ENERGY,
  TRANSFORM_GRID,
  TRANSFORM_HBOND,
  TRANSFORM_IMAGE,
  TRANSFORM_MATRIX,
  TRANSFORM_PRINCIPAL,
  TRANSFORM_PRNLEV,
  TRANSFORM_PROJECTION,
  TRANSFORM_PUCKER,
  TRANSFORM_RADIAL,
  TRANSFORM_RADIUSOFGYRATION,
  TRANSFORM_RANDOMIZEIONS,
  TRANSFORM_REFERENCE,
  TRANSFORM_RMS,
  TRANSFORM_RUNNINGAVERAGE,
  TRANSFORM_SCALE,
  TRANSFORM_SECONDARYSTRUCT,
  TRANSFORM_STRIP,
  TRANSFORM_SMARTIMAGE,
  TRANSFORM_SOLVENT,
  TRANSFORM_TRAJIN,
  TRANSFORM_TRAJOUT,
  TRANSFORM_TRANSFORM,
  TRANSFORM_TRANSLATE,
  TRANSFORM_TRUNCOCT,
  TRANSFORM_TEST,
  TRANSFORM_UNWRAP,
  TRANSFORM_VECTOR,
  TRANSFORM_CORRELATION,
  TRANSFORM_WATERSHELL,
  TRANSFORM_2DRMS
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

#ifdef __cplusplus
}
#endif
#endif
