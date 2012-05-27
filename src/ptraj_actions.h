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

// ---------- Ptraj actions ----------------------------------------------------
// possible ptraj actions
typedef enum _actionType {
  TRANSFORM_NOOP, TRANSFORM_ANALYZE, TRANSFORM_TRUNCOCT 
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

// ---------- Functions --------------------------------------------------------
//int transformTruncOct(actionInformation *, double *, double *, double *, double *, int);

#ifdef __cplusplus
}
#endif
#endif
