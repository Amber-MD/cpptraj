#ifndef INC_PTRAJ_SCALAR_H
#define INC_PTRAJ_SCALAR_H

#include "ptraj_state.h" // ptrajState: scalarInfo, transformMatrixInfo
#include "ptraj_evec.h"  // modesInfo: transformVectorInfo
#include "ptraj_stack.h" // stackType
  /*
   *  Information for the stackStack, i.e. saving various 
   *  scalar values for later processing
   */

typedef enum _scalarMode { SCALAR_NULL,
                           SCALAR_DISTANCE,
                           SCALAR_ANGLE,
                           SCALAR_TORSION,
                           SCALAR_PUCKER,
                           SCALAR_RMS
} scalarMode;

typedef enum _scalarType {
  SCALAR_TYPE_UNDEFINED,
  SCALAR_TYPE_ALPHA,
  SCALAR_TYPE_BETA,
  SCALAR_TYPE_DELTA,
  SCALAR_TYPE_GAMMA,
  SCALAR_TYPE_EPSILON,
  SCALAR_TYPE_ZETA,
  SCALAR_TYPE_PUCKER,
  SCALAR_TYPE_CHI,
  SCALAR_TYPE_H1P,
  SCALAR_TYPE_C2P,
  SCALAR_TYPE_PHI,
  SCALAR_TYPE_PSI,
  SCALAR_TYPE_PCHI,
  SCALAR_TYPE_HBOND,
  SCALAR_TYPE_NOE
} scalarType;

typedef struct _scalarInfo {
  char *name;
  char *filename;
  scalarMode mode;
  scalarType type;
  int atom1, atom2, atom3, atom4, atom5;
  int *mask1, *mask2, *mask3, *mask4, *mask5;
  int totalFrames;
  int frame;
  double mean;
  double stddev;
  double max;
  double min;
  double bound;
  double boundh;
  ptrajState *state;
  double *value;
  double *cos;
  double *sin;
  void *results;
  //actionInformation *action;
} scalarInfo;

#define INITIALIZE_scalarInfo(_p_) \
  _p_->name      = NULL; \
  _p_->filename  = NULL; \
  _p_->mode      = SCALAR_NULL; \
  _p_->type      = SCALAR_TYPE_UNDEFINED; \
  _p_->atom1 = -1; _p_->atom2 = -1; _p_->atom3 = -1; _p_->atom4 = -1; _p_->atom5 = -1; \
  _p_->mask1=NULL; _p_->mask2=NULL; _p_->mask3=NULL; _p_->mask4=NULL; _p_->mask5=NULL; \
  _p_->totalFrames = 0;  \
  _p_->frame     = 0;    \
  _p_->mean      = 0.0; \
  _p_->stddev     = 0.0; \
  _p_->max       = 0.0; \
  _p_->min       = 0.0; \
  _p_->bound     = 0.0;    \
  _p_->boundh    = 0.0;    \
  _p_->state     = NULL; \
  _p_->value     = NULL; \
  _p_->cos      = NULL; \
  _p_->sin      = NULL; \
  _p_->results   = NULL
//  _p_->action    = NULL; 
// -------------------------------------
scalarInfo *ptrajCopyScalar(scalarInfo **scalarinp);
scalarInfo *scalarStackGetName(stackType **, char *);

#endif
