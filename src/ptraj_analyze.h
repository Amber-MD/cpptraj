#ifndef INC_PTRAJ_ANALYZE_H
#define INC_PTRAJ_ANALYZE_H
#ifdef __cplusplus
extern "C" {
#endif
// This file was adapted from ptraj/actions.h for cpptraj
// Dan Roe, Rutgers, November 2011
// INCLUDES
#include "ptraj_stack.h" // stackType
#include "ptraj_state.h" // ptrajState

// ---------- ptraj analyze routines -------------------------------------------
// Possible analysis types - HBOND disabled
typedef enum _analyzeType { 
  ANALYZE_NOOP, ANALYZE_CORRELATIONCOEFFICIENT, ANALYZE_STATISTICS
} analyzeType;

// ANALYZE FUNCTION TYPE DEFINITION
typedef int (*analyzeFunction)(void *, stackType *, int);

// ANALYZE INFORMATION, passed in as first arg to analyze functions
typedef struct _analyzeInformation { 
  analyzeFunction fxn; 
  analyzeType type;
  ptrajState *state;
  int iarg1;
  int iarg2;
  int iarg3;
  int iarg4;
  double darg1;
  double darg2;
  double darg3;
  double darg4;
  void *carg1;
  void *carg2;
  void *carg3;
  void *carg4;
} analyzeInformation;

#define INITIALIZE_analyzeInformation(_p_) \
  _p_->fxn = NULL;            \
  _p_->type = ANALYZE_NOOP;   \
  _p_->state = NULL;          \
  _p_->iarg1 = 0;             \
  _p_->iarg2 = 0;             \
  _p_->iarg3 = 0;             \
  _p_->iarg4 = 0;             \
  _p_->darg1 = 0.0;           \
  _p_->darg2 = 0.0;           \
  _p_->darg3 = 0.0;           \
  _p_->darg4 = 0.0;           \
  _p_->carg1 = NULL;          \
  _p_->carg2 = NULL;          \
  _p_->carg3 = NULL;          \
  _p_->carg4 = NULL

// ---------- Functions --------------------------------------------------------
int analyzeCorrelationCoefficient(analyzeInformation *, stackType *, int);
int analyzeStatistics(analyzeInformation *, stackType *, int);

#ifdef __cplusplus
}
#endif
#endif
