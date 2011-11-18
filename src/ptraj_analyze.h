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
  ANALYZE_NOOP, ANALYZE_CORRELATIONCOEFFICIENT, ANALYZE_CRANKSHAFT,
  ANALYZE_HBOND, ANALYZE_MATRIX, ANALYZE_MODES, ANALYZE_SET,
  ANALYZE_STATISTICS, ANALYZE_TIMECORR, ANALYZE_TEST
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
int analyzeTest(analyzeInformation *, stackType *, int);
int analyzeCorrelationCoefficient(analyzeInformation *, stackType *, int);
int analyzeCrankshaft(analyzeInformation *, stackType *, int);
//int analyzeHBond(analyzeInformation *, stackType *, int);
int analyzeModes(analyzeInformation *, stackType *, int);
int analyzeSet(analyzeInformation *, stackType *, int);
int analyzeStatistics(analyzeInformation *, stackType *, int);
#ifndef NO_PTRAJ_ANALYZE
int analyzeMatrix(analyzeInformation *, stackType *, int);
int analyzeTimecorr(analyzeInformation *, stackType *, int);
#endif
int analyzeTest(analyzeInformation *, stackType *, int);

  
// ---------- LOCALLY VISIBLE FUNCTION PROTOTYPES ------------------------------

#ifdef HPUX
  void dsaupd(int *, char *, int *, char *, int *, double *, double*, 
#elif defined UNICOS
  void DSAUPD(int *, char *, int *, char *, int *, double *, double*, 
#else
  void dsaupd_(int *, char *, int *, char *, int *, double *, double*, 
#endif
               int *, double *, int *, int *, int *, double *, double *,
               int *, int *);

#ifdef HPUX
  void dseupd(int *, char *, int *, double *, double *, int *, double *,
#elif defined UNICOS
  void DSEUPD(int *, char *, int *, double *, double *, int *, double *,
#else
  void dseupd_(int *, char *, int *, double *, double *, int *, double *,
#endif
               char *, int *, char *, int*, double *, double *,
               int *, double *, int *, int *, int *, double *, double *,
               int *, int *);

#ifdef HPUX
  void thermo(int *, int *, int *, double *, double *, double *, 
#elif defined UNICOS
  void THERMO(int *, int *, int *, double *, double *, double *, 
#else
  void thermo_(int *, int *, int *, double *, double *, double *, 
#endif
               double *, double *, double *, double *,
               double *, double *);

#ifdef __cplusplus
}
#endif
#endif
