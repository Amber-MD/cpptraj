#ifdef __cplusplus
extern "C" {
#endif
// This file was adapted from ptraj/actions.h for cpptraj
// Dan Roe, Rutgers, November 2011

#include "ptraj_stack.h" // stackType
#include "ptraj_state.h" // ptrajState

/*
 *  This is the header file for analyze.c which contains the basic structures
 *  necessary for analyzing data acculated by the various actions in ptraj()
 *
 * (1) EXTERNALLY VISIBLE DEFINITIONS
 * (2) GLOBAL VARIABLES
 * (3) EXTERNALLY VISIBLE FUNCTION PROTOTYPES
 * (4) LOCAL STRUCTURES
 * (5) LOCALLY VISIBLE FUNCTION PROTOTYPES
 */


/*
 * (1) EXTERNALLY VISIBLE DEFINITIONS
 */


  /*
   *  ANALYZE INFORMATION
   */

typedef enum _analyzeType { 
  ANALYZE_NOOP, 
  ANALYZE_CORRELATIONCOEFFICIENT,
  ANALYZE_CRANKSHAFT,
  ANALYZE_HBOND,
  ANALYZE_MATRIX,
  ANALYZE_MODES,
  ANALYZE_SET,
  ANALYZE_STATISTICS,
  ANALYZE_TIMECORR,
  ANALYZE_TEST
} analyzeType;


//#  ifdef __STDC__
typedef int (*analyzeFunction)(void *, stackType *, int);
//#  else
//typedef int (*analyzeFunction)();
//#  endif

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


/*
 * (2) GLOBAL VARIABLES
 */


static char *pucker_ss[10] = 
  { "C3'-endo", "C4'-exo ", "O4'-endo", "C1'-exo ", "C2'-endo",
    "C3'-exo ", "C4'-endo", "O4'-exo ", "C1'-endo", "C2'-exo " };


static char *torsion_ss_2D[6][6] = 
  { {"g+ g+", "g+ a+", "g+ t", "g+ a-", "g+ g-", "g+ c" },
    {"a+ g+", "a+ a+", "a+ t", "a+ a-", "a+ g-", "a+ c" },
    {"t  g+", "t  a+", "t  t", "t  a-", "t  g-", "t  c" },
    {"a- g+", "a- a+", "a- t", "a- a-", "a- g-", "a- c" },
    {"g- g+", "g- a+", "g- t", "g- a-", "g- g-", "g- c" },
    {"c  g+", "c  a+", "c  t", "c  a-", "c  g-", "c  c" } 
  };

static char *torsion_ss[6] =
  { "g+     ", "a+     ", "t      ", "a-     ", "g-     ", "c      " };

static double torsion_offset[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 180.0 };

static char *distance_ss_2D[6][6] = 
  { {"< 2, < 2", "< 2, 2-3", "< 2, 3-4", "< 2, 4-5", "< 2, 5-6", "< 2, > 6" },
    {"2-3, < 2", "2-3, 2-3", "2-3, 3-4", "2-3, 4-5", "2-3, 5-6", "2-3, > 6" },
    {"3-4, < 2", "3-4, 2-3", "3-4, 3-4", "3-4, 4-5", "3-4, 5-6", "3-4, > 6" },
    {"4-5, < 2", "4-5, 2-3", "4-5, 3-4", "4-5, 4-5", "4-5, 5-6", "4-5, > 6" },
    {"5-6, < 2", "5-6, 2-3", "5-6, 3-4", "5-6, 4-5", "5-6, 5-6", "5-6, > 6" },
    {"> 6, < 2", "> 6, 2-3", "> 6, 3-4", "> 6, 4-5", "> 6, 5-6", "> 6, > 6" }
  };

static char *distance_ss[6] = {" < 2.5 ", "2.5-3.5", "3.5-4.5", "4.5-5.5", "5.5-6.5", " > 6.5 " };




/*
 * (3) EXTERNALLY VISIBLE FUNCTION PROTOTYPES
 */


#ifndef ANALYZE_MODULE


#  ifdef __STDC__

extern int analyzeTest(analyzeInformation *, stackType *, int);
extern int analyzeCorrelationCoefficient(analyzeInformation *, stackType *, int);
extern int analyzeCrankshaft(analyzeInformation *, stackType *, int);
extern int analyzeHBond(analyzeInformation *, stackType *, int);
extern int analyzeMatrix(analyzeInformation *, stackType *, int);
extern int analyzeModes(analyzeInformation *, stackType *, int);
extern int analyzeSet(analyzeInformation *, stackType *, int);
extern int analyzeStatistics(analyzeInformation *, stackType *, int);
extern int analyzeTimecorr(analyzeInformation *, stackType *, int);
extern int analyzeTest(analyzeInformation *, stackType *, int);

#  else  /* __STDC__ */

extern int analyzeTest();
extern int analyzeCorrelationCoefficient();
extern int analyzeCrankshaft();
extern int analyzeHBond();
extern int analyzeMatrix();
extern int analyzeModes();
extern int analyzeSet();
extern int analyzeStatistics();
extern int analyzeTimecorr();
extern int analyzeTest();

#  endif /* __STDC__ */


#endif /* ANALYZE_MODULE */



/*
 * (4) LOCAL STRUCTURES
 */

typedef enum _unaryOperator { UNARY_NOOP,     UNARY_NORMALIZE, 
			      UNARY_COSINE,   UNARY_SINE,
			      UNARY_ARCCOS,   UNARY_ARCSIN,
			      UNARY_SHIFT 
} unaryOperator;

typedef enum _binaryOperator { BINARY_NOOP,   BINARY_PLUS,
			       BINARY_MINUS,  BINARY_TIMES,
			       BINARY_DIVIDE
} binaryOperator;

typedef enum _analyzeModesType {
  ANALYZEMODES_UNKNOWN,
  ANALYZEMODES_FLUCT,
  ANALYZEMODES_DISPL,
  ANALYZEMODES_CORR
} analyzeModesType;

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

typedef enum _analyzeTimecorrMode {
  ATCM_UNKNOWN,
  ATCM_AUTO,
  ATCM_CROSS
} analyzeTimecorrMode;

typedef enum _analyzeTimecorrType {
  ATCT_UNKNOWN,
  ATCT_IRED,
  ATCT_NORMAL
} analyzeTimecorrType;

typedef struct _timecorrResults {
  analyzeTimecorrMode mode;
  analyzeTimecorrType type;
  int dplr;
  int norm;
  int drct;
  int relax;
  double tstep;
  double tcorr;
  double distnh;
  double freq;
  char *filename;
  char *noeFilename;
  int ndata;
  double *table;
  double *data1;
  double *data2;
  double *cf;
  double *cfinf;
  double *p2cf;
  double *rcf;
} timecorrResults;

#define INITIALIZE_timecorrResults(_p_) \
  _p_->mode     = ATCM_UNKNOWN;         \
  _p_->type     = ATCT_UNKNOWN;         \
  _p_->dplr     = 0;                    \
  _p_->norm     = 0;                    \
  _p_->drct     = 0;                    \
  _p_->relax    = 0;                    \
  _p_->tstep    = 1.0;                  \
  _p_->tcorr    = 10000.0;              \
  _p_->distnh   = 1.02;                 \
  _p_->freq     = -1.0;                 \
  _p_->filename = NULL;                 \
  _p_->noeFilename = NULL;              \
  _p_->ndata    = 0;                    \
  _p_->table    = NULL;                 \
  _p_->data1    = NULL;                 \
  _p_->data2    = NULL;                 \
  _p_->cf       = NULL;                 \
  _p_->cfinf    = NULL;                 \
  _p_->p2cf     = NULL;                 \
  _p_->rcf      = NULL;

  
/*
 * (5) LOCALLY VISIBLE FUNCTION PROTOTYPES
 */

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
