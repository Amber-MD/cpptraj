#ifndef INC_PTRAJ_EVEC_H
#define INC_PTRAJ_EVEC_H
// This file was adapted from ptraj/evec.h for cpptraj
// Dan Roe, Rutgers, November 2011
/*
 *  header file for evec.c (Gohlke)
 */

/*
 *  possible mode types
 */
typedef enum _modesType{
  MT_UNKNOWN,
  MT_DIST,
  MT_COVAR,
  MT_MWCOVAR,
  MT_DISTCOVAR,
  MT_CORREL,
  MT_IDEA,
  MT_IRED
} modesType;

/*
 *  possible "sources" for modes
 */
typedef enum _modesSource{
  MS_UNKNOWN,
  MS_STACK,
  MS_FILE
} modesSource;

/*
 *  information relating to modes
 *
 *  eigenvectors are stored in evec as:
 *  [evec(1,1,x),evec(1,1,y),evec(1,1,z), ..., evec(1,n,x),
 *   ...,
 *   evec(n,1,x), ...,                         evec(n,n,x)]
 */

typedef struct _modesInfo {
  char *name;
  modesType type;
  modesSource source;
  int navgelem;
  double *avg;
  int nvect;
  int nvectelem;
  double *freq;
  double *evec;
} modesInfo;

#define INITIALIZE_modesInfo(_p_) \
  _p_->name      = NULL;          \
  _p_->type      = MT_UNKNOWN;    \
  _p_->source    = MS_UNKNOWN;    \
  _p_->navgelem  = 0;             \
  _p_->avg       = NULL;          \
  _p_->nvect     = 0;             \
  _p_->nvectelem = 0;             \
  _p_->freq      = NULL;          \
  _p_->evec      = NULL;

//#ifndef EVEC_MODULE
//#  ifdef __STDC__

int readEvecFile(FILE *, int, int, modesInfo *);

//#  else

//extern int readEvecFile();

//#  endif
//#endif

#endif
