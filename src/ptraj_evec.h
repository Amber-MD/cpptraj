/*  _______________________________________________________________________
 *
 *                        RDPARM/PTRAJ: 2008
 *  _______________________________________________________________________
 *
 *  This file is part of rdparm/ptraj.
 *
 *  rdparm/ptraj is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  rdparm/ptraj is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You can receive a copy of the GNU General Public License from
 *  http://www.gnu.org or by writing to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *  ________________________________________________________________________
 *
 *  CVS tracking:
 *
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/evec.h,v 10.0 2008/04/15 23:24:11 case Exp $
 *
 *  Revision: $Revision: 10.0 $
 *  Date: $Date: 2008/04/15 23:24:11 $
 *  Last checked in by $Author: case $
 *  ________________________________________________________________________
 *
 *
 *  CONTACT INFO: To learn who the code developers are, who to contact for
 *  more information, and to know what version of the code this is, refer
 *  to the CVS information and the include files (contributors.h && version.h)
 *
 */

/*  ________________________________________________________________________
 */



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

#ifndef EVEC_MODULE
#  ifdef __STDC__

extern int readEvecFile(FILE *, int, int, modesInfo *);

#  else

extern int readEvecFile();

#  endif
#endif
