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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/evec.c,v 10.0 2008/04/15 23:24:11 case Exp $
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


#include <stdio.h>
#include <string.h>
#include <math.h>

#define EVEC_MODULE
#include "ptraj_evec.h"
#include "ptraj_common.h"

/*
 *  This source file contains routines for I/O of eigenvector files
 *    (Holger Gohlke, Scripps)
 */

/*-----------------------------------------------------------------------*/

   int
readEvecFile(FILE *fp, int ibeg, int iend, modesInfo *modinfo){

  char buffer[BUFFER_SIZE];
  int i, j, nlines, nent, nrest;
  int nno, nvec, ncoords;
  float tmpval;
  double *avg, *freq, *evec;
  char *tmpbuf;

  /*
   * Read title line 
   */
  if(fgets(buffer, BUFFER_SIZE, fp) == NULL){
    fprintf(stderr,
            "WARNING in ptraj(), readEvecFile: error while reading title\n");
    return 1;
  }
  modinfo->type = (strstr(buffer, "DISTCOVAR") != NULL ? MT_DISTCOVAR :
                   (strstr(buffer, "MWCOVAR") != NULL ? MT_MWCOVAR :
                    (strstr(buffer, "COVAR") != NULL ? MT_COVAR :
                     (strstr(buffer, "DIST") != NULL ? MT_DIST :
                      (strstr(buffer, "CORREL") != NULL ? MT_CORREL :
                       (strstr(buffer, "IDEA") != NULL ? MT_IDEA :
                        (strstr(buffer, "IRED") != NULL ? MT_IRED : MT_UNKNOWN
                        )
                       )
                      )
                     )
                    )
                   )
                  );

  /*
   * For compatibility with quasih and nmode output
   */
  if(modinfo->type == MT_UNKNOWN){
    fprintf(stderr,
            "FYI: No evec type found in %s, assuming it is MWCOVAR\n", modinfo->name);
    modinfo->type = MT_MWCOVAR;
  }

  /* 
   * Read number of coords for avg and evec
   */
  if(fgets(buffer, BUFFER_SIZE, fp) == NULL){
    fprintf(stderr,
            "WARNING in ptraj(), readEvecFile: error while reading number of atoms\n");
    return 1;
  }
  if (sscanf(buffer, "%i", &(modinfo->navgelem)) != 1) {
    fprintf(stderr,
	    "WARNING in ptraj(), readEvecFile: sscanf on coords failed\n");
    return 1;
  }
 switch (sscanf(buffer,"%i %i", &(modinfo->navgelem), &(modinfo->nvectelem))) {
 case 0:
   fprintf(stderr,
    "WARNING in ptraj(), readEvecFile: sscanf on coords failed\n");
   return 1;
 case 1: /* assume the first was read in */
   fprintf(stderr,
    "FYI: No value for nvectelem found in %s, assuming it is navgelem\n", modinfo->name);
   modinfo->nvectelem = modinfo->navgelem;
   break;
 }

  /*
   * Allocate memory for avg, freq, evec
   */
  avg  = (double *) safe_malloc(sizeof(double) * modinfo->navgelem);
  freq = (double *) safe_malloc(sizeof(double) * (iend - ibeg + 1));
  evec = (double *) safe_malloc(sizeof(double) * (iend - ibeg + 1) * modinfo->nvectelem);

  /* 
   * Read average coordinates 
   */
  ncoords = modinfo->navgelem;
  nlines = (int) (ncoords / 7);
  if(ncoords > 0 && ncoords % 7 != 0)
    nlines++;

  nent = 0;
  for(i = 0; i < nlines; i++){
    if(fgets(buffer, BUFFER_SIZE, fp) == NULL){
      fprintf(stderr,
              "WARNING in ptraj(), readEvecFile: error while reading avg coords\n");
      return 1;
    }

    if(strchr(buffer, '*') != NULL ){
      fprintf(stderr,
              "WARNING in ptraj(), readEvecFile: avg coords out of bounds (i.e. ****'s)\n");
      return 1;
    }
     
    tmpbuf = buffer;
    nrest = i < nlines - 1 ? 7 : ncoords - nent;
    for(j = 0; j < nrest; j++){
      if(sscanf(tmpbuf, "%f", &tmpval) != 1){
        fprintf(stderr,
                "WARNING in ptraj(), readEvecFile: error while scanning avg coords\n");
        return 1;
      }
      else{
        tmpbuf += 11;
        avg[nent++] = tmpval;
      }
    }
  }
  modinfo->avg = avg;

  /* 
   * Read eigen vectors
   */
  ncoords = modinfo->nvectelem;
  nlines = (int) (ncoords / 7);
  if(ncoords > 0 && ncoords % 7 != 0)
    nlines++;

  nvec = 0;
  if(fgets(buffer, BUFFER_SIZE, fp) == NULL){
    fprintf(stderr,
            "WARNING in ptraj(), readEvecFile: error while reading eigen vector stars\n");
    return 1;
  }

  while(strstr(buffer, "****") != NULL && feof(fp) == 0){
    /*
     *  Read number and freq
     */
    if(fgets(buffer, BUFFER_SIZE, fp) == NULL){
      fprintf(stderr,
              "WARNING in ptraj(), readEvecFile: error while reading number and freq\n");
      return 1;
    }
    if(sscanf(buffer, "%i%f", &nno, &tmpval) != 2){
      fprintf(stderr,
              "WARNING in ptraj(), readEvecFile: error while scanning number and freq\n");
      return 1;
    }
    if(nno >= ibeg && nno <= iend)
      freq[nvec] = tmpval;
    else if(nno > iend)
      break;

    /*
     *  Read coords
     */
    nent = 0;
    for(i = 0; i < nlines; i++){
      if(fgets(buffer, BUFFER_SIZE, fp) == NULL){
        fprintf(stderr,
                "WARNING in ptraj(), readEvecFile: error while reading evec coords\n");
        return 1;
      }

      if(nno >= ibeg && nno <= iend){
        if(strchr(buffer, '*') != NULL ){
          fprintf(stderr,
                  "WARNING in ptraj(), readEvecFile: evec coords out of bounds (i.e. ****'s)\n");
          return 1;
        }
     
        tmpbuf = buffer;
        nrest = i < nlines - 1 ? 7 : ncoords - nent;
        for(j = 0; j < nrest; j++){
          if(sscanf(tmpbuf, "%f", &tmpval) != 1){
            fprintf(stderr,
                    "WARNING in ptraj(), readEvecFile: error while scanning evec coords\n");
            return 1;
          }
          else{
            tmpbuf += 11;
            evec[ncoords * nvec + nent] = tmpval;
            nent++;
          }
        }
      }
    }

    if(nno >= ibeg && nno <= iend)
      nvec++;

    /*
     *  Read next star line
     */
    if(fgets(buffer, BUFFER_SIZE, fp) == NULL && feof(fp) == 0){
      fprintf(stderr,
              "WARNING in ptraj(), readEvecFile: error while reading eigen vector stars\n");
      return 1;
    }

  } /* end while */
  modinfo->nvect = nvec;
  modinfo->freq = freq;
  modinfo->evec = evec;

  return 0;
}
