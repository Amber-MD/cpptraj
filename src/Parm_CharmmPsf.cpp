// Parm_CharmmPsf.cpp
#include <cstring> // memcpy, strcpy
#include <cstdio> // sscanf
#include "Parm_CharmmPsf.h"
#include "CpptrajStdio.h"

// CharmmPsfParmFile::ReadParm()
/** Open the Charmm PSF file specified by filename and set up topology data.
  * Mask selection requires natom, nres, names, resnames, resnums.
  */
int CharmmPsfParmFile::ReadParm(AmberParm &parmOut, CpptrajFile &parmfile) {
  char buffer[256],tag[256];
  NAME psfname, psfresname;
  int bondatoms[8];
  int currResnum;
  int psfresnum;
  int psfattype;
  int nbond,nlines;

  mprintf("    Reading Charmm PSF file %s as topology file.\n",parmOut.parmName);
  currResnum=-1;
  memset(buffer,' ',256);
  memset(tag,' ',256);
  tag[0]='\0';

  // Read the first line, should contain PSF...
  if (parmfile.IO->Gets(buffer,256)) return 1;
  // Sanity check
  if (buffer[0]!='P' || buffer[1]!='S' || buffer[2]!='F') {
    mprinterr("Error: ReadParmPSF(): Could not read Charmm PSF file.\n");
    return 1;
  }
  // Advance to <natom> !NATOM
  while (strncmp(tag,"!NATOM",6)!=0) {
    if (parmfile.IO->Gets(buffer,256)) return 1;
    sscanf(buffer,"%i %s",&parmOut.natom,tag);
  }
  mprintf("\tPSF: !NATOM tag found, natom=%i\n",parmOut.natom);
  // If no atoms, probably issue with PSF file
  if (parmOut.natom<=0) {
    mprintf("Error: No atoms in PSF file.\n");
    return 1;
  }

  // Allocate memory for atom name, charge, mass.
  parmOut.names=(NAME*)    new NAME[ parmOut.natom ];
  parmOut.mass=(double*)   new double[ parmOut.natom ];
  parmOut.charge=(double*) new double[ parmOut.natom ];

  // Read the next natom lines
  for (int atom=0; atom < parmOut.natom; atom++) {
    if (parmfile.IO->Gets(buffer,256) ) {
      mprinterr("Error: ReadParmPSF(): Reading atom %i\n",atom+1);
      return 1;
    }
    // Detect and remove trailing newline
    //bufferLen = strlen(buffer);
    //if (buffer[bufferLen-1] == '\n') buffer[bufferLen-1]='\0';
    // Read line
    // ATOM# SEGID RES# RES ATNAME ATTYPE CHRG MASS (REST OF COLUMNS ARE LIKELY FOR CMAP AND CHEQ)
    sscanf(buffer,"%*8i %*4s %i %4s %4s %4i %14lf %14lf",&psfresnum,psfresname,psfname,
           &psfattype,parmOut.charge+atom,parmOut.mass+atom);
    // Ensure name has 4 chars
    PadWithSpaces( psfname );
    strcpy(parmOut.names[atom],psfname);
    // If this residue number is different than the last, allocate mem for new res
    if (currResnum!=psfresnum) {
        NAME *temprname = new NAME[ parmOut.nres+1 ];
        memcpy(temprname, parmOut.resnames, parmOut.nres * sizeof(NAME));
        delete[] parmOut.resnames;
        parmOut.resnames = temprname;
        // Ensure resname has 4 chars
        PadWithSpaces( psfresname );
        strcpy(parmOut.resnames[parmOut.nres],psfresname);
        if (debug>3) mprintf("        PSFRes %i [%s]\n",parmOut.nres,parmOut.resnames[parmOut.nres]);
        int *temprnum = new int[ parmOut.nres+1 ];
        memcpy(temprnum, parmOut.resnums, parmOut.nres * sizeof(int));
        delete[] parmOut.resnums;
        parmOut.resnums = temprnum;
        parmOut.resnums[parmOut.nres]=atom; 
        currResnum=psfresnum;
        parmOut.nres++;
    }
    // Clear the buffer
    memset(buffer,' ',256);
  } // END loop over atoms 

  // Advance to <nbond> !NBOND
  while (strncmp(tag,"!NBOND",6)!=0) {
    if (parmfile.IO->Gets(buffer,256)) return 1;
    sscanf(buffer,"%i %s",&nbond,tag);
  }
  nlines = nbond / 4;
  if ( (nbond % 4) != 0) nlines++;
  for (int bondline=0; bondline < nlines; bondline++) {
    if (parmfile.IO->Gets(buffer,256) ) {
      mprinterr("Error: ReadParmPSF(): Reading bond line %i\n",bondline+1);
      return 1;
    }
    // Each line has 4 pairs of atom numbers
    int nbondsread = sscanf(buffer,"%i %i %i %i %i %i %i %i",bondatoms,bondatoms+1,
                            bondatoms+2,bondatoms+3, bondatoms+4,bondatoms+5,
                            bondatoms+6,bondatoms+7);
    // NOTE: Charmm atom nums start from 1
    for (int bondidx=0; bondidx < nbondsread; bondidx+=2)
      parmOut.AddBond(bondatoms[bondidx]-1,bondatoms[bondidx+1]-1,-1);
  }
  //mprintf("DEBUG: Charmm PSF last line after bond read:\n");
  //mprintf("\t[%s]\n",buffer);
  mprintf("\t%i bonds to hydrogen.\n\t%i bonds to non-hydrogen.\n",parmOut.NbondsWithH,parmOut.NbondsWithoutH);
    
  //boxType = NOBOX;

  //if (debug>0) 
    mprintf("    PSF contains %i atoms, %i residues, %i molecules.\n",
            parmOut.natom,parmOut.nres,parmOut.molecules);

  return 0;
}

