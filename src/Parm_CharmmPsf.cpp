// Parm_CharmmPsf.cpp
#include <cstring> // memcpy, strcpy
#include <cstdio> // sscanf
#include "Parm_CharmmPsf.h"
#include "CpptrajStdio.h"

bool Parm_CharmmPsf::ID_ParmFormat() {
  // Assumes already set up
  if (OpenFile()) return false;
  if (IO->Gets(buffer_, BUF_SIZE_)) return false;
  if (buffer_[0]=='P' && buffer_[1]=='S' && buffer_[2]=='F') {
    CloseFile();
    return true;
  }
  CloseFile();
  return false;
}

// Parm_CharmmPsf::ReadParm()
/** Open the Charmm PSF file specified by filename and set up topology data.
  * Mask selection requires natom, nres, names, resnames, resnums.
  */
int Parm_CharmmPsf::ReadParm(Topology &parmOut) {
  const size_t TAGSIZE = 10; 
  char tag[TAGSIZE];

  if (OpenFile()) return 1;

  mprintf("    Reading Charmm PSF file %s as topology file.\n",parmOut.c_str());
  memset(buffer_,' ',BUF_SIZE_);
  memset(tag,' ',TAGSIZE);
  tag[0]='\0';

  // Read the first line, should contain PSF...
  if (IO->Gets(buffer_,BUF_SIZE_)) return 1;

  // TODO: Assign title
  std::string psftitle;
  parmOut.SetParmName( psftitle, BaseFileStr() );

  // Advance to <natom> !NATOM
  int natom = 0;
  while (strncmp(tag,"!NATOM",6)!=0) {
    if (IO->Gets(buffer_,BUF_SIZE_)) return 1;
    sscanf(buffer_,"%i %10s",&natom,tag);
  }
  mprintf("\tPSF: !NATOM tag found, natom=%i\n", natom);
  // If no atoms, probably issue with PSF file
  if (natom<=0) {
    mprinterr("Error: No atoms in PSF file.\n");
    return 1;
  }

  // Allocate memory for atom name, charge, mass.
  //parmOut.names=(NAME*)    new NAME[ parmOut.natom ];
  //parmOut.mass=(double*)   new double[ parmOut.natom ];
  //parmOut.charge=(double*) new double[ parmOut.natom ];

  // Read the next natom lines
  int psfresnum;
  char psfresname[6];
  char psfname[6];
  char psftype[6];
  double psfcharge;
  double psfmass;
  for (int atom=0; atom < natom; atom++) {
    if (IO->Gets(buffer_,BUF_SIZE_) ) {
      mprinterr("Error: ReadParmPSF(): Reading atom %i\n",atom+1);
      return 1;
    }
    // Detect and remove trailing newline
    //bufferLen = strlen(buffer);
    //if (buffer[bufferLen-1] == '\n') buffer[bufferLen-1]='\0';
    // Read line
    // ATOM# SEGID RES# RES ATNAME ATTYPE CHRG MASS (REST OF COLUMNS ARE LIKELY FOR CMAP AND CHEQ)
    sscanf(buffer_,"%*i %*s %i %s %s %s %lf %lf",&psfresnum, psfresname, 
           psfname, psftype, &psfcharge, &psfmass);
    parmOut.AddAtom( Atom( psfname, psfcharge, 0, psfmass, 0, psftype, 0, 0, psfresnum),
                     Residue(psfresnum, psfresname), NULL );
    // Clear the buffer
    //memset(buffer,' ',256);
  } // END loop over atoms 

  // Advance to <nbond> !NBOND
  int nbond = 0;
  int bondatoms[8];
  while (strncmp(tag,"!NBOND",6)!=0) {
    if (IO->Gets(buffer_,BUF_SIZE_)) return 1;
    sscanf(buffer_,"%i %10s",&nbond,tag);
  }
  int nlines = nbond / 4;
  if ( (nbond % 4) != 0) nlines++;
  for (int bondline=0; bondline < nlines; bondline++) {
    if (IO->Gets(buffer_,BUF_SIZE_) ) {
      mprinterr("Error: ReadParmPSF(): Reading bond line %i\n",bondline+1);
      return 1;
    }
    // Each line has 4 pairs of atom numbers
    int nbondsread = sscanf(buffer_,"%i %i %i %i %i %i %i %i",bondatoms,bondatoms+1,
                            bondatoms+2,bondatoms+3, bondatoms+4,bondatoms+5,
                            bondatoms+6,bondatoms+7);
    // NOTE: Charmm atom nums start from 1
    for (int bondidx=0; bondidx < nbondsread; bondidx+=2)
      parmOut.AddBond(bondatoms[bondidx]-1, bondatoms[bondidx+1]-1);
  }
  //mprintf("DEBUG: Charmm PSF last line after bond read:\n");
  //mprintf("\t[%s]\n",buffer);
  /*mprintf("\t%i bonds to hydrogen.\n\t%i bonds to non-hydrogen.\n",
          parmOut.NbondsWithH,parmOut.NbondsWithoutH);*/
    
  //boxType = NOBOX;

  //if (debug_>0) 
    mprintf("    PSF contains %i atoms, %i residues.\n",
            parmOut.Natom(), parmOut.Nres());

  CloseFile();

  return 0;
}

