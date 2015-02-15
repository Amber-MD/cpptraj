// Parm_CharmmPsf.cpp
#include <cstdio> // sscanf
#include <cstring> // strncmp
#include "Parm_CharmmPsf.h"
#include "CpptrajStdio.h"

bool Parm_CharmmPsf::ID_ParmFormat(CpptrajFile& fileIn) {
  // Assumes already set up
  if (fileIn.OpenFile()) return false;
  std::string nextLine = fileIn.GetLine();
  if (nextLine.empty()) return false;
  bool isPSF = ( nextLine.compare(0, 3, "PSF") == 0 );
  fileIn.CloseFile();
  return isPSF;
}

int Parm_CharmmPsf::FindTag(char* tag, const char* target, int tgtsize, CpptrajFile& infile) {
  int nval = 0;
  while (strncmp(tag,target,tgtsize)!=0) {
    const char* buffer = infile.NextLine();
    if ( buffer == 0 ) return 0;
    sscanf(buffer,"%i %10s",&nval,tag);
  }
  return nval;
}

// Parm_CharmmPsf::ReadParm()
/** Open the Charmm PSF file specified by filename and set up topology data.
  * Mask selection requires natom, nres, names, resnames, resnums.
  */
int Parm_CharmmPsf::ReadParm(std::string const& fname, Topology &parmOut) {
  const size_t TAGSIZE = 10; 
  char tag[TAGSIZE];
  tag[0]='\0';

  CpptrajFile infile;
  if (infile.OpenRead(fname)) return 1;
  mprintf("    Reading Charmm PSF file %s as topology file.\n",infile.Filename().base());
  // Read the first line, should contain PSF...
  const char* buffer = 0;
  if ( (buffer=infile.NextLine()) == 0 ) return 1;
  // Advance to <ntitle> !NTITLE
  int ntitle = FindTag(tag, "!NTITLE", 7, infile); 
  // Only read in 1st title
  std::string psftitle;
  if (ntitle > 0) {
    buffer = infile.NextLine();
    psftitle.assign( buffer );
  }
  parmOut.SetParmName( psftitle, infile.Filename() );
  // Advance to <natom> !NATOM
  int natom = FindTag(tag, "!NATOM", 6, infile);
  mprintf("\tPSF: !NATOM tag found, natom=%i\n", natom);
  // If no atoms, probably issue with PSF file
  if (natom < 1) {
    mprinterr("Error: No atoms in PSF file.\n");
    return 1;
  }
  // Read the next natom lines
  int psfresnum = 0;
  char psfresname[6];
  char psfname[6];
  char psftype[6];
  double psfcharge;
  double psfmass;
  for (int atom=0; atom < natom; atom++) {
    if ( (buffer=infile.NextLine()) == 0 ) {
      mprinterr("Error: ReadParmPSF(): Reading atom %i\n",atom+1);
      return 1;
    }
    // Read line
    // ATOM# SEGID RES# RES ATNAME ATTYPE CHRG MASS (REST OF COLUMNS ARE LIKELY FOR CMAP AND CHEQ)
    sscanf(buffer,"%*i %*s %i %s %s %s %lf %lf",&psfresnum, psfresname, 
           psfname, psftype, &psfcharge, &psfmass);
    parmOut.AddTopAtom( Atom( psfname, psfcharge, psfmass, psftype), psfresnum, psfresname, 0 );
  } // END loop over atoms 
  // Advance to <nbond> !NBOND
  int bondatoms[8];
  int nbond = FindTag(tag, "!NBOND", 6, infile);
  if (nbond > 0) {
    mprintf("\tPSF: !NBOND tag found, nbond=%i\n", nbond);
    int nlines = nbond / 4;
    if ( (nbond % 4) != 0) nlines++;
    for (int bondline=0; bondline < nlines; bondline++) {
      if ( (buffer=infile.NextLine()) == 0 ) {
        mprinterr("Error: ReadParmPSF(): Reading bond line %i\n",bondline+1);
        return 1;
      }
      // Each line has 4 pairs of atom numbers
      int nbondsread = sscanf(buffer,"%i %i %i %i %i %i %i %i",bondatoms,bondatoms+1,
                              bondatoms+2,bondatoms+3, bondatoms+4,bondatoms+5,
                              bondatoms+6,bondatoms+7);
      // NOTE: Charmm atom nums start from 1
      for (int bondidx=0; bondidx < nbondsread; bondidx+=2)
        parmOut.AddBond(bondatoms[bondidx]-1, bondatoms[bondidx+1]-1);
    }
  } else
    mprintf("Warning: PSF has no bonds.\n");
  // Advance to <nangles> !NTHETA
  int nangle = FindTag(tag, "!NTHETA", 7, infile);
  if (nangle > 0) {
    mprintf("\tPSF: !NTHETA tag found, nangle=%i\n", nangle);
    int nlines = nangle / 3;
    if ( (nangle % 3) != 0) nlines++;
    for (int angleline=0; angleline < nlines; angleline++) {
      if ( (buffer=infile.NextLine()) == 0) {
        mprinterr("Error: Reading angle line %i\n", angleline+1);
        return 1;
      }
      // Each line has 3 groups of atom numbers
      int nanglesread = sscanf(buffer,"%i %i %i %i %i %i %i %i %i",bondatoms,bondatoms+1,
                              bondatoms+2,bondatoms+3, bondatoms+4,bondatoms+5,
                              bondatoms+6,bondatoms+7, bondatoms+8);
      for (int angleidx=0; angleidx < nanglesread; angleidx+=3) {
        parmOut.AddAngle( bondatoms[angleidx  ]-1,
                          bondatoms[angleidx+1]-1,
                          bondatoms[angleidx+2]-1 );
      }
    }
  } else
    mprintf("Warning: PSF has no angles.\n");
  // Advance to <ndihedrals> !NPHI
  int ndihedral = FindTag(tag, "!NPHI", 5, infile);
  if (ndihedral > 0) {
    mprintf("\tPSF: !NPHI tag found, ndihedral=%i\n", ndihedral);
    int nlines = ndihedral / 2;
    if ( (ndihedral % 2) != 0) nlines++;
    for (int dihline = 0; dihline < nlines; dihline++) {
      if ( (buffer=infile.NextLine()) == 0) {
        mprinterr("Error: Reading dihedral line %i\n", dihline+1);
        return 1;
      }
      // Each line has 2 groups of atom numbers
      int ndihread = sscanf(buffer,"%i %i %i %i %i %i %i %i",bondatoms,bondatoms+1,
                              bondatoms+2,bondatoms+3, bondatoms+4,bondatoms+5,
                              bondatoms+6,bondatoms+7);
      for (int dihidx=0; dihidx < ndihread; dihidx += 4) {
        parmOut.AddDihedral( bondatoms[dihidx  ]-1,
                             bondatoms[dihidx+1]-1,
                             bondatoms[dihidx+2]-1,
                             bondatoms[dihidx+3]-1 );
      }
    }
  } else
    mprintf("Warning: PSF has no dihedrals.\n");
  mprintf("\tPSF contains %i atoms, %i residues.\n", parmOut.Natom(), parmOut.Nres());

  infile.CloseFile();

  return 0;
}

int Parm_CharmmPsf::WriteParm(std::string const& fname, Topology const& parm) {
  // TODO: CMAP etc info
  CpptrajFile outfile;
  if (outfile.OpenWrite(fname)) return 1;
  // Write PSF
  outfile.Printf("PSF\n\n");
  // Write title
  std::string titleOut = parm.ParmName();
  titleOut.resize(78);
  outfile.Printf("%8i !NTITLE\n* %-78s\n\n", 1, titleOut.c_str());
  // Write NATOM section
  outfile.Printf("%8i !NATOM\n", parm.Natom());
  unsigned int idx = 1;
  // Make fake segment ids for now.
  char segid[2];
  segid[0] = 'A';
  segid[1] = '\0';
  mprintf("Warning: Assigning single letter segment IDs.\n");
  int currentMol = 0;
  bool inSolvent = false;
  for (Topology::atom_iterator atom = parm.begin(); atom != parm.end(); ++atom, ++idx) {
    int resnum = atom->ResNum();
    if (atom->MolNum() != currentMol) {
      if (!inSolvent) {
        inSolvent = parm.Mol(atom->MolNum()).IsSolvent();
        currentMol = atom->MolNum();
        segid[0]++;
      } else
        inSolvent = parm.Mol(atom->MolNum()).IsSolvent();
    }
    // TODO: Print type name for xplor-like PSF
    outfile.Printf("%8i %-4s %4i %4s %4s %4i %14.6G %14.6g %8i\n", idx, segid,
                   parm.Res(resnum).OriginalResNum(), parm.Res(resnum).c_str(),
                   atom->c_str(), atom->TypeIndex()+1, atom->Charge(),
                   atom->Mass(), 0);
  }
  outfile.Printf("\n");
  // Write NBOND section
  outfile.Printf("%8u !NBOND: bonds\n", parm.Bonds().size() + parm.BondsH().size());
  idx = 1;
  for (BondArray::const_iterator bond = parm.BondsH().begin();
                                 bond != parm.BondsH().end(); ++bond, ++idx)
  {
    outfile.Printf("%8i%8i", bond->A1()+1, bond->A2()+1);
    if ((idx % 4)==0) outfile.Printf("\n"); 
  }
  for (BondArray::const_iterator bond = parm.Bonds().begin();
                                 bond != parm.Bonds().end(); ++bond, ++idx)
  {
    outfile.Printf("%8i%8i", bond->A1()+1, bond->A2()+1);
    if ((idx % 4)==0) outfile.Printf("\n"); 
  }
  if ((idx % 4)!=0) outfile.Printf("\n");
  outfile.Printf("\n");

  outfile.CloseFile();
  return 0;
}
