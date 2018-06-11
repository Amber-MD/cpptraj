// Parm_CharmmPsf.cpp
#include <cstdio> // sscanf
#include <cstring> // strncmp
#include <cctype> // isdigit
#include "Parm_CharmmPsf.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include "Mol.h" // UniqueCount()
#include "CharmmParamFile.h"

// Parm_CharmmPsf::ID_ParmFormat()
bool Parm_CharmmPsf::ID_ParmFormat(CpptrajFile& fileIn) {
  // Assumes already set up
  if (fileIn.OpenFile()) return false;
  std::string nextLine = fileIn.GetLine();
  if (nextLine.empty()) return false;
  bool isPSF = ( nextLine.compare(0, 3, "PSF") == 0 );
  fileIn.CloseFile();
  return isPSF;
}

// Parm_CharmmPsf::ReadHelp()
void Parm_CharmmPsf::ReadHelp() {
  mprintf("\tparam <file> : Read CHARMM parameters from given file. Can do multiple times.\n");
}

// Parm_CharmmPsf::processReadArgs()
int Parm_CharmmPsf::processReadArgs(ArgList& argIn) {
  int err = 0;
  // Read CHARMM parameters if specified.
  std::string parFileName = argIn.GetStringKey("param");
  while (!parFileName.empty()) {
    CharmmParamFile infile;
    mprintf("\tReading CHARMM parameters from '%s'\n", parFileName.c_str());
    err += infile.ReadParams(params_, parFileName, debug_ );
    parFileName = argIn.GetStringKey("param");
  }
  return err;
}

// Parm_CharmmPsf::FindTag()
int Parm_CharmmPsf::FindTag(char* tag, const char* target, CpptrajFile& infile) {
  int nval = 0;
  int tgtsize = strlen( target );
  while (strncmp(tag,target,tgtsize)!=0) {
    const char* buffer = infile.NextLine();
    if ( buffer == 0 ) return 0;
    sscanf(buffer,"%i %10s",&nval,tag);
  }
  return nval;
}

//  Parm_CharmmPsf::ReadDihedrals()
int Parm_CharmmPsf::ReadDihedrals(CpptrajFile& infile, int ndihedral, const char* typestr, Topology& parmOut) const
{
    bool found;
    int bondatoms[8];
    const char* buffer = 0; 
    int nlines = ndihedral / 2;
    if ( (ndihedral % 2) != 0) nlines++;
    for (int dihline = 0; dihline < nlines; dihline++) {
      if ( (buffer=infile.NextLine()) == 0) {
        mprinterr("Error: Reading %s line %i\n", typestr, dihline+1);
        return 1;
      }
      // Each line has 2 groups of 4 atom numbers
      int ndihread = sscanf(buffer,"%i %i %i %i %i %i %i %i",bondatoms,bondatoms+1,
                              bondatoms+2,bondatoms+3, bondatoms+4,bondatoms+5,
                              bondatoms+6,bondatoms+7);
      if (params_.DP().empty())
        for (int dihidx=0; dihidx < ndihread; dihidx += 4)
          // TODO: Determine end dihedrals
          parmOut.AddDihedral( DihedralType(bondatoms[dihidx  ]-1,
                                            bondatoms[dihidx+1]-1,
                                            bondatoms[dihidx+2]-1,
                                            bondatoms[dihidx+3]-1,
                                            DihedralType::NORMAL) );
      else {
        for (int dihidx=0; dihidx < ndihread; dihidx += 4) {
          int a1 = bondatoms[dihidx]-1;
          int a2 = bondatoms[dihidx+1]-1;
          int a3 = bondatoms[dihidx+2]-1;
          int a4 = bondatoms[dihidx+3]-1;
          DihedralType dih(a1, a2, a3, a4, DihedralType::NORMAL);
          AtomTypeHolder types(4);
          types.AddName( parmOut[a1].Type() );
          types.AddName( parmOut[a2].Type() );
          types.AddName( parmOut[a3].Type() );
          types.AddName( parmOut[a4].Type() );
          DihedralParmType dpt;
          if (typestr[0] == 'd')
            dpt = params_.DP().FindParam( types, found );
          else
            dpt = params_.IP().FindParam( types, found );
          if (found) {
            if (typestr[0] == 'd')
              parmOut.AddDihedral( dih, dpt );
            else
              parmOut.AddCharmmImproper( dih, dpt );
          } else {
            mprintf("Warning: Parameters not found for %s %s - %s - %s - %s\n", typestr, parmOut.AtomMaskName(a1).c_str(), parmOut.AtomMaskName(a2).c_str(), parmOut.AtomMaskName(a3).c_str(), parmOut.AtomMaskName(a4).c_str());
            if (typestr[0] == 'd')
              parmOut.AddDihedral( dih );
            else
              parmOut.AddCharmmImproper( dih );
          }
        }
      }
    }
  return 0;
}

// Parm_CharmmPsf::ReadParm()
/** Open the Charmm PSF file specified by filename and set up topology data.
  * Mask selection requires natom, nres, names, resnames, resnums.
  */
int Parm_CharmmPsf::ReadParm(FileName const& fname, Topology &parmOut) {
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
  int ntitle = FindTag(tag, "!NTITLE", infile); 
  // Only read in 1st title. Skip any asterisks.
  std::string psftitle;
  if (ntitle > 0) {
    buffer = infile.NextLine();
    const char* ptr = buffer;
    while (*ptr != '\0' && (*ptr == ' ' || *ptr == '*')) ++ptr;
    psftitle.assign( ptr );
  }
  parmOut.SetParmName( NoTrailingWhitespace(psftitle), infile.Filename() );
  // Advance to <natom> !NATOM
  int natom = FindTag(tag, "!NATOM", infile);
  if (debug_>0) mprintf("\tPSF: !NATOM tag found, natom=%i\n", natom);
  // If no atoms, probably issue with PSF file
  if (natom < 1) {
    mprinterr("Error: No atoms in PSF file.\n");
    return 1;
  }
  bool found; // Used when assigning parameters
  // DEBUG
  //params_.Debug();
  // Read the next natom lines
  int psfresnum = 0;
  char psfresname[9];
  char psfname[9];
  char psftype[9];
  char segmentID[9];
  double psfcharge;
  double psfmass;
  typedef std::vector<std::string> Sarray;
  // TODO AtomTypeArray should eventually be in Topology
  AtomTypeArray atomTypes;
  atomTypes.SetDebug( 1 );
  Sarray SegIDs;
  for (int atom=0; atom < natom; atom++) {
    if ( (buffer=infile.NextLine()) == 0 ) {
      mprinterr("Error: ReadParmPSF(): Reading atom %i\n",atom+1);
      return 1;
    }
    // Read line
    // ATOM# SEGID RES# RES ATNAME ATTYPE CHRG MASS (REST OF COLUMNS ARE LIKELY FOR CMAP AND CHEQ)
    sscanf(buffer,"%*i %s %i %s %s %s %lf %lf", segmentID, &psfresnum, psfresname, 
           psfname, psftype, &psfcharge, &psfmass);
    // Search for segment ID
    int idx = -1;
    for (int i = 0; i != (int)SegIDs.size(); i++)
      if (SegIDs[i].compare( segmentID )==0) {
        idx = i;
        break;
      }
    if (idx == -1) {
      idx = (int)SegIDs.size();
      SegIDs.push_back( segmentID );
      if (debug_>0) mprintf("DEBUG: New segment ID %i '%s'\n", idx, SegIDs.back().c_str());
    }
    // TODO the type index stuff should be in Topology
    int typeidx = atomTypes.CheckForAtomType( NameType(psftype), AtomType(psfmass) );
    Atom chmAtom( psfname, psfcharge, psfmass, psftype );
    chmAtom.SetTypeIndex( typeidx );
    parmOut.AddTopAtom( chmAtom, Residue( psfresname, psfresnum, idx) );
  } // END loop over atoms 
  // Advance to <nbond> !NBOND
  int bondatoms[9];
  int nbond = FindTag(tag, "!NBOND", infile);
  if (nbond > 0) {
    if (debug_>0) mprintf("\tPSF: !NBOND tag found, nbond=%i\n", nbond);
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
      if (params_.BP().empty()) {
        for (int bondidx=0; bondidx < nbondsread; bondidx+=2)
          parmOut.AddBond(bondatoms[bondidx]-1, bondatoms[bondidx+1]-1);
      } else {
        for (int bondidx = 0; bondidx < nbondsread; bondidx += 2) {
          int a1 = bondatoms[bondidx]-1;
          int a2 = bondatoms[bondidx+1]-1;
          AtomTypeHolder types(2);
          types.AddName( parmOut[a1].Type() );
          types.AddName( parmOut[a2].Type() );
          BondParmType bpt = params_.BP().FindParam( types, found );
          if (found)
            parmOut.AddBond( a1, a2, bpt );
          else {
            mprintf("Warning: Parameters not found for bond %s - %s\n", parmOut.AtomMaskName(a1).c_str(), parmOut.AtomMaskName(a2).c_str());
            parmOut.AddBond( a1, a2 );
          }
        }
      }
    }
  } else
    mprintf("Warning: PSF has no bonds.\n");
  // Advance to <nangles> !NTHETA
  int nangle = FindTag(tag, "!NTHETA", infile);
  if (nangle > 0) {
    if (debug_>0) mprintf("\tPSF: !NTHETA tag found, nangle=%i\n", nangle);
    int nlines = nangle / 3;
    if ( (nangle % 3) != 0) nlines++;
    for (int angleline=0; angleline < nlines; angleline++) {
      if ( (buffer=infile.NextLine()) == 0) {
        mprinterr("Error: Reading angle line %i\n", angleline+1);
        return 1;
      }
      // Each line has 3 groups of 3 atom numbers
      int nanglesread = sscanf(buffer,"%i %i %i %i %i %i %i %i %i",bondatoms,bondatoms+1,
                              bondatoms+2,bondatoms+3, bondatoms+4,bondatoms+5,
                              bondatoms+6,bondatoms+7, bondatoms+8);
      if (params_.AP().empty()) {
        for (int angleidx=0; angleidx < nanglesread; angleidx += 3)
          parmOut.AddAngle( bondatoms[angleidx  ]-1,
                            bondatoms[angleidx+1]-1,
                            bondatoms[angleidx+2]-1 );
      } else {
        for (int angleidx=0; angleidx < nanglesread; angleidx += 3) {
          int a1 = bondatoms[angleidx]-1;
          int a2 = bondatoms[angleidx+1]-1;
          int a3 = bondatoms[angleidx+2]-1;
          AtomTypeHolder types(3);
          types.AddName( parmOut[a1].Type() );
          types.AddName( parmOut[a2].Type() );
          types.AddName( parmOut[a3].Type() );
          AngleParmType apt = params_.AP().FindParam( types, found );
          if (found)
            parmOut.AddAngle( a1, a2, a3, apt );
          else {
            mprintf("Warning: Parameters not found for angle %s - %s - %s\n", parmOut.AtomMaskName(a1).c_str(), parmOut.AtomMaskName(a2).c_str(), parmOut.AtomMaskName(a3).c_str());
            parmOut.AddAngle( a1, a2, a3 );
          }
        }
      }
    }
  } else
    mprintf("Warning: PSF has no angles.\n");
  // Advance to <ndihedrals> !NPHI
  int ndihedral = FindTag(tag, "!NPHI", infile);
  if (ndihedral > 0) {
    if (debug_>0) mprintf("\tPSF: !NPHI tag found, ndihedral=%i\n", ndihedral);
    if (ReadDihedrals(infile, ndihedral, "dihedral", parmOut)) return 1;
  } else
    mprintf("Warning: PSF has no dihedrals.\n");
  // Advance to <nimpropers> !NIMPHI
  int nimproper = FindTag(tag, "!NIMPHI", infile);
  if (nimproper > 0) {
    if (debug_ > 0) mprintf("\tPSF: !NIMPHI tag found, nimproper=%i\n", nimproper);
    if (ReadDihedrals(infile, nimproper, "improper", parmOut)) return 1;
  } else
    mprintf("Warning: PSF has no impropers.\n");
  mprintf("\tPSF contains %i atoms, %i residues.\n", parmOut.Natom(), parmOut.Nres());

  infile.CloseFile();

  // Add nonbonded parameters
  if (params_.HasLJparams()) {
    int ntypes = (int)atomTypes.Size();
    parmOut.SetNonbond().SetupLJforNtypes( ntypes );
    mprintf("\tAtom Types:\n");
    for (AtomTypeArray::const_iterator it = atomTypes.begin(); it != atomTypes.end(); ++it)
    {
      int idx1 = it->second;
      int idx0 = params_.AT().AtomTypeIndex( it->first );
      if (idx0 < 0)
        mprintf("Warning: No LJ parameters for type '%s'\n", *(it->first));
      else {
        atomTypes.UpdateType(idx1).SetRadius( params_.AT()[idx0].Radius() );
        atomTypes.UpdateType(idx1).SetDepth( params_.AT()[idx0].Depth() );
      }
      mprintf("\t\t%3i '%s' mass=%10.4f radius=%10.4f depth=%10.4f\n",
              idx1, *(it->first),
              atomTypes[idx1].Mass(),
              atomTypes[idx1].Radius(),
              atomTypes[idx1].Depth());
    }
    mprintf("\tAdding Lennard-Jones parameters using Lorentz-Berthelot combining rules.\n");
    for (AtomTypeArray::const_iterator it1 = atomTypes.begin(); it1 != atomTypes.end(); ++it1)
    {
      int type1 = it1->second;
      for (AtomTypeArray::const_iterator it2 = it1; it2 != atomTypes.end(); it2++)
      {
        int type2 = it2->second;
        NonbondType LJ = atomTypes[type1].Combine_LB( atomTypes[type2] );
        mprintf("\t%3i - %3i : Ri=%10.4f Ei=%10.4f Rj=%10.4f Ej=%10.4f A=%10.4f B=%10.4f\n",
                type1, type2,
                atomTypes[type1].Radius(), atomTypes[type1].Depth(),
                atomTypes[type2].Radius(), atomTypes[type2].Depth(),
                LJ.A(), LJ.B());
        parmOut.SetNonbond().AddLJterm(type1, type2, LJ);
      }
    }
  }

  return 0;
}

static int FindMolType(int molNum, Mol::Marray const& mols) {
  for (Mol::Marray::const_iterator mol = mols.begin(); mol != mols.end(); ++mol)
    for (Mol::Iarray::const_iterator it = mol->idxs_.begin();
                                     it != mol->idxs_.end(); ++it)
      if (*it == molNum) return (mol - mols.begin());
  return -1;
}

int Parm_CharmmPsf::WriteParm(FileName const& fname, Topology const& parm) {
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
  // Make segment ids based on molecule type for now.
  Mol::Marray mols = Mol::UniqueCount(parm);
  mprintf("Warning: Assigning segment IDs based on molecule type.\n");
  int currentMol = 0;
  int currentMtype = FindMolType(currentMol, mols);
  const char* segid = mols[currentMtype].name_.c_str();
  Mol::Iarray::const_iterator mit = mols[currentMtype].idxs_.begin();
//  bool inSolvent = false;
  for (Topology::atom_iterator atom = parm.begin(); atom != parm.end(); ++atom, ++idx) {
    int resnum = atom->ResNum();
    if (atom->MolNum() != currentMol) {
      currentMol = atom->MolNum();
      ++mit;
      if (mit == mols[currentMtype].idxs_.end() || *mit != currentMol) {
        currentMtype = FindMolType(currentMol, mols);
        segid = mols[currentMtype].name_.c_str();
      }
    }
    // TODO: Print type name for xplor-like PSF
    int typeindex = atom->TypeIndex() + 1;
    // If type begins with digit, assume charmm numbers were read as
    // type. Currently Amber types all begin with letters.
    if (isdigit(atom->Type()[0]))
      typeindex = convertToInteger( *(atom->Type()) );
    // ATOM# SEGID RES# RES ATNAME ATTYPE CHRG MASS (REST OF COLUMNS ARE LIKELY FOR CMAP AND CHEQ)
    outfile.Printf("%8i %-4s %-4i %-4s %-4s %4i %14.6G %9g  %10i\n", idx, segid,
                   parm.Res(resnum).OriginalResNum(), parm.Res(resnum).c_str(),
                   atom->c_str(), typeindex, atom->Charge(),
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
  // Write NTHETA section
  outfile.Printf("%8u !NTHETA: angles\n", parm.Angles().size() + parm.AnglesH().size());
  idx = 1;
  for (AngleArray::const_iterator ang = parm.AnglesH().begin();
                                  ang != parm.AnglesH().end(); ++ang, ++idx)
  {
    outfile.Printf("%8i%8i%8i", ang->A1()+1, ang->A2()+1, ang->A3()+1);
    if ((idx % 3)==0) outfile.Printf("\n");
  }
  for (AngleArray::const_iterator ang = parm.Angles().begin();
                                  ang != parm.Angles().end(); ++ang, ++idx)
  {
    outfile.Printf("%8i%8i%8i", ang->A1()+1, ang->A2()+1, ang->A3()+1);
    if ((idx % 3)==0) outfile.Printf("\n");
  }
  if ((idx % 3)==0) outfile.Printf("\n");
  outfile.Printf("\n");
  // Write out NPHI section
  outfile.Printf("%8u !NPHI: dihedrals\n", parm.Dihedrals().size() + parm.DihedralsH().size());
  idx = 1;
  for (DihedralArray::const_iterator dih = parm.DihedralsH().begin();
                                     dih != parm.DihedralsH().end(); ++dih, ++idx)
  {
    outfile.Printf("%8i%8i%8i%8i", dih->A1()+1, dih->A2()+1, dih->A3()+1, dih->A4()+1);
    if ((idx % 2)==0) outfile.Printf("\n");
  }
  for (DihedralArray::const_iterator dih = parm.Dihedrals().begin();
                                     dih != parm.Dihedrals().end(); ++dih, ++idx)
  {
    outfile.Printf("%8i%8i%8i%8i", dih->A1()+1, dih->A2()+1, dih->A3()+1, dih->A4()+1);
    if ((idx % 2)==0) outfile.Printf("\n");
  }
  if ((idx % 2)==0) outfile.Printf("\n");
  outfile.Printf("\n");

  outfile.CloseFile();
  return 0;
}
