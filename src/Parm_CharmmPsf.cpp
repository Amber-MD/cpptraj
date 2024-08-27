// Parm_CharmmPsf.cpp
#include <cstdio> // sscanf
#include <cstdlib> // atoi
#include <cstring> // strncmp
#include <cctype> // isdigit
#include "Parm_CharmmPsf.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include "Mol.h" // UniqueCount()
#include "CharmmParamFile.h"
#include "BufferedLine.h"

/// CONSTRUCTOR
Parm_CharmmPsf::Parm_CharmmPsf() :
  extfmt_(false),
  cheq_(false),
  xplor_(true)
{}

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
int Parm_CharmmPsf::FindTag(char* tag, const char* target, BufferedLine& infile) {
  int nval = 0;
  int tgtsize = strlen( target );
  while (strncmp(tag,target,tgtsize)!=0) {
    const char* buffer = infile.Line();
    if ( buffer == 0 ) return 0;
    sscanf(buffer,"%i %10s",&nval,tag);
  }
  return nval;
}

//  Parm_CharmmPsf::ReadDihedrals()
int Parm_CharmmPsf::ReadDihedrals(BufferedLine& infile, int ndihedral, const char* typestr, Topology& parmOut) const
{
    bool found;
    int bondatoms[8];
    const char* buffer = 0; 
    int nlines = ndihedral / 2;
    if ( (ndihedral % 2) != 0) nlines++;
    for (int dihline = 0; dihline < nlines; dihline++) {
      if ( (buffer=infile.Line()) == 0) {
        mprinterr("Error: Reading %s line %i\n", typestr, dihline+1);
        return 1;
      }
      // Each line has 2 groups of 4 atom numbers
      int ndihread = sscanf(buffer,"%i %i %i %i %i %i %i %i",bondatoms,bondatoms+1,
                              bondatoms+2,bondatoms+3, bondatoms+4,bondatoms+5,
                              bondatoms+6,bondatoms+7);
      if (params_.DP().empty()) {
        for (int dihidx=0; dihidx < ndihread; dihidx += 4) {
          // TODO: Determine end dihedrals
          if (typestr[0] == 'd')
            parmOut.AddDihedral( DihedralType(bondatoms[dihidx  ]-1,
                                              bondatoms[dihidx+1]-1,
                                              bondatoms[dihidx+2]-1,
                                              bondatoms[dihidx+3]-1,
                                              DihedralType::NORMAL) );
          else
            parmOut.AddCharmmImproper( DihedralType(bondatoms[dihidx  ]-1,
                                              bondatoms[dihidx+1]-1,
                                              bondatoms[dihidx+2]-1,
                                              bondatoms[dihidx+3]-1,
                                              DihedralType::IMPROPER), -1 );
        }
      } else {
        for (int dihidx=0; dihidx < ndihread; dihidx += 4) {
          int a1 = bondatoms[dihidx]-1;
          int a2 = bondatoms[dihidx+1]-1;
          int a3 = bondatoms[dihidx+2]-1;
          int a4 = bondatoms[dihidx+3]-1;
          DihedralType dih(a1, a2, a3, a4, DihedralType::NORMAL);
          TypeNameHolder types(4);
          types.AddName( parmOut[a1].Type() );
          types.AddName( parmOut[a2].Type() );
          types.AddName( parmOut[a3].Type() );
          types.AddName( parmOut[a4].Type() );
          if (typestr[0] == 'd') {
            // Dihedral. Can have multiple multiplicities.
            DihedralParmArray dpa = params_.DP().FindParam( types, found );
            if (found) {
              for (DihedralParmArray::const_iterator dp = dpa.begin(); dp != dpa.end(); ++dp)
                parmOut.AddDihedral( dih, *dp );
            }
          } else {
            // Charmm Improper. Expect only one paramter per type.
            DihedralParmType ipt = params_.IP().FindParam( types, found );
            parmOut.AddCharmmImproper( dih, ipt );
          } 
          if (!found) {
            mprintf("Warning: Parameters not found for %s %s - %s - %s - %s\n", typestr, parmOut.AtomMaskName(a1).c_str(), parmOut.AtomMaskName(a2).c_str(), parmOut.AtomMaskName(a3).c_str(), parmOut.AtomMaskName(a4).c_str());
            if (typestr[0] == 'd')
              parmOut.AddDihedral( dih );
            else
              parmOut.AddCharmmImproper( dih );
          }
        } // END loop over number of dihedrals read
      } // END if dihedral params present
    } // END loop over lines
  return 0;
}

// ----- Charmm Lone Pair ------------------------
class Parm_CharmmPsf::LonePair {
  public:
    //LonePair() : numSupport_(0), idx_(-1), dist_(0), ang_(0), dih_(0) {}
    LonePair(int n, int i, const char* t, double dist, double ang, double dih) :
      nat_(n), idx_(i), type_(t), dist_(dist), ang_(ang), dih_(dih) {}

    int Nat() const { return nat_; }
    int Idx() const { return idx_; }
    const char* lptype() const { return type_.c_str(); }
    double Dist() const { return dist_; }
    double Ang() const { return ang_; }
    double Dih() const { return dih_; }

    void Print() const { mprintf(" %i %i %s %f %f %f", Nat(), Idx()+1, lptype(), Dist(), Ang(), Dih()); }
  private:
    int nat_;          ///< Number of support atoms
    int idx_;          ///< Index of support atoms in support (host) array
    std::string type_; ///< Lone pair type TODO make char?
    double dist_;      ///< Distance in angstroms
    double ang_;       ///< Angle in degrees
    double dih_;       ///< Dihedral in degrees
};
// -----------------------------------------------

/** Read lone pairs. */
int Parm_CharmmPsf::ReadLonePairs(BufferedLine& infile, int numlp, int numlph, Topology& parmOut)
const
{
  std::vector<LonePair> LPs;
  LPs.reserve( numlp );
  // Read the lone pairs
  const char* buffer = 0;
  char type[9];
  type[8] = '\0';
  int num = 0;
  int idx = -1;
  double dist = 0;
  double ang = 0;
  double dih = 0;
  for (int ilp = 0; ilp != numlp; ilp++) {
    if ( (buffer = infile.Line()) == 0) {
      mprinterr("Error: Reading lone pair %i on line %i\n", ilp+1, infile.LineNumber());
      return 1;
    }
    //mprintf("DEBUG: LP: %s\n", buffer);
    num = 0;
    idx = -1;
    dist = 0;
    ang = 0;
    dih = 0;
    type[0] = '\0';
    int nscan = sscanf(buffer, "%i %i %s %lf %lf %lf", &num, &idx, type, &dist, &ang, &dih);
    if (nscan < 4) {
      mprinterr("Error: Only read %i values for lone pair line, expected at least 4.\n", nscan);
      mprinterr("Error: Line %i: %s\n", infile.LineNumber(), buffer);
      return 1;
    }
    LPs.push_back( LonePair(num, idx-1, type, dist, ang, dih) );
  }
  // Read the lone pair atoms
  int lpat[8];
  int nlines = numlph / 8;
  if ( (numlph%8) > 0 )
    nlines++;
  std::vector<int> LPatoms;
  LPatoms.reserve( numlph );
  for (int ilp = 0; ilp != nlines; ilp++) {
    if ( (buffer = infile.Line()) == 0) {
      mprinterr("Error: Reading lone pair atoms on line %i\n", infile.LineNumber());
      return 1;
    }
    int nscan = sscanf(buffer, "%i %i %i %i %i %i %i %i",
                       lpat, lpat+1, lpat+2, lpat+3, lpat+4, lpat+5, lpat+6, lpat+7);
    for (int ii = 0; ii != nscan; ii++)
      LPatoms.push_back( lpat[ii]-1 );
  }
  if ( (unsigned int)numlph != LPatoms.size() ) {
    mprinterr("Error: Number of lone pair atoms read (%zu) != expected (%i)\n", LPatoms.size(), numlph);
    return 1;
  }
  // Add bonds for lone pairs TODO angles and dihedrals?
  for (std::vector<LonePair>::const_iterator it = LPs.begin(); it != LPs.end(); ++it) {
    // Sanity check
    int lpAtomIdx = LPatoms[it->Idx()];
    int bondedAtomIdx = LPatoms[it->Idx()+1];
    if ( parmOut[lpAtomIdx].Element() != Atom::EXTRAPT ) {
      mprintf("Warning: PSF defines atom %s as lone pair, but it does not appear to be one.\n",
               parmOut.AtomMaskName( lpAtomIdx ).c_str());
    }
    if (params_.BP().empty()) {
      // Just add the bond to lone pair
      parmOut.AddBond( lpAtomIdx, bondedAtomIdx );
    } else {
      // Add the bond and distance as a parameter FIXME is Rk=0 ok?
      // FIXME what does a negative distance mean?
      double req = it->Dist();
      if (req < 0) req = -req;
      parmOut.AddBond( lpAtomIdx, bondedAtomIdx, BondParmType(0, req) );
    }
  }
  // DEBUG - Print out lone pair information
  if (debug_ > 0) {
    for (std::vector<LonePair>::const_iterator it = LPs.begin(); it != LPs.end(); ++it) {
      mprintf("DEBUG: LP:");
      it->Print();
      mprintf(" Atoms:");
      int jj = it->Idx();
      // Nat() is the number of support atoms; +1 for the lone pair atom itself
      for (int ii = 0; ii <= it->Nat(); ii++, jj++)
        mprintf(" %i", LPatoms[jj]+1);
      mprintf("\n");
    }
  }

  return 0;
}

const unsigned int Parm_CharmmPsf::ChmStrMax_ = 9;

/** Extract residue number and alternatively the insertion code. */
int Parm_CharmmPsf::ParseResID(char& psficode, const char* psfresid)
{
  char buf[ChmStrMax_];
  int bidx = -1;
  const char* ptr = psfresid;
  // Parse out residue number
  int resnum = 0;
  while ( isdigit( *ptr ) && bidx < (int)ChmStrMax_ )
  {
    buf[++bidx] = *ptr;
    ++ptr;
  }
  bidx++;
  buf[bidx] = '\0';
  // Sanity check
  if (bidx < 1)
    mprintf("Warning: PSF residue ID does not begin with a digit: '%s'\n", psfresid);
  else
    resnum = atoi( buf );
  // Check for icode
  if (*ptr != '\0')
    psficode = *ptr;
  else
    psficode = ' ';
  return resnum;
}

// Parm_CharmmPsf::ReadParm()
/** Open the Charmm PSF file specified by filename and set up topology data.
  * Mask selection requires natom, nres, names, resnames, resnums.
  */
int Parm_CharmmPsf::ReadParm(FileName const& fname, Topology &parmOut) {
  const size_t TAGSIZE = 10; 
  char tag[TAGSIZE];
  tag[0]='\0';

  BufferedLine infile;
  if (infile.OpenFileRead(fname)) return 1;
  mprintf("    Reading Charmm PSF file %s as topology file.\n",infile.Filename().base());
  // Read the first line, should contain PSF...
  const char* buffer = 0;
  if ( (buffer=infile.Line()) == 0 ) return 1;
  // Advance to <ntitle> !NTITLE
  int ntitle = FindTag(tag, "!NTITLE", infile); 
  // Only read in 1st title. Skip any asterisks.
  std::string psftitle;
  if (ntitle > 0) {
    buffer = infile.Line();
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
  if (debug_ > 0) params_.Debug();
  // Read the next natom lines
  char psfresid[ChmStrMax_];
  char psfresname[ChmStrMax_];
  char psfname[ChmStrMax_];
  char psftype[ChmStrMax_];
  char segmentID[ChmStrMax_];
  char psficode;
  double psfcharge;
  double psfmass;
//  typedef std::vector<std::string> Sarray;
  ParmHolder<AtomType>& atomTypes = params_.AT();
//  Sarray SegIDs;
  bool firstLine = true;
  enum PsfFormatType { T_CHARMM = 0, T_VMD };
  PsfFormatType psfFormatType = T_CHARMM;
  for (int atom=0; atom < natom; atom++) {
    if ( (buffer=infile.Line()) == 0 ) {
      mprinterr("Error: ReadParmPSF(): Reading atom %i\n",atom+1);
      return 1;
    }
    // Check first line
    if (firstLine) {
      int ntokens = infile.TokenizeLine(" \t");
      //mprintf("DEBUG: first line: ntokens= %i\n", ntokens);
      // A CHARMM psf will have 9-11 columns:
      // AtomID, SegID, ResID, ResName, AtomName, AtomType, Charge, Mass, Constrained[, Polarizability, TholeScaleFactor]
      if (ntokens < 9) {
        mprintf("Warning: PSF has non-standard format; atoms line contains less than 11 columns (%i).\n", ntokens);
        if (ntokens == 8) {
          mprintf("Warning: 8 columns in atoms line; assuming VMD or similarly formatted PSF.\n");
          psfFormatType = T_VMD;
        } else {
          mprinterr("Error: Unrecognized number of columns in atoms line (%i); PSF has bad format.\n", ntokens);
          return 1;
        }
      }
      firstLine = false;
    }
    // Read line
    if (psfFormatType == T_CHARMM) {
      // ATOM# SEGID RESID RES ATNAME ATTYPE CHRG MASS CONST POL THOLE 
      int nread = sscanf(buffer,"%*i %s %s %s %s %s %lf %lf", segmentID, psfresid, psfresname, 
                         psfname, psftype, &psfcharge, &psfmass);
      //mprintf("DEBUG: Read %i columns.\n", nread);
      if (nread != 7) {
        mprintf("Warning: During read of PSF atoms, expected to read 7 columns, got %i (line %i)\n", nread, infile.LineNumber());
      }
    } else if (psfFormatType == T_VMD) {
      // ATOM# RESID RES ATNAME ATTYPE CHRG MASS CONST
      int nread = sscanf(buffer,"%*i %s %s %s %s %lf %lf", psfresid, psfresname,
                         psfname, psftype, &psfcharge, &psfmass);
      segmentID[0] = '\0';
      if (nread != 6) {
        mprintf("Warning: During read of PSF atoms, expected to read 6 columns, got %i (line %i)\n", nread, infile.LineNumber());
      }
    }
    // Determine if this is a Drude particle
    Atom::AtomicElementType psfElt = Atom::UNKNOWN_ELEMENT;
    if (psftype[0] == 'D') {
      //mprintf("DEBUG: Potential Drude particle: %s %s\n", psfname, psftype);
      psfElt = Atom::DRUDE;
    }
    // Extract residue number and alternatively insertion code.
    int psfresnum = ParseResID(psficode, psfresid);
    //mprintf("DEBUG: resnum %10i  icode %c\n", psfresnum, psficode);
    // Search for segment ID
    /*int idx = -1;
    if (segmentID[0] != '\0') {
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
    }*/
    atomTypes.AddParm( TypeNameHolder(NameType(psftype)), AtomType(psfmass), false );
    Atom chmAtom( psfname, psfcharge, psfmass, psftype, psfElt );
    parmOut.AddTopAtom( chmAtom, Residue(psfresname, psfresnum, ' ', std::string(segmentID)) );
  } // END loop over atoms 
  // Advance to <nbond> !NBOND
  int bondatoms[9];
  int nbond = FindTag(tag, "!NBOND", infile);
  if (nbond > 0) {
    if (debug_>0) mprintf("\tPSF: !NBOND tag found, nbond=%i\n", nbond);
    int nlines = nbond / 4;
    if ( (nbond % 4) != 0) nlines++;
    for (int bondline=0; bondline < nlines; bondline++) {
      if ( (buffer=infile.Line()) == 0 ) {
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
          TypeNameHolder types(2);
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
      if ( (buffer=infile.Line()) == 0) {
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
          TypeNameHolder types(3);
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
  // SKIPPING NDON, NACC, NNB, NGRP
  // Advance to <# lone pairs> <# lone pair hosts> NUMLP NUMLPH
  int numlp = -1;
  int numlph = -1;
  while (strncmp(tag, "!NUMLP", 6) !=0) {
    const char* buffer = infile.Line();
    if ( buffer == 0 ) break;
    sscanf(buffer,"%i %i %10s", &numlp, &numlph, tag);
  }
  if (numlp > -1) {
    mprintf("DEBUG: PSF contains %i lone pairs, %i lone pair hosts.\n", numlp, numlph);
    if (ReadLonePairs(infile, numlp, numlph, parmOut)) return 1;
  }

  mprintf("\tPSF contains %i atoms, %i residues.\n", parmOut.Natom(), parmOut.Nres());

  infile.CloseFile();

  // Add nonbonded parameters
  if (params_.HasLJparams()) {
    parmOut.AssignNonbondParams( atomTypes, params_.NB() );
  }

  return 0;
}

// =============================================================================
// Parm_CharmmPsf::WriteHelp()
void Parm_CharmmPsf::WriteHelp() {
  mprintf("\toldpsf : Write atom type indices instead of type names (not recommended).\n"
          "\text    : Use extended format.\n");
}

// Parm_CharmmPsf::processWriteArgs()
int Parm_CharmmPsf::processWriteArgs(ArgList& argIn) {
  if (argIn.hasKey("oldpsf")) xplor_ = false;
  if (argIn.hasKey("ext")) extfmt_ = true;
  return 0;
}

//  FindMolType()
static int FindMolType(int molNum, Mol::Marray const& mols) {
  for (Mol::Marray::const_iterator mol = mols.begin(); mol != mols.end(); ++mol)
    for (Mol::Iarray::const_iterator it = mol->idxs_.begin();
                                     it != mol->idxs_.end(); ++it)
      if (*it == molNum) return (mol - mols.begin());
  return -1;
}

//  Parm_CharmmPsf::WriteSectionHeader()
void Parm_CharmmPsf::WriteSectionHeader(CpptrajFile& outfile, const char* title, int ival)
const
{
  if (extfmt_)
    outfile.Printf("%10i %s\n", ival, title);
  else
    outfile.Printf("%8i %s\n", ival, title);
}

// WriteDihedrals() // FIXME sort out impropers
static inline void WriteDihedrals(
  DihedralArray const& dihedrals,
  unsigned int& idx,
  CpptrajFile& outfile,
  const char* dihfmt)
{
  for (DihedralArray::const_iterator dih = dihedrals.begin(); dih != dihedrals.end(); ++dih, ++idx)
  {
    outfile.Printf(dihfmt, dih->A1()+1, dih->A2()+1, dih->A3()+1, dih->A4()+1);
    if ((idx % 2)==0) outfile.Printf("\n");
  }
}

// Parm_CharmmPsf::WriteParm()
int Parm_CharmmPsf::WriteParm(FileName const& fname, Topology const& parm) {
  // TODO: CMAP etc info
  CpptrajFile outfile;
  if (outfile.OpenWrite(fname)) return 1;
  // Write PSF header
  std::string header("PSF");
  if (extfmt_) header.append(" EXT");
  if (cheq_)   header.append(" CHEQ");
  if (xplor_)  header.append(" XPLOR");
  outfile.Printf("%s\n\n", header.c_str());
  // Write title FIXME spread out long title over multiple lines
  WriteSectionHeader(outfile, "!NTITLE", 1);
  std::string titleOut = parm.ParmName();
  if (titleOut.size() > 78)
    mprintf("Warning: PSF title will be truncated.\n");
  titleOut.resize(78);
  outfile.Printf("* %-78s\n\n", titleOut.c_str());
  // Write NATOM section
  WriteSectionHeader(outfile, "!NATOM", parm.Natom());
  unsigned int idx = 1;
  // Make segment ids based on molecule type for now.
  Mol::Marray mols = Mol::UniqueCount(parm);
  mprintf("Warning: Unknown segment IDs will be assigned based on molecule type.\n");
  int currentMol = 0;
  int currentMtype = FindMolType(currentMol, mols);
  const char* segid = mols[currentMtype].name_.c_str();
  Mol::Iarray::const_iterator mit = mols[currentMtype].idxs_.begin();
//  bool inSolvent = false;
  // Output format.
  // ATOM# SEGID RES# RES ATNAME ATTYPE CHRG MASS IMOVE
  // Remaining columns (CHEQ)
  //   XPLOR  & DRUDE  : ALPHADP THOLEI
  //   XPLOR  & !DRUDE : ECH     EHA
  //   !XPLOR & DRUDE  : ALPHADP THOLEI ISDRUDE
  //   !XPLOR & !DRUDE : ECH     EHA
  // Where ECH is electronegativity for atoms and EHA is hardness for atoms.
  const char* atmfmt = 0;
  // segIdMax is used to ensure ChainID will not overflow
  unsigned int segIdMax;
  if (extfmt_) {
    segIdMax = 9;
    if (xplor_)
      atmfmt = "%10i %-8s %-8i %-8s %-8s %-6s %14.6G%14.6G%8i\n";
    else
      atmfmt = "%10i %-8s %-8i %-8s %-8s %4s %14.6G%14.6G%8i\n";
  } else {
    segIdMax = 5;
    if (xplor_)
      atmfmt = "%8i %-4s %-4i %-4s %-4s %-4s %14.6G%14.6G%8i\n";
    else
      atmfmt = "%8i %-4s %-4i %-4s %-4s %4s %14.6G%14.6G%8i\n";
  }
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
    // If there is a chain ID, use that instead.
    if (parm.Res(resnum).HasChainID() && parm.Res(resnum).ChainID().size() < segIdMax)
      segid = parm.Res(resnum).chainID();
    // Figure out how atom type is being printed.
    std::string psftype;
    if (xplor_)
      psftype.assign( *(atom->Type()) );
    else
      psftype = integerToString( atom->TypeIndex() + 1 );
    // Figure out what IMOVE should be.
    // IMOVE : 1 = fixed atom, 0 = mobile, -1 = a lonepair (may move, no DoF).
    int imove;
    if (atom->Element() == Atom::EXTRAPT) // TODO is this OK?
      imove = -1;
    else
      imove = 0;
    // Write atom line
    outfile.Printf(atmfmt, idx, segid,
                   parm.Res(resnum).OriginalResNum(), parm.Res(resnum).c_str(),
                   atom->c_str(), psftype.c_str(), atom->Charge(),
                   atom->Mass(), imove);
  }
  outfile.Printf("\n");
  // Write NBOND section
  WriteSectionHeader(outfile, "!NBOND: bonds", parm.Bonds().size() + parm.BondsH().size());
  const char* bndfmt = 0;
  if (extfmt_)
    bndfmt = "%10i%10i";
  else
    bndfmt = "%8i%8i";
  idx = 1;
  for (BondArray::const_iterator bond = parm.BondsH().begin();
                                 bond != parm.BondsH().end(); ++bond, ++idx)
  {
    outfile.Printf(bndfmt, bond->A1()+1, bond->A2()+1);
    if ((idx % 4)==0) outfile.Printf("\n"); 
  }
  for (BondArray::const_iterator bond = parm.Bonds().begin();
                                 bond != parm.Bonds().end(); ++bond, ++idx)
  {
    outfile.Printf(bndfmt, bond->A1()+1, bond->A2()+1);
    if ((idx % 4)==0) outfile.Printf("\n"); 
  }
  if ((idx % 4)!=0) outfile.Printf("\n");
  outfile.Printf("\n");
  // Write NTHETA section
  WriteSectionHeader(outfile, "!NTHETA: angles", parm.Angles().size() + parm.AnglesH().size());
  const char* angfmt = 0;
  if (extfmt_)
    angfmt = "%10i%10i%10i";
  else
    angfmt = "%8i%8i%8i";
  idx = 1;
  for (AngleArray::const_iterator ang = parm.AnglesH().begin();
                                  ang != parm.AnglesH().end(); ++ang, ++idx)
  {
    outfile.Printf(angfmt, ang->A1()+1, ang->A2()+1, ang->A3()+1);
    if ((idx % 3)==0) outfile.Printf("\n");
  }
  for (AngleArray::const_iterator ang = parm.Angles().begin();
                                  ang != parm.Angles().end(); ++ang, ++idx)
  {
    outfile.Printf(angfmt, ang->A1()+1, ang->A2()+1, ang->A3()+1);
    if ((idx % 3)==0) outfile.Printf("\n");
  }
  if ((idx % 3)==0) outfile.Printf("\n");
  outfile.Printf("\n");
  // Write out NPHI section
  WriteSectionHeader(outfile, "!NPHI: dihedrals", parm.Dihedrals().size() + parm.DihedralsH().size());
  const char* dihfmt = 0;
  if (extfmt_)
    dihfmt = "%10i%10i%10i%10i";
  else
    dihfmt = "%8i%8i%8i%8i";
  idx = 1;
  WriteDihedrals(parm.DihedralsH(), idx, outfile, dihfmt);
  WriteDihedrals(parm.Dihedrals(), idx, outfile, dihfmt);
  if ((idx % 2)==0) outfile.Printf("\n");
  outfile.Printf("\n");
  // Write out NIMPHI section
  WriteSectionHeader(outfile, "!NIMPHI: impropers", parm.Chamber().Impropers().size());
  idx = 1;
  WriteDihedrals(parm.Chamber().Impropers(), idx, outfile, dihfmt);
  if ((idx % 2)==0) outfile.Printf("\n");
  outfile.Printf("\n");

  outfile.CloseFile();
  return 0;
}
