#include <cstdlib> // atoi, atof
#include <cstring> // strncmp
#include "PDBfile.h"

/// PDB record types
// NOTE: Must correspond with PDB_RECTYPE
const char* PDBfile::PDB_RECNAME[] = { "ATOM  ", "HETATM", "TER   " };

// PDBfile::IsPDBkeyword()
bool PDBfile::IsPDBkeyword(std::string const& recname) {
  // Title Section
  if (recname.compare(0,6,"HEADER")==0) return true;
  if (recname.compare(0,6,"SOURCE")==0) return true;
  if (recname.compare(0,6,"AUTHOR")==0) return true;
  if (recname.compare(0,6,"OBSLTE")==0) return true;
  if (recname.compare(0,6,"KEYWDS")==0) return true;
  if (recname.compare(0,6,"REVDAT")==0) return true;
  if (recname.compare(0,6,"TITLE ")==0) return true;
  if (recname.compare(0,6,"EXPDTA")==0) return true;
  if (recname.compare(0,6,"SPRSDE")==0) return true;
  if (recname.compare(0,6,"SPLT  ")==0) return true;
  if (recname.compare(0,6,"NUMMDL")==0) return true;
  if (recname.compare(0,6,"JRNL  ")==0) return true;
  if (recname.compare(0,6,"CAVEAT")==0) return true;
  if (recname.compare(0,6,"MDLTYP")==0) return true;
  if (recname.compare(0,6,"REMARK")==0) return true;
  if (recname.compare(0,6,"COMPND")==0) return true;
  // Primary Structure Section
  if (recname.compare(0,5,"DBREF" )==0) return true; // DBREF[|1|2]
  if (recname.compare(0,6,"SEQADV")==0) return true;
  if (recname.compare(0,6,"MODRES")==0) return true;
  if (recname.compare(0,6,"SEQRES")==0) return true;
  // Coordinate Section
  if (recname.compare(0,6,"MODEL ")==0) return true;
  if (recname.compare(0,6,"ATOM  ")==0) return true;
  if (recname.compare(0,6,"HETATM")==0) return true;
  if (recname.compare(0,6,"TER   ")==0) return true;
  // Crystallographic and Coordinate Transformation Section 
  if (recname.compare(0,6,"CRYST1")==0) return true;
  if (recname.compare(0,5,"SCALE" )==0) return true; // SCALEn
  if (recname.compare(0,5,"ORIGX" )==0) return true; // ORIGXn
  if (recname.compare(0,5,"MTRIX" )==0) return true; // MTRIXn
  return false;
}

// PDBfile::ID_PDB()
bool PDBfile::ID_PDB(CpptrajFile& fileIn) {
  // NOTE: ASSUME FILE SET UP FOR READ
  if (fileIn.OpenFile()) return false;
  std::string line1 = fileIn.GetLine();
  std::string line2 = fileIn.GetLine();
  fileIn.CloseFile();
  if (!IsPDBkeyword( line1 )) return false;
  if (!IsPDBkeyword( line2 )) return false;
  return true;
}

// PDBfile::IsPDBatomKeyword()
bool PDBfile::IsPDBatomKeyword() {
  if (strncmp(linebuffer_,"ATOM  ",6)==0) return true;
  if (strncmp(linebuffer_,"HETATM",6)==0) return true;
  return false;
}

bool PDBfile::IsPDB_TER() {
  if (linebuffer_[0]=='T' && linebuffer_[1]=='E' && linebuffer_[2]=='R')
    return true;
  return false;
}

bool PDBfile::IsPDB_END() {
  if (linebuffer_[0]=='E' && linebuffer_[1]=='N' && linebuffer_[2]=='D')
    return true;
  return false;
}

// PDBfile::pdb_Atom()
Atom PDBfile::pdb_Atom() {
  // Atom number (6-11)
  // Atom name (12-16)
  char savechar = linebuffer_[16];
  linebuffer_[16] = '\0';
  NameType aname(linebuffer_+12);
  // Replace asterisks with single quotes
  aname.ReplaceAsterisk();
  linebuffer_[16] = savechar;
  return Atom(aname);
}

// PDBfile::pdb_Residue()
NameType PDBfile::pdb_Residue(int& current_res) {
  // Res name (16-20)
  char savechar = linebuffer_[20];
  linebuffer_[20] = '\0';
  NameType resname(linebuffer_+16);
  // Replace asterisks with single quotes
  resname.ReplaceAsterisk();
  linebuffer_[20] = savechar;
  // Res num (22-27)
  savechar = linebuffer_[27];
  linebuffer_[27] = '\0';
  current_res = atoi( linebuffer_+22 );
  linebuffer_[27] = savechar;
  return resname;
}

// PDBfile::pdb_XYZ()
// NOTE: Should check for null Xout
void PDBfile::pdb_XYZ(double *Xout) {
  // X coord (30-38)
  char savechar = linebuffer_[38];
  linebuffer_[38] = '\0';
  Xout[0] = atof( linebuffer_+30 );
  linebuffer_[38] = savechar;
  // Y coord (38-46)
  savechar = linebuffer_[46];
  linebuffer_[46] = '\0';
  Xout[1] = atof( linebuffer_+38 );
  linebuffer_[46] = savechar;
  // Z coord (46-54)
  savechar = linebuffer_[54];
  linebuffer_[54] = '\0';
  Xout[2] = atof( linebuffer_+46 );
  linebuffer_[54] = savechar;
}

void PDBfile::WriteHET(int res, double x, double y, double z) {
  WriteRec(HETATM, anum_++, "XX", "XXX", ' ', res, x, y, z);
}

void PDBfile::WriteATOM(int res, double x, double y, double z, 
                        const char* resnameIn, double Occ)
{
  WriteRec(ATOM, anum_++, "XX", resnameIn, ' ',
           res, x, y, z, (float)Occ, 0, "", 0, false);
}

void PDBfile::WriteATOM(const char* anameIn, int res, double x, double y, double z, 
                        const char* resnameIn, double Occ)
{
  WriteRec(ATOM, anum_++, anameIn, resnameIn, ' ',
           res, x, y, z, (float)Occ, 0, "", 0, false);
}

void PDBfile::WriteTER(int anum, NameType const& resnameIn, char chain, int resnum)
{
  WriteRec(TER, anum, "", resnameIn, chain, resnum, 0, 0, 0, 0, 0, "", 0, false);
}

void PDBfile::WriteRec(PDB_RECTYPE Record, int anum, NameType const& name,
                       NameType const& resnameIn, char chain, int resnum,
                       double X, double Y, double Z)
{
  WriteRec(Record, anum, name, resnameIn, chain, resnum, X, Y, Z, 0, 0, "", 0, false);
}

/// Write out an ATOM or HETATM record
/** \return the number of characters written */
void PDBfile::WriteRec(PDB_RECTYPE Record, int anum, NameType const& name,
                       NameType const& resnameIn, char chain, int resnum,
                       double X, double Y, double Z, float Occ, float B, 
                       const char* Elt, int charge, bool highPrecision) 
{
  char resName[5], atomName[5];

  resName[4]='\0';
  atomName[4]='\0';
  // Residue number in PDB format can only be 4 digits wide
  while (resnum>9999) resnum-=9999;
  // Atom number in PDB format can only be 5 digits wide
  while (anum>99999) anum-=99999;
  // Residue names in PDB format are 3 chars long starting at column 18. 
  // However in Amber residues are 4 characters long, usually with a space
  // at the end. If this is the case remove the space so that the residue name
  // lines up properly.
  resName[0] = resnameIn[0];
  resName[1] = resnameIn[1];
  resName[2] = resnameIn[2];
  if (resnameIn[3]!=' ')
    resName[3] = resnameIn[3];
  else
    resName[3] = '\0';
  // Atom names in PDB format start from col 14 when <= 3 chars, 13 when 4 chars.
  if (name[3]!=' ') { // 4 chars
    atomName[0] = name[0];
    atomName[1] = name[1];
    atomName[2] = name[2];
    atomName[3] = name[3];
  } else {            // <= 3 chars
    atomName[0] = ' ';
    atomName[1] = name[0];
    atomName[2] = name[1];
    atomName[3] = name[2];
  }

  Printf("%-6s%5i %-4s%4s %c%4i",PDB_RECNAME[Record], anum, atomName, 
               resName, chain, resnum);
  if (Record == TER) 
    Printf("\n");
  else if (highPrecision)
    Printf("    %8.3f%8.3f%8.3f%8.4f%8.4f      %2s%2s\n", X, Y, Z, Occ, B, Elt, "");
  else
    Printf("    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n", X, Y, Z, Occ, B, Elt, "");
}

/* Additional Values:
 *  Chain ID : buffer[21]
 *  Occupancy: buffer[54] to buffer[59]
 *  B-factor : buffer[60] to buffer[65]
 *  10 chars between Bfactor and element: buffer[66] to buffer[75]
 *  Element  : buffer[76] and buffer[77]
 *  Charge   : buffer[78] and buffer[79]
 */
