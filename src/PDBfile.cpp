#include <cstdio>  // sscanf
#include <cstdlib> // atoi, atof
#include <cstring> // strncmp
#include "PDBfile.h"
#include "CpptrajStdio.h"

/// PDB record types
// NOTE: Must correspond with PDB_RECTYPE
const char* PDBfile::PDB_RECNAME[] = { "ATOM  ", "HETATM", "CRYST1", "TER   ",
  "END   ", "ANISOU", "EndRec", "CONECT", 0 };

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
  if (recname.compare(0,6,"ANISOU")==0) return true;
  if (recname.compare(0,3,"TER"   )==0) return true; // To recognize blank TER cards.
  if (recname.compare(0,6,"HETATM")==0) return true;
  // Crystallographic and Coordinate Transformation Section 
  if (recname.compare(0,6,"CRYST1")==0) return true;
  if (recname.compare(0,5,"SCALE" )==0) return true; // SCALEn
  if (recname.compare(0,5,"ORIGX" )==0) return true; // ORIGXn
  if (recname.compare(0,5,"MTRIX" )==0) return true; // MTRIXn
  if (recname.compare(0,9,"USER  MOD")==0) return true; // reduce
  // Bookkeeping Section
  if (recname.compare(0,6,"MASTER")==0) return true;
  if (recname.compare(0,3,"END"   )==0) return true;
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
  if (!line2.empty() && !IsPDBkeyword( line2 )) return false;
  return true;
}

// PDBfile::NextRecord()
PDBfile::PDB_RECTYPE PDBfile::NextRecord() {
  if (NextLine() == 0) { // No more records to read
    recType_ = END_OF_FILE;
    return END_OF_FILE;
  }
  recType_ = UNKNOWN;
  // Try to ID the current record.
  if (strncmp(linebuffer_,"ATOM  ",6)==0 ||
      strncmp(linebuffer_,"HETATM",6)==0   )
    // Just return ATOM for ATOM/HETATM when reading.
    recType_ = ATOM;
  else if (strncmp(linebuffer_,"CONECT",6)==0)
    recType_ = CONECT;
  else if (strncmp(linebuffer_,"CRYST1",6)==0)
    recType_ = CRYST1;
  else if (linebuffer_[0]=='T' && linebuffer_[1]=='E' && linebuffer_[2]=='R')
    recType_ = TER;
  else if (linebuffer_[0]=='E' && linebuffer_[1]=='N' && linebuffer_[2]=='D')
    recType_ = END;
  return recType_;
}

Atom PDBfile::pdb_Atom(char& altLoc, int& atnum) {
  // ATOM or HETATM keyword.
  // Check line length before any modification.
  size_t lineLength = strlen( linebuffer_ );
  // Atom number (6-11)
  altLoc = linebuffer_[11];
  linebuffer_[11] = '\0';
  atnum = atoi(linebuffer_+6);
  linebuffer_[11] = altLoc;
  // Atom name (12-16), alt location indicator (16)
  // Replace asterisks in name with single quotes.
  altLoc = linebuffer_[16];
  linebuffer_[16] = '\0';
  NameType aname(linebuffer_+12);
  aname.ReplaceAsterisk();
  linebuffer_[16] = altLoc;
  // Element (76-77), Protect against broken PDB files (lines too short).
  char eltString[2]; eltString[0] = ' '; eltString[1] = ' ';
  if (lineLength > 77) {
    eltString[0] = linebuffer_[76];
    eltString[1] = linebuffer_[77];
  } else if (!lineLengthWarning_) {
    mprintf("Warning: PDB line length is short (%zu chars, expected 80).\n", lineLength);
    lineLengthWarning_ = true;
  }
  // NOTE: Additional values:
  //       10 chars between Bfactor and element: buffer[66] to buffer[75]
  //       Charge: buffer[78] and buffer[79]
  return Atom(aname, eltString);
}

Residue PDBfile::pdb_Residue() {
  // Res name (17-20), Replace asterisks with single quotes.
  char savechar = linebuffer_[20];
  linebuffer_[20] = '\0';
  NameType resName(linebuffer_+17);
  linebuffer_[20] = savechar;
  resName.ReplaceAsterisk();
  // Chain ID (21)
  // Res num (22-26), insertion code (26)
  char icode = linebuffer_[26];
  linebuffer_[26] = '\0';
  int resnum = atoi( linebuffer_+22 );
  linebuffer_[26] = icode;
  return Residue( resName, resnum, icode, linebuffer_[21] );
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

// PDBfile::pdb_OccupancyAndBfactor()
void PDBfile::pdb_OccupancyAndBfactor(float& occ, float& bfac) {
  // Occupancy (54-60)
  char savechar = linebuffer_[60];
  linebuffer_[60] = '\0';
  occ = atof(linebuffer_ + 54);
  linebuffer_[60] = savechar;
  // B-factor (60-66)
  savechar = linebuffer_[66];
  linebuffer_[66] = '\0';
  bfac = atof(linebuffer_ + 60);
  linebuffer_[66] = savechar;
}
  
/** Read charge and radius from PQR file (where occupancy and B-factor would be
  * in a PDB). Use sscanf() since these columns could have different widths.
  * Could fail if reading a PDB with values > 99.99 in B-factor column.
  */
void PDBfile::pdb_ChargeAndRadius(float& charge, float& radius) {
  sscanf(linebuffer_+54, "%f %f", &charge, &radius);
}

void PDBfile::pdb_Box(double* box) {
  // CRYST1 keyword. RECORD A B C ALPHA BETA GAMMA SGROUP Z
  unsigned int lb_size = strlen(linebuffer_);
  if (lb_size < 54) {
    mprintf("Warning: Malformed CRYST1 record. Skipping.\n");
    return;
  }
  // A=6-15 B=15-24 C=24-33
  unsigned int lb = 6;
  for (unsigned int ib = 0; ib != 3; ib++, lb += 9) {
    unsigned int end = lb + 9;
    char savechar = linebuffer_[end];
    linebuffer_[end] = '\0';
    box[ib] = atof( linebuffer_ + lb );
    linebuffer_[end] = savechar;
  }
  // alpha=33-40 beta=40-47 gamma=47-54
  for (unsigned int ib = 3; ib != 6; ib++, lb += 7) {
    unsigned int end = lb + 7;
    char savechar = linebuffer_[end];
    linebuffer_[end] = '\0';
    box[ib] = atof( linebuffer_ + lb );
    linebuffer_[end] = savechar;
  }
  mprintf("\tRead CRYST1 info from PDB: a=%g b=%g c=%g alpha=%g beta=%g gamma=%g\n",
          box[0], box[1], box[2], box[3], box[4], box[5]);
  // Warn if the box looks strange.
  if (box[0] == 1.0 && box[1] == 1.0 && box[2] == 1.0)
    mprintf("Warning: PDB cell lengths are all 1.0 Ang.;"
            " this usually indicates an invalid box.\n");
}

int PDBfile::pdb_Bonds(int* bnd) {
  unsigned int lb_size = strlen(linebuffer_);
  int Nscan = 0;
  for (unsigned int lb = 6; lb < lb_size; lb += 5) {
    // Check if first char is newline or final char is blank.
    if (linebuffer_[lb] == '\n' || linebuffer_[lb + 4] == ' ') break;
    // Safety valve
    if (Nscan == 5) {
      mprintf("Warning: CONECT record has more than 4 bonds. Only using first 4 bonds.\n");
      break;
    }
    unsigned int end = lb + 5;
    char savechar = linebuffer_[end];
    linebuffer_[end] = '\0';
    bnd[Nscan++] = atof( linebuffer_ + lb );
    linebuffer_[end] = savechar;
  }
  if (Nscan < 2)
    mprintf("Warning: Malformed CONECT record: %s", linebuffer_);
  //mprintf("DEBUG: CONECT: Atom record %i to", bnd[0]);
  //for (int i = 1; i < Nscan; i++)
  //  mprintf(" %i", bnd[i]);
  //mprintf("\n");
  return Nscan;
}

// -----------------------------------------------------------------------------
// PDBfile::WriteRecordHeader()
void PDBfile::WriteRecordHeader(PDB_RECTYPE Record, int anum, NameType const& name,
                                char altLoc, NameType const& resnameIn, char chain, 
                                int resnum, char icode, const char* Elt)
{
  char resName[5], atomName[5];

  resName[4]='\0';
  atomName[4]='\0';
  // Residue number in PDB format can only be 4 digits wide
  if (resnum > 9999) resnum = resnum % 10000;
  // Atom number in PDB format can only be 5 digits wide
  if (anum > 99999) anum = anum % 100000;
  // Residue names in PDB format are 3 chars long, right-justified, starting
  // at column 18, while the alternate location indicator is column 17. 
  // However in Amber residues can be 4 characters long; in this case overwrite
  // the alternate location indicator.
  resName[0] = altLoc;
  resName[1] = ' ';
  resName[2] = ' '; // TODO set location 3 as well?
  const char* ptr = *resnameIn;
  while (*ptr != ' ' && *ptr != '\0') ++ptr;
  int rn_size = (int)(ptr - *resnameIn);
  // Protect against residue names larger than 4 chars.
  if (rn_size > 4) rn_size = 4;
  int rn_idx = 3;
  for (int i = rn_size - 1; i > -1; i--, rn_idx--)
    resName[rn_idx] = resnameIn[i];
  // Determine size in characters of element name if given.
  int eNameChars = 0;
  if (Elt != 0) eNameChars = strlen( Elt );
  // For atoms with element names of 1 character, names in PDB format start
  // from col 14 when <= 3 chars, 13 when 4 chars. Atoms with element names of
  // 2 characters start from col 13.
  if (eNameChars == 2 || name[3] != ' ') { // 4 chars or 2 char elt name
    atomName[0] = name[0];
    atomName[1] = name[1];
    atomName[2] = name[2];
    atomName[3] = name[3];
  } else {            // <= 3 chars or 1 char elt name
    atomName[0] = ' ';
    atomName[1] = name[0];
    atomName[2] = name[1];
    atomName[3] = name[2];
  }

  Printf("%-6s%5i %-4s%4s %c%4i%c",PDB_RECNAME[Record], anum, atomName,
               resName, chain, resnum, icode);
  if (Record == TER) Printf("\n");
}

// PDBfile::WriteHET()
void PDBfile::WriteHET(int res, double x, double y, double z) {
  WriteCoord(HETATM, anum_++, "XX", ' ', "XXX", ' ', 
             res, ' ', x, y, z, 0.0, 0.0, "", 0, false);
}

// PDBfile::WriteATOM()
void PDBfile::WriteATOM(int res, double x, double y, double z, 
                        const char* resnameIn, double Occ)
{
  WriteCoord(ATOM, anum_++, "XX", ' ', resnameIn, ' ',
             res, ' ', x, y, z, (float)Occ, 0.0, "", 0, false);
}

// PDBfile::WriteATOM()
void PDBfile::WriteATOM(const char* anameIn, int res, double x, double y, double z, 
                        const char* resnameIn, double Occ)
{
  WriteCoord(ATOM, anum_++, anameIn, ' ', resnameIn, ' ',
             res, ' ', x, y, z, (float)Occ, 0.0, "", 0, false);
}

// PDBfile::WriteCoord()
void PDBfile::WriteCoord(PDB_RECTYPE Record, int anum, NameType const& name,
                         NameType const& resnameIn, char chain, int resnum,
                         double X, double Y, double Z)
{
  WriteCoord(Record, anum, name, ' ', resnameIn, chain, 
             resnum, ' ', X, Y, Z, 0.0, 0.0, "", 0, false);
}

void PDBfile::WriteCoord(PDB_RECTYPE Record, int anum, NameType const& aname,
                         NameType const& rname, int resnum,
                         double X, double Y, double Z,
                         float Occ, float Bfac, const char* Elt, int charge)
{
  WriteCoord(Record, anum, aname, ' ', rname, ' ', resnum, ' ', X, Y, Z, 
             Occ, Bfac, Elt, charge, false);
}

// PDBfile::WriteCoord()
void PDBfile::WriteCoord(PDB_RECTYPE Record, int anum, NameType const& name,
                         char altLoc,
                         NameType const& resnameIn, char chain, int resnum,
                         char icode,
                         double X, double Y, double Z, float Occ, float B, 
                         const char* Elt, int charge, bool highPrecision) 
{
  WriteRecordHeader(Record, anum, name, altLoc, resnameIn,  chain, resnum, icode, Elt);
  if (highPrecision)
    Printf("   %8.3f%8.3f%8.3f%8.4f%8.4f      %2s%2s\n", X, Y, Z, Occ, B, Elt, "");
  else
    Printf("   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n", X, Y, Z, Occ, B, Elt, "");
}

// PDBfile::WriteANISOU()
void PDBfile::WriteANISOU(int anum, NameType const& name, 
                          NameType const& resnameIn, char chain, int resnum,
                          int u11, int u22, int u33, int u12, int u13, int u23,
                          const char* Elt, int charge)
{ // TODO icode, altLoc
  WriteRecordHeader(ANISOU, anum, name, ' ', resnameIn, chain, resnum, ' ', Elt);
  Printf(" %7i%7i%7i%7i%7i%7i      %2s%2i\n", u11, u22, u33, 
         u12, u13, u23, Elt, charge);
}

// PDBfile::WriteTITLE()
void PDBfile::WriteTITLE(std::string const& titleIn) {
  std::string titleOut;
  titleOut.reserve(70);
  int lineNum = 1;
  for (std::string::const_iterator t = titleIn.begin(); t != titleIn.end(); ++t) {
    if (titleOut.empty()) { // Start of a new line
      if (lineNum > 1) // Need to write continuation
        Printf("TITLE   %2i ", lineNum);
      else
        Printf("TITLE      ");
    }
    titleOut.push_back( *t );
    if (titleOut.size() == 69) {
      // Flush line
      Printf("%-69s\n", titleOut.c_str());
      lineNum++;
      titleOut.clear();
    }
  }
  // Flush any unwritten chars.
  if (!titleOut.empty())
    Printf("%-69s\n", titleOut.c_str());
}

/** Expect x, y, z, alpha, beta, gamma */
void PDBfile::WriteCRYST1(const double* box, const char* space_group) {
  if (box==0) return;
  // RECROD A B C ALPHA BETA GAMMA SGROUP Z
  Printf("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4i\n",
         box[0], box[1], box[2], box[3], box[4], box[5], space_group, 1);
}

/** Write MODEL record: 1-6 MODEL, 11-14 model serial # */
void PDBfile::WriteMODEL(int model) {
  // Since num frames could be large, do not format the integer with width - OK?
  Printf("MODEL     %i\n", model);
}

/** Write CONECT record: 1-6 CONECT, 7-11 serial #, 12-16 serial # ... */
void PDBfile::WriteCONECT(int atnum, std::vector<int> const& atrec, Atom const& atomIn) {
  if (atomIn.Nbonds() < 1) return;
  // PDB V3.3 spec: target-atom serial #s carried on these records also occur in increasing order.
  Atom atom(atomIn);
  atom.SortBonds();
  int nbond = 0;
  while (nbond < atom.Nbonds()) {
    if ((nbond % 4) == 0)
      Printf("CONECT%5i", atnum);
    Printf("%5i", atrec[atom.Bond(nbond++)]);
    if ((nbond % 4) == 0 || nbond == atom.Nbonds())
      Printf("\n");
  }
}

/** This version is primarily intended for writing CONECT for disulfides. */
void PDBfile::WriteCONECT(int atnum1, int atnum2) {
  Printf("CONECT%5i%5i\n", atnum1, atnum2);
}

void PDBfile::WriteSSBOND(int num, SSBOND const& ss, float distIn) {
  // TODO: SymOp
  //Printf("SSBOND %3i %3s %c %4i%c   %3s %c %4i%c                       %6i %6i %5.2f\n", num,
  Printf("SSBOND %3i %3s %c %4i%c   %3s %c %4i%c                       %6s %6s %5.2f\n", num,
         ss.name1(), ss.Chain1(), ss.Rnum1(), ss.Icode1(),
         ss.name2(), ss.Chain2(), ss.Rnum2(), ss.Icode2(), "", "", distIn);
}

void PDBfile::WriteENDMDL() { Printf("ENDMDL\n"); }

void PDBfile::WriteEND()    { Printf("END   \n"); }
// -----------------------------------------------------------------------------
PDBfile::SSBOND::SSBOND() :
  idx1_(-1), idx2_(-1), rnum1_(-1), rnum2_(-1), 
  chain1_(' '), chain2_(' '), icode1_(' '), icode2_(' ')
{
  std::fill(name1_, name1_+4, '\0');
  std::fill(name2_, name1_+4, '\0');
}

PDBfile::SSBOND::SSBOND(int idx1, int idx2, Residue const& r1, Residue const& r2) :
  idx1_(  idx1),                idx2_(  idx2),
  rnum1_( r1.OriginalResNum()), rnum2_( r2.OriginalResNum()),
  chain1_(r1.ChainID()),        chain2_(r2.ChainID()),
  icode1_(r1.Icode()),          icode2_(r2.Icode())
{
  std::copy(r1.c_str(), r1.c_str()+3, name1_);
  std::copy(r2.c_str(), r2.c_str()+3, name2_);
}

PDBfile::SSBOND::SSBOND(SSBOND const& rhs) :
  idx1_(  rhs.idx1_),   idx2_(  rhs.idx2_),
  rnum1_( rhs.rnum1_),  rnum2_( rhs.rnum2_),
  chain1_(rhs.chain1_), chain2_(rhs.chain2_),
  icode1_(rhs.icode1_), icode2_(rhs.icode2_)
{
  std::copy(rhs.name1_, rhs.name1_+3, name1_);
  std::copy(rhs.name2_, rhs.name2_+3, name2_);
}

PDBfile::SSBOND PDBfile::SSBOND::operator=(SSBOND const& rhs) {
  if (this != &rhs) {
    idx1_ = rhs.idx1_;
    idx2_ = rhs.idx2_;
    rnum1_ = rhs.rnum1_;
    rnum2_ = rhs.rnum2_;
    chain1_ = rhs.chain1_;
    chain2_ = rhs.chain2_;
    icode1_ = rhs.icode1_;
    icode2_ = rhs.icode2_;
    std::copy(rhs.name1_, rhs.name1_+3, name1_);
    std::copy(rhs.name2_, rhs.name2_+3, name2_);
  }
  return *this;
}
