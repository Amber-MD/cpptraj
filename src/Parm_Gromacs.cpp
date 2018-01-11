#include <cstdlib> // atof, getenv
#include <algorithm> // std::remove
#include "Parm_Gromacs.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include "ArgList.h"

bool Parm_Gromacs::LineContainsKey(std::string const& line, std::string const& directive)
const
{
  return (line.compare(0, directive.size(), directive)==0);
}

int Parm_Gromacs::LineContainsDirective(std::string const& line, std::string const& directive,
                                        std::string& dir_val)
{
  directive_err_ = 0;
  if (LineContainsKey(line, directive)) {
    // Remove quotes/newline from directive name.
    dir_val = line.substr( directive.size() );
    std::string::iterator pend = std::remove( dir_val.begin(), dir_val.end(), '"' );
    pend = std::remove( dir_val.begin(), pend, '\n' );
    size_t newsize = pend - dir_val.begin();
    dir_val.resize( newsize );
    if (dir_val.empty()) {
      mprinterr("Error: Malformed %sin '%s'\n", directive.c_str(), line.c_str());
      directive_err_ = 1;
    }
    return 1;
  }
  return 0;
}

int Parm_Gromacs::Defined(std::string const& dir_val) const {
  for (Sarray::const_iterator val = defines_.begin(); val != defines_.end(); ++val)
    if (dir_val == *val) return 1;
  return 0;
}

/** Seek to #else or #endif */
int Parm_Gromacs::AdvanceToElse( BufferedLine& infile ) const {
  //mprintf("DBG: Advancing to else or endif.\n");
  const char* ptr = infile.Line();
  while (ptr != 0) {
    if (*ptr == '#') {
      std::string gmx_line(ptr);
      if (LineContainsKey(gmx_line, "#else") ||
          LineContainsKey(gmx_line, "#endif"))
        return 0;
    }
    ptr = infile.Line();
  }
  mprinterr("Error: Missing #else or #endif\n");
  return 1;
}

Parm_Gromacs::KeyType Parm_Gromacs::FindKey(std::string const& line) const
{
  // Conversion to arglist will skip leading bracket and whitespace
  ArgList keyArg(line, "[ ]");
  if (keyArg.Nargs() < 1) return G_UNKNOWN_KEY;
  if      ( LineContainsKey(keyArg[0], "moleculetype"  ) ) return G_MOLECULE_TYPE;
  else if ( LineContainsKey(keyArg[0], "atoms"         ) ) return G_ATOMS;
  else if ( LineContainsKey(keyArg[0], "bonds"         ) ) return G_BONDS;
  else if ( LineContainsKey(keyArg[0], "system"        ) ) return G_SYSTEM;
  else if ( LineContainsKey(keyArg[0], "molecules"     ) ) return G_MOLECULES;
  else if ( LineContainsKey(keyArg[0], "settles"       ) ) return G_SETTLES;
  else if ( LineContainsKey(keyArg[0], "virtual_sites3") ) return G_VIRTUAL_SITES3;
  return G_UNKNOWN_KEY;
}

const char* Parm_Gromacs::SEP = " \t";

int Parm_Gromacs::ReadAtomsSection( BufferedLine& infile ) {
  // ; <#> <type> <res#> <resname> <atomname> <cgnr> <charge> <mass>
  // It appears that the labels for the [ atoms ] section can vary, so
  // do not rely on them. Expect at least 7 columns, mass may be
  // ommitted and defined elsewhere.
  if (gmx_molecules_.empty()) {
    mprinterr("Error: Encountered [ atoms ] before [ moleculetype ]\n");
    return 1;
  }
  if (debug_ > 0)
    mprintf("DEBUG: Reading atoms for molecule %s\n", gmx_molecules_.back().Mname());
  AtomArray& MolAtoms = gmx_molecules_.back().atoms_;
  if (!MolAtoms.empty())
    mprintf("Warning: Encountered second [ atoms ] section before [ moleculetype ]\n");
  NameType aname, atype, rname;
  double mass = 0.0, chrg = 0.0;
  int rnum = 0;
  int currentCols = infile.TokenizeLine(SEP);
  if ( currentCols < 7 ) {
    mprinterr("Error: Line %i: Expected at least 7 columns for [ atoms ], got %i\n",
              infile.LineNumber(), currentCols);
    return 1;
  }
  const char* ptr = infile.CurrentLine(); 
  while (ptr != 0 && currentCols > 6) {
    if (*ptr != ';') {
      for (int col = 0; col < currentCols; col++) {
        if      (col == 1) atype = infile.NextToken();
        else if (col == 2) rnum  = atoi(infile.NextToken());
        else if (col == 3) rname = infile.NextToken();
        else if (col == 4) aname = infile.NextToken();
        else if (col == 6) chrg  = atof(infile.NextToken());
        else if (col == 7) mass  = atof(infile.NextToken());
        else infile.NextToken(); // Blank read.
      }
      if (currentCols == 7)
        MolAtoms.push_back( gmx_atom(aname, atype, rname, chrg, -1.0, rnum) );
      else // currentCols == 8
        MolAtoms.push_back( gmx_atom(aname, atype, rname, chrg, mass, rnum) );
    } 
    ptr = infile.Line();
    currentCols = infile.TokenizeLine(SEP);
  }
  if (debug_ > 0)
    mprintf("DEBUG: Molecule %s contains %zu atoms.\n", gmx_molecules_.back().Mname(),
            MolAtoms.size());
  return 0;
}

int Parm_Gromacs::ReadSettles(BufferedLine& infile) {
  // ; OW    funct   doh     dhh
  // The two hydrogens are OW + 1 and OW + 2
   if (gmx_molecules_.empty()) {
    mprinterr("Error: Encountered [ settles ] before [ moleculetype ]\n");
    return 1;
  }
  // Bond OW-H1, OW-H2, H1-H2
  BondArray& MolBonds = gmx_molecules_.back().bonds_;
  // FIXME: Check if already defined?
  if (infile.TokenizeLine(SEP) < 1) return 1;
  int OW = atoi(infile.NextToken()) - 1; // Internal atom #s start from 0
  int H1 = OW + 1;
  int H2 = OW + 2;
  MolBonds.push_back( OW ); MolBonds.push_back( H1 );
  MolBonds.push_back( OW ); MolBonds.push_back( H2 );
  MolBonds.push_back( H1 ); MolBonds.push_back( H2 );
  if (debug_ > 0)
    mprintf("DEBUG: Processed [ settles ], bonds %i-%i, %i-%i, %i-%i\n",
            OW+1, H1+1, OW+1, H2+1, H1+1, H2+1);
  return 0;
}

int Parm_Gromacs::ReadVsite3(BufferedLine& infile) {
  // ; Vsite from_{i j k}                funct   a               b
   if (gmx_molecules_.empty()) {
    mprinterr("Error: Encountered [ virtual_sites3 ] before [ moleculetype ]\n");
    return 1;
  }
  BondArray& MolBonds = gmx_molecules_.back().bonds_;
  // FIXME: Check if already defined?
  int currentCols = infile.TokenizeLine(SEP);
  if (currentCols != 7) {
    mprinterr("Error: Malformed [ virtual_sites3 ]\n");
    return 1;
  }
  const char* ptr = infile.CurrentLine();
  while (ptr != 0 && currentCols == 7) {
    int site = atoi( infile.NextToken() ) - 1;
    int a0   = atoi( infile.NextToken() ) - 1;
    infile.NextToken(); // a1
    infile.NextToken(); // a2
    int func = atoi( infile.NextToken() );
    if (func != 1) {
      mprinterr("Error: Only virtual_site3 function 1 supported.\n");
      return 1;
    }
    // TODO A B
    // Bond site and a0, angle a1-a0-a2?
    MolBonds.push_back( site ); MolBonds.push_back( a0 );
    ptr = infile.Line();
    currentCols = infile.TokenizeLine(SEP);
  }
  if (debug_ > 0)
    mprintf("DEBUG: Processed [ virtual_sites3 ]\n");
  return 0;
}

int Parm_Gromacs::ReadBondsSection(BufferedLine& infile) {
  // ; i     j       funct   length  force_constant
  // Gromacs bond lengths are in nm, bond energy in kJ/mol*nm^2
  if (gmx_molecules_.empty()) {
    mprinterr("Error: Encountered [ bonds ] before [ moleculetype ]\n");
    return 1;
  }
  if (debug_ > 0)
    mprintf("DEBUG: Reading bonds for molecule %s\n", gmx_molecules_.back().Mname());
  BondArray& MolBonds = gmx_molecules_.back().bonds_;
  if (!MolBonds.empty())
    mprintf("Warning: Encountered second [ bonds ] section before [ moleculetype ]\n");
  int currentCols = infile.TokenizeLine(SEP);
  if (currentCols < 2) {
    mprinterr("Error: Empty [ bonds ] section.\n");
    return 1;
  }
  const char* ptr = infile.CurrentLine();
  while (ptr != 0 && currentCols > 1) {
    // Only need first two columns. Internal atom #s start from 0.
    MolBonds.push_back( atoi(infile.NextToken()) - 1 );
    MolBonds.push_back( atoi(infile.NextToken()) - 1 );
    ptr = infile.Line();
    currentCols = infile.TokenizeLine(SEP);
  }
  if (debug_ > 0)
    mprintf("DEBUG: Processed [ bonds ], %zi bonds.\n", MolBonds.size() / 2);
  return 0;
}

int Parm_Gromacs::ReadMolsSection(BufferedLine& infile) {
  int currentCols = infile.TokenizeLine(SEP);
  if (currentCols != 2) {
    mprinterr("Error: [ molecules ]: Line %i, expected 2 entries (<name> <count>)\n",
              infile.LineNumber());
    return 1;
  }
  const char* ptr = infile.CurrentLine();
  while (ptr != 0 && currentCols == 2) {
    mols_.push_back( std::string(infile.NextToken()) );
    nums_.push_back( atoi(infile.NextToken()) );
    ptr = infile.Line();
    if (ptr == 0) break;
    currentCols = infile.TokenizeLine(SEP);
  }
  if (debug_ > 0)
    mprintf("DEBUG: Processed [ molecules ], %zu mols.\n", mols_.size());
  return 0;
}

// Parm_Gromacs::ReadGmxFile()
int Parm_Gromacs::ReadGmxFile(FileName const& fnameIn) {
  if (fnameIn.empty()) {
    mprinterr("Error: No file name.\n");
    return 1;
  }
  numOpen_++;
  if (numOpen_ == 100) {
    mprinterr("Error: Gromacs topology opening too many files %i; possible infinite recursion.\n",
              numOpen_);
    return 1;
  }
  // First see if absolute path exists. Then look in GMXDATA/top
  FileName fname(fnameIn);
  if (!File::Exists(fname)) {
    if (!currentWorkDir_.empty()) {
      // Check w.r.t. current working directory.
      fname = FileName( currentWorkDir_ + fnameIn.Full() );
    }
    if (!File::Exists( fname )) {
      // Look in GMXDATA/top
      const char* env = getenv("GMXDATA");
      if (env != 0)
        fname = FileName( std::string(env) + "/top/" + fnameIn.Full() );
    }
  }
  BufferedLine infile;
  if (debug_ > 0)
    mprintf("DEBUG: Opening GMX file '%s' (%i)\n", fname.full(), numOpen_);
  if (infile.OpenFileRead( fname )) {
    if (numOpen_ == 1) // Do not print errors for #includes
      mprinterr("Error: Could not open '%s'\n", fname.full());
    return 1;
  }
  currentWorkDir_ = fname.DirPrefix();
  if (infileName_.empty())
    infileName_ = infile.Filename();
  bool elseActive = false; // When true, seek to #endif when #else encountered.
  std::string directive_val;
  const char* ptr = infile.Line(); // First line
  KeyType currentKey = G_UNKNOWN_KEY;
  while (ptr != 0) {
    //mprintf("[%-8i] (%1i) %s\n", infile.LineNumber(), (int)currentKey, ptr);
    // -------------------------------------------
    if ( ptr[0] == '#' ) {
      // Pre-processor directive
      std::string gmx_line( ptr );
      if ( LineContainsDirective(gmx_line, "#include ", directive_val) ) {
        //mprintf("DBG: Processing #include %s\n", directive_val.c_str());
        if (directive_err_) return 1;
        if (ReadGmxFile( directive_val ))
          mprintf("Warning: Could not process '#include' file '%s', skipping.\n",
                  directive_val.c_str());
        numOpen_--;
      } else if ( LineContainsDirective(gmx_line, "#define ", directive_val) ) {
        //mprintf("DBG: Processing #define %s\n", directive_val.c_str());
        if (directive_err_) return 1;
        defines_.push_back( directive_val );
      } else if ( LineContainsDirective(gmx_line, "#ifdef ", directive_val) ) {
        //mprintf("DBG: Processing #ifdef %s\n", directive_val.c_str());
        if (directive_err_) return 1;
        if (Defined(directive_val)) {
          //mprintf("DBG: Define '%s' is active\n", directive_val.c_str());
          elseActive = true;
        } else {
          //mprintf("DBG: Define '%s' is inactive\n", directive_val.c_str());
          if (AdvanceToElse( infile )) return 1;
        }
      } else if ( LineContainsDirective(gmx_line, "#ifndef ", directive_val) ) {
        //mprintf("DBG: Processing #ifndef %s\n", directive_val.c_str());
        if (directive_err_) return 1;
        if (!Defined(directive_val))
          elseActive = true;
        else {
          if (AdvanceToElse( infile )) return 1;
        }
      } else if ( LineContainsKey(gmx_line, "#endif") ) {
        //mprintf("DBG: Processing #endif\n");
        elseActive = false;
      } else if ( LineContainsKey(gmx_line, "#else") ) {
        //mprintf("DBG: Processing #else\n");
        if (!elseActive) {
          mprinterr("Error: #else encountered without #ifdef or #ifndef\n");
          return 1;
        }
        // Seek to endif
        while (ptr != 0) {
          if (*ptr == '#') {
            gmx_line.assign( ptr );
            if (LineContainsKey(gmx_line, "#endif")) break;
          }
          ptr = infile.Line();
        }
        if (ptr == 0) {
          mprinterr("Error: Missing #endif\n");
          return 1;
        } 
      } else
        mprintf("Warning: Unsupported directive: %s\n", gmx_line.c_str());
      if (directive_err_) return 1;
    // -------------------------------------------
    } else if ( ptr[0] == '[' ) {
      // Bracket keyword
      std::string gmx_line( ptr );
      // Determine which keyword this is
      currentKey = FindKey( gmx_line );
      // Warn for unknown key
      if (debug_ > 0 && currentKey == G_UNKNOWN_KEY)
        mprintf("Warning: Skipping section %s\n", gmx_line.c_str());
    // -------------------------------------------
    } else if ( ptr[0] != ';' && ptr[0] != '\n' && ptr[0] != '\r') {
      // Section Read
      int err = 0;
      // -------------------------------------------------
      if (currentKey == G_MOLECULE_TYPE) { // New molecule
        // ;  molname nrexcl
        if (infile.TokenizeLine(SEP) < 2) {
          mprinterr("Error: After [ moleculetype ] expected name, nrexcl.\n");
          err = 1;
        } else {
          gmx_molecules_.push_back( gmx_mol(std::string(infile.NextToken())) );
          if (debug_ > 0)
            mprintf("DEBUG: Processed [ moleculetype ] %s\n", gmx_molecules_.back().Mname());
        }
      // -------------------------------------------------
      } else if (currentKey == G_ATOMS) { // Atoms for current molecule
        err = ReadAtomsSection(infile);
      // -------------------------------------------------
      } else if (currentKey == G_BONDS) { // Bonds for current molecule
        err = ReadBondsSection(infile);
      // -------------------------------------------------
      } else if (currentKey == G_SETTLES) {
        err = ReadSettles(infile);
      // -------------------------------------------------
      } else if (currentKey == G_VIRTUAL_SITES3) {
        err = ReadVsite3(infile);
      // -------------------------------------------------
      } else if (currentKey == G_SYSTEM) { // Title
        title_.assign( ptr );
        RemoveTrailingWhitespace( title_ );
        if (debug_ > 0)
          mprintf("DEBUG: Processed [ system ]\n");
      // -------------------------------------------------
      } else if (currentKey == G_MOLECULES) { // System layout
        err = ReadMolsSection(infile);
      }
      if (err != 0) return 1;
      currentKey = G_UNKNOWN_KEY;
    } // End section read 
    // Get next line
    ptr = infile.Line();
  }
  return 0;
}
        
// Parm_Gromacs::ReadParm()
int Parm_Gromacs::ReadParm(FileName const& fname, Topology &TopIn) {
  mprintf("Warning: Currently only basic topology info (no parameters) read"
          " from gromacs topologies.\n");
  // Reads topology and #included files, sets up gmx_molXXX arrays.
  if (ReadGmxFile(fname)) return 1;
  // Set title/filename
  TopIn.SetParmName( title_, infileName_ );
  int resoffset = 0;
  int atomoffset = 0;
  // Set up <count> of each <molecule>
  for (unsigned int m = 0; m != mols_.size(); m++) {
    mprintf("\t%i instances of molecule %s\n", nums_[m], mols_[m].c_str());
    // Find molecule
    int tgtmol = -1;
    for (unsigned int n = 0; n != gmx_molecules_.size(); n++)
      if (gmx_molecules_[n].mname_ == mols_[m]) {
        tgtmol = (int)n;
        break;
      }
    if (tgtmol == -1) {
      mprinterr("Error: Molecule %s is not defined in gromacs topology.\n", mols_[m].c_str());
      return 1;
    }
    AtomArray const& Mol = gmx_molecules_[tgtmol].atoms_;
    BondArray const& Bonds = gmx_molecules_[tgtmol].bonds_;
    for (int molcount = 0; molcount != nums_[m]; molcount++) {
      for (AtomArray::const_iterator atom = Mol.begin(); atom != Mol.end(); ++atom)
      {
        if (atom->mass_ > -1.0)
          TopIn.AddTopAtom( Atom( atom->aname_, atom->charge_, atom->mass_, atom->atype_),
                            Residue(atom->rname_, atom->rnum_ + resoffset, ' ', ' ') );
        else
          TopIn.AddTopAtom( Atom( atom->aname_, atom->atype_, atom->charge_ ),
                            Residue(atom->rname_, atom->rnum_ + resoffset, ' ', ' ') );
      }
      for (BondArray::const_iterator bond = Bonds.begin(); bond != Bonds.end(); bond += 2)
        TopIn.AddBond( *bond + atomoffset , *(bond+1) + atomoffset );
      resoffset = TopIn.Nres();
      atomoffset = TopIn.Natom();
    }
  }
  return 0;
}

// Parm_Gromacs::ID_ParmFormat()
bool Parm_Gromacs::ID_ParmFormat(CpptrajFile& fileIn) {
  // Assumes already set up for READ
  if (fileIn.OpenFile()) return false;
  // Advance to first non-blank / non-comment line
  const char* ptr = fileIn.NextLine();
  while (ptr != 0 && (*ptr == ' ' || *ptr == ';' || *ptr == '\n' || *ptr == '\r'))
    ptr = fileIn.NextLine();
  bool is_gmx = false;
  if (ptr != 0) {
    // Look for an #include directive or gromacs-style directive
    std::string gmx_line( ptr );
    if      ( gmx_line.compare(0, 9,"#include "       )==0 ) is_gmx = true;
    else if ( gmx_line.compare(0,10,"[ system ]"      )==0 ) is_gmx = true;
    else if ( gmx_line.compare(0,16,"[ moleculetype ]")==0 ) is_gmx = true;
    else if ( gmx_line.compare(0,12,"[ defaults ]"    )==0 ) is_gmx = true;
    else if ( gmx_line.compare(0,13,"[ molecules ]"   )==0 ) is_gmx = true;
    else if ( gmx_line.compare(0, 9,"[ atoms ]"       )==0 ) is_gmx = true;
  }
  fileIn.CloseFile();
  return is_gmx;
}
