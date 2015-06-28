#include <cstdlib> // atof
#include <algorithm> // std::remove
#include "Parm_Gromacs.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"

// Parm_Gromacs::ReadGmxFile()
int Parm_Gromacs::ReadGmxFile(std::string const& fname) {
  numOpen_++;
  if (numOpen_ == 100) {
    mprinterr("Error: Gromacs topology opening too many files %i; possible infinite recursion.\n",
              numOpen_);
    return 1;
  }
  const char* SEP = " \t";
  BufferedLine infile;
  if (debug_ > 0)
    mprintf("DEBUG: Opening GMX file '%s' (%i)\n", fname.c_str(), numOpen_);
  if (infile.OpenFileRead( fname )) {
    if (numOpen_ == 1) // Do not print errors for #includes
      mprinterr("Error: Could not open '%s'\n", fname.c_str());
    return 1;
  }
  if (infileName_.empty())
    infileName_ = infile.Filename(); 
  const char* ptr = infile.Line();
  while (ptr != 0) {
    if ( ptr[0] == '#' ) {
      // Directive
      std::string gmx_line( ptr );
      if ( gmx_line.compare(0, 9,"#include "       )==0 ) {
        std::string inc_fname = gmx_line.substr( 9 );
        // Remove any quotes/newline
        std::string::iterator pend = std::remove( inc_fname.begin(), inc_fname.end(), '"' );
        pend = std::remove( inc_fname.begin(), pend, '\n' );
        size_t newsize = pend - inc_fname.begin();
        inc_fname.resize( newsize );
        if (ReadGmxFile( inc_fname ))
          mprintf("Warning: Could not process '#include' directive for '%s', skipping.\n",
                  fname.c_str());
        numOpen_--;
      } /*else if ( gmx_line.compare(0, 8, "#define ") == 0 ) {
        ArgLine define( gmx_line, SEP );
        if (define.Nargs() > 1)
          defines_.push_back( define[1] );
      } else if ( gmx_line.compare(0, 7, "#ifdef ") == 0 ) {
        ArgLine ifdef( gmx_line, SEP );
        if (define.Nargs() < 2) {
          mprinterr("Error: Malformed ifdef '%s'\n", gmx_line.c_str());
          return 1;
        }
      } else if ( gmx_line.compare(0, 8, "#ifndef ") == 0) {

      } */else
        mprintf("Warning: Unsupported directive: %s\n", gmx_line.c_str());
    } else if ( ptr[0] == '[' ) {
      // Bracket keyword
      std::string gmx_line( ptr );
      if ( gmx_line.compare(0,16,"[ moleculetype ]")==0 ) {
        // New molecule. Read next line. Continue advancing past comments.
        if ( (ptr = infile.Line()) == 0 ) return 1; // Read past '; name  nrexcl'
        // If line was comment read next line.
        if (*ptr == ';') {
          if ( (ptr = infile.Line()) == 0 ) return 1; // Read name, nrexcl
        }
        if (infile.TokenizeLine(SEP) < 2) {
          mprinterr("Error: After [ moleculetype ] expected name, nrexcl.\n");
          return 1;
        }
        gmx_molecules_.push_back( gmx_mol(std::string(infile.NextToken())) );
      } else if ( gmx_line.compare(0, 9,"[ atoms ]"       )==0 ) {
        // Atoms for current molecule
        if (gmx_molecules_.empty()) {
          mprinterr("Error: Encountered [ atoms ] before [ moleculetype ]\n");
          return 1;
        }
        if (debug_ > 0)
          mprintf("DEBUG: Reading atoms for molecule %s\n", gmx_molecules_.back().Mname());
        AtomArray& MolAtoms = gmx_molecules_.back().atoms_;
        if (!MolAtoms.empty())
          mprintf("Warning: Encountered second [ atoms ] section before [ moleculetype ]\n");
        // Read header.
        // ; <#> <type> <res#> <resname> <atomname> <cgnr> <charge> <mass>
        // It appears that the labels for the [ atoms ] section can vary, so
        // do not rely on them. Expect at least 7 columns, mass may be
        // ommitted and defined elsewhere.
        if ( (ptr = infile.Line()) == 0 ) return 1;
        // If line was comment read next line.
        if (*ptr == ';') {
          if ( (ptr = infile.Line()) == 0 ) return 1; // Read first atom 
        }
        bool readAtoms = true;
        NameType aname, atype, rname;
        double mass = 0.0, chrg = 0.0;
        int rnum = 0;
        while (readAtoms) {
          if (ptr == 0)
            readAtoms = false;
          else if (*ptr != ';') {
            int currentcols = infile.TokenizeLine(SEP);
            if ( currentcols < 2 ) {
              // Assume blank line, done reading atoms.
              readAtoms = false;
              break;
            } else if ( currentcols < 7 ) {
              mprinterr("Error: Line %i: Expected at least 7 columns for [ atoms ], got %i\n",
                        infile.LineNumber(), currentcols);
              return 1;
            }
            for (int col = 0; col < currentcols; col++) {
              if      (col == 1) atype = infile.NextToken();
              else if (col == 2) rnum  = atoi(infile.NextToken());
              else if (col == 3) rname = infile.NextToken();
              else if (col == 4) aname = infile.NextToken();
              else if (col == 6) chrg  = atof(infile.NextToken());
              else if (col == 7) mass  = atof(infile.NextToken());
              else infile.NextToken(); // Blank read.
            }
            if (currentcols == 7)
              MolAtoms.push_back( gmx_atom(aname, atype, rname, chrg, -1.0, rnum) );
            else // currentcols == 8
              MolAtoms.push_back( gmx_atom(aname, atype, rname, chrg, mass, rnum) );
          } 
          ptr = infile.Line();
        }
        if (debug_ > 0)
          mprintf("DEBUG: Molecule %s contains %zu atoms.\n", gmx_molecules_.back().Mname(),
                  MolAtoms.size());
      } else if ( gmx_line.compare(0,  9,"[ bonds ]"        )==0 ) {
        // Bonds for current molecule
        if (gmx_molecules_.empty()) {
          mprinterr("Error: Encountered [ bonds ] before [ moleculetype ]\n");
          return 1;
        }
        if (debug_ > 0)
          mprintf("DEBUG: Reading bonds for molecule %s\n", gmx_molecules_.back().Mname());
        BondArray& MolBonds = gmx_molecules_.back().bonds_;
        if (!MolBonds.empty())
          mprintf("Warning: Encountered second [ bonds ] section before [ moleculetype ]\n");
        // Read header line
        // ; i     j       funct   length  force_constant
        if ( (ptr = infile.Line()) == 0 ) return 1;
        // If line was comment read next line.
        if (*ptr == ';') {
          if ( (ptr = infile.Line()) == 0 ) return 1; // Read first bond 
        }
        // Gromacs bond lengths are in nm, bond energy in kJ/mol*nm^2
        bool readBonds = true;
        while (readBonds) {
          if (ptr == 0)
            readBonds = false;
          else {
            int currentcols = infile.TokenizeLine(SEP);
            if ( currentcols < 2 ) {
              // Assume blank line, done reading atoms.
              readBonds = false;
              break;
            }
            // Only need first two columns. Internal atom #s start from 0.
            MolBonds.push_back( atoi(infile.NextToken()) - 1 );
            MolBonds.push_back( atoi(infile.NextToken()) - 1 );
          }
          ptr = infile.Line();
        }
      } else if ( gmx_line.compare(0, 10,"[ system ]"       )==0 ) {
        // Title.
        if ( (ptr = infile.Line()) == 0 ) return 1;
        // If line was comment read next line.
        if (*ptr == ';') {
          if ( (ptr = infile.Line()) == 0 ) return 1; // Read first bond 
        }
        title_.assign( ptr );
        // Remove any quotes/newline
        std::string::iterator pend = std::remove( title_.begin(), title_.end(), '\n' );
        size_t newsize = pend - title_.begin();
        title_.resize( newsize );
      } else if ( gmx_line.compare(0, 13,"[ molecules ]"       )==0 ) {
        // System layout
        ptr = infile.Line();
        while (ptr != 0 && ptr[0] != ' ') { // TODO: Allow blank to end?
          if (*ptr != ';') {
            if (infile.TokenizeLine(SEP) != 2) {
              mprinterr("Error: [ molecules ]: Line %i, expected 2 entries (<name> <count>)\n",
                        infile.LineNumber());
              return 1;
            }
            mols_.push_back( std::string(infile.NextToken()) );
            nums_.push_back( atoi(infile.NextToken()) );
          }
          ptr = infile.Line(); 
        }
      }
    } // End bracket '[' read
    // Get next line
    ptr = infile.Line();
  }
  return 0;
}
        
// Parm_Gromacs::ReadParm()
int Parm_Gromacs::ReadParm(std::string const& fname, Topology &TopIn) {
  mprintf("Warning: Currently only basic topology info read from gromacs topologies.\n");
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
                            atom->rnum_ + resoffset, atom->rname_, 0 );
        else
          TopIn.AddTopAtom( Atom( atom->aname_, atom->atype_, atom->charge_ ),
                            atom->rnum_ + resoffset, atom->rname_, 0 );
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
