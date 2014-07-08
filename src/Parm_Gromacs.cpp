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
  mprintf("DEBUG: Opening GMX file '%s' (%i)\n", fname.c_str(), numOpen_);
  if (infile.OpenFileRead( fname )) {
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
      }
    } else if ( ptr[0] == '[' ) {
      // Bracket keyword
      std::string gmx_line( ptr );
      if ( gmx_line.compare(0,16,"[ moleculetype ]")==0 ) {
        // New molecule
        if ( (ptr = infile.Line()) == 0 ) return 1; // Read past '; name  nrexcl'
        if ( (ptr = infile.Line()) == 0 ) return 1; // Read name, nrexcl
        if (infile.TokenizeLine(SEP) < 2) {
          mprinterr("Error: After [ moleculetype ] expected name, nrexcl.\n");
          return 1;
        }
        gmx_molnames_.push_back( std::string(infile.NextToken()) );
        gmx_molecules_.push_back( AtomArray() );
      } else if ( gmx_line.compare(0, 9,"[ atoms ]"       )==0 ) {
        // Atoms for current molecule
        if (gmx_molecules_.empty()) {
          mprinterr("Error: Encountered [ atoms ] before [ moleculetype ]\n");
          return 1;
        }
        mprintf("DEBUG: Reading atoms for molecule %s\n", gmx_molnames_.back().c_str());
        if (!gmx_molecules_.back().empty())
          mprintf("Warning: Encountered second [ atoms ] section before [ moleculetype ]\n");
        // Read header.
        // ; <#> <type> <res#> <resname> <atomname> <cgnr> <charge> <mass>
        // It appears that the labels for the [ atoms ] section can vary, so
        // do not rely on them. Expect at least 7 columns, mass may be
        // ommitted and defined elsewhere.
        if ( (ptr = infile.Line()) == 0 ) return 1;
        bool readAtoms = true;
        NameType aname, atype, rname;
        double mass = 0.0, chrg = 0.0;
        int rnum = 0;
        while (readAtoms) {
          ptr = infile.Line();
          if (ptr == 0)
            readAtoms = false;
          else {
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
              gmx_molecules_.back().push_back( gmx_atom(aname, atype, rname, chrg, -1.0, rnum) );
            else // currentcols == 8
              gmx_molecules_.back().push_back( gmx_atom(aname, atype, rname, chrg, mass, rnum) );
          } 
        }
        mprintf("DEBUG: Molecule %s contains %zu atoms.\n", gmx_molnames_.back().c_str(),
                gmx_molecules_.back().size());
      } else if ( gmx_line.compare(0, 10,"[ system ]"       )==0 ) {
        // Title.
        ptr = infile.Line();
        if (ptr == 0) return 1;
        title_.assign( ptr );
      } else if ( gmx_line.compare(0, 13,"[ molecules ]"       )==0 ) {
        // System layout
        ptr = infile.Line();
        while (ptr != 0 && ptr[0] != ' ') { // TODO: Allow blank to end?
          if (infile.TokenizeLine(SEP) != 2) {
            mprinterr("Error: [ molecules ]: Line %i, expected 2 entries (<name> <count>)\n",
                      infile.LineNumber());
            return 1;
          }
          mols_.push_back( std::string(infile.NextToken()) );
          nums_.push_back( atoi(infile.NextToken()) );
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
  // Reads topology and #included files, sets up gmx_molXXX arrays.
  if (ReadGmxFile(fname)) return 1;
  // Set title/filename
  TopIn.SetParmName( title_, infileName_ );
  int resoffset = 0;
  // Set up <count> of each <molecule>
  for (unsigned int m = 0; m != mols_.size(); m++) {
    mprintf("\t%i instances of molecule %s\n", nums_[m], mols_[m].c_str());
    // Find molecule
    int tgtmol = -1;
    for (unsigned int n = 0; n != gmx_molnames_.size(); n++)
      if (gmx_molnames_[n] == mols_[m]) {
        tgtmol = (int)n;
        break;
      }
    if (tgtmol == -1) {
      mprinterr("Error: Molecule %s is not defined in gromacs topology.\n", mols_[m].c_str());
      return 1;
    }
    AtomArray const& Mol = gmx_molecules_[tgtmol];
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
      resoffset = TopIn.Nres();
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
  while (ptr != 0 && (ptr[0] == ' ' || ptr[0] == ';'))
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
