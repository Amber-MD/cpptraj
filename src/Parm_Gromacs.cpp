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
  mprintf("DEBUG: Opening GMX file '%s'\n", fname.c_str());
  if (infile.OpenFileRead( fname )) return 1;
  const char* ptr = infile.Line();
  while (ptr != 0) {
    if ( ptr[0] == '#' ) {
      // Directive
      std::string gmx_line( ptr );
      if ( gmx_line.compare(0, 9,"#include "       )==0 ) {
        std::string inc_fname = gmx_line.substr( 9 );
        // Remove any quotes
        std::string::iterator pend = std::remove( inc_fname.begin(), inc_fname.end(), '"' );
        size_t newsize = pend - inc_fname.begin();
        inc_fname.resize( newsize );
        if (ReadGmxFile( inc_fname )) return 1;
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
        if (!gmx_molecules_.back().empty())
          mprintf("Warning: Encountered second [ atoms ] section before [ moleculetype ]\n");
        // ;   nr    type   resnr  residu    atom    cgnr  charge
        if ( (ptr = infile.Line()) == 0 ) return 1;
        // How many columns to expect?
        int ncols = infile.TokenizeLine(SEP);
        // Which columns to expect.
        int namecol = -1;
        int typecol = -1;
        int chrgcol = -1;
        int masscol = -1;
        for (int col = 0; col < ncols; col++) {
          std::string col_label( infile.NextToken() );
          if      (col_label == "type"  ) typecol = col;
          else if (col_label == "atom"  ) namecol = col;
          else if (col_label == "charge") chrgcol = col;
          else if (col_label == "mass"  ) masscol = col;
        }
        // Require at least name/type
        if (namecol == -1 || typecol == -1) {
          mprinterr("Error: In [ atoms ], could not find either type or atom columns.\n");
          return 1;
        }
        bool readAtoms = true;
        NameType aname, atype;
        double mass = 0.0, chrg = 0.0;
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
            } else if ( currentcols != ncols ) {
              mprinterr("Error: Number of columns changes at %i (expected %i)\n",
                        infile.LineNumber(), ncols);
              return 1;
            }
            for (int col = 0; col < ncols; col++) {
              if      (col == typecol) atype = infile.NextToken();
              else if (col == namecol) aname = infile.NextToken();
              else if (col == chrgcol) chrg  = atof(infile.NextToken());
              else if (col == masscol) mass  = atof(infile.NextToken());
              else infile.NextToken(); // Blank read.
            }
            if (masscol == -1)
              gmx_molecules_.back().push_back( Atom(aname, atype, chrg) );
            else
              gmx_molecules_.back().push_back( Atom(aname, chrg, mass, atype) );
          } 
        }
        mprintf("DEBUG: Molecule %s contains %zu atoms.\n", gmx_molnames_.back().c_str(),
                gmx_molecules_.back().size());
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
