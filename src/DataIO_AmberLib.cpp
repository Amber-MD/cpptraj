#include "DataIO_AmberLib.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "AssociatedData_Connect.h"

/// CONSTRUCTOR
DataIO_AmberLib::DataIO_AmberLib() :
  allowOverwrite_(false)
{

}

// DataIO_AmberLib::ID_DataFormat()
bool DataIO_AmberLib::ID_DataFormat(CpptrajFile& infile)
{
  if (infile.OpenFile()) return false;
  std::string line = infile.GetLine();
  infile.CloseFile();
  bool isLib = (line == "!!index array str");

  return isLib;
}

// DataIO_AmberLib::ReadHelp()
void DataIO_AmberLib::ReadHelp()
{
  mprintf("\tallowoverwrite : Allow existing sets to be overwritten.\n");
}

// DataIO_AmberLib::processReadArgs()
int DataIO_AmberLib::processReadArgs(ArgList& argIn)
{
  allowOverwrite_ = argIn.hasKey("allowoverwrite");
  return 0;
}

// DataIO_AmberLib::ReadData()
int DataIO_AmberLib::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  ClearAddedByMe();
  BufferedLine infile;
  if (infile.OpenFileRead( fname )) {
    mprinterr("Error: Could not open Amber lib file '%s' for reading.\n", fname.full());
    return 1;
  }

  // Read first line
  std::string line = infile.GetLine();
  if (line != "!!index array str") {
    mprinterr("Error: Expected first line to be '!!index array str', got '%s'\n", line.c_str());
    return 1;
  }
  typedef std::vector<std::string> Sarray;
  Sarray UnitNames;
  // Read units
  line = infile.GetLine();
  while (!line.empty() && line[0] != '!')
  {
    // Using ArgList here is a hacky way to get rid of the quotes
    ArgList tmparg(line);
    UnitNames.push_back( tmparg[0] );
    line = infile.GetLine();
  }
  // Now should be at first unit
  unsigned int units_read = 0;
  if (debug_ > 0) {
    mprintf("DEBUG: Units:");
    for (Sarray::const_iterator it = UnitNames.begin(); it != UnitNames.end(); ++it)
      mprintf(" %s", it->c_str());
    mprintf("\n");
  }
  while (units_read < UnitNames.size())
  {
    int is_unit = -1;
    // Which unit are we reading: !entry.<name>.unit.atoms most likely
    ArgList tmparg( line, " !." );
    if (tmparg.Nargs() >= 3) {
      if (tmparg[2] == "unit")
        is_unit = 1;
      else if (tmparg[2] == "parm")
        is_unit = 0;
    }
    std::string currentName = tmparg.GetStringKey("entry");
    if (debug_ > 0) {
       if (is_unit == 1)
         mprintf("DEBUG: Reading unit %s\n", currentName.c_str());
       else if (is_unit == 0)
         mprintf("DEBUG: Reading parameter set %s\n", currentName.c_str());
       else
         mprintf("Warning: Unknown '!entry' : %s\n", tmparg[2].c_str());
    }
    // Ensure this is in the list of entries
    bool found = false;
    for (Sarray::const_iterator it = UnitNames.begin(); it != UnitNames.end(); ++it) {
      if (*it == currentName) {
        found = true;
        break;
      }
    }
    if (!found) {
      mprinterr("Error: Entry %s was not found in the initial list of names in lib file.\n",
                currentName.c_str());
      return 1;
    }

    DataSet_Coords* ds = 0;
    if (is_unit == 1) {
      MetaData meta(dsname, currentName);
      DataSet* set = dsl.CheckForSet( meta );
      if (set == 0) {
        ds = (DataSet_Coords*)dsl.AddSet( DataSet::COORDS, meta );
      } else if (allowOverwrite_) {
        mprintf("Warning: Overwriting existing set %s\n", set->legend());
        dsl.RemoveSet( set );
        ds = (DataSet_Coords*)dsl.AddSet( DataSet::COORDS, meta );
      } else {
        mprinterr("Error: Set %s already exists and 'allowoverwrite' not specified.\n");
        return 1;
      }
      if (ds == 0) {
        mprinterr("Error: Could not create data set for unit %s\n", currentName.c_str());
        return 1;
      }
    }
    if (read_unit( ds, infile, line, currentName, (is_unit==1) )) {
      mprinterr("Error: Reading unit '%s'\n", currentName.c_str());
      return 1;
    }
    AddedByMe( ds );
    units_read++;
  }
  mprintf("\tRead %zu units from Amber OFF file %s.\n", UnitNames.size(), fname.base());

  return 0;
}

/** Strings corresponding to SectionType */
const char* DataIO_AmberLib::sectionStr_[] = {
  "atoms", "atomspertinfo", "boundbox", "childsequence", "connect",
  "connectivity", "hierarchy", "name", "positions", "residueconnect",
  "residues", "residuesPdbSequenceNumber", "solventcap", "velocities", "unknown" };

/** \return Section type from entry line. */
DataIO_AmberLib::SectionType DataIO_AmberLib::id_section(std::string const& line,
                                                         std::string const& unitName)
{
  //mprintf("DEBUG: id_section %s : %s\n", unitName.c_str(), line.c_str());
  std::string entry = "!entry." + unitName + ".unit.";
  for (int idx = 0; idx < (int)UNKNOWN_SECTION; idx++) {
    std::string sectionName = entry + std::string(sectionStr_[idx]);
    std::size_t found = line.find_first_of(" ");
    if (found == std::string::npos) {
      mprinterr("Error: Malformed entry line: %s\n", line.c_str());
      break;
    }
    if (line.compare(0, found, sectionName) == 0)
      return (SectionType)idx;
  }
  return UNKNOWN_SECTION;
}

/// Remove quotes from string
static inline std::string noquotes(const char* ptrIn) {
 std::string out;
 for (const char* ptr = ptrIn; *ptr != '\0'; ++ptr)
   if (*ptr != '"')
     out += *ptr;
  return out;
}

/** Read atoms line */
int DataIO_AmberLib::read_atoms(Topology& topOut, std::string const& line, std::string const& unitName) {
  // Format: "Atom name" "Type" "Type index (unused)" "resnum" "flags" "sequence" "element" "charge"
  char aname[16];
  char atype[16];
  int typex;
  int resx;
  int flags;
  int seq;
  int elt;
  double charge;

  if (sscanf(line.c_str(), "%s %s %i %i %i %i %i %lf",
             aname, atype, &typex, &resx, &flags, &seq, &elt, &charge) != 8)
  {
    mprinterr("Error: Expected 8 columns for atoms table line: %s\n", line.c_str());
    return 1;
  }
  // Sanity check
  //if (seq-1 != topOut.Natom()) {
  //  mprinterr("Error: For unit %s expected sequence %i, got %i\n", unitName.c_str(), topOut.Natom()+1, seq);
  //  return 1;
  //}
  Atom atm;
  atm.SetName( NameType(noquotes(aname)) );
  atm.SetTypeName( NameType(noquotes(atype)) );
  atm.DetermineElement( elt );
  atm.SetMassFromElement();
  atm.SetCharge( charge );
  // We dont know the actual residue name yet
  Residue res( "TMP", resx, ' ', "" );
  topOut.AddTopAtom( atm, res );
  return 0;
}

/** Read positions line. */
int DataIO_AmberLib::read_positions(std::vector<Vec3>& positions, std::string const& line)
{
  double x, y, z;
  if (sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &z) != 3) {
    mprinterr("Error: Expected 3 columns for positions line: %s\n", line.c_str());
    return 1;
  }
  positions.push_back( Vec3(x, y, z) );
  return 0;
}

/** Read bonds line */
int DataIO_AmberLib::read_bonds(Topology& topOut, std::string const& line) {
  int at0, at1, flags;
  if (sscanf(line.c_str(), "%i %i %i", &at0, &at1, &flags) != 3) {
    mprinterr("Error: Expected 3 columns for connectivity line: %s\n", line.c_str());
    return 1;
  }
  topOut.AddBond( at0-1, at1-1 );
  return 0;
}

/** Read connections line */
int DataIO_AmberLib::read_connect(AssociatedData_Connect& ConnectAtoms, std::string const& line) const {
  int connectAtom = -1;
  if (sscanf(line.c_str(), "%i", &connectAtom) != 1) {
    mprinterr("Error: Expected 1 column for connect line: %s\n", line.c_str());
    return 1;
  }
  // Amber lib atoms start from 1
  if (connectAtom < 1) {
    // This may be legit when only one connect atom
    if (debug_ > 0)
      mprintf("Warning: Atom index < 1 in connect line: %s\n", line.c_str());
  }
  ConnectAtoms.AddConnectAtom( connectAtom-1 );
  return 0;
}

/** Hold information for a unit entry. */
class unitEntry {
  public:
    /// CONSTRUCTOR - Take unit name
    unitEntry(std::string const& unitName) : ridx_(0) {
      top_.SetParmName(unitName, FileName());
      boundbox_.reserve(5);
    }
    /// \return Reference to topology
    Topology& Top() { return top_; }
    /// \return Reference to Positions
    std::vector<Vec3>& Positions() { return positions_; }
    /// \return Reference to connect atoms
    AssociatedData_Connect& ConnectAtoms() { return ConnectAtoms_; }
    /// Add residue
    int AddResname( std::string const& resname ) {
      if (ridx_ >= top_.Nres()) {
        mprinterr("Error: Too many residues in residues section, or residues section before atom table.\n");
        return 1;
      }
      top_.SetRes( ridx_++ ).SetName( resname );
      return 0;
    }
    /// Add box value
    void AddBoxVal(double d) { boundbox_.push_back( d ); }
    /// Finish setting up the unit
    int FinishUnit(DataSet_Coords* crd, int debug_) {
      // Set up topology; determine molecules, but no residue renumber
      top_.CommonSetup(true, false);
      if (debug_ > 1) top_.Summary();
      Frame frm;
      frm.SetupFrameV( top_.Atoms(), CoordinateInfo() );
      frm.ClearAtoms();
      for (std::vector<Vec3>::const_iterator it = positions_.begin(); it != positions_.end(); ++it)
        frm.AddVec3( *it );
      // Set up box if needed
      if (!boundbox_.empty()) {
        if (boundbox_.size() != 5) {
          mprinterr("Error: Expected 5 values for boundbox entry for %s, got %zu\n", top_.c_str(), boundbox_.size());
          return 1;
        }
        // hasbox, angle, x, y, z
        if (boundbox_[0] > 0) {
          double boxcrd[6];
          boxcrd[0] = boundbox_[2];
          boxcrd[1] = boundbox_[3];
          boxcrd[2] = boundbox_[4];
          // determine angles
          double beta = boundbox_[1];
          if (beta < 3.15) {
            // Assume the angle is given in radians
            beta = beta * Constants::RADDEG;
            //mprintf("DEBUG: Converted beta: %f\n", beta);
            // Fix imperfectly converted orthorhombic
            double deltaBeta = 90.0 - beta;
            if (deltaBeta < 0) deltaBeta = -deltaBeta;
            if (deltaBeta < 0.0001) {
              beta = 90.0;
              //mprintf("DEBUG: Fixed beta.\n");
            }
          }
          boxcrd[3] = beta;
          boxcrd[4] = beta;
          boxcrd[5] = beta;
          if ( Box::IsTruncOct( beta ) ) {
            // Use trunc oct angle from Box; higher precision
            boxcrd[3] = Box::TruncatedOctAngle();
            boxcrd[4] = boxcrd[3];
            boxcrd[5] = boxcrd[3];
          } else if (beta == 60.0) {
            // Special case - rhombic dodecahedron
            boxcrd[3] = 60.0;
            boxcrd[4] = 90.0;
            boxcrd[5] = 60.0;
          }
          frm.ModifyBox().SetupFromXyzAbg( boxcrd );
          if (debug_ > 0)
            frm.BoxCrd().PrintInfo(); // DEBUG
        } else {
          frm.ModifyBox().SetNoBox();
        }
      }
      crd->CoordsSetup(top_, frm.CoordsInfo());
      crd->AddFrame( frm );
      crd->AssociateData( &ConnectAtoms_ );
      if (debug_ > 1) ConnectAtoms_.Ainfo();
      if (debug_ > 0) mprintf("\n");
      return 0;
    }
  private:
    AssociatedData_Connect ConnectAtoms_;
    Topology top_;
    std::vector<Vec3> positions_;
    int ridx_;
    std::vector<double> boundbox_; ///< hold values from box entry: hasbox, angle, x, y, zi
};

/** Read a unit from OFF file. It is expected that the next line from
  * infile is the first entry in the unit.atoms table.
  */
int DataIO_AmberLib::read_unit(DataSet_Coords* crd,
                               BufferedLine& infile, std::string& Line,
                               std::string const& unitName,
                               bool isUnit)
const
{
  SectionType currentSection = id_section( Line, unitName );
  //if (currentSection == UNKNOWN_SECTION) {
  //  mprinterr("Error: Could not ID first section for entry '%s': %s\n", unitName.c_str(), Line.c_str());
  //  return 1;
  //}
  //if (debug_ > 1) mprintf("DEBUG: First section is %s\n", sectionStr_[currentSection]);

  //Frame frm;
  unitEntry currentUnit( unitName );

  bool readUnit = true;
  while (readUnit) {
    const char* lineptr = infile.Line();
    if (lineptr == 0) {
      readUnit = false;
      break;
    }
    // Advance past any leading whitespace
    while (*lineptr != 0 && (*lineptr == ' ' || *lineptr == '\t'))
      ++lineptr;
    if (*lineptr == 0)
      Line.clear();
    else
      Line.assign(lineptr);
    if (!Line.empty()) {
      if (debug_ > 2) mprintf("DEBUG: Section '%s' Line: '%s'\n", sectionStr_[currentSection], Line.c_str());
      //ArgList cols( Line, " \t\n" );
      if (Line[0] == '!') {
        // See if we are at another unit
        ArgList tmparg( Line, ". " );
        if (tmparg[1] != unitName) {
          //mprintf("DEBUG: New entry %s detected.\n", tmparg[1].c_str());
          readUnit = false;
          break;
        }
/*        if (tmparg[2] == "parm" && tmparg[4] == "table") {
          // This is a PARM table, not a unit.
          mprintf("Parm section %s detected.\n", tmparg[3].c_str());
          readUnit = false;
          break;
        }
        if (tmparg[2] == "unit" && tmparg[3] == "atoms" && tmparg[4] == "table") {
          readUnit = false;
          break;
        }*/
        currentSection = id_section( Line, unitName );
        if (currentSection == UNKNOWN_SECTION) {
          mprintf("Warning: Could not ID section: %s\n", Line.c_str());
        } else {
          if (debug_ > 1) mprintf("DEBUG: Section is %s\n", sectionStr_[currentSection]);
        }
      } else if (currentSection == ATOMTABLE) {
        if (read_atoms(currentUnit.Top(), Line, unitName)) return 1;
      } else if (currentSection == CONNECTIVITY) {
        if (read_bonds(currentUnit.Top(), Line)) return 1;
      } else if (currentSection == POSITIONS) {
        if (read_positions(currentUnit.Positions(), Line)) return 1;
      } else if (currentSection == CONNECT) {
        if (read_connect(currentUnit.ConnectAtoms(), Line)) return 1;
      } else if (currentSection == RESIDUES) {
        // Rely on ArgList to remove quotes
        ArgList tmpArg( Line );
        if (tmpArg.Nargs() < 1) {
          mprinterr("Error: Could not read residue from residues section.\n");
          mprinterr("Error: Line: %s\n", Line.c_str());
          return 1;
        }
        if (currentUnit.AddResname( tmpArg[0] )) return 1;
      } else if (currentSection == BOUNDBOX) {
        // hasbox, angle, x, y, z
        double bval = 0;
        sscanf( Line.c_str(), "%lf", &bval );
        currentUnit.AddBoxVal( bval );
      }
      // else {
      //  mprintf("Warning: Unhandled section for %s: %s\n", unitName.c_str(), sectionStr_[currentSection]);
      //}
    }
  }

  if (isUnit) {
    if (currentUnit.FinishUnit(crd, debug_)) return 1;
  }

  return 0;
}

// DataIO_AmberLib::WriteHelp()
void DataIO_AmberLib::WriteHelp()
{

}

// DataIO_AmberLib::processWriteArgs()
int DataIO_AmberLib::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_AmberLib::WriteData()
int DataIO_AmberLib::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}
