#include "DataIO_AmberLib.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "AssociatedData_Connect.h"

/// CONSTRUCTOR
DataIO_AmberLib::DataIO_AmberLib()
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

}

// DataIO_AmberLib::processReadArgs()
int DataIO_AmberLib::processReadArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_AmberLib::ReadData()
int DataIO_AmberLib::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
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
  if (debug_ > 0) mprintf("DEBUG: Units:\n");
  for (Sarray::const_iterator it = UnitNames.begin(); it != UnitNames.end(); ++it)
  {
    if (debug_ > 0) mprintf("DEBUG: Reading unit %s\n", it->c_str());
    std::string entryColumn = "!entry." + *it + ".unit.atoms";
    ArgList tmparg( line );
    if (tmparg.Nargs() < 1 || tmparg[0] != entryColumn) {
      mprinterr("Error: Expected '%s', got '%s'\n", entryColumn.c_str(), line.c_str());
      return 1;
    }
    DataSet_Coords* ds = (DataSet_Coords*)dsl.AddSet( DataSet::COORDS, MetaData(dsname, *it) );
    if (ds == 0) {
      mprinterr("Error: Could not create data set for unit %s\n", it->c_str());
      return 1;
    }
    if (read_unit( ds, infile, line, *it )) {
      mprinterr("Error: Reading unit '%s'\n", it->c_str());
      return 1;
    }
  }
  mprintf("\tRead %zu units from Amber OFF file.\n", UnitNames.size());

  return 0;
}

/** Strings corresponding to SectionType */
const char* DataIO_AmberLib::sectionStr_[] = {
  "atoms", "atomspertinfo", "boundbox", "childsequence", "connect",
  "connectivity", "hierarchy", "name", "positions", "residueconnect",
  "residues", "solventcap", "velocities", 0 };

/** \return Section type from entry line. */
DataIO_AmberLib::SectionType DataIO_AmberLib::id_section(std::string const& line,
                                                         std::string const& unitName)
{
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
  if (seq-1 != topOut.Natom()) {
    mprinterr("Error: For unit %s expected sequence %i, got %i\n", unitName.c_str(), topOut.Natom()+1, seq);
    return 1;
  }
  Atom atm;
  atm.SetName( NameType(noquotes(aname)) );
  atm.SetTypeName( NameType(noquotes(atype)) );
  atm.DetermineElement( elt );
  atm.SetMassFromElement();
  atm.SetCharge( charge );
  Residue res( unitName, resx, ' ', ' ' );
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
int DataIO_AmberLib::read_connect(AssociatedData_Connect& ConnectAtoms, std::string const& line) {
  int connectAtom = -1;
  if (sscanf(line.c_str(), "%i", &connectAtom) != 1) {
    mprinterr("Error: Expected 1 column for connect line: %s\n", line.c_str());
    return 1;
  }
  // Amber lib atoms start from 1
  if (connectAtom < 1) {
    mprintf("Warning: Atom index < 1 in connect line: %s\n", line.c_str());
  }
  ConnectAtoms.AddConnectAtom( connectAtom-1 );
  return 0;
}

/** Read a unit from OFF file. It is expected that the next line from
  * infile is the first entry in the unit.atoms table.
  */
int DataIO_AmberLib::read_unit(DataSet_Coords* crd,
                               BufferedLine& infile, std::string& Line,
                               std::string const& unitName)
const
{
  AssociatedData_Connect ConnectAtoms;
  SectionType currentSection = id_section( Line, unitName );
  if (currentSection == UNKNOWN_SECTION) {
    mprinterr("Error: Could not ID first section: %s\n", Line.c_str());
    return 1;
  }
  if (debug_ > 1) mprintf("DEBUG: First section is %s\n", sectionStr_[currentSection]);

  Topology top;
  top.SetParmName( unitName, FileName() );
  std::vector<Vec3> positions;
  Frame frm;

  bool readUnit = true;
  while (readUnit) {
    const char* lineptr = infile.Line();
    if (lineptr == 0) {
      readUnit = false;
      break;
    }
    Line.assign(lineptr);
    if (!Line.empty()) {
      if (debug_ > 2) mprintf("DEBUG: Line: %s\n", Line.c_str());
      //ArgList cols( Line, " \t\n" );
      if (Line[0] == '!') {
        // See if we are at another unit
        ArgList tmparg( Line, ". " );
        if (tmparg[2] == "unit" && tmparg[3] == "atoms" && tmparg[4] == "table") {
          readUnit = false;
          break;
        }
        currentSection = id_section( Line, unitName );
        if (currentSection == UNKNOWN_SECTION) {
          mprintf("Warning: Could not ID section: %s\n", Line.c_str());
        } else {
          if (debug_ > 1) mprintf("DEBUG: Section is %s\n", sectionStr_[currentSection]);
        }
      } else if (currentSection == ATOMTABLE) {
        if (read_atoms(top, Line, unitName)) return 1;
      } else if (currentSection == CONNECTIVITY) {
        if (read_bonds(top, Line)) return 1;
      } else if (currentSection == POSITIONS) {
        if (read_positions(positions, Line)) return 1;
      } else if (currentSection == CONNECT) {
        if (read_connect(ConnectAtoms, Line)) return 1;
      }
    }
  }
  top.CommonSetup();
  if (debug_ > 1) top.Summary();
  frm.SetupFrameV( top.Atoms(), CoordinateInfo() );
  frm.ClearAtoms();
  for (std::vector<Vec3>::const_iterator it = positions.begin(); it != positions.end(); ++it)
    frm.AddVec3( *it );

  crd->CoordsSetup(top, frm.CoordsInfo());
  crd->AddFrame( frm );
  crd->AssociateData( &ConnectAtoms );
  if (debug_ > 1) ConnectAtoms.Ainfo();
  mprintf("\n");
  
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
