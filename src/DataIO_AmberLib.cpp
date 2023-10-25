#include "DataIO_AmberLib.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"

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
  mprintf("DEBUG: Units:\n");
  for (Sarray::const_iterator it = UnitNames.begin(); it != UnitNames.end(); ++it)
  {
    mprintf("DEBUG: Reading unit %s\n", it->c_str());
    std::string entryColumn = "!entry." + *it + ".unit.atoms";
    ArgList tmparg( line );
    if (tmparg.Nargs() < 1 || tmparg[0] != entryColumn) {
      mprinterr("Error: Expected '%s', got '%s'\n", entryColumn.c_str(), line.c_str());
      return 1;
    }
    read_unit( infile, line, *it );
  }


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

/** Read a unit from OFF file. It is expected that the next line from
  * infile is the first entry in the unit.atoms table.
  */
int DataIO_AmberLib::read_unit(BufferedLine& infile, std::string& Line,
                               std::string const& unitName)
const
{
  // Format: "Atom name" "Type" "Type index (unused)" "resnum" "flags" "sequence" "element" "charge"
  char aname[16];
  char atype[16];
  int typex;
  int resx;
  int flags;
  int seq;
  int elt;
  double charge;
  
  SectionType currentSection = id_section( Line, unitName );
  if (currentSection == UNKNOWN_SECTION) {
    mprinterr("Error: Could not ID first section: %s\n", Line.c_str());
    return 1;
  }
  mprintf("DEBUG: First section is %s\n", sectionStr_[currentSection]);
  bool readUnit = true;
  while (readUnit) {
    const char* lineptr = infile.Line();
    if (lineptr == 0) {
      readUnit = false;
      break;
    }
    Line.assign(lineptr);
    if (!Line.empty()) {
      mprintf("DEBUG: Line: %s\n", Line.c_str());
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
        } else
          mprintf("DEBUG: Section is %s\n", sectionStr_[currentSection]);
      }
    }
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
