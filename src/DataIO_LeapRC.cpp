#include "DataIO_LeapRC.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "DataIO_AmberFF.h"
#include "DataIO_AmberFrcmod.h"
#include "DataIO_AmberLib.h"
#include "DataIO_AmberPrep.h"
#include "DataIO_Coords.h"
#include "DataSet_LeapOpts.h"
#include "DataSet_NameMap.h"
#include "DataSet_Parameters.h"
#include "DataSet_PdbResMap.h"
#include "Exec_Build.h"
#include "Parm_Amber.h"
#include "StringRoutines.h" // ToLower
#include "Trajout_Single.h"
#include "Parm/AssignParams.h"
#include "Structure/Creator.h"
#include <cstdlib> //getenv
#include <cstring> // strncmp

// -----------------------------------------------------------------------------
// LeapEltHybrid class
/// CONSTRUCTOR
DataIO_LeapRC::LeapEltHybrid::LeapEltHybrid() : hybrid_(AtomType::UNKNOWN_HYBRIDIZATION) {
  elt_[0]=' ';
  elt_[1]=' ';
  elt_[2]='\0';
}

/// COPY CONSTRUCTOR
DataIO_LeapRC::LeapEltHybrid::LeapEltHybrid(LeapEltHybrid const& rhs) : hybrid_(rhs.hybrid_) {
  elt_[0] = rhs.elt_[0];
  elt_[1] = rhs.elt_[1];
}

/// ASSIGNMENT
DataIO_LeapRC::LeapEltHybrid& DataIO_LeapRC::LeapEltHybrid::operator=(const LeapEltHybrid& rhs) {
  if (&rhs == this) return *this;
  hybrid_ = rhs.hybrid_;
  elt_[0] = rhs.elt_[0];
  elt_[1] = rhs.elt_[1];
  return *this;
}

/** Print a list of supported leap commands to stdout. */
void DataIO_LeapRC::PrintSupportedLeapCommands() {
  mprintf("  Supported leap commands: addatomtypes, addpath, addpdbatommap, addpdbresmap,\n"
          "                           loadamberparams, loadamberprep, loadmol2, loadoff,\n"
          "                           loadpdb, quit, saveamberparm, set, source\n");
}

/// \return true if not equal.
bool DataIO_LeapRC::LeapEltHybrid::operator!=(const LeapEltHybrid& rhs) const {
  return ( hybrid_ != rhs.hybrid_ ||
           elt_[0] != rhs.elt_[0] ||
           elt_[1] != rhs.elt_[1] );
}

/// Set from hybridization string, element string
void DataIO_LeapRC::LeapEltHybrid::SetEltHybrid(std::string const& eltStr, std::string const& hybStr) {
  // Hybridization
  if (hybStr == "sp3")
    hybrid_ = AtomType::SP3;
  else if (hybStr == "sp2")
    hybrid_ = AtomType::SP2;
  else if (hybStr == "sp")
    hybrid_ = AtomType::SP;
  else {
    mprintf("Warning: Unknown hybridization in addAtomTypes entry %s\n", hybStr.c_str());
    hybrid_ = AtomType::UNKNOWN_HYBRIDIZATION;
  }
  elt_[0] = ' ';
  elt_[1] = ' ';
  elt_[2] = '\0';
  if (eltStr.empty()) {
    // Assume extra point
    elt_[0] = 'X'; elt_[1] = 'P';
  } else if (eltStr.size() == 1) {
    elt_[0] = eltStr[0];
  } else {
    // 2 or more characters
    elt_[0] = eltStr[0];
    elt_[1] = eltStr[1];
    if (eltStr.size() > 2)
      mprintf("Warning: Element string '%s' has more than 2 chars. Only using first 2 chars.\n", eltStr.c_str());
  }
}

/// \return Atom type hybridization
AtomType::HybridizationType DataIO_LeapRC::LeapEltHybrid::AtypeHybridization() const { return hybrid_; }

/// \return Atom type element string
const char* DataIO_LeapRC::LeapEltHybrid::AtypeElementStr() const { return elt_; }

// -----------------------------------------------------------------------------
/// CONSTRUCTOR
DataIO_LeapRC::DataIO_LeapRC() :
  leapopts_(0),
  pdbResidueMap_(0)
{}

/** Track already loaded parm files. */
DataIO_LeapRC::Sarray DataIO_LeapRC::paramFiles_ = Sarray();

/** Track already loaded lib/prep files. */
DataIO_LeapRC::Sarray DataIO_LeapRC::libFiles_ = Sarray();

// DataIO_LeapRC::ID_DataFormat()
bool DataIO_LeapRC::ID_DataFormat(CpptrajFile& infile)
{
  if (infile.OpenFile()) return false;
  bool isLeaprc = false;
  // Scan the first 5 non-blank non-comment lines
  int nLinesScanned = 0;
  while (nLinesScanned < 5 && !isLeaprc) {
    const char* ptr = infile.NextLine();
    if (ptr == 0) break;
    if (ptr[0] == '\0' || ptr[0] == '#') continue;
    nLinesScanned++;
    // LEaP commands are case-insensitive
    ArgList line(ToLower(std::string(ptr)), " \t");
    if (line.Nargs() > 0) {
      if (line[0] == "logfile")
        isLeaprc = true;
      else if (line[0] == "source")
        isLeaprc = true;
      else if (line[0] == "addatomtypes")
        isLeaprc = true;
      // TODO more commands
    }
  }
  infile.CloseFile();
  return isLeaprc;
}

// DataIO_LeapRC::ReadHelp()
void DataIO_LeapRC::ReadHelp()
{

}

// DataIO_LeapRC::processReadArgs()
int DataIO_LeapRC::processReadArgs(ArgList& argIn)
{

  return 0;
}

/** Check if file already loaded. */
bool DataIO_LeapRC::check_already_loaded(Sarray const& files, std::string const& filename) {
  for (Sarray::const_iterator it = files.begin(); it != files.end(); ++it)
    if (*it == filename) return true;
  return false;
}

/** First look for filename, then look for AMBERHOME/dir/filename. */
std::string DataIO_LeapRC::find_path(std::string const& filename,
                                     std::string const& dir)
const
{
  if (File::Exists( filename )) return filename;
  // Check searchPaths_ next
  for (Sarray::const_iterator it = searchPaths_.begin(); it != searchPaths_.end(); ++it) {
    std::string spath = *it + "/" + filename;
    if (File::Exists( spath )) return spath;
  }
  // Check AMBERHOME last
  if (amberhome_.empty()) {
    mprinterr("Error: '%s' not found.\n", filename.c_str());
    return filename;
  }
  std::string amberpath = amberhome_ + dir + filename;
  if (File::Exists( amberpath )) return amberpath;
  mprinterr("Error: '%s' not found.\n", amberpath.c_str());
  return amberpath;
}

/** LEaP addPath command. */
int DataIO_LeapRC::AddPath(std::string const& path) {
  mprintf("%s added to file search path.\n", path.c_str());
  searchPaths_.push_back( path );
  return 0;
}

/** Detect whether a file is an Amber main FF file or frcmod.
  * \return 1 if frcmod, 0 if FF, -1 if an error occurred.
  */
int DataIO_LeapRC::is_frcmod(std::string const& filename, bool& hasMass, bool& hasNonBond) {
  hasMass = false;
  hasNonBond = false;
  BufferedLine infile;
  if (infile.OpenFileRead( filename )) {
    mprinterr("Error: Could not open '%s' to determine if Amber FF/frcmod.\n", filename.c_str());
    return -1;
  }
  int stat = 0;
  const char* ptr = infile.Line();
  while (ptr != 0) {
    if ( strncmp( ptr, "MASS", 4 ) == 0 ) {
        stat = 1;
        hasMass = true;
    } else if ( strncmp( ptr, "BOND", 4 ) == 0 ) {
      stat = 1;
    } else if ( strncmp( ptr, "ANGL", 4 ) == 0 ) {
      stat = 1;
    } else if ( strncmp( ptr, "DIHE", 4 ) == 0 ) {
      stat = 1;
    } else if ( strncmp( ptr, "IMPR", 4 ) == 0 ) {
      stat = 1;
    } else if ( strncmp( ptr, "HBON", 4 ) == 0 ) {
      stat = 1;
    } else if ( strncmp( ptr, "NONB", 4 ) == 0 ) {
        stat = 1;
        hasNonBond = true;
    }
    ptr = infile.Line();
  }
  return stat;
}

/** LEaP loadAmberParams command. */
int DataIO_LeapRC::LoadAmberParams(std::string const& filename, DataSetList& dsl,
                                   std::string const& dsname,
                                   AtypeEltHybridPairMap const& atomHybridizations)
const
{
  DataSet* paramSet = 0;
  // Detect whether we have a frcmod file or not. Need to scan entire file.
  bool hasMass = false;
  bool hasNonBond = false;
  std::string full_path = find_path( filename, "parm/" );
  int stat = is_frcmod( full_path, hasMass, hasNonBond );
  if (stat == 1) {
    // Amber frcmod file
    mprintf("\tLoading force field modifications from '%s'\n", filename.c_str());
    if ( hasMass != hasNonBond ) {
      mprinterr("Error: Modified force field files must have both MASS and NONBOND entries or neither.\n"
                "Error: Could not load parameter set from %s.\n", filename.c_str());
      return 1;
    }
    DataIO_AmberFrcmod infile;
    infile.SetDebug( debug_ );
    if (infile.ReadData( full_path, dsl, dsname)) {
      mprinterr("Error: Could not load force field modifications from '%s'\n", filename.c_str());
      return 1;
    }
    if (infile.Nadded() != 1) {
      mprinterr("Internal Error: DataIO_LeapRC::LoadAmberParams(): Expected only 1 parameter set added, got %u\n", infile.Nadded());
      return 1;
    }
    paramSet = infile.added_back();
  } else if (stat == 0) {
    // Amber FF file
    if (check_already_loaded(paramFiles_, filename)) {
      mprintf("Warning: Force field %s has already been loaded, skipping.\n", filename.c_str());
      return 0;
    } else {
      mprintf("\tLoading force field from '%s'\n", filename.c_str());
      DataIO_AmberFF infile;
      if (infile.ReadData( full_path, dsl, dsname)) {
        mprinterr("Error: Could not load force field from '%s'\n", filename.c_str());
        return 1;
      }
      paramFiles_.push_back( filename );
      if (infile.Nadded() != 1) {
        mprinterr("Internal Error: DataIO_LeapRC::LoadAmberParams(): Expected only 1 parameter set added, got %u\n", infile.Nadded());
        return 1;
      }
      paramSet = infile.added_back();
    }
  } else {
    mprinterr("Error: loadamberparams failed for '%s'\n", filename.c_str());
    return 1;
  }
  if (paramSet == 0) {
    mprinterr("Internal Error: DataIO_LeapRC::LoadAmberParams(): Parameter set is null.\n");
    return 1;
  }
  // Update hybridizations for parameter atom types
  //for (DataIO::set_iterator ds = paramDSL.begin(); ds != paramDSL.end(); ++ds)
  //{
    if ( paramSet->Type() == DataSet::PARAMETERS ) {
      DataSet_Parameters& param = static_cast<DataSet_Parameters&>( *paramSet );
      mprintf("\tUpdating atom hybridizations and elements in set %s\n", param.legend());
      for (Cpptraj::Parm::ParmHolder<AtomType>::iterator it = param.AT().begin();
                                          it != param.AT().end(); ++it)
      {
        AtypeEltHybridPairMap::const_iterator ah = atomHybridizations.find( it->first[0] );
        if (ah == atomHybridizations.end())
          mprintf("Warning: No element/hybridization set for atom type '%s'\n", *(it->first[0]));
        else {
          it->second.SetHybridization( ah->second.AtypeHybridization() );
          it->second.SetEltStr( ah->second.AtypeElementStr() );
        }
      }
    } else {
      mprinterr("Internal Error: DataIO_LeapRC::LoadAmberParams(): Set %s is not parameter set.\n",
                paramSet->legend());
    }
  //}
  return 0;
}

/** LEaP loadOff command. */
int DataIO_LeapRC::LoadOFF(std::string const& filename, DataSetList& dsl, std::string const& dsname, DSarray& units) const {
  if (check_already_loaded(libFiles_, filename)) {
    mprintf("Warning: Library %s has already been loaded, skipping.\n", filename.c_str());
  } else {
    DataIO_AmberLib infile;
    infile.SetDebug( debug_ );
    // Allow lib to overwrite e.g. something from previous prep
    ArgList tmpArgs("allowoverwrite");
    infile.processReadArgs(tmpArgs);
    if (infile.ReadData( find_path(filename, "lib/"), dsl, dsname)) {
      mprinterr("Error: Could not load library file '%s'\n", filename.c_str());
      return 1;
    }
    if (infile.Nadded() < 1) {
      mprinterr("Internal Error: DataIO_LeapRC::LoadOFF(): No unit sets added.\n");
      return 1;
    }
    for (DataIO::set_iterator it = infile.added_begin(); it != infile.added_end(); ++it)
      units.push_back( *it );
    libFiles_.push_back( filename );
  }
  return 0;
}

/** LEaP loadAmberPrep command. */
int DataIO_LeapRC::LoadAmberPrep(std::string const& filename, DataSetList& dsl, std::string const& dsname, DSarray& units) const {
  if (check_already_loaded(libFiles_, filename)) {
    mprintf("Warning: Prep file %s has already been loaded, skipping.\n", filename.c_str());
  } else {
    DataIO_AmberPrep infile;
    infile.SetDebug( debug_ );
    if (infile.ReadData( find_path(filename, "prep/"), dsl, dsname)) {
      mprinterr("Error: Could not load prep file '%s'\n", filename.c_str());
      return 1;
    }
    if (infile.Nadded() < 1) {
      mprinterr("Internal Error: DataIO_LeapRC::LoadAmberPrep(): No unit sets added.\n");
      return 1;
    }
    for (DataIO::set_iterator it = infile.added_begin(); it != infile.added_end(); ++it)
      units.push_back( *it );
    libFiles_.push_back( filename );
  }
  return 0;
}

/** LEaP addAtomTypes command. */
int DataIO_LeapRC::AddAtomTypes(AtypeEltHybridPairMap& atomHybridizations, BufferedLine& infile)
const
{
  int bracketCount = 0;
  // Count duplicate entries
  Sarray duplicateEntries;
  Sarray overwrittenEntries;
  // First line should contain the command
  const char* line = infile.CurrentLine();
  while (line != 0) {
    // Process the line
    std::string tmp;
    for (const char* ptr = line; *ptr != '\0'; ++ptr)
    {
      if (*ptr == '#') {
        // Comment - skip everything else
        break;
      } else if (*ptr == '{')
        bracketCount++;
      else if (*ptr == '}') {
        bracketCount--;
        if (bracketCount == 1) {
          if (debug_ > 0) mprintf("DEBUG: addAtomTypes: %s\n", tmp.c_str());
          ArgList aline( tmp );
          // Some entries (like LP and EP) are not required to have elements.
          // Set the hybridization index to 1 or 2.
          LeapEltHybrid eltHybrid;
          if (aline.Nargs() == 3)
            eltHybrid.SetEltHybrid( aline[1], aline[2] );
          else if (aline.Nargs() == 2)
            eltHybrid.SetEltHybrid( "", aline[1] );
          else {
            mprinterr("Error: Malformed addAtomTypes entry %s\n", tmp.c_str());
            return 1;
          }
          NameType atype(aline[0]);
          AtypeEltHybridPairMap::iterator it = atomHybridizations.lower_bound( atype );
          if (it == atomHybridizations.end() || it->first != atype) {
            it = atomHybridizations.insert( it, AtypeEltHybridPairType(atype, eltHybrid) );
          } else {
            //mprintf("Warning: Duplicate entry for '%s' in addAtomTypes.", *atype);
            if (it->second != eltHybrid) {
              //mprintf(" Overwriting.\n");
              overwrittenEntries.push_back( atype.Truncated() );
              if (debug_ > 0) mprintf("Warning: Line is %s\n", tmp.c_str());
              it->second = eltHybrid;
            } else
              duplicateEntries.push_back( atype.Truncated() );
              //mprintf("\n");
          }
          tmp.clear();
        }
      } else {
        if (bracketCount == 2)
          tmp += *ptr;
      }
    }
    if (bracketCount < 0) {
      mprinterr("Error: Too many close brackets '}' in addAtomTypes command.\n");
      return 1;
    } else if (bracketCount == 0) {
      break;
    } //else {
      //mprintf("DEBUG: addAtomTypes: %s\n", tmp.c_str());
    //}
    line = infile.Line();
  }
  if (bracketCount != 0) {
    mprinterr("Error: Not enough close brackets '}' in addAtomTypes command.\n");
    return 1;
  }
  mprintf("\tRead %zu atom hybridizations.\n", atomHybridizations.size());
  if (!duplicateEntries.empty()) {
    mprintf("Warning: %zu duplicate entries in addAtomTypes; this is expected if other leaprc files have been loaded.\n",
            duplicateEntries.size());
    //if (debug_ > 0) {
      mprintf("Warning: Duplicates:");
      for (Sarray::const_iterator it = duplicateEntries.begin(); it != duplicateEntries.end(); ++it)
        mprintf(" %s", it->c_str());
      mprintf("\n");
    //}
  }
  if (!overwrittenEntries.empty()) {
    mprintf("Warning: %zu entries were overwritten in addAtomTypes:", overwrittenEntries.size());
    for (Sarray::const_iterator it = overwrittenEntries.begin(); it != overwrittenEntries.end(); ++it)
      mprintf(" %s", it->c_str());
    mprintf("\n");
  }
  return 0;
}

/** LEaP addPdbResMap command. */
int DataIO_LeapRC::AddPdbResMap(BufferedLine& infile)
const
{
  int bracketCount = 0;
  // First line should contain the command
  const char* line = infile.CurrentLine();
  while (line != 0) {
    // Process the line
    std::string tmp;
    for (const char* ptr = line; *ptr != '\0'; ++ptr)
    {
      if (*ptr == '#') {
        // Comment - skip everything else
        break;
      } else if (*ptr == '{')
        bracketCount++;
      else if (*ptr == '}') {
        bracketCount--;
        if (bracketCount == 1) {
          if (debug_ > 0) mprintf("DEBUG: addPdbResMap: %s\n", tmp.c_str());
          ArgList aline( tmp );
          // 3 tokens: terminal type (0=beg 1=end), PDB name, unit name
          if (aline.Nargs() < 2 || aline.Nargs() > 3) {
            mprinterr("Error: Malformed entry in addPdbResMap: %s\n", tmp.c_str());
            return 1;
          }
          Cpptraj::Structure::TerminalType termType = Cpptraj::Structure::NON_TERMINAL;
          int pdbidx = 0;
          int unitidx = 1;
          if (aline.Nargs() == 3) {
            if (aline[0] == "0")
              termType = Cpptraj::Structure::BEG_TERMINAL;
            else if (aline[0] == "1")
              termType = Cpptraj::Structure::END_TERMINAL;
            else
              mprintf("Warning: Unrecognized terminal type in addPdbResMap: %s\n", aline[0].c_str());
            pdbidx = 1;
            unitidx = 2;
          }
          ((DataSet_PdbResMap*)pdbResidueMap_)->AddPdbResMap( Cpptraj::PdbResMapType(aline[unitidx], aline[pdbidx], termType) );
          tmp.clear();
        }
      } else {
        if (bracketCount == 2)
          tmp += *ptr;
      }
    }
    if (bracketCount < 0) {
      mprinterr("Error: Too many close brackets '}' in addPdbResMap command.\n");
      return 1;
    } else if (bracketCount == 0) {
      break;
    } //else {
    //mprintf("DEBUG: END OF LINE: addPdbResMap: %s\n", tmp.c_str());
    //}
    line = infile.Line();
  }
  if (bracketCount != 0) {
    mprinterr("Error: Not enough close brackets '}' in addPdbResMap command.\n");
    return 1;
  }
  return 0;
}
/** LEaP addPdbAtomMap command. */
int DataIO_LeapRC::AddPdbAtomMap(std::string const& dsname, DataSetList& DSL, BufferedLine& infile)
const
{
  MetaData meta(dsname, "atommap");
  DataSet* ds = DSL.CheckForSet(meta);
  if (ds == 0) {
    ds = DSL.AddSet(DataSet::NAMEMAP, meta);
    if (ds == 0) return 1;
  }
  DataSet_NameMap& namemap = static_cast<DataSet_NameMap&>( *ds );
  if (debug_ > 0)
    mprintf("DEBUG: Name map set: %s\n", namemap.legend());
  int bracketCount = 0;
  // First line should contain the command
  const char* line = infile.CurrentLine();
  while (line != 0) {
    // Process the line
    std::string tmp;
    for (const char* ptr = line; *ptr != '\0'; ++ptr)
    {
      if (*ptr == '#') {
        // Comment - skip everything else
        break;
      } else if (*ptr == '{')
        bracketCount++;
      else if (*ptr == '}') {
        bracketCount--;
        if (bracketCount == 1) {
          if (debug_ > 0)
            mprintf("DEBUG: addPdbAtomMap: %s\n", tmp.c_str());
          ArgList aline( tmp );
          // 2 tokens: Old name, new name
          if (aline.Nargs() != 2) {
            mprinterr("Error: Malformed entry in addPdbAtomMap: %s\n", tmp.c_str());
            return 1;
          }
          if (debug_ > 0) mprintf("DEBUG: old= %s  new= %s\n", aline[0].c_str(), aline[1].c_str());
          namemap.AddNameMap( aline[0], aline[1] );
          tmp.clear();
        }
      } else {
        if (bracketCount == 2)
          tmp += *ptr;
      }
    }
    if (bracketCount < 0) {
      mprinterr("Error: Too many close brackets '}' in addPdbAtomMap command.\n");
      return 1;
    } else if (bracketCount == 0) {
      break;
    } //else {
    //mprintf("DEBUG: END OF LINE: addPdbResMap: %s\n", tmp.c_str());
    //}
    line = infile.Line();
  }
  if (bracketCount != 0) {
    mprinterr("Error: Not enough close brackets '}' in addPdbAtomMap command.\n");
    return 1;
  }
  return 0;
}

/** Load mol2 as a COORDS set */
int DataIO_LeapRC::LoadMol2(ArgList const& argIn, DataSetList& dsl) const {
  ArgList args( argIn.ArgLineStr(), " =" );
  //mprintf("DEBUG: LoadMol2\n");
  //args.PrintList();
  // Should be at least 3 args: NAME loadmol2 FILE
  if (args.Nargs() < 3) {
    mprinterr("Error: Expected at least <NAME> = loadmol2 <FILE>, got: %s\n", argIn.ArgLine());
    return 1;
  }
  DataIO_Coords coordsIn;
  coordsIn.SetDebug( debug_ );
  std::string mol2path = find_path( args[2], "" );
  if (mol2path.empty()) {
    mprinterr("Error: Could not find mol2 %s\n", args[2].c_str());
    return 1;
  }
  if (coordsIn.ReadData( mol2path, dsl, args[0] )) {
    mprinterr("Error: Could not load structure from '%s' into '%s'\n",
              args[2].c_str(), args[0].c_str());
    return 1;
  }
  mprintf("\tLoaded file '%s' into '%s'\n",
          args[2].c_str(), args[0].c_str());

  // Assume we want the remaining args to be passed to the builder.
  args.MarkArg(0);
  args.MarkArg(1);
  args.MarkArg(2);
  ArgList tmparg = args.RemainingArgs();

  Exec_Build build;
  Exec::RetType ret = build.BuildStructure( coordsIn.added_back(), dsl, debug_, "", "default_name", tmparg.hasKey("prepareforleap") );
  if (ret == CpptrajState::ERR) {
    mprinterr("Error: Build of '%s' failed.\n", args[2].c_str());
    return 1;
  }

  return 0;
}

/** Load PDB and build it. */
int DataIO_LeapRC::LoadPDB(ArgList const& argIn, DataSetList& dsl) const {
  ArgList args( argIn.ArgLineStr(), " =" );
  // Should be at least 3 args: NAME loadpdb FILE
  if (args.Nargs() < 3) {
    mprinterr("Error: Expected at least <NAME> = loadpdb <FILE>, got: %s\n", argIn.ArgLine());
    return 1;
  }
  DataIO_Coords coordsIn;
  coordsIn.SetDebug( debug_ );
//  DataSetList tmpdsl;
  std::string pdbpath = find_path( args[2], "" );
  if (pdbpath.empty()) {
    mprinterr("Error: Could not find PDB %s\n", args[2].c_str());
    return 1;
  }
  if (coordsIn.ReadData( pdbpath, dsl, args[0] )) {
    mprinterr("Error: Could not load structure from '%s' into '%s'\n",
              args[2].c_str(), args[0].c_str());
    return 1;
  }
  if (coordsIn.Nadded() != 1) {
    mprinterr("Internal Error: DataIO_LeapRC::LoadPDB(): Expected 1 COORDS set loaded from PDB, got %u\n", coordsIn.Nadded());
    return 1;
  }
  mprintf("\tLoaded file '%s' into '%s'\n", args[2].c_str(), args[0].c_str());

  // Assume we want the remaining args to be passed to the builder.
  args.MarkArg(0);
  args.MarkArg(1);
  args.MarkArg(2);
  ArgList tmparg = args.RemainingArgs();

  Exec_Build build;
  Exec::RetType ret = build.BuildStructure( coordsIn.added_back(), dsl, debug_, "", "default_name", tmparg.hasKey("prepareforleap") );
  if (ret == CpptrajState::ERR) {
    mprinterr("Error: Build of '%s' failed.\n", args[2].c_str());
    return 1;
  }

  return 0;
}

/** Leap options data set default name. */
const char* DataIO_LeapRC::LEAPOPTSNAME_ = "_LEAP_OPTIONS_";

/** Leap pdb residue map data set default name. */
const char* DataIO_LeapRC::PDBRESMAPNAME_ = "_LEAP_PDBRESMAP_";

/** LEaP 'set' command. */
int DataIO_LeapRC::LeapSet(ArgList const& argIn, DataSetList& dsl) const {
  if (leapopts_ == 0) {
    mprinterr("Internal Error: DataIO_LeapRC::LeapSet(): LeapSet() called before ReadData().\n");
    return 1;
  }
  DataSet_LeapOpts& OPTS = static_cast<DataSet_LeapOpts&>( *leapopts_ );

  if (argIn.Nargs() > 3 && ToLower(argIn[2]) == "name") {
    // set <unit> name <name>
    // Get the unit
    DataSet* ds = getUnit( argIn[1], dsl );
    if (ds == 0) {
      return 1;
    }
    DataSet_Coords& crd = static_cast<DataSet_Coords&>( *ds );
    crd.TopPtr()->SetParmTitle( argIn[3] );
  } else if (argIn.Nargs() > 3 && ToLower(argIn[1]) == "default") {
    // set default <parameter> <value>
    if ( ToLower(argIn[2]) == "pbradii" ) {
      if ( OPTS.SetGbRadii( ToLower(argIn[3]) ) ) {
        mprinterr("Error: Could not set default pbradii %s\n", argIn[3].c_str());
        return 1;
      }
    } else if ( ToLower(argIn[2]) == "cmap" ) {
      if ( ToLower(argIn[3]) == "on" ) {
        mprintf("Warning: LEaP 'set default cmap on' command ignored; CMAP terms will always be added if present.\n");
      } else if ( ToLower(argIn[3]) == "off" ) {
        mprinterr("Error: LEaP 'set default cmap off' is not supported; CMAP terms will always be added if present.\n");
        return 1;
      } else
        mprintf("Warning: Skipping unrecognized 'set default cmap' mode: %s\n", argIn[3].c_str());
    } else {
      mprintf("Warning: Unhandled 'set default' command: %s\n", argIn.ArgLine());
    }
  } else {
    mprintf("Warning: Unhandled 'set' command: %s\n", argIn.ArgLine());
  }

  return 0;
}

/** Get COORDS unit from data set list */
DataSet* DataIO_LeapRC::getUnit(std::string const& unitName, DataSetList const& dsl) const {
  //DataSet* ds = dsl.CheckForSet( MetaData(unitName) );
  DataSet* ds = dsl.GetDataSet( unitName );
  if (ds == 0) {
    mprinterr("Error: Set '%s' not found.\n", unitName.c_str());
    return 0;
  }
  if (ds->Group() != DataSet::COORDINATES) {
    mprinterr("Error: Set '%s' is not COORDINATES.\n", ds->legend());
    return 0;
  }
  mprintf("\tCOORDINATES set %s\n", ds->legend());
  return ds;
}

/** Save specified unit to a topology and restart file. */
int DataIO_LeapRC::SaveAmberParm(std::string const& unitName, ArgList& line, DataSetList& dsl)
const
{
  // saveamberparm <unit> <topfile> <restartfile>
  if (unitName.empty()) {
    mprinterr("Error: 'saveamberparm' : Unit name missing.\n");
    return 1;
  }
  std::string topName = line.GetStringNext();
  if (topName.empty()) {
    mprinterr("Error: 'saveamberparm' : Topology name missing.\n");
    return 1;
  }
  std::string crdName = line.GetStringNext();
  if (crdName.empty()) {
    mprinterr("Error: 'saveamberparm' : Coordinates name missing.\n");
    return 1;
  }
  mprintf("\tSaving unit '%s' to topology '%s' and coordinates '%s'\n",
          unitName.c_str(), topName.c_str(), crdName.c_str());
  // Get the unit
  DataSet* ds = getUnit( unitName, dsl );
//  DataSet* ds = dsl.CheckForSet( MetaData(unitName) );
  if (ds == 0) {
//    mprinterr("Error: Set '%s' not found.\n", unitName.c_str());
    return 1;
  }
//  if (ds->Group() != DataSet::COORDINATES) {
//    mprinterr("Error: Set '%s' is not COORDINATES.\n", ds->legend());
//    return 1;
//  }
//  mprintf("\tCOORDINATES set %s\n", ds->legend());
  if (ds->Size() < 1) {
    mprinterr("Error: '%s' is empty.\n", ds->legend());
    return 1;
  }
  DataSet_Coords& crd = static_cast<DataSet_Coords&>( *ds ); // FIXME this really should be const

  // Give the Topology parameters 
  ArgList buildarg = line.RemainingArgs();
  // Default GB radii
  DataSet_LeapOpts& OPTS = static_cast<DataSet_LeapOpts&>( *leapopts_ );
  Cpptraj::Parm::GB_Params gbradii;
  if (gbradii.Init_GB_Radii(OPTS.PbRadii())) return 1;
  // Get parameters
  Cpptraj::Structure::Creator creator;
  if (creator.InitCreator(buildarg, dsl, debug_)) {
    return 1;
  }
  if (!creator.HasMainParmSet()) {
    mprinterr("Error: No parameter sets.\n");
    return 1;
  }
  // Assign parameters
  Cpptraj::Parm::AssignParams AP;
  AP.SetDebug( debug_ );
  AP.SetVerbose( debug_ );
  Topology& topOut = *(crd.TopPtr());
  int nmissing = 0;
  if ( AP.AssignParameters( topOut, *(creator.MainParmSetPtr()), nmissing ) ) {
    mprinterr("Error: Could not assign parameters for '%s'.\n", topOut.c_str());
    return 1;
  }
  if (nmissing != 0) {
    mprinterr("Error: Missing %i parameters for '%s'\n", nmissing, topOut.c_str());
    return 1;
  }
  // Assign GB parameters
  gbradii.GB_Info();
  if (gbradii.Assign_GB_Radii(topOut)) {
    mprinterr("Error: Could not assign GB parameters for '%s'\n", topOut.c_str());
    return 1;
  }
 
  Parm_Amber parmOut;
  parmOut.SetDebug( debug_ );
  if (parmOut.WriteParm( topName, crd.Top() )) {
    mprinterr("Error: Write of '%s' failed.\n", topName.c_str());
    return 1;
  }

  if (crd.Size() > 1)
    mprintf("Warning: '%s' has more than 1 frame. Only using first frame.\n");

  Trajout_Single trajOut;
  ArgList tmpargs;
  if (trajOut.PrepareTrajWrite( crdName, tmpargs, dsl, crd.TopPtr(), crd.CoordsInfo(),
                                1, TrajectoryFile::UNKNOWN_TRAJ ))
  {
    mprinterr("Error: Could not set up '%s' for write.\n", crdName.c_str());
    return 1;
  }
  trajOut.PrintInfo(0);
  Frame currentFrame = crd.AllocateFrame();
  crd.GetFrame( 0, currentFrame );
  if ( trajOut.WriteSingle( 0, currentFrame )) {
    mprinterr("Error: Could not write frame to '%s'\n", crdName.c_str());
    return 1;
  }
  trajOut.EndTraj();

  return 0;
}

/// Move sets from paramDSL to dsl
/*static inline int addSetsToList(DataSetList& dsl, DataSetList& paramDSL)
{
  // Add data sets to the main data set list
  for (DataSetList::const_iterator ds = paramDSL.begin(); ds != paramDSL.end(); ++ds) {
    DataSet* dtmp = dsl.CheckForSet( (*ds)->Meta() );
    if (dtmp != 0) {
      mprinterr("Error: Set '%s' already exists.\n", (*ds)->legend());
      return 1;
    }
    dsl.AddSet( *ds );
  }
  paramDSL.SetHasCopies( true );
  return 0;
}*/

/** Read (source) a leaprc (input) file. */
int DataIO_LeapRC::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  atomHybridizations_.clear();
  units_.clear();
  // Get/allocate leap options
  leapopts_ = dsl.FindSetOfType( std::string(LEAPOPTSNAME_), DataSet::LEAPOPTS );
  if ( leapopts_ == 0)
    leapopts_ = dsl.AddSet( DataSet::LEAPOPTS, MetaData( std::string(LEAPOPTSNAME_) ) );
  if (leapopts_ == 0) {
    mprinterr("Internal Error: DataIO_LeapRC::ReadData(): Could not allocate leap options data set.\n");
    return 1;
  }
  if (debug_ > 0)
    mprintf("DEBUG: Leap options set: %s\n", leapopts_->legend());
  // Get/allocate pdb residue map
  pdbResidueMap_ = dsl.FindSetOfType( std::string(PDBRESMAPNAME_), DataSet::PDBRESMAP );
  if (pdbResidueMap_ == 0)
    pdbResidueMap_ = dsl.AddSet( DataSet::PDBRESMAP, MetaData( std::string(PDBRESMAPNAME_) ) );
  if (pdbResidueMap_ == 0) {
    mprinterr("Internal Error: DataIO_LeapRC::ReadData(): Could not allocate leap PDB residue map data set.\n");
    return 1;
  }
  if (debug_ > 0)
    mprintf("DEBUG: Leap PDB residue map set: %s\n", pdbResidueMap_->legend());

  // First, need to determine where the Amber FF files are
  const char* env = getenv("AMBERHOME");
  if (env != 0)
    amberhome_ = std::string(env) + "/dat/leap/";
  else {
    mprintf("Warning: AMBERHOME is not set. Determining FF file location based on leaprc file.\n");
    // Try to guess based on where the leaprc file is
    FileName leapcmddir( fname.DirPrefix_NoSlash() );
    if (leapcmddir.Base() == "oldff") {
      FileName leapcmddir2( leapcmddir.DirPrefix_NoSlash() );
      amberhome_ = leapcmddir2.DirPrefix();
    } else
      amberhome_ = leapcmddir.DirPrefix();
  }
  if (!amberhome_.empty())
    mprintf("\tForce field files located in '%s'\n", amberhome_.c_str());

  if (Source(fname, dsl, dsname)) {
    mprinterr("Error: Could not read LEaP input from '%s'\n", fname.full());
    return 1;
  }

  // DEBUG
  if (debug_ > 0)
    ((DataSet_PdbResMap*)pdbResidueMap_)->PrintPdbResMap();
  return 0;
}

/** Copy rhs to lhs */
int DataIO_LeapRC::unitAlias(std::string const& lhs, std::string const& rhs,
                             DataSetList& dsl, std::string const& dsname)
const
{
  if (debug_ > 0)
    mprintf("DEBUG: %s = %s\n", lhs.c_str(), rhs.c_str());
  // Find the unit to make a copy of
  DataSet* ds0 = dsl.CheckForSet( MetaData(dsname, rhs) );
  if (ds0 == 0) {
    // Its possible the unit in question is loaded by a previous command.
    ds0 = dsl.GetDataSet( "*[" + rhs + "]" );
    if (ds0 != 0) {
      mprintf("Info: Using unit '%s' from previously loaded set '%s'\n",
              rhs.c_str(), ds0->Meta().Name().c_str());
    }
  }
  if (ds0 == 0) {
    //mprinterr("Error: Could not find unit '%s' to copy to '%s'\n", rhs.c_str(), lhs.c_str());
    //return 1;
    // NOTE: Make only a warning to replicate LEAP behavior
    mprintf("Warning: Could not find unit '%s' to copy to '%s' - skipping.\n", rhs.c_str(), lhs.c_str());
    return 0;
  }
  DataSet_Coords& crd0 = static_cast<DataSet_Coords&>( *ds0 );
  // Allocate copy
  DataSet* ds1 = dsl.AddSet( DataSet::COORDS, MetaData(dsname, lhs) );
  if (ds1 == 0) {
    mprinterr("Error: Could not allocate unit '%s' for '%s'\n", lhs.c_str(), rhs.c_str());
    return 1;
  }
  DataSet_Coords& crd1 = static_cast<DataSet_Coords&>( *ds1 );
  if (crd1.CoordsSetup( crd0.Top(), crd0.CoordsInfo() )) {
    mprinterr("Error: Could not set up unit '%s' for '%s'\n", lhs.c_str(), rhs.c_str());
    return 1;
  }
  crd1.Allocate( DataSet::SizeArray(1, 1) );
  // Copy
  Frame tmpFrm = crd0.AllocateFrame();
  crd0.GetFrame(0, tmpFrm);
  crd1.SetCRD(0, tmpFrm );
  // Copy associated data
  crd1.CopyAssociatedDataFrom( crd0 );
  if (debug_ > 0)
    mprintf("DEBUG: Created unit set %s\n", crd1.legend());
  return 0;
}

/** Execute leap source command */
int DataIO_LeapRC::Source(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  BufferedLine infile;
  if (infile.OpenFileRead( find_path(fname.Full(), "cmd/") )) {
    mprinterr("Error: Could not open leaprc file '%s'\n", fname.full());
    return 1;
  }
  mprintf("\tReading LEaP input from '%s'\n", fname.base());
  enum LeapCmdType { LOADAMBERPARAMS = 0, LOADOFF, LOADAMBERPREP, ADDATOMTYPES,
                     ADDPDBRESMAP, ADDPDBATOMMAP, LOADMOL2, LOADPDB, SOURCE,
                     QUIT, SAVEAMBERPARM, ADDPATH, SET, UNKNOWN_CMD };
  int err = 0;
  const char* ptr = infile.Line();
  // FIXME need to convert to all lowercase for matching commands; leap allows
  //       mixed case
  typedef std::vector<std::string> Sarray;
  Sarray aliasPairs;
  while (ptr != 0) {
    // Use lineBuf to store all non-comment characters
    std::string lineBuf;
    if (ptr[0] != '\0' && ptr[0] != '#') {
      // Note if this line contains an equals sign.
      bool has_equals = false;
      for (const char* p = ptr; *p != '\0'; ++p) {
        if (*p == '=') {
          has_equals = true;
          //break;
        } else if (*p == '#' && p != ptr && *(p-1) != '\\') {
          // Ignore comments
          break;
        }
        lineBuf += *p;
      }

      //ArgList line( ptr, " =\t" );
      ArgList line( lineBuf, " =\t" );
      if (debug_ > 0)
        mprintf("\tLEAP> %s\n", ptr);
      //line.PrintDebug(); // DEBUG

      LeapCmdType leapcmd = UNKNOWN_CMD;
      // Look through all args, lowercase, for recognized commands.
      int pos = -1;
      for (int arg = 0; arg != line.Nargs(); arg++) {
        std::string argStr = ToLower( line[arg] );
        if      (argStr == "loadamberparams") { pos = arg; leapcmd = LOADAMBERPARAMS; break; }
        else if (argStr == "loadoff"        ) { pos = arg; leapcmd = LOADOFF; break; }
        else if (argStr == "loadamberprep"  ) { pos = arg; leapcmd = LOADAMBERPREP; break; }
        else if (argStr == "addatomtypes"   ) { pos = arg; leapcmd = ADDATOMTYPES; break; }
        else if (argStr == "addpdbresmap"   ) { pos = arg; leapcmd = ADDPDBRESMAP; break; }
        else if (argStr == "addpdbatommap"  ) { pos = arg; leapcmd = ADDPDBATOMMAP; break; }
        else if (argStr == "loadmol2"       ) { pos = arg; leapcmd = LOADMOL2; break; }
        else if (argStr == "loadpdb"        ) { pos = arg; leapcmd = LOADPDB; break; }
        else if (argStr == "source"         ) { pos = arg; leapcmd = SOURCE; break; }
        else if (argStr == "quit"           ) { pos = arg; leapcmd = QUIT; break; }
        else if (argStr == "saveamberparm"  ) { pos = arg; leapcmd = SAVEAMBERPARM; break; }
        else if (argStr == "addpath"        ) { pos = arg; leapcmd = ADDPATH; break; }
        else if (argStr == "set"            ) { pos = arg; leapcmd = SET; break; }
      }

      err = 0;
      if (leapcmd == LOADAMBERPARAMS)
        err = LoadAmberParams(line.GetStringKey(line[pos]), dsl, dsname, atomHybridizations_ );
      else if (leapcmd == SAVEAMBERPARM)
        err = SaveAmberParm(line.GetStringKey(line[pos]), line, dsl);
      else if (leapcmd == LOADOFF)
        err = LoadOFF( line.GetStringKey(line[pos]), dsl, dsname, units_ );
      else if (leapcmd == LOADAMBERPREP)
        err = LoadAmberPrep( line.GetStringKey(line[pos]), dsl, dsname, units_ );
      else if (leapcmd == ADDATOMTYPES)
        err = AddAtomTypes(atomHybridizations_, infile);
      else if (leapcmd == ADDPDBRESMAP) {
        err = AddPdbResMap(infile);
      } else if (leapcmd == ADDPDBATOMMAP)
        err = AddPdbAtomMap(dsname, dsl, infile);
      else if (leapcmd == LOADMOL2)
        err = LoadMol2(line, dsl);
      else if (leapcmd == LOADPDB)
        err = LoadPDB(line, dsl);
      else if (leapcmd == SOURCE) {
        std::string fname1 = line.GetStringKey("source");
        if (fname1.empty()) {
          mprinterr("Error: No filename given for 'source'\n");
          return 1;
        } else if (fname1 == fname.Full()) {
          mprinterr("Error: File '%s' attempting to source itself '%s'\n", 
                    fname.full(), fname1.c_str());
          return 1;
        }
        err = Source(fname1, dsl, dsname);
      } else if (leapcmd == QUIT) {
        // Do not read any more.
        mprintf("\tEncountered 'quit' in leaprc file, not reading any more.\n");
        break;
      } else if (leapcmd == ADDPATH) {
        // Add path to directories to search for files specified by other commands.
        err = AddPath( line.GetStringKey(line[pos]) );
      } else if (leapcmd == SET) {
        err = LeapSet(line, dsl);
      } else {
        // Unrecognized so far. See if this is a unit alias (interpret as 'alias = unit')
        if (has_equals && line.Nargs() == 2) {
          //err = unitAlias( line[0], line[1], dsl, dsname );
          aliasPairs.push_back( line[0] );
          aliasPairs.push_back( line[1] );
        } else {
          mprintf("Warning: Skipping unhandled LEaP command line: %s\n", ptr);
        }
      }
    }
    if (err != 0) break;
    ptr = infile.Line();
  }
  infile.CloseFile();

  // Do any unit aliases
  for (unsigned int idx = 0; idx < aliasPairs.size(); idx += 2)
  {
    if ( unitAlias( aliasPairs[idx], aliasPairs[idx+1], dsl, dsname ) ) {
      mprinterr("Error: Creating unit alias %s = %s\n", aliasPairs[idx].c_str(), aliasPairs[idx+1].c_str());
      err++;
    }
  }

  // Add data sets to the main data set list
  //if (addSetsToList(dsl, paramDSL)) return err+1;

  //if (addSetsToList(dsl, unitDSL)) return err+1;

  return err;
}

// DataIO_LeapRC::WriteHelp()
void DataIO_LeapRC::WriteHelp()
{

}

// DataIO_LeapRC::processWriteArgs()
int DataIO_LeapRC::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_LeapRC::WriteData()
int DataIO_LeapRC::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}
