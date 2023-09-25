#include "DataIO_AmberPrep.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "DataSet_Topology.h"
#include "StringRoutines.h"
#include "Structure/Zmatrix.h"
#include "Structure/InternalCoords.h"

/// CONSTRUCTOR
DataIO_AmberPrep::DataIO_AmberPrep()
{

}

// DataIO_AmberPrep::ID_DataFormat()
bool DataIO_AmberPrep::ID_DataFormat(CpptrajFile& infile)
{

  return false;
}

// DataIO_AmberPrep::ReadHelp()
void DataIO_AmberPrep::ReadHelp()
{

}

// DataIO_AmberPrep::processReadArgs()
int DataIO_AmberPrep::processReadArgs(ArgList& argIn)
{

  return 0;
}

static inline int CheckLine(const char* line) {
  if (line==0) {
    mprinterr("Error: Unexpected end of prep file.\n");
    return 1;
  }
  return 0;
}

// DataIO_AmberPrep::ReadData()
int DataIO_AmberPrep::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  BufferedLine infile;

  if (infile.OpenFileRead(fname)) {
    mprinterr("Error: Could not open '%s'\n", fname.full());
    return 1;
  }
  // 1 - Control for data base generation
  // IDBGEN, IREST, ITYPF
  // Format (3I)
  const char* line = infile.Line();
  if (CheckLine(line)) return 1;
  // 2 - Name of the data base file. Blank if not data base gen.
  // NAMDBF
  // Format (A80)
  line = infile.Line();
  if (CheckLine(line)) return 1;
  // 3 - Title
  // Descriptive header for the residue
  line = infile.Line();
  if (CheckLine(line)) return 1;
  mprintf("DEBUG: Prep title: '%s'\n", line);
  // 4 - Name of output file if an individual res file is being generated
  // NAMF
  // Format (A80)
  line = infile.Line();
  if (CheckLine(line)) return 1;
  // 5 - NAMRES, INTX, KFORM
  // NAMERES - a unique name for the residue 
  // INTX - Flag for the type of coords to be saved for LINK module. INT - internal, XYZ - cart.
  // KFORM - Format of output for individual residue files 0 - formatted, 1 - binary
  // FORMAT (2A, I)
  line = infile.Line();
  if (CheckLine(line)) return 1;
  ArgList args( line );
  if (args.Nargs() != 3) {
    mprinterr("Error: Expected NAMRES, INTX, KFORM: %s\n", line);
    return 1;
  }
  std::string resName = args.GetStringNext();
  std::string coordFlag = args.GetStringNext();
  int kform = args.getNextInteger(-1);
  mprintf("DEBUG: %s %s %i\n", resName.c_str(), coordFlag.c_str(), kform);
  // 6 - IFIXC , IOMIT , ISYMDU , IPOS
  // FORMAT (4A)
  // IFIXC      Flag for the type of input geometry of the residue(s)
  //    'CORRECT' The geometry is input as internal coordinates with
  //              correct order according to the tree structure.
  //              NOTE: the tree structure types ('M', 'S', etc) and order
  //              must be defined correctly: NA(I), NB(I), and NC(I) on card
  //              8 are always ignored.
  //    'CHANGE'  It is input as cartesian coordinates or part cartesian
  //              and part internal.  Cartesians should precede internals
  //              to ensure that the resulting coordinates are correct.
  //              Coordinates need not be in correct order, since each
  //              is labeled with its atom number. NOTE: NA(I), NB(I), and
  //              NC(I) on card 8 must be omitted for cartesian coordinates
  //              with this option.
  //
  // IOMIT      Flag for the omission of dummy atoms
  //    'OMIT'    dummy atoms will be deleted after generating all the
  //              information (this is used for all but the first residue
  //              in the system)
  //    'NOMIT'   they will not be deleted (dummy atoms are retained for
  //              the first residue of the system.  others are omitted)
  //
  // ISYMDU     Symbol for the dummy atoms.  The symbol must be
  //            be unique.  It is preferable to use 'DU' for it.
  //
  // IPOS       Flag for the position of dummy atoms to be deleted
  //    'ALL'     all the dummy atoms will be deleted
  //    'BEG'     only the beginning dummy atoms will be deleted
  line = infile.Line();
  if (CheckLine(line)) return 1;
  args.SetList( line, " \t" );
  if (args.Nargs() != 4) {
    mprinterr("Error: Expected IFIXC, IOMIT, ISYMDU, IPOS: %s\n", line);
    return 1;
  }
  std::string IFIXC  = args.GetStringNext();
  std::string IOMIT  = args.GetStringNext();
  std::string ISYMDU = args.GetStringNext();
  std::string IPOS   = args.GetStringNext();
  mprintf("DEBUG: %s %s %s %s\n", IFIXC.c_str(), IOMIT.c_str(), ISYMDU.c_str(), IPOS.c_str());
  if (IFIXC != "CORRECT") {
    mprinterr("Error: IFIXC is not 'CORRECT' (%s=)\n", IFIXC.c_str());
    return 1;
  }
  // 7 - CUT
  // The cutoff distance for loop closing bonds which
  // cannot be defined by the tree structure.  Any pair of
  // atoms within this distance is assumed to be bonded.
  line = infile.Line();
  if (CheckLine(line)) return 1;
  //     0 1         2         3        4     5     6     7    8        9      10
  // 8 - I IGRAPH(I) ISYMBL(I) ITREE(I) NA(I) NB(I) NC(I) R(I) THETA(I) PHI(I) CHG(I) [I = 1, NATOM]
  // FORMAT(I,3A,3I,4F)
  // Terminated by a blank line
  using namespace Cpptraj::Structure;
  Zmatrix zmatrix;
  // Topology
  DataSet* ds = dsl.AddSet( DataSet::TOPOLOGY, MetaData(dsname, "top") );
  if (ds == 0) {
    mprinterr("Error: Could not create topology for prep.\n");
    return 1;
  }
  Topology& top = ((DataSet_Topology*)ds)->ModifyTop();
  // Residue
  Residue res(resName, 1, ' ', ' ');
  // Loop over entries
  line = infile.Line();
  if (CheckLine(line)) return 1;
  while (line[0] != '\0') {
    mprintf("DEBUG: '%s'\n", line);
    args.SetList( line, " \t" );
    if (args.Nargs() != 11) {
      mprinterr("Error: Expected 11 columns in prep line, got %i\n", args.Nargs());
      return 1;
    }
    double charge = convertToDouble(args[10]);
    // Add top atom
    if (args[2] == "EP") {
      // Special case for extra points
      Atom atm(args[1], "XP");
      atm.SetTypeName(args[2]);
      atm.SetCharge(charge);
      top.AddTopAtom(atm, res);
    } else
      top.AddTopAtom( Atom(args[1], args[2], charge), res );
    // Add zmatrix entry. In prep file, indices start from 1.
    int atJ = convertToInteger(args[4]) - 1;
    int atK = convertToInteger(args[5]) - 1;
    int atL = convertToInteger(args[6]) - 1;
    double dist  = convertToDouble(args[7]);
    double theta = convertToDouble(args[8]);
    double phi   = convertToDouble(args[9]);
    zmatrix.AddIC( InternalCoords(atJ, atK, atL, dist, theta, phi) );
    line = infile.Line();
    if (line == 0) break;
  }
  // Ignore everything else for now
  infile.CloseFile();
  // Set up topology
  top.CommonSetup(true, false);
  top.Summary();
  zmatrix.print();

  return 0;
}

// DataIO_AmberPrep::WriteHelp()
void DataIO_AmberPrep::WriteHelp()
{

}

// DataIO_AmberPrep::processWriteArgs()
int DataIO_AmberPrep::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_AmberPrep::WriteData()
int DataIO_AmberPrep::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}
