#include "DataIO_AmberPrep.h"
#include "CpptrajStdio.h"
#include "BondSearch.h"
#include "BufferedLine.h"
#include "DataSet_Coords.h"
#include "StringRoutines.h"
#include "Structure/Zmatrix.h"

/// CONSTRUCTOR
DataIO_AmberPrep::DataIO_AmberPrep() :
  removeDummyAtoms_(true)
{ }

// DataIO_AmberPrep::ID_DataFormat()
bool DataIO_AmberPrep::ID_DataFormat(CpptrajFile& infile)
{

  return false;
}

// DataIO_AmberPrep::ReadHelp()
void DataIO_AmberPrep::ReadHelp()
{
  mprintf("\tkeepdummyatoms : Keep dummy atoms instead of removing them.\n");
}

// DataIO_AmberPrep::processReadArgs()
int DataIO_AmberPrep::processReadArgs(ArgList& argIn)
{
  removeDummyAtoms_ = !argIn.hasKey("keepdummyatoms");
  return 0;
}

/// Used to check for premature EOF
static inline int CheckLine(const char* line) {
  if (line==0) {
    mprinterr("Error: Unexpected end of prep file.\n");
    return 1;
  }
  return 0;
}

/** Read CHARGE section. */
int DataIO_AmberPrep::readCHARGE(BufferedLine& infile) const {
  const char* line = infile.Line();
  while (line != 0 && line[0] != '\0') {
    if (debug_ > 0) mprintf("DEBUG: CHARGE: %s\n", line);
    line = infile.Line();
  }
  return 0;
}

/** Read LOOP section. */
int DataIO_AmberPrep::readLOOP(BufferedLine& infile, AtPairArray& BondPairs) const {
  const char* line = infile.Line();
  while (line != 0 && line[0] != '\0') {
    if (debug_ > 0) mprintf("DEBUG: LOOP: %s\n", line);
    ArgList atpair(line, " \t");
    if (atpair.Nargs() != 2) {
      mprinterr("Error: Expected 2 atom names, got %i: %s\n", atpair.Nargs(), line);
      return 1;
    }
    BondPairs.push_back(AtPairType(atpair[0], atpair[1]));
    line = infile.Line();
  }
  return 0;
}

/** Read IMPROPER section. */
int DataIO_AmberPrep::readIMPROPER(BufferedLine& infile) const {
  const char* line = infile.Line();
  while (line != 0 && line[0] != '\0') {
    if (debug_ > 0) mprintf("DEBUG: IMPROPER: %s\n", line);
    line = infile.Line();
  }
  return 0;
}

/** Read 1 or more PREP residue sections (cards 3-10).
  * \return 1 if error, -1 if STOP/EOF, 0 if more residues to be read.
  */
int DataIO_AmberPrep::readAmberPrep(BufferedLine& infile, DataSetList& dsl, std::string const& dsname)
const
{
  // 3 - Title
  // Descriptive header for the residue
  const char* line = infile.Line();
  // If EOF here, assume missing STOP
  if (line == 0) {
    mprintf("Warning: Prep file '%s' missing STOP.\n", infile.Filename().base());
    return -1;
  }
  // Check for STOP
  if (std::string(line) == "STOP")
    return -1;
  if (debug_ > 0)
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
  if (debug_ > 0)
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
  if (debug_ > 0)
    mprintf("DEBUG: %s %s %s %s\n", IFIXC.c_str(), IOMIT.c_str(), ISYMDU.c_str(), IPOS.c_str());
  if (IFIXC != "CORRECT") {
    mprinterr("Error: IFIXC is not 'CORRECT' (%s=); only internal coordinates currently supported.\n",
              IFIXC.c_str());
    return 1;
  }
  // 7 - CUT
  // The cutoff distance for loop closing bonds which
  // cannot be defined by the tree structure.  Any pair of
  // atoms within this distance is assumed to be bonded.
  // TODO implement bond search
  line = infile.Line();
  if (CheckLine(line)) return 1;
  double bondCutoff = 0;
  if (line[0] != '\0')
    bondCutoff = convertToDouble( std::string(line) );
  if (bondCutoff > 0.0)
    mprintf("Warning: Non-zero cutoff in prep file (%g); will search for additional\n"
            "Warning:  bonds based on distances.\n", bondCutoff);
  //     0 1         2         3        4     5     6     7    8        9      10
  // 8 - I IGRAPH(I) ISYMBL(I) ITREE(I) NA(I) NB(I) NC(I) R(I) THETA(I) PHI(I) CHG(I) [I = 1, NATOM]
  // FORMAT(I,3A,3I,4F)
  // Terminated by a blank line
  using namespace Cpptraj::Structure;
  Zmatrix zmatrix;
  zmatrix.SetDebug( debug_ );
  // Topology
  //DataSet* ds = dsl.AddSet( DataSet::TOPOLOGY, MetaData(dsname, "top") );
  //if (ds == 0) {
  //  mprinterr("Error: Could not create topology for prep.\n");
  //  return 1;
  //}
  //Topology& top = ((DataSet_Topology*)ds)->ModifyTop();
  Topology top;
  top.SetParmName( resName, infile.Filename() );
  // Residue
  Residue res(resName, 1, ' ', ' ');
  // Loop over entries
  line = infile.Line();
  if (CheckLine(line)) return 1;
  int atIdx = 0;
  while (line[0] != '\0') {
    if (debug_ > 1)
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
    int atI = convertToInteger(args[0]) - 1;
    if (atI != atIdx) {
      mprinterr("Error: Expected index %i, got %i\n", atIdx + 1, atI + 1);
      return 1;
    }
    int atJ = convertToInteger(args[4]) - 1;
    int atK = convertToInteger(args[5]) - 1;
    int atL = convertToInteger(args[6]) - 1;
    double dist  = convertToDouble(args[7]);
    double theta = convertToDouble(args[8]);
    double phi   = convertToDouble(args[9]);
    //if (args[2] == ISYMDU) {
      if (zmatrix.AddIC( InternalCoords(atIdx, atJ, atK, atL, dist, theta, phi) ))
        return 1;
    //} else
    //  zmatrix.AddIC( InternalCoords(atIdx, atJ, atK, atL, dist, theta, phi) );
    atIdx++;
    line = infile.Line();
    if (line == 0) break;
  }
  // 9 - Read additional information about the residue.
  // CHARGE - Read additional partial atomic charges. FORMAT (5f)
  // LOOP - Read explicit loop closing bonds
  // IMPROPER - Read improper torsion angles. 
  AtPairArray BondPairs;
  while (line != 0) {
    if (line[0] != '\0') {
      std::string lineStr(line);
      if (lineStr == "DONE") {

        break;
      } else if (lineStr == "CHARGE") {
        readCHARGE(infile);
      } else if (lineStr == "LOOP") {
        if (readLOOP(infile, BondPairs)) {
          mprinterr("Error: Reading LOOP section.\n");
          return 1;
        }
      } else if (lineStr == "IMPROPER") {
        readIMPROPER(infile);
      } 
    }
    line = infile.Line();
  }
  // -----------------------------------
  // Add bonds to topology
  for (Zmatrix::const_iterator it = zmatrix.begin(); it != zmatrix.end(); ++it)
    if (it->AtJ() != InternalCoords::NO_ATOM)
      top.AddBond( it - zmatrix.begin(), it->AtJ() );
  for (AtPairArray::const_iterator it = BondPairs.begin(); it != BondPairs.end(); ++it) {
    int idx1 = top.FindAtomInResidue(0, it->first);
    if (idx1 < 0) {
      mprinterr("Error: Could not find LOOP atom '%s'\n", it->first.c_str());
      return 1;
    }
    int idx2 = top.FindAtomInResidue(0, it->second);
    if (idx2 < 0) {
      mprinterr("Error: Could not find LOOP atom '%s'\n", it->second.c_str());
      return 1;
    }
    top.AddBond( idx1, idx2 );
  }
  // Frame
  Frame frm( top.Natom() );
  // Internal coords to cart
  if (zmatrix.SetToFrame( frm )) {
    mprinterr("Error: IC to Cartesian coords failed.\n");
    return 1;
  }
  // Bond search
  if (bondCutoff > 0) {
    BondSearch bondSearch;
    bondSearch.FindBonds( top, BondSearch::SEARCH_REGULAR, frm, 0.2, debug_ );
  }
  // Set up topology
  top.CommonSetup(true, false);
  if (debug_ > 0)
    top.Summary();
  if (debug_ > 1)
    zmatrix.print();
  // Create COORDS set
  //ds = dsl.AddSet( DataSet::REF_FRAME, MetaData(dsname, "crd") );
  DataSet* ds = dsl.AddSet( DataSet::REF_FRAME, MetaData(dsname, resName) );
  if (ds == 0) {
    mprinterr("Error: Could not create coordinates for prep.\n");
    return 1;
  }
  DataSet_Coords* CRD = static_cast<DataSet_Coords*>( ds );
  
  // Delete dummy atoms if needed
  if (removeDummyAtoms_) {
    AtomMask keepAtoms;
    if (keepAtoms.SetMaskString( "!@%" + ISYMDU )) {
      mprinterr("Error: Could not set up mask string to remove dummy atoms.\n");
      return 1;
    }
    if (top.SetupIntegerMask( keepAtoms )) {
      mprinterr("Error: Could not set up mask '%s' to remove dummy atoms.\n", keepAtoms.MaskString());
      return 1;
    }
    if (debug_ > 0) keepAtoms.MaskInfo();
    if (keepAtoms.Nselected() == top.Natom()) {
      mprintf("\tNo dummy atoms to remove.\n");
    } else {
      Topology* newTop = top.modifyStateByMask( keepAtoms );
      if (newTop == 0) {
        mprinterr("Error: Could not remove dummy atoms.\n");
        return 1;
      }
      Frame newFrame( frm, keepAtoms );
      top = *newTop;
      delete newTop;
      frm = newFrame;
      if (debug_ > 0) top.Summary();
    }
  }
  // DEBUG - back convert
/*  mprintf("DEBUG: BACK CONVERT %s\n", resName.c_str());
  mprinterr("DEBUG: BACK CONVERT %s\n", resName.c_str());
  Zmatrix tempZ;
  tempZ.SetDebug(10);
  //tempZ.SetSeedPositions(frm, top, 5, 0, 1);
  int tempErr = tempZ.SetFromFrame( frm, top, 0 );
  if (tempErr != 0) {
    mprintf("DEBUG Error: Back converting %s\n", resName.c_str());
    mprinterr("DEBUG Error: Back converting %s\n", resName.c_str());
    return 1;
  }
  tempZ.print();
  Frame tmpFrame(top.Natom());
  if (tempErr == 0) {
    tempZ.SetToFrame(tmpFrame);
    for (int ii = 0; ii < tmpFrame.Natom(); ii++)
      tmpFrame.printAtomCoord(ii);
  }*/
  // Output Set up frame set
  if (CRD->CoordsSetup(top, CoordinateInfo())) {
    mprinterr("Error: Could not set up COORDS set for prep.\n");
    return 1;
  }
  CRD->SetCRD(0, frm);

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

  bool readPrep = true;
  while (readPrep) {
    int errStat = readAmberPrep(infile, dsl, dsname);
    if (errStat == 1) {
      mprinterr("Error: Could not read residue(s) from prep file.\n");
      return 1;
    } else if (errStat == -1) {
      readPrep = false;
    }
  }

  infile.CloseFile();

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
