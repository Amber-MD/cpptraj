#include "Exec_AddMissingRes.h"
#include "Exec_AddMissingRes_Pres.h"
#include "BufferedLine.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords_CRD.h"
#include "ParmFile.h"
#include "Trajin_Single.h"
#include <cstring>
#include <algorithm>
#include <list>

/** Get missing residues from PDB, organize them into "gaps", i.e.
  * contiguous sequences.
  */
int Exec_AddMissingRes::FindGaps(Garray& Gaps, CpptrajFile& outfile, std::string const& pdbname)
const
{
  BufferedLine infile;
  if (infile.OpenFileRead( pdbname )) {
    mprinterr("Error: Could not open '%s' for reading.\n", pdbname.c_str());
    return 1;
  }
  const char* linePtr = infile.Line();
  int inMissing = 0;
  int nmissing = 0;
  while (linePtr != 0) {
    if (strncmp(linePtr, "REMARK", 6) == 0)
    {
      ArgList line(linePtr);
      if (line.Nargs() > 2) {
        if (inMissing == 0) {
          // MISSING section not yet encountered.
          if (line[0] == "REMARK" && line[1] == "465" && line[2] == "MISSING" && line[3] == "RESIDUES") {
            inMissing = 1;
          }
        } else if (inMissing == 1) {
          // In MISSING, looking for start of missing residues
          if (line[0] == "REMARK" && line[2] == "M") {
            inMissing = 2;
            nmissing = 0;
          }
        } else if (inMissing == 2) {
          // Reading MISSING residues
          if (line[1] != "465") {
            //mprinterr("END REACHED.\n"); // DEBUG
            break; 
          } else {
            // This is a missing residue
            nmissing++;
            //           11111111112222222
            // 012345678901234567890123456
            // REMARK 465   M RES C SSSEQI
            std::string const& currentname = line[2];
            char currentchain = linePtr[19];
            // Need to be able to parse out insertion code
            char currenticode = linePtr[26];
            int currentres;
            if (currenticode == ' ') {
              currentres = atoi(line[4].c_str());
            } else {
              char numbuf[6];
              std::copy(linePtr+21, linePtr+26, numbuf);
              numbuf[5] = '\0';
              currentres = atoi(numbuf);
            }
            mprintf("DEBUG: Missing residue %s %i icode= %c chain= %c\n",
                    currentname.c_str(), currentres, currenticode, currentchain);
            Pres thisRes(currentname, currentres, currenticode, currentchain);
            // Is this the first "gap"?
            if (Gaps.empty())
              Gaps.push_back( ResArray(1, thisRes) );
            else {
              ResArray& currentGap = Gaps.back();
              if ( currentres - currentGap.back().Onum() > 1 ||
                   currentchain != currentGap.back().Chain() )
              {
                // Starting a new "gap"
                Gaps.push_back( ResArray(1, thisRes) );
              } else {
                // Add to existing "gap"
                currentGap.push_back( thisRes );
              }
            }
          } // END missing residue
        } // END inMissing == 2
      } // END nargs > 2
    } // END REMARK
    linePtr = infile.Line();
  } // END while linePtr != 0

  // Printout
  for (Garray::const_iterator it = Gaps.begin(); it != Gaps.end(); ++it) {
    outfile.Printf("  Gap %c %4s %6i to %4s %6i %6zu\n",
                   it->front().Chain(),
                   it->front().Name().c_str(), it->front().Onum(),
                   it->back().Name().c_str(), it->back().Onum(),
                   it->size());
    // Print residues
    unsigned int col = 1;
    for (ResArray::const_iterator res = it->begin(); res != it->end(); ++res) {
      outfile.Printf("%c", Residue::ConvertResName(res->Name()));
      col++;
      if (col > 80) {
        outfile.Printf("\n");
        col = 1;
      }
    }
    if (col > 1)
      outfile.Printf("\n");
  }
  outfile.Printf("%i missing residues.\n", nmissing);
  if (Gaps.empty()) {
    mprintf("Warning: No gaps found.\n");
  }
  return 0;
}

/** Try to add in missing residues.
  * \param dataOut Output COORDS set with missing residues added in.
  * \param topIn Input Topology that is missing residues.
  * \param frameIn Input Frame that is missing residues.
  * \param Gaps Array containing info on missing residues.
  */
int Exec_AddMissingRes::AddMissingResidues(DataSet_Coords_CRD* dataOut,
                                           Topology const& topIn,
                                           Frame const& frameIn,
                                           Garray const& Gaps)
const
{
  // Use a list to preserve original topology order and make insertion easy.
  typedef std::list<Pres> ResList;
  ResList AllResidues;
  // First create a list that contains all existing residues.
  for (int rnum = 0; rnum != topIn.Nres(); ++rnum) {
    Residue const& Res = topIn.Res(rnum);
    AllResidues.push_back( Pres(Res.Name().Truncated(), Res.OriginalResNum(), rnum,
                                Res.Icode(), Res.ChainID()) );
  }
  // Sanity check
  if (AllResidues.empty()) {
    mprinterr("Error: No residues in input PDB.\n");
    return 1;
  }
  // Next, loop over gaps, add missing residues.
  ResList::iterator resPtr = AllResidues.begin();
  for (Garray::const_iterator gap = Gaps.begin(); gap != Gaps.end(); ++gap)
  {
    Pres const& gapRes0 = gap->front();
    mprintf("\tAttempting to insert gap %c %s %i to %s %i:\n", gapRes0.Chain(),
            gapRes0.Name().c_str(), gapRes0.Onum(),
            gap->back().Name().c_str(), gap->back().Onum());
    // Search until we find the correct chain
    while (resPtr->Chain() != gapRes0.Chain() && resPtr != AllResidues.end()) ++resPtr;
    if (resPtr == AllResidues.end()) {
      mprinterr("Error: Chain %c not found\n", gapRes0.Chain());
      return 1;
    }
    mprintf("\t  Chain %c found: %s %i\n", gapRes0.Chain(),
            resPtr->Name().c_str(), resPtr->Onum());
  }

  return 0;
}

// Exec_AddMissingRes::Help()
void Exec_AddMissingRes::Help() const
{
  mprintf("\tpdbname <pdbname> name <setname> [out <filename>]\n"
          "\t[parmargs <parm args>] [trajargs <trajin args>]\n"
          "\t[pdbout <pdb>] [nminsteps <nmin>] [noopt]\n");
}

// Exec_AddMissingRes::Execute()
Exec::RetType Exec_AddMissingRes::Execute(CpptrajState& State, ArgList& argIn)
{
  debug_ = State.Debug();
  if (argIn.hasKey("usenewmin")) {
    mprintf("Warning: usenewmin is deprecated.");
  }
  std::string pdbname = argIn.GetStringKey("pdbname");
  if (pdbname.empty()) {
    mprinterr("Error: provide PDB name.\n");
    return CpptrajState::ERR;
  }
  mprintf("\tPDB name: %s\n", pdbname.c_str());
  CpptrajFile* outfile = State.DFL().AddCpptrajFile(argIn.GetStringKey("out"),
                                                    "AddMissingRes", DataFileList::TEXT, true);
  if (outfile==0) {
    mprinterr("Internal Error: Unable to allocate 'out' file.\n");
    return CpptrajState::ERR;
  }
  mprintf("\tOutput file: %s\n", outfile->Filename().full());
  nMinSteps_ = argIn.getKeyInt("nminsteps", 1000);
  mprintf("\t# minimization steps: %i\n", nMinSteps_);
  optimize_ = !argIn.hasKey("noopt");
  if (optimize_)
    mprintf("\tWill attempt to optimize missing coordinates.\n");
  else
    mprintf("\tWill not attempt to optimize missing coordinates.\n");
  // Arg lists
  ArgList parmArgs;
  std::string parmArgStr = argIn.GetStringKey("parmargs");
  if (!parmArgStr.empty()) {
    parmArgs.SetList(parmArgStr, ",");
    mprintf("\tParm args: %s\n", parmArgStr.c_str());
  }
  ArgList trajArgs;
  std::string trajArgStr = argIn.GetStringKey("trajargs");
  if (!trajArgStr.empty()) {
    trajArgs.SetList(trajArgStr, ",");
    mprintf("\tTraj args: %s\n", trajArgStr.c_str());
  }
  std::string dsname = argIn.GetStringKey("name");
  if (dsname.empty())  {
    mprinterr("Error: Output set name must be specified with 'name'.\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords_CRD* dataOut = (DataSet_Coords_CRD*)State.DSL().AddSet(DataSet::COORDS, dsname);
  if (dataOut == 0) {
    mprinterr("Error: Unable to allocate output coords data set.\n");
    return CpptrajState::ERR;
  }
  mprintf("\tOutput set: %s\n", dataOut->legend());

  // Find missing residues/gaps in the PDB
  Garray Gaps;
  if (FindGaps(Gaps, *outfile, pdbname)) {
    mprinterr("Error: Finding missing residues failed.\n");
    return CpptrajState::ERR;
  }
  mprintf("\tThere are %zu gaps in the PDB.\n", Gaps.size());

  // Read in topology
  ParmFile parmIn;
  Topology topIn;
  if (parmIn.ReadTopology(topIn, pdbname, parmArgs, debug_)) {
    mprinterr("Error: Read of topology from PDB failed.\n");
    return CpptrajState::ERR;
  }
  topIn.Summary();

  // Set up input trajectory
  Trajin_Single trajIn;
  if (trajIn.SetupTrajRead(pdbname, trajArgs, &topIn)) {
    mprinterr("Error: Setup of PDB for coordinates read failed.\n");
    return CpptrajState::ERR;
  }
  trajIn.PrintInfo(1);
  // Create input frame
  Frame frameIn;
  frameIn.SetupFrameV(topIn.Atoms(), trajIn.TrajCoordInfo());
  // Read input
  if (trajIn.BeginTraj()) {
    mprinterr("Error: Opening PDB for coordinates read failed.\n");
    return CpptrajState::ERR;
  }
  trajIn.GetNextFrame(frameIn);
  trajIn.EndTraj();

  // Try to add in missing residues
  if (AddMissingResidues(dataOut, topIn, frameIn, Gaps)) {
    mprinterr("Error: Attempt to add missing residues failed.\n");
    return CpptrajState::ERR;
  }

  return CpptrajState::OK;
}
