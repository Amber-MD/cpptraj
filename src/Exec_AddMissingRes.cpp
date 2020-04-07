#include "Exec_AddMissingRes.h"
#include "BufferedLine.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords_CRD.h"
#include "ParmFile.h"
#include "StringRoutines.h"
#include "Trajin_Single.h"
#include <cstdlib> // atoi
#include <cstring> //strncmp

/** Get gap info from PDB */
int Exec_AddMissingRes::FindGaps(Garray& Gaps, CpptrajFile& outfile, std::string const& pdbname)
const
{
  BufferedLine infile;
  if (infile.OpenFileRead( pdbname )) return 1; 
  const char* linePtr = infile.Line();
  int inMissing = 0;
  int nmissing = 0;
  bool firstGap = true;
  bool atTheEnd = false;
  std::string lastname, lastchain;
  int lastres = 0;
  std::string currentname, currentchain;
  int currentres;
  while (linePtr != 0) {
    if (strncmp(linePtr, "REMARK", 6) == 0)
    {
      ArgList line(linePtr);
      if (line.Nargs() > 2) {
        if (inMissing == 0) {
          // MISSING section not yet encountered.
          if (line[0] == "REMARK" && line[2] == "MISSING") {
            inMissing = 1;
          }
        } else if (inMissing == 1) {
          // In MISSING, looking for start of missing residues
          if (line[0] == "REMARK" && line[2] == "M") {
            inMissing = 2;
            nmissing = 0;
            firstGap = true;
            atTheEnd = false;
          }
        } else if (inMissing == 2) {
          // Reading MISSING residues
          if (line[1] != "465") {
            atTheEnd = true;
          } else
            nmissing++;
          if (firstGap) {
            // The very first gap TODO check Gaps.empty
            Gaps.push_back( Gap(line[2], atoi(line[4].c_str()), line[3]) );
            lastname = Gaps.back().LastName();
            lastres = Gaps.back().StartRes();
            lastchain = Gaps.back().Chain();
            firstGap = false;
          } else {
            currentname = line[2];
            currentres = atoi(line[4].c_str());
            currentchain = line[3];
            if (atTheEnd || currentres - lastres > 1 || currentchain != lastchain) {
              // New sequence starting or end. Finish current.
              Gaps.back().SetStopRes(lastres);
              /*
              mprintf("  Gap %c %4s %6i to %4s %6i %6u\n",
                      Gaps.back().Chain(),
                      Gaps.back().FirstName().c_str(), Gaps.back().StartRes(),
                      Gaps.back().LastName().c_str(), Gaps.back().StopRes(),
                      Gaps.back().Nres());*/
              if (atTheEnd) {
                break;
              }
              Gaps.push_back( Gap(currentres, currentchain) );
            }
            // Continue the current sequence
            Gaps.back().AddGapRes(currentname);
            lastname = currentname;
            lastres = currentres;
            lastchain = currentchain;
          }
        } // END inMissing == 2
      } // END nargs > 2
    } // END REMARK
    linePtr = infile.Line();
  } // END while linePtr != 0

  // Printout
  for (Garray::const_iterator it = Gaps.begin(); it != Gaps.end(); ++it) {
    outfile.Printf("  Gap %c %4s %6i to %4s %6i %6u\n",
                   it->Chain(),
                   it->FirstName().c_str(), it->StartRes(),
                   it->LastName().c_str(), it->StopRes(),
                   it->Nres());
    // Print residues
    unsigned int col = 1;
    for (Gap::name_iterator name = it->nameBegin(); name != it->nameEnd(); ++name) {
      outfile.Printf("%c", Residue::ConvertResName(*name));
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

/** Try to add in missing residues. */
int Exec_AddMissingRes::AddMissingResidues(DataSet_Coords_CRD* dataOut,
                                           Topology const& topIn,
                                           Frame const& coordsIn,
                                           Garray const& Gaps)
{
  // Loop over gaps
  for (Garray::const_iterator gap = Gaps.begin(); gap != Gaps.end(); ++gap)
  {
    mprintf("\tGap %c %i to %i\n", gap->Chain(), gap->StartRes(), gap->StopRes());
    // Find start res connector
    std::string maskStr0("::" + std::string(1,gap->Chain()) + "&:;" + integerToString(gap->StartRes()-1));
    // Find stop res connector
    std::string maskStr1("::" + std::string(1,gap->Chain()) + "&:;" + integerToString(gap->StopRes()+1));

    mprintf("\t  Mask0=[%s] Mask1=[%s]\n", maskStr0.c_str(), maskStr1.c_str());
    
  }

  return 0;
}

// Exec_AddMissingRes::Help()
void Exec_AddMissingRes::Help() const
{
  mprintf("\tpdbname <pdbname> name <setname> [out <filename>]\n"
          "\t[parmargs <parm args>] [trajargs <trajin args>]\n");
}

// Exec_AddMissingRes::Execute()
Exec::RetType Exec_AddMissingRes::Execute(CpptrajState& State, ArgList& argIn)
{
  std::string pdbname = argIn.GetStringKey("pdbname");
  if (pdbname.empty()) {
    mprinterr("Error: provide PDB name.\n");
    return CpptrajState::ERR;
  }
  mprintf("\tPDB name: %s\n", pdbname.c_str());
  CpptrajFile* outfile = State.DFL().AddCpptrajFile(argIn.GetStringKey("out"),
                                                    "AddMissingRes", DataFileList::TEXT, true);
  if (outfile==0) return CpptrajState::ERR;
  mprintf("\tOutput file: %s\n", outfile->Filename().full());
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
  if (dataOut == 0) return CpptrajState::ERR;
  mprintf("\tOutput set: %s\n", dataOut->legend());

  // Find missing residues/gaps in the PDB
  Garray Gaps;
  if (FindGaps(Gaps, *outfile, pdbname))
    return CpptrajState::ERR;
  mprintf("\tThere are %zu gaps in the PDB.\n", Gaps.size());

  // Read in topology
  ParmFile parmIn;
  Topology topIn;
  if (parmIn.ReadTopology(topIn, pdbname, parmArgs, State.Debug()))
    return CpptrajState::ERR;
  topIn.Summary();

  // Set up input trajectory
  Trajin_Single trajIn;
  if (trajIn.SetupTrajRead(pdbname, trajArgs, &topIn)) return CpptrajState::ERR;
  trajIn.PrintInfo(1);
  // Create input frame
  Frame frameIn;
  frameIn.SetupFrameV(topIn.Atoms(), trajIn.TrajCoordInfo());
  // Read input
  if (trajIn.BeginTraj()) return CpptrajState::ERR;
  trajIn.GetNextFrame(frameIn);
  trajIn.EndTraj();

  // Try to add in missing residues
  if (AddMissingResidues(dataOut, topIn, frameIn, Gaps)) return CpptrajState::ERR;

  return CpptrajState::OK;
}
