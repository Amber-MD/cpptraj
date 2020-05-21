#include "Exec_AddMissingRes.h"
#include "Exec_AddMissingRes_Pres.h"
#include "BufferedLine.h"
#include "CharMask.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords_CRD.h"
#include "ParmFile.h"
#include "Trajin_Single.h"
#include "Trajout_Single.h"
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

/** Write topology and frame to pdb. */
int Exec_AddMissingRes::WriteStructure(std::string const& fname, Topology* newTop, Frame const& newFrame,
                                       TrajectoryFile::TrajFormatType typeOut)
const
{
  Trajout_Single trajOut;
  if (trajOut.InitTrajWrite(fname, ArgList("pdbter"), DataSetList(), typeOut))
    return 1;
  if (trajOut.SetupTrajWrite(newTop, CoordinateInfo(), 1))
    return 1;
  if (trajOut.WriteSingle(0, newFrame)) return 1;
  trajOut.EndTraj();
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
    Pres const& gapRes1 = gap->back();
    mprintf("\tAttempting to insert gap %c %s %i to %s %i:\n", gapRes0.Chain(),
            gapRes0.Name().c_str(), gapRes0.Onum(),
            gapRes1.Name().c_str(), gapRes1.Onum());
    // Search until we find the correct chain
    while (resPtr != AllResidues.end() && resPtr->Chain() != gapRes0.Chain())
      ++resPtr;
    if (resPtr == AllResidues.end()) {
      mprinterr("Error: Chain %c not found\n", gapRes0.Chain());
      return 1;
    }
    mprintf("\t  Chain %c found: %s %i\n", gapRes0.Chain(),
            resPtr->Name().c_str(), resPtr->Onum());
    // Search until we find an appropriate resnum to insert after
    int resDelta0 = gapRes0.Onum() - resPtr->Onum();
    int resDelta1 = gapRes1.Onum() - resPtr->Onum();
    while (abs(resDelta0) > 1 && abs(resDelta1) > 1)
    {
      ++resPtr;
      if (resPtr == AllResidues.end()) {
        mprinterr("Error: Could not find appropriate # to insert %i in Chain %c.\n",
                  gapRes0.Onum(), gapRes0.Chain());
        return 1;
      }
      resDelta0 = gapRes0.Onum() - resPtr->Onum();
      resDelta1 = gapRes1.Onum() - resPtr->Onum();
    }
    //while (resPtr != AllResidues.end()) {
    //        mprintf("\t    Res %4s %6i delta0= %6i delta1= %6i\n",
    //          resPtr->Name().c_str(), resPtr->Onum(), resDelta0, resDelta1);
    //  resPtr++;
    //}
    // resPtr = AllResidues.begin(); // DEBUG
    mprintf("\t  Res found: %s %i (delta0= %i delta1= %i)\n",
            resPtr->Name().c_str(), resPtr->Onum(),
            resDelta0, resDelta1);
    // Determine if this gap should come before or after resPtr
    // TODO handle insertion codes
    if (resDelta0 == 1) {
      mprintf("\t    Gap comes AFTER residue.\n");
      resPtr++;
      AllResidues.insert(resPtr, gap->begin(), gap->end());
    } else if (resDelta1 == -1) {
      mprintf("\t    Gap comes BEFORE residue.\n");
      AllResidues.insert(resPtr, gap->begin(), gap->end());
    } else {
      mprintf("\t    UNHANDLED.\n");
    }
  }
  // Print all residues
  // Count the number of present atoms and missing residues.
  int nAtomsPresent = 0;
  int nResMissing = 0;
  for (ResList::const_iterator it = AllResidues.begin(); it != AllResidues.end(); ++it)
  {
    mprintf("\t%4s %6i %6i %c %c\n", it->Name().c_str(), it->Onum(), it->Tnum(),
            it->Icode(), it->Chain());
    if (it->Tnum() < 0)
      nResMissing++;
    else
      nAtomsPresent += topIn.Res(it->Tnum()).NumAtoms();
  }
  mprintf("\t%i atoms present, %i residues missing.\n", nAtomsPresent, nResMissing);
  // Create new Frame
  Frame newFrame(nAtomsPresent + nResMissing);
  newFrame.ClearAtoms();
  // Create a new topology with all residues. For missing residues, create a CA atom.
  Topology newTop;
  // Zero coord for new CA atoms
  Vec3 Zero(0.0);
  // Topology for CA atoms
  Topology CAtop;
  // Frame for CA atoms
  Frame CAframe;
  // Mask for missing CA atoms
  CharMask CAmissing;
  // Loop over all residues
  for (ResList::const_iterator it = AllResidues.begin(); it != AllResidues.end(); ++it)
  {
    if (it->Tnum() < 0) {
      // This was a missing residue
      newTop.AddTopAtom( Atom("CA", "C "),
                         Residue(it->Name(), it->Onum(), it->Icode(), it->Chain()) );
      newFrame.AddVec3( Zero );
      // CA top
      CAtop.AddTopAtom( Atom("CA", "CA", 0),
                        Residue(it->Name(), it->Onum(), it->Icode(), it->Chain()) );
      CAframe.AddVec3( Zero );
      CAmissing.AddAtom(true);
    } else {
      // This residue is present in the original PDB
      Residue const& topres = topIn.Res(it->Tnum());
      Residue newres(topres.Name(), topres.OriginalResNum(), topres.Icode(), topres.ChainID());
      int caidx = -1;
      // Calculate the center of the residue as we go in case we need it
      Vec3 vcenter(0.0);
      for (int at = topres.FirstAtom(); at < topres.LastAtom(); at++) {
        if (topIn[at].Name() == "CA") caidx = at;
        newTop.AddTopAtom( Atom(topIn[at].Name(), topIn[at].ElementName()), newres );
        //frameIn.printAtomCoord(at);
        const double* txyz = frameIn.XYZ(at);
        newFrame.AddXYZ( txyz );
        vcenter[0] += txyz[0];
        vcenter[1] += txyz[1];
        vcenter[2] += txyz[2];
        //newFrame.printAtomCoord(newTop->Natom()-1);
      }
      // CA top
      if (caidx == -1) {
        mprintf("Warning: No CA atom found for residue %s\n",
                topIn.TruncResNameNum(it->Tnum()).c_str());
        // Use the center of the residue
        vcenter /= (double)topres.NumAtoms();
        mprintf("Warning: Using center: %g %g %g\n", vcenter[0], vcenter[1], vcenter[2]);
        CAtop.AddTopAtom( Atom("CA", "C"), newres );
        CAframe.AddVec3( vcenter );
        CAmissing.AddAtom(false);
      } else {
        CAtop.AddTopAtom( Atom(topIn[caidx].Name(), topIn[caidx].ElementName()), newres );
        CAframe.AddXYZ( frameIn.XYZ(caidx) );
        CAmissing.AddAtom(false);
      }
    }
  } // END loop over all residues

  // Try to determine which frames are terminal so pdbter works since
  // we wont be able to determine molecules by bonds.
  for (int ridx = 0; ridx != newTop.Nres(); ridx++)
  {
    if (ridx + 1 == newTop.Nres() ||
        newTop.Res(ridx).ChainID() != newTop.Res(ridx+1).ChainID())
      newTop.SetRes(ridx).SetTerminal(true);
  }
  // Finish new top and write
  newTop.SetParmName("newpdb", "temp.pdb");
  newTop.CommonSetup( false ); // No molecule search
  newTop.Summary();
  if (WriteStructure("temp.pdb", &newTop, newFrame, TrajectoryFile::PDBFILE)) {
    mprinterr("Error: Write of temp.pdb failed.\n");
    return 1;
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
