#include "Exec_AddMissingRes.h"
#include "Exec_AddMissingRes_Pres.h"
#include "BufferedLine.h"
#include "CharMask.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords_CRD.h"
#include "Minimize_SteepestDescent.h"
#include "ParmFile.h"
#include "PotentialFunction.h"
#include "Trajin_Single.h"
#include "Trajout_Single.h"
#include <cmath>
#include <cstring>
#include <algorithm>

/** CONSTRUCTOR. */
Exec_AddMissingRes::Exec_AddMissingRes() :
  Exec(GENERAL),
  debug_(0),
  nMinSteps_(0),
  optimize_(true)
{
  SetHidden(true);
}

/** Search for closest residue to the given residue in terms
  * of sequence.
*/
Exec_AddMissingRes::Rlist::iterator Exec_AddMissingRes::FindClosestRes(Rlist& ResList, Residue const& tgtRes,
                                                                       int& lowestAbsDist)
const
{
  mprintf("DEBUG: Searching for closest residue to %s %i %c %c\n",
          *(tgtRes.Name()), tgtRes.OriginalResNum(), tgtRes.Icode(), tgtRes.ChainId());
  Rlist::iterator pos = ResList.begin();
  Rlist::iterator closest = pos;
  lowestAbsDist = tgtRes.AbsResDist( *pos );
  ++pos;
  while (pos != ResList.end()) {
    int absDist = tgtRes.AbsResDist( *pos );
    // As soon as we are farther away, stop looking
    if (absDist > lowestAbsDist) break;
    if (absDist < lowestAbsDist) {
      closest = pos;
      lowestAbsDist = absDist;
    }
    ++pos;
  }
  return closest;
}

/** Get a complete sequence of residues from the PDB, including
  * missing residues.
  */
int Exec_AddMissingRes::GetSequenceFromPDB(Rlist &ResList, std::string const& pdbname)
const
{
  Rlist missingResidues;
  BufferedLine infile;
  if (infile.OpenFileRead( pdbname )) {
    mprinterr("Error: Could not open '%s' for reading.\n", pdbname.c_str());
    return 1;
  }
  // Loop over lines from PDB
  int inMissing = 0;
  int nmissing = 0;
  const char* linePtr = infile.Line();
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
            //mprintf("DEBUG: Missing residue %s %i icode= %c chain= %c\n",
            //        currentname.c_str(), currentres, currenticode, currentchain);
            missingResidues.push_back( Residue(currentname, currentres, currenticode, currentchain) );
          } // END missing residue
        } // END inMissing == 2
      } // END nargs > 2
    } // END REMARK
    linePtr = infile.Line();
  } // END while linePtr != 0
  infile.CloseFile();

  // DEBUG
  for (Rlist::const_iterator it = missingResidues.begin(); it != missingResidues.end(); ++it)
    mprintf("DEBUG: Missing residue %s %i icode= %c chain= %c\n",
            *(it->Name()), it->OriginalResNum(), it->Icode(), it->ChainId());

  // Read existing residues from the PDB
  ParmFile pdbIn;
  Topology topIn;
  if (pdbIn.ReadTopology(topIn, pdbname, ArgList(), debug_)) {
    mprinterr("Error: Read of topology from PDB failed.\n");
    return CpptrajState::ERR;
  }
  topIn.Summary();

  // Put existing residues into a list
  ResList.clear();
  for (Topology::res_iterator res = topIn.ResStart(); res != topIn.ResEnd(); ++res)
    ResList.push_back( Residue(res->Name(), res->OriginalResNum(), res->Icode(), res->ChainId()) );

  // Figure out where to insert the missing residues.
  while (!missingResidues.empty()) {
    Rlist::iterator gapStart = missingResidues.begin();
    // Search for gap end
    Rlist::iterator gapEnd = gapStart;
    bool searchGap = true;
    while ( searchGap ) {
      Rlist::iterator currentRes = gapEnd;
      ++gapEnd;
      if (gapEnd == missingResidues.end()) break;
      int resDelta = gapEnd->OriginalResNum() - currentRes->OriginalResNum();
      if ( gapEnd->ChainId() != currentRes->ChainId() )
        searchGap = false;
      else if ( resDelta == 0 ) {
        int icodeDelta = (int)gapEnd->Icode() - (int)currentRes->Icode();
        if (icodeDelta == 0) {
          mprinterr("Error: Missing residues have same res # %i and insertion codes %c %c\n",
                    currentRes->OriginalResNum(), currentRes->Icode(), gapEnd->Icode());
          return 1;
        }
        if (icodeDelta < 0) icodeDelta = -icodeDelta;
        if (icodeDelta > 1)
        searchGap = false;
      } else if (resDelta != 1) {
        searchGap = false;
      }
    }
    Rlist::iterator gapFinalRes = gapEnd;
    --gapFinalRes;
    mprintf("Gap:\n");
    int gapLength = 0;
    for (Rlist::iterator it = gapStart; it != gapEnd; ++it) {
      mprintf("\tRes %s %i icode= %c chain= %c\n",
              *(it->Name()), it->OriginalResNum(), it->Icode(), it->ChainId());
      gapLength++;
    }
    mprintf("\t  Gap length= %i\n", gapLength);
    // Find closest res to start
    int startDist = -1, endDist = -1;
    Rlist::iterator closestStart = FindClosestRes(ResList, *gapStart, startDist);
    mprintf("\t  Closest res to gap start is %s %i %c %c (%i)\n",
            *(closestStart->Name()), closestStart->OriginalResNum(), closestStart->Icode(), closestStart->ChainId(), startDist);
    Rlist::iterator closestEnd = ResList.end();
    if (gapStart != gapFinalRes) {
      closestEnd = FindClosestRes(ResList, *gapFinalRes, endDist);
      mprintf("\t  Closest res to gap end is %s %i %c %c (%i)\n",
              *(closestEnd->Name()), closestEnd->OriginalResNum(), closestEnd->Icode(), closestEnd->ChainId(), endDist);
    }
    if (startDist > 1 && endDist > 1) {
      mprinterr("Error: Cannot find place to attach gap start or end.\n");
      return 1;
    } else if (startDist == 1 && endDist == 1) {
      // Somewhere in the middle of the chain.
      mprintf("\t  Attach in the middle of a chain.\n");
    } else if (startDist == 1) {
      if (endDist == -1) {
        if (gapLength != 1) {
          mprinterr("Internal Error: No gap end found but gap length is greater than 1.\n");
          return 1;
        }
        mprintf("\t  Attach single residue in middle of a chain.\n");
      } else
        mprintf("\t  Attach C-terminal.\n");
    } else if (endDist == 1) {
      mprintf("\t  Attach N-terminal.\n");
    } else {
      // SANITY CHECK
      mprinterr("Internal Error: Unexpected gap distances (start=%i, end=%i)\n", startDist, endDist);
      return 1;
    }
    
    missingResidues.erase(gapStart, gapEnd);
    gapStart = gapEnd;
  }
 
  


  return 0;
}

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

// -----------------------------------------------------------------------------
/** Try to minimize using steepest descent. */
int Exec_AddMissingRes::Minimize(Topology const& topIn, Frame& frameIn, CharMask const& maskIn)
const
{
  double min_tol = 1.0E-5;

  PotentialFunction potential;
  potential.AddTerm( PotentialTerm::BOND );
  potential.AddTerm( PotentialTerm::SIMPLE_LJ_Q );
  Minimize_SteepestDescent SD;

  if (potential.SetupPotential( topIn, maskIn )) return 1;

  if (SD.SetupMin("min.nc", min_tol, 1.0, nMinSteps_)) return 1;

  // Add force info to frame.
  frameIn.AddForces(Frame::Darray(frameIn.size(), 0.0));

  CpptrajFile outfile;
  outfile.OpenWrite("");
  if (SD.RunMin(potential, frameIn, outfile)) return 1;

  return 0;
}

// -----------------------------------------------------------------------------
/** Calculate pseudo force vector at given index. Only take into account
  * atoms within a certain cutoff, i.e. the local environment.
  */
int Exec_AddMissingRes::CalcFvecAtIdx(Vec3& vecOut, Vec3& XYZ0, int tgtidx, 
                                      Topology const& CAtop, Frame const& CAframe,
                                      CharMask const& isMissing)
{
  double cut2 = 100.0; // 10 ang
  vecOut = Vec3(0.0);
  XYZ0 = Vec3( CAframe.XYZ(tgtidx) );
  double e_total = 0.0;
  for (int idx = 0; idx != CAtop.Nres(); idx++)
  {
    if (idx != tgtidx && !isMissing.AtomInCharMask(idx))
    {
      const double* XYZ1 = CAframe.XYZ( idx );
      double rx = XYZ0[0] - XYZ1[0];
      double ry = XYZ0[1] - XYZ1[1];
      double rz = XYZ0[2] - XYZ1[2];
      double rij2 = rx*rx + ry*ry + rz*rz;
      if (rij2 > 0 && rij2 < cut2) {
        double rij = sqrt( rij2 );
        // COULOMB
        double qiqj = .01; // Give each atom charge of .1
        double e_elec = 1.0 * (qiqj / rij); // 1.0 is electrostatic constant, not really needed
        e_total += e_elec;
        // COULOMB force
        double felec = e_elec / rij; // kes * (qiqj / r) * (1/r)
        vecOut[0] += rx * felec;
        vecOut[1] += ry * felec;
        vecOut[2] += rz * felec;
      } //else {
        //mprinterr("Error: Atom clash between CA %i and %i\n", tgtidx+1, idx+1);
        //return 1;
      //}
    }
  }
  vecOut.Normalize();
  return 0;
}

/** \return A vector containing residues from start up to and including end. */
static inline std::vector<int> ResiduesToSearch(int startRes, int endRes) {
  std::vector<int> residuesToSearch;
  if (startRes < endRes) {
    for (int i = startRes; i <= endRes; i++)
      residuesToSearch.push_back( i );
  } else {
    for (int i = startRes; i >= endRes; i--)
      residuesToSearch.push_back( i );
  }
  return residuesToSearch;
}

/** Calculate a guiding force connecting XYZ0 and XYZ1 */
static inline void CalcGuideForce(Vec3 const& XYZ0, Vec3 const& XYZ1, double maxDist,
                                  double Rk, Vec3& fvec0, Vec3& fvec1)
{
  // Vector from 0 to 1
  Vec3 v01 = XYZ1 - XYZ0;
  // Distance
  double r2 = v01.Magnitude2();
  double r01 = sqrt( r2 );
  // Only apply the guiding force above maxDist 
  if (r01 > maxDist) {
    // Normalize
    v01.Normalize();
    // Augment
    v01 *= Rk;
    v01.Print("guide");
    // Add
    fvec0 += v01;
    fvec1 -= v01;
  }
}

/** Try to generate linear coords beteween idx0 and idx1. */
void Exec_AddMissingRes::GenerateLinearGapCoords(int idx0, int idx1, Frame& frm)
{
  Vec3 vec0( frm.XYZ(idx0) );
  Vec3 vec1( frm.XYZ(idx1) );
  vec0.Print("vec0");
  vec1.Print("vec1");
  Vec3 V10 = vec1 - vec0;
  int nsteps = idx1 - idx0;
  if (nsteps < 1) {
    mprinterr("Internal Error: GenerateLinearGapCoords: Invalid steps from %i to %i (%i)\n",
              idx0, idx1, nsteps);
    return;
  }
  mprintf("DEBUG: Generating %i steps from %i to %i\n", nsteps, idx0+1, idx1+1);
  Vec3 delta = V10 / (double)nsteps;
  double* Xptr = frm.xAddress() + ((idx0+1)*3);
  for (int i = 1; i < nsteps; i++, Xptr += 3)
  {
    Vec3 xyz = vec0 + (delta * (double)i);
    xyz.Print("xyz");
    Xptr[0] = xyz[0];
    Xptr[1] = xyz[1];
    Xptr[2] = xyz[2];
  }
}
 
/** Search for coords from anchor0 to anchor1, start at start, end at end. */
int Exec_AddMissingRes::CoordSearchGap(int anchor0, int anchor1, int startRes, int endRes, 
                                     Topology const& CAtop, CharMask& isMissing, Frame& CAframe)
const
{
  Iarray residues = ResiduesToSearch(startRes, endRes);
  if (residues.size() < 2) {
    GenerateLinearGapCoords(anchor0, anchor1, CAframe);
    return 0;
  }
  // First calculate the force vector at anchor0
  mprintf("Anchor Residue 0: %i\n", anchor0+1);
  Vec3 XYZ0, Vec0;
  CalcFvecAtIdx(Vec0, XYZ0, anchor0, CAtop, CAframe, isMissing);
  XYZ0.Print("Anchor 0 coords");
  Vec0.Print("anchor 0 vec");
  // Calculate the force vector at anchor1
  mprintf("Anchor Residue 1: %i\n", anchor1+1);
  Vec3 XYZ1, Vec1;
  CalcFvecAtIdx(Vec1, XYZ1, anchor1, CAtop, CAframe, isMissing);
  XYZ1.Print("Anchor 1 coords");
  Vec1.Print("anchor 1 vec");
  // Guide vector
  double guidefac = 3.8;
  double guidek = 1.0;
  CalcGuideForce(XYZ0, XYZ1, guidefac, guidek, Vec0, Vec1);

  // Determine the halfway index
  int halfIdx = residues.size() / 2; 
  // N-terminal
  Iarray fromAnchor0 = ResiduesToSearch(residues.front(), residues[halfIdx-1]);
  // C-terminal
  Iarray fromAnchor1 = ResiduesToSearch(residues.back(), residues[halfIdx]);

  mprintf("DEBUG: Generating Gap residues:\n");
  mprintf("\tFrom %i:", anchor0+1);
  for (Iarray::const_iterator it = fromAnchor0.begin(); it != fromAnchor0.end(); ++it)
    mprintf(" %i", *it + 1);
  mprintf("\n");
  mprintf("\tFrom %i:", anchor1+1);
  for (Iarray::const_iterator it = fromAnchor1.begin(); it != fromAnchor1.end(); ++it)
    mprintf(" %i", *it + 1);
  mprintf("\n");

  // Loop over fragment to generate coords for
  double fac = 2.0;
  Iarray::const_iterator a0 = fromAnchor0.begin();
  Iarray::const_iterator a1 = fromAnchor1.begin();
  while (a0 != fromAnchor0.end() || a1 != fromAnchor1.end()) {
    if (a0 != fromAnchor0.end()) {
      double* Xptr = CAframe.xAddress() + (*a0 * 3);
      Vec3 xyz = XYZ0 + (Vec0 * fac);
      mprintf("  %i %12.4f %12.4f %12.4f\n", *a0+1, xyz[0], xyz[1], xyz[2]);
      Xptr[0] = xyz[0];
      Xptr[1] = xyz[1];
      Xptr[2] = xyz[2];
      // Mark as not missing
      isMissing.SelectAtom(*a0, false);
      // Update the anchor
      CalcFvecAtIdx(Vec0, XYZ0, *a0, CAtop, CAframe, isMissing);
      ++a0;
    }
    if (a1 != fromAnchor1.end()) {
      double* Xptr = CAframe.xAddress() + (*a1 * 3);
      Vec3 xyz = XYZ1 + (Vec1 * fac);
      mprintf("  %i %12.4f %12.4f %12.4f\n", *a1+1, xyz[0], xyz[1], xyz[2]);
      Xptr[0] = xyz[0];
      Xptr[1] = xyz[1];
      Xptr[2] = xyz[2];
      // Mark as not missing
      isMissing.SelectAtom(*a1, false);
      // Update the anchor
      CalcFvecAtIdx(Vec1, XYZ1, *a1, CAtop, CAframe, isMissing);
      ++a1;
    }
    CalcGuideForce(XYZ0, XYZ1, guidefac, guidek, Vec0, Vec1);
  }

  return 0;
}

/** Search for coords using anchorRes as an anchor, start at start, end at end. */
int Exec_AddMissingRes::CoordSearchTerminal(int anchorRes, int startRes, int endRes,
                                     Topology const& CAtop, CharMask& isMissing, Frame& CAframe)
const
{
  // First calculate the force vector at the anchorRes
  mprintf("Anchor Residue %i\n", anchorRes+1);
  Vec3 XYZ0, anchorVec;
  CalcFvecAtIdx(anchorVec, XYZ0, anchorRes, CAtop, CAframe, isMissing);
  XYZ0.Print("Anchor coords");
  anchorVec.Print("DEBUG: anchorVec");
  // Determine the direction
  Iarray residuesToSearch = ResiduesToSearch(startRes, endRes);
  // Loop over fragment to generate coords for
  double fac = 2.0;
  mprintf("DEBUG: Generating linear fragment extending from %i for indices %i to %i (%zu)\n",
          anchorRes+1, startRes+1, endRes+1, residuesToSearch.size());
  for (Iarray::const_iterator it = residuesToSearch.begin();
                              it != residuesToSearch.end(); ++it)
  {
    double* Xptr = CAframe.xAddress() + (*it * 3);
    Vec3 xyz = XYZ0 + (anchorVec * fac);
    mprintf("  %i %12.4f %12.4f %12.4f\n", *it+1, xyz[0], xyz[1], xyz[2]);
    Xptr[0] = xyz[0];
    Xptr[1] = xyz[1];
    Xptr[2] = xyz[2];
    //CAframe.printAtomCoord( *it );
    // Mark as not missing
    isMissing.SelectAtom(*it, false);
    // Update the anchor
    CalcFvecAtIdx(anchorVec, XYZ0, *it, CAtop, CAframe, isMissing);
    
  }

  return 0;
}

/** Generate coords using an energy search.
  * \param newTop The new topology containing all existing atoms and CA atoms for missing residues.
  * \param newFrame The new frame containins coordinates for all existing atoms and 0 for missing residues.
  * \param CAtop Only CA atoms from the new topology.
  * \param CAframe Coords for CAtop
  * \param CAmissing True if CA was missing, false otherwise.
  */ 
int Exec_AddMissingRes::AssignCoordsBySearch(Topology const& newTop, Frame const& newFrame,
                                             Topology const& CAtop, Frame& CAframe,
                                             CharMask const& CAmissing)
const
{
  CharMask isMissing = CAmissing;
  int ridx = 0;
  int gapStart = -1;
  char chain = ' ';
  while (ridx < newTop.Nres())
  {
    // Determine the next gap
    if (gapStart == -1) {
      // Not in a gap. Check for gap start.
      if (CAmissing.AtomInCharMask(ridx)) {
        //mprintf("Start gap %i\n", ridx);
        gapStart = ridx;
        chain = CAtop.Res(ridx).ChainId();
        // Find gap end
        int gapEnd = gapStart + 1;
        for (; gapEnd <= newTop.Nres(); gapEnd++) {
          if (gapEnd == newTop.Nres() ||
              !CAmissing.AtomInCharMask(gapEnd) ||
              CAtop.Res(gapEnd).ChainId() != chain)
          {
            // The real end of the gap is the previous res
            gapEnd = gapEnd - 1;
            mprintf("\tNew Gap: %c %i to %i\n", chain, gapStart+1, gapEnd+1);
            // -----------------------------------
            // Is there a previous residue
            int prev_res = gapStart - 1;
            if (prev_res < 0 ||
                newTop.Res(prev_res).ChainId() != newTop.Res(gapStart).ChainId())
              prev_res = -1;
            // Is there a next residue
            int next_res = gapEnd + 1;
            if (next_res == newTop.Nres() ||
                newTop.Res(next_res).ChainId() != newTop.Res(gapEnd).ChainId())
              next_res = -1;
            mprintf("\t  Prev res %i  Next res %i\n", prev_res + 1, next_res + 1);
            if (prev_res == -1 && next_res == -1) {
              mprinterr("Error: Gap is unconnected.\n");
              return 1;
            }
            if (prev_res > -1 && next_res > -1) {
              CoordSearchGap(prev_res, next_res, gapStart, gapEnd, CAtop, isMissing, CAframe); 
            } else if (prev_res == -1) {
              // N-terminal
              CoordSearchTerminal(next_res, gapEnd, gapStart, CAtop, isMissing, CAframe);
            } else if (next_res == -1) {
              // C-terminal
              CoordSearchTerminal(prev_res, gapStart, gapEnd, CAtop, isMissing, CAframe);
            }
            // -----------------------------------
            gapStart = -1;
            ridx = gapEnd + 1;
            break;
          }
        }
      } else
        ridx++;
    }
  }
  return 0;
}

// -----------------------------------------------------------------------------
/** Write topology and frame to specified format. */
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
                                Res.Icode(), Res.ChainId()) );
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
  // Map original residue index to new residue index
  Iarray originalIdxToNew;
  // Map original atom number to new atom number
  Iarray originalAtToNew;
  // Loop over all residues
  int newResNum = 0;
  int newAtNum = 0;
  for (ResList::const_iterator it = AllResidues.begin(); it != AllResidues.end(); ++it, ++newResNum)
  {
    if (it->Tnum() < 0) {
      // This was a missing residue
      newTop.AddTopAtom( Atom("CA", "C "),
                         Residue(it->Name(), it->Onum(), it->Icode(), it->Chain()) );
      newAtNum++;
      newFrame.AddVec3( Zero );
      // CA top
      CAtop.AddTopAtom( Atom("CA", "CA", 0),
                        Residue(it->Name(), it->Onum(), it->Icode(), it->Chain()) );
      CAframe.AddVec3( Zero );
      CAmissing.AddAtom(true);
    } else {
      originalIdxToNew.push_back( newResNum );
      // This residue is present in the original PDB
      Residue const& topres = topIn.Res(it->Tnum());
      Residue newres(topres.Name(), topres.OriginalResNum(), topres.Icode(), topres.ChainId());
      int caidx = -1;
      // Calculate the center of the residue as we go in case we need it
      Vec3 vcenter(0.0);
      for (int at = topres.FirstAtom(); at < topres.LastAtom(); at++) {
        if (topIn[at].Name() == "CA") caidx = at;
        newTop.AddTopAtom( Atom(topIn[at].Name(), topIn[at].ElementName()), newres );
        originalAtToNew.push_back( newAtNum++ );
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
        newTop.Res(ridx).ChainId() != newTop.Res(ridx+1).ChainId())
      newTop.SetRes(ridx).SetTerminal(true);
  }
  // Add original bonds to new topology.
  for (BondArray::const_iterator bnd = topIn.Bonds().begin();
                                 bnd != topIn.Bonds().end(); ++bnd)
    newTop.AddBond( originalAtToNew[bnd->A1()],
                    originalAtToNew[bnd->A2()] );
  // Finish new top and write
  newTop.SetParmName("newpdb", "temp.pdb");
  newTop.CommonSetup( false ); // No molecule search
  newTop.Summary();
  if (WriteStructure("temp.pdb", &newTop, newFrame, TrajectoryFile::PDBFILE)) {
    mprinterr("Error: Write of temp.pdb failed.\n");
    return 1;
  }

  // Determine pseudo-bonds for CA topology
  BondParmType CAbond(300.0, 3.8);
  int ridx = 0; // Residue index for new/CA topology
  for (ResList::const_iterator it = AllResidues.begin(); it != AllResidues.end(); ++it, ++ridx)
  {
    if (it->Tnum() < 0) {
      // This was a missing residue. Bond it to the previous/next residue in same chain.
      int pidx = ridx - 1;
      if (pidx > 0 && CAtop.Res(pidx).ChainId() == CAtop.Res(ridx).ChainId())
        CAtop.AddBond(pidx, ridx, CAbond);
      int nidx = ridx + 1;
      if (nidx < CAtop.Nres() && CAtop.Res(ridx).ChainId() == CAtop.Res(nidx).ChainId())
        CAtop.AddBond(ridx, nidx, CAbond);
    } else {
      // Check which residues this was originally bonded to
      //mprintf("New index %i original index %i originally bonded to", ridx, it->Tnum());
      Residue const& ores = topIn.Res(it->Tnum());
      //std::set<int> bondedRes;
      for (int at = ores.FirstAtom(); at != ores.LastAtom(); ++at)
      {
        for (Atom::bond_iterator bat = topIn[at].bondbegin();
                                 bat != topIn[at].bondend(); ++bat)
        {
          if (topIn[*bat].ResNum() != it->Tnum()) {
            //mprintf(" %i", topIn[*bat].ResNum());
            CAtop.AddBond(ridx, originalIdxToNew[topIn[*bat].ResNum()]);
          }
        }
      }
      //mprintf("\n");
    }
  }
  // Add pseudo parameters for the "CA" atom type (0)
  LJparmType CAtype(3.8, 10.0);
  NonbondType AB = CAtype.Combine_LB( CAtype );
  CAtop.SetNonbond().SetupLJforNtypes(1);
  CAtop.SetNonbond().AddLJterm(0, 0, AB);
  mprintf("DEBUG: LJ radius= %g\n", CAtop.GetVDWradius(0));
  // Final setup
  CAtop.SetParmName("capdb", "temp.ca.mol2");
  CAtop.CommonSetup(true, 2); // molecule search, exclude bonds
  CAtop.Summary();

  // Write CA top
  if (WriteStructure("temp.ca.mol2", &CAtop, CAframe, TrajectoryFile::MOL2FILE)) {
    mprinterr("Error: Write of temp.ca.mol2 failed.\n");
    return 1;
  }

  if (optimize_) {
    // Try to assign new coords for the missing residues.
    if (AssignCoordsBySearch(newTop, newFrame, CAtop, CAframe, CAmissing)) {
      mprinterr("Error: Could not assign coords by search.\n");
      return 1;
    }
    // Minimize
    if (Minimize(CAtop, CAframe, CAmissing)) {
      mprinterr("Error: Minimization of CA atoms failed.\n");
      return 1;
    }
    // Transfer final CA coords for missing residues to newFrame
    for (int idx = 0; idx != CAtop.Nres(); idx++) {
      if (CAmissing.AtomInCharMask(idx)) {
        const double* XYZ = CAframe.XYZ(idx);
        Residue const& newres = newTop.Res(idx);
        int newat = newres.FirstAtom();
        double* Xptr = newFrame.xAddress() + (3*newat);
        std::copy(XYZ, XYZ+3, Xptr);
      }
    }
  }

  // Write final CA top
  if (WriteStructure("temp.ca.final.mol2", &CAtop, CAframe, TrajectoryFile::MOL2FILE)) {
    mprinterr("Error: Write of temp.ca.final.mol2 failed.\n");
    return 1;
  }

  // Set output coords
  dataOut->CoordsSetup( newTop, CoordinateInfo() );
  dataOut->AddFrame( newFrame );
  
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
  Rlist FullResSequence;
  std::string pdbseq = argIn.GetStringKey("pdbseq");
  if (!pdbseq.empty()) {
    if (GetSequenceFromPDB(FullResSequence, pdbseq))
      return CpptrajState::ERR;
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
