#include "Exec_Graft.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords.h"
#include <algorithm> // std::copy

// Exec_Graft::Help()
void Exec_Graft::Help() const
{
  mprintf("\tsrc <source COORDS> [srcframe <#>] [srcfitmask <mask>] [srcmask <mask>]\n"
          "\ttgt <target COORDS> [tgtframe <#>] [tgtfitmask <mask>] [tgtmask <mask>]\n"
          "\tname <output COORDS> [bond <tgt>,<src> ...]\n"
          "  Graft coordinates from source to coordinates in target.\n");
}

/** Update indices in the given array to what the new indices will
  * be after processing with the given mask.
  * \return 1 if an index is not in the mask, 0 if all indices updated successfully.
  */
static int UpdateIndices(std::vector<int>& Idxs, AtomMask const& maskIn, int offset)
{
  for (std::vector<int>::iterator it = Idxs.begin();
                                  it != Idxs.end(); ++it)
  {
    // Search for index in mask
    int newidx = 0;
    for (; newidx != maskIn.Nselected(); newidx++)
    {
      if (*it == maskIn[newidx]) {
        *it = newidx + offset;
        break;
      }
    }
    if (newidx == maskIn.Nselected()) {
      mprinterr("Error: Bonded index is in a removed section.\n");
      return 1;
    }
  }
  return 0;
}

// Exec_Graft::Execute()
Exec::RetType Exec_Graft::Execute(CpptrajState& State, ArgList& argIn)
{
  typedef std::vector<int> Iarray;
  // Get source coords
  std::string kw = argIn.GetStringKey("src");
  if (kw.empty()) {
    mprinterr("Error: Source COORDS must be specified with 'src'.\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords* srcCoords = (DataSet_Coords*)State.DSL().FindSetOfGroup(kw, DataSet::COORDINATES);
  if (srcCoords == 0) {
    mprinterr("Error: Source COORDS %s not found.\n", kw.c_str());
    return CpptrajState::ERR;
  }
  Frame srcFrame = srcCoords->AllocateFrame();
  srcCoords->GetFrame(argIn.getKeyInt("srcframe", 1)-1, srcFrame);
  // Get target coords
  kw = argIn.GetStringKey("tgt");
  if (kw.empty()) {
    mprinterr("Error: Target COORDS must be specified with 'tgt'.\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords* tgtCoords = (DataSet_Coords*)State.DSL().FindSetOfGroup(kw, DataSet::COORDINATES);
  if (tgtCoords == 0) {
    mprinterr("Error: Target COORDS %s not found.\n", kw.c_str());
    return CpptrajState::ERR;
  }
  Frame tgtFrame = tgtCoords->AllocateFrame();
  tgtCoords->GetFrame(argIn.getKeyInt("tgtframe", 1)-1, tgtFrame);
  // Create output coords
  kw = argIn.GetStringKey("name");
  if (kw.empty()) {
    mprinterr("Error: Output COORDS must be specified with 'name'.\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords* outCoords = (DataSet_Coords*)State.DSL().AddSet(DataSet::COORDS, MetaData(kw));
  if (outCoords == 0) {
    mprinterr("Error: Output COORDS %s could not be created.\n", kw.c_str());
    return CpptrajState::ERR;
  }
  // Get atoms to keep from source.
  AtomMask srcMask;
  if (srcMask.SetMaskString( argIn.GetStringKey("srcmask") ))
    return CpptrajState::ERR;
  if (srcCoords->Top().SetupIntegerMask(srcMask))
    return CpptrajState::ERR;
  // Get atoms to keep from target.
  AtomMask tgtMask;
  if (tgtMask.SetMaskString( argIn.GetStringKey("tgtmask") ))
    return CpptrajState::ERR;
  if (tgtCoords->Top().SetupIntegerMask(tgtMask))
    return CpptrajState::ERR;
  // Get bonds
  Iarray tgtBondAtoms, srcBondAtoms;
  kw = argIn.GetStringKey("bond");
  while (!kw.empty()) {
    ArgList bndarg(kw, ",");
    if (bndarg.Nargs() != 2) {
      mprinterr("Error: Expected 2 atom masks for 'bond' (target, source).\n");
      return CpptrajState::ERR;
    }
    AtomMask tb, sb;
    if (tb.SetMaskString(bndarg[0])) return CpptrajState::ERR;
    if (sb.SetMaskString(bndarg[1])) return CpptrajState::ERR;
    if (tgtCoords->Top().SetupIntegerMask(tb)) return CpptrajState::ERR;
    if (srcCoords->Top().SetupIntegerMask(sb)) return CpptrajState::ERR;
    if (tb.Nselected() != 1) {
      mprinterr("Error: 'bond' target mask does not select only 1 atom.\n");
      return CpptrajState::ERR;
    }
    if (sb.Nselected() != 1) {
      mprinterr("Error: 'bond' source mask does not select only 1 atom.\n");
      return CpptrajState::ERR;
    }
    tgtBondAtoms.push_back( tb[0] );
    srcBondAtoms.push_back( sb[0] );
    mprintf("\tWill bond target %s (%i) to source %s (%i)\n",
            tb.MaskString(), tgtBondAtoms.back()+1,
            sb.MaskString(), srcBondAtoms.back()+1);
    kw = argIn.GetStringKey("bond");
  }
  // Update the bond indices for the new topologies
  if (UpdateIndices(tgtBondAtoms, tgtMask, 0))
    return CpptrajState::ERR;
  if (UpdateIndices(srcBondAtoms, srcMask, tgtMask.Nselected()))
    return CpptrajState::ERR;
  mprintf("\tUpdated bond indices:\n");
  for (unsigned int ii = 0; ii != tgtBondAtoms.size(); ii++)
    mprintf("\t  tgt= %i  src= %i\n", tgtBondAtoms[ii]+1, srcBondAtoms[ii]+1);
  // Get atoms from source to fit on target, and atoms from target
  // for source to fit on.
  AtomMask srcFitMask, tgtFitMask;
  std::string srcfitstr = argIn.GetStringKey("srcfitmask");
  std::string tgtfitstr = argIn.GetStringKey("tgtfitmask");
  bool doRmsFit = false;
  if (!srcfitstr.empty() || !tgtfitstr.empty()) {
    doRmsFit = true;
    // If either is empty, fill with the other.
    if (srcfitstr.empty())
      srcfitstr = tgtfitstr;
    else if (tgtfitstr.empty())
      tgtfitstr = srcfitstr;
    if (srcFitMask.SetMaskString( srcfitstr ))
      return CpptrajState::ERR;
    if (tgtFitMask.SetMaskString( tgtfitstr ))
      return CpptrajState::ERR;
    // Set up the masks.
    if (srcCoords->Top().SetupIntegerMask(srcFitMask))
      return CpptrajState::ERR;
    if (tgtCoords->Top().SetupIntegerMask(tgtFitMask))
      return CpptrajState::ERR;
  }
  
  // Info
  mprintf("\tSource coords   : %s\n", srcCoords->legend());
  mprintf("\tTarget coords   : %s\n", tgtCoords->legend());
  mprintf("\tOutput coords   : %s\n", outCoords->legend());
  mprintf("\tSource mask     :");
  srcMask.BriefMaskInfo();
  mprintf("\n");
  mprintf("\tTarget mask     :");
  tgtMask.BriefMaskInfo();
  mprintf("\n");
  if (doRmsFit) {
    mprintf(  "\tSource fit mask :");
    srcFitMask.BriefMaskInfo();
    mprintf("\n\tTarget fit mask :");
    tgtFitMask.BriefMaskInfo();
    mprintf("\n");
    if (srcFitMask.Nselected() != tgtFitMask.Nselected()) {
      mprinterr("Error: RMS-fit requires same # of atoms selected in source and target.\n");
      return CpptrajState::ERR;
    }
    // Source gets RMS fit to target (reference)
    Frame srcFitFrame;
    srcFitFrame.SetupFrameFromMask(srcFitMask, srcCoords->Top().Atoms());
    srcFitFrame.SetCoordinates(srcFrame, srcFitMask);
    Frame tgtFitFrame;
    tgtFitFrame.SetupFrameFromMask(tgtFitMask, tgtCoords->Top().Atoms());
    tgtFitFrame.SetCoordinates(tgtFrame, tgtFitMask);
    Vec3 refTrans = tgtFitFrame.CenterOnOrigin(false);
    Matrix_3x3 Rot;
    Vec3 Trans;
    srcFitFrame.RMSD_CenteredRef( tgtFitFrame, Rot, Trans, false );
    srcFrame.Trans_Rot_Trans( Trans, Rot, refTrans );
  }

  // Modify source if needed.
  Topology* srcTopPtr = srcCoords->TopPtr();
  Frame*    srcFrmPtr = &srcFrame;
  if (srcMask.Nselected() != srcCoords->Top().Natom()) {
    srcTopPtr = srcCoords->Top().modifyStateByMask( srcMask );
    if (srcTopPtr == 0) {
      mprinterr("Error: Could not modify source topology.\n");
      return CpptrajState::ERR;
    }
    srcFrmPtr = new Frame();
    srcFrmPtr->SetupFrameV(srcTopPtr->Atoms(), srcCoords->CoordsInfo());
    srcFrmPtr->SetFrame(srcFrame, srcMask);
  }

  // Modify target if needed.
  Topology* tgtTopPtr = tgtCoords->TopPtr();
  Frame*    tgtFrmPtr = &tgtFrame;
  if (tgtMask.Nselected() != tgtCoords->Top().Natom()) {
    tgtTopPtr = tgtCoords->Top().modifyStateByMask( tgtMask );
    if (tgtTopPtr == 0) {
      mprinterr("Error: Could not modify target topology.\n");
      return CpptrajState::ERR;
    }
    tgtFrmPtr = new Frame();
    tgtFrmPtr->SetupFrameV(tgtTopPtr->Atoms(), tgtCoords->CoordsInfo());
    tgtFrmPtr->SetFrame(tgtFrame, tgtMask);
  }

  // Combine topologies. Use target box info.
  Topology combinedTop;
  combinedTop.SetDebug( State.Debug() );
  combinedTop.SetParmName( outCoords->Meta().Name(), FileName() );
  combinedTop.AppendTop( *tgtTopPtr );
  combinedTop.AppendTop( *srcTopPtr );
  // Add any bonds
  for (unsigned int ii = 0; ii != tgtBondAtoms.size(); ii++)
    combinedTop.AddBond( tgtBondAtoms[ii], srcBondAtoms[ii] );
  combinedTop.SetParmBox( tgtFrmPtr->BoxCrd() );
  combinedTop.Brief("Grafted parm:");

  // Output coords.
  // Only coords+box for now.
  CoordinateInfo outInfo(tgtFrmPtr->BoxCrd(), false, false, false);
  if (outCoords->CoordsSetup(combinedTop, outInfo)) return CpptrajState::ERR;
  Frame CombinedFrame = outCoords->AllocateFrame();
  std::copy(tgtFrmPtr->xAddress(), tgtFrmPtr->xAddress()+tgtFrmPtr->size(), CombinedFrame.xAddress());
  std::copy(srcFrmPtr->xAddress(), srcFrmPtr->xAddress()+srcFrmPtr->size(), CombinedFrame.xAddress()+tgtFrmPtr->size());
  CombinedFrame.SetBox( tgtFrmPtr->BoxCrd() );
  outCoords->AddFrame( CombinedFrame );

  // Free memory if needed
  if (srcTopPtr != srcCoords->TopPtr()) {
    delete srcTopPtr;
    delete srcFrmPtr;
  }
  if (tgtTopPtr != tgtCoords->TopPtr()) {
    delete tgtTopPtr;
    delete tgtFrmPtr;
  }

  return CpptrajState::OK;
}
