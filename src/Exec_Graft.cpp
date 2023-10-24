#include "Exec_Graft.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords.h"
#include "Structure/Builder.h"
#include <algorithm> // std::copy

using namespace Cpptraj::Structure;
// Exec_Graft::Help()
void Exec_Graft::Help() const
{
  mprintf("\tsrc <source COORDS> [srcframe <#>] [srcfitmask <mask>] [srcmask <mask>]\n"
          "\t[srccharge <charge>\n"
          "\ttgt <target COORDS> [tgtframe <#>] [tgtfitmask <mask>] [tgtmask <mask>]\n"
          "\t[tgtcharge <charge>\n"
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

/** Redistribute charge on atoms in topology to match a target charge. */
int Exec_Graft::redistribute_charge(Topology& topIn, double charge) {
  //mprintf("DEBUG: Redistribute charge for %s, total charge = %g\n", topIn.c_str(), charge);
  double pcharge = 0;
  double ncharge = 0;
  for (int iat = 0; iat != topIn.Natom(); iat++) {
    if (topIn[iat].Charge() > 0)
      pcharge += topIn[iat].Charge();
    else if (topIn[iat].Charge() < 0)
      ncharge += topIn[iat].Charge();
  }
  //if (fabs(pcharge) < Constants::SMALL)
  bool PchargeZero = false;
  if (pcharge == 0) {
    mprintf("\tTotal positive charge is 0.0\n");
    PchargeZero = true;
  }
  bool NchargeZero = false;
  //if (fabs(ncharge) < Constants::SMALL)
  if (ncharge == 0) {
    mprintf("\tTotal negative charge is 0.0\n");
    NchargeZero = true;
  }
  if (!PchargeZero && !NchargeZero) {
    //double total_charge = 0;
    for (int iat = 0; iat != topIn.Natom(); iat++) {
      double delta = topIn[iat].Charge() * (charge - pcharge - ncharge) / (pcharge - ncharge);
      if (topIn[iat].Charge() >= 0) {
        topIn.SetAtom(iat).SetCharge( topIn[iat].Charge() + delta );
      } else {
        topIn.SetAtom(iat).SetCharge( topIn[iat].Charge() - delta );
      }
      //total_charge += topIn[iat].Charge();
    }
    //mprintf("DEBUG: Total charge after redistribute: %g\n", total_charge);
  }
  return 0;
}

// Exec_Graft::Execute()
Exec::RetType Exec_Graft::Execute(CpptrajState& State, ArgList& argIn)
{
  if (argIn.hasKey("ic"))
    return graft_ic(State,argIn);
  else
    return graft_rms(State, argIn);
}

/** Get COORDS set. */
DataSet_Coords* Exec_Graft::get_crd(ArgList& argIn, DataSetList const& DSL,
                                    const char* key, const char* desc, Frame& srcFrame, const char* frameKey)
{
  std::string kw = argIn.GetStringKey(key);
  if (kw.empty()) {
    mprinterr("Error: %s must be specified with '%s'.\n", desc, key);
    return 0;
  }
  DataSet_Coords* srcCoords = (DataSet_Coords*)DSL.FindSetOfGroup(kw, DataSet::COORDINATES);
  if (srcCoords == 0) {
    mprinterr("Error: %s %s not found.\n", desc, kw.c_str());
    return 0;
  }
  srcFrame = srcCoords->AllocateFrame();
  srcCoords->GetFrame(argIn.getKeyInt(frameKey, 1)-1, srcFrame);
  return srcCoords;
}

/** Modify given topology and frame by specified mask. Ensure the bonding
  * atom will be present in the stripped topology.
  * \return topology modified by mask.
  */
Topology* Exec_Graft::modify_top(Topology const& topIn, AtomMask const& mask, Frame& srcFrame)
                                 
{
  mprintf("\tAtoms to keep from '%s' : ", topIn.c_str());
  mask.BriefMaskInfo();
  mprintf("\n");

  Topology* srcTopPtr = topIn.modifyStateByMask( mask );
  if (srcTopPtr == 0) {
    mprinterr("Error: Could not modify topology %s.\n", topIn.c_str());
    return 0;
  }
  Frame* srcFrmPtr = new Frame();
  srcFrmPtr->SetupFrameV(srcTopPtr->Atoms(), srcFrame.CoordsInfo());
  srcFrmPtr->SetFrame(srcFrame, mask);
  srcFrame = *srcFrmPtr;
  delete srcFrmPtr;
  return srcTopPtr;
}

/** Select bond atom index. */
int Exec_Graft::select_bond_idx(std::string const& bond0maskstr, Topology const& mol0Top) {
  // Select bond atom indices
  AtomMask bondmask0;
  if (bondmask0.SetMaskString( bond0maskstr )) return -1;
  if (mol0Top.SetupIntegerMask( bondmask0 )) return -1;
  if (bondmask0.None()) {
    mprinterr("Error: Bond mask '%s' selects no atoms in topology '%s'\n", bondmask0.MaskString(), mol0Top.c_str());
    return -1;
  }
  if (bondmask0.Nselected() > 1) {
    mprinterr("Error: Bond mask '%s' selects more than 1 atom in topology '%s'\n", bondmask0.MaskString(), mol0Top.c_str());
    return -1;
  }
  return bondmask0[0];
}

/** Graft using internal coordinates to build the final structure. */
Exec::RetType Exec_Graft::graft_ic(CpptrajState& State, ArgList& argIn)
const
{
  // Source (fragment)
  Frame mol1frm;
  DataSet_Coords* mol1crd = get_crd(argIn, State.DSL(), "src", "Source COORDS", mol1frm, "srcframe");
  if (mol1crd == 0) return CpptrajState::ERR;
  // Target (base)
  Frame mol0frm;
  DataSet_Coords* mol0crd = get_crd(argIn, State.DSL(), "tgt", "Target COORDS", mol0frm, "tgtframe");
  if (mol0crd == 0) return CpptrajState::ERR;
  // Create output coords
  std::string kw = argIn.GetStringKey("name");
  if (kw.empty()) {
    mprinterr("Error: Output COORDS must be specified with 'name'.\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords* outCoords = (DataSet_Coords*)State.DSL().AddSet(DataSet::COORDS, MetaData(kw));
  if (outCoords == 0) {
    mprinterr("Error: Output COORDS %s could not be created.\n", kw.c_str());
    return CpptrajState::ERR;
  }
 
  // Get atoms to bond
  std::string bondargstr = argIn.GetStringKey("bond");
  if (bondargstr.empty()) {
    mprinterr("Error: No 'bond' keyword specified.\n");
    return CpptrajState::ERR;
  }
  ArgList bondarg( bondargstr, "," );
  if (bondarg.Nargs() != 2) {
    mprinterr("Error: Expected 2 comma-separated masks for 'bond' keyword, got %i\n", bondarg.Nargs());
    return CpptrajState::ERR;
  }
  std::string const& tgtbondmask = bondarg[0];
  std::string const& srcbondmask = bondarg[1];

  // Get atoms to keep from source.
  AtomMask mol1Mask;
  if (mol1Mask.SetMaskString( argIn.GetStringKey("srcmask") ))
    return CpptrajState::ERR;
  if (mol1crd->Top().SetupIntegerMask(mol1Mask))
    return CpptrajState::ERR;
  // Get atoms to keep from target.
  AtomMask mol0Mask;
  if (mol0Mask.SetMaskString( argIn.GetStringKey("tgtmask") ))
    return CpptrajState::ERR;
  if (mol0crd->Top().SetupIntegerMask(mol0Mask))
    return CpptrajState::ERR;

  // Modify the target (mol0) topology, update bond atom indices.
  bool newMol0Top = false;
  Topology* mol0Top = 0;
  if (mol0Mask.Nselected() == mol0crd->Top().Natom()) {
    mol0Top = mol0crd->TopPtr();
  } else {
    newMol0Top = true;
    mol0Top = modify_top(mol0crd->Top(), mol0Mask, mol0frm);
    if (mol0Top == 0) return CpptrajState::ERR; // FIXME need to free mol0Top memory
  }
  // Modify the source (mol1) topology, update bond atom indices.
  bool newMol1Top = false;
  Topology* mol1Top = 0;
  if (mol1Mask.Nselected() == mol1crd->Top().Natom()) {
    mol1Top = mol1crd->TopPtr();
  } else {
    newMol1Top = true;
    mol1Top = modify_top(mol1crd->Top(), mol1Mask, mol1frm);
    if (mol1Top == 0) return CpptrajState::ERR; // FIXME need to free mol1Top memory
  }

  // Select bond atom indices
  int bondat0 = select_bond_idx(tgtbondmask, *mol0Top);
  if (bondat0 < 0) {
    mprinterr("Error: Could not select target bond atom '%s'\n", tgtbondmask.c_str());
    return CpptrajState::ERR;
  }
  int bondat1 = select_bond_idx(srcbondmask, *mol1Top);
  if (bondat1 < 0) {
    mprinterr("Error: Could not select source bond atom '%s'\n", srcbondmask.c_str());
    return CpptrajState::ERR;
  }

  // Combine topologies.
  Topology combinedTop;
  combinedTop.SetDebug( State.Debug() );
  combinedTop.SetParmName( outCoords->Meta().Name(), FileName() );
  combinedTop.AppendTop( *mol0Top );

  Frame CombinedFrame = mol0frm;
  Builder builder;
  if (builder.Combine( combinedTop, CombinedFrame, *mol1Top, mol1frm, bondat0, bondat1 )) {
    mprinterr("Error: Fragment combine failed.\n");
    return CpptrajState::ERR;
  }

  // Add topology to output COORDS set
  if (outCoords->CoordsSetup(combinedTop, CombinedFrame.CoordsInfo())) return CpptrajState::ERR; // FIXME free molXTop memory

  // Add frame to the output data set
  outCoords->AddFrame( CombinedFrame );

  if (newMol0Top) delete mol0Top;
  if (newMol1Top) delete mol1Top;

  return CpptrajState::OK;
}

/** Graft with RMS-fitting. */
Exec::RetType Exec_Graft::graft_rms(CpptrajState& State, ArgList& argIn)
const
{
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
  bool hasSrcCharge = argIn.Contains("srccharge");
  double srccharge = argIn.getKeyDouble("srccharge", 0);
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
  bool hasTgtCharge = argIn.Contains("tgtcharge");
  double tgtcharge = argIn.getKeyDouble("tgtcharge", 0);
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
  if (hasSrcCharge && (srcMask.Nselected() != srcCoords->Top().Natom()))
    mprintf("\tAdjusting source charge to %g\n", srccharge);
  if (hasTgtCharge && (tgtMask.Nselected() != tgtCoords->Top().Natom()))
    mprintf("\tAdjusting target charge to %g\n", tgtcharge);
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
    // Modify charges if needed
    if (hasSrcCharge) {
      if (redistribute_charge(*srcTopPtr, srccharge)) {
        mprinterr("Error: Redistribute src charge failed.\n");
        return CpptrajState::ERR;
      }
    }
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
    // Modify charges if needed
    if (hasTgtCharge) {
      if (redistribute_charge(*tgtTopPtr, tgtcharge)) {
        mprinterr("Error: Redistribute tgt charge failed.\n");
        return CpptrajState::ERR;
      }
    }
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
  // Regenerate the molecule info FIXME should Topology just do this?
  if (combinedTop.DetermineMolecules()) return CpptrajState::ERR;
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
