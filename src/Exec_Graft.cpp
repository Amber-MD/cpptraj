#include "Exec_Graft.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords.h"
#include "Structure/Builder.h"
#include <algorithm> // std::copy

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Exec_Graft::Exec_Graft() :
  Exec(COORDS),
  debug_(0),
  newMol0Top_(0),
  newMol1Top_(0)
{}

/** DESTRUCTOR */
Exec_Graft::~Exec_Graft() {
  if (newMol0Top_ != 0) delete newMol0Top_;
  if (newMol1Top_ != 0) delete newMol1Top_;
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

/** Modify given topology and frame by specified mask.
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

// Exec_Graft::Help()
void Exec_Graft::Help() const
{
  mprintf("\tsrc <source COORDS> [srcframe <#>] [srcmask <mask> [srccharge <charge>]]\n"
          "\ttgt <target COORDS> [tgtframe <#>] [tgtmask <mask> [tgtcharge <charge>]]\n"
          "\t{ic | [srcfitmask <mask>] [tgtfitmask <mask>]}\n"
          "\tname <output COORDS> [bond <tgt>,<src> ...]\n"
          "  Graft coordinates from source to coordinates in target.\n"
          "  If 'ic' is specified use internal coordinates to link the coordinates,\n"
          "  otherwise rely on rms-fitting. If 'ic' is specified, exactly 1 bond\n"
          "  must be specified.\n");
}

// Exec_Graft::Execute()
Exec::RetType Exec_Graft::Execute(CpptrajState& State, ArgList& argIn)
{
  debug_ = State.Debug();
  bool use_ic = argIn.hasKey("ic");

  // Source (1, fragment)
  Frame mol1frm;
  DataSet_Coords* mol1crd = get_crd(argIn, State.DSL(), "src", "Source COORDS", mol1frm, "srcframe");
  if (mol1crd == 0) return CpptrajState::ERR;
  // Target (0, base)
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
  Sarray bond0ArgStrings;
  Sarray bond1ArgStrings;
  std::string bondargstr = argIn.GetStringKey("bond");
  while (!bondargstr.empty()) {
    ArgList bondarg( bondargstr, "," );
    if (bondarg.Nargs() != 2) {
      mprinterr("Error: Expected 2 comma-separated masks for 'bond' keyword, got %i\n", bondarg.Nargs());
      return CpptrajState::ERR;
    }
    bond0ArgStrings.push_back(bondarg[0]);
    bond1ArgStrings.push_back(bondarg[1]);
    bondargstr = argIn.GetStringKey("bond");
  }

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

  // Determine if charges will be redistributed
  bool hasSrcCharge = argIn.Contains("srccharge");
  double srccharge = argIn.getKeyDouble("srccharge", 0);
  bool hasTgtCharge = argIn.Contains("tgtcharge");
  double tgtcharge = argIn.getKeyDouble("tgtcharge", 0);

  // Get atoms from source to fit on target, and atoms from target
  // for source to fit on.
  AtomMask srcFitMask, tgtFitMask;
  bool doRmsFit = false;
  if (!use_ic) { // TODO allow with IC?
    std::string srcfitstr = argIn.GetStringKey("srcfitmask");
    std::string tgtfitstr = argIn.GetStringKey("tgtfitmask");
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
      if (mol1crd->Top().SetupIntegerMask(srcFitMask))
        return CpptrajState::ERR;
      if (mol0crd->Top().SetupIntegerMask(tgtFitMask))
        return CpptrajState::ERR;
    }
  }

  // Info
  mprintf("\tSource coords   : %s\n", mol1crd->legend());
  mprintf("\tTarget coords   : %s\n", mol0crd->legend());
  mprintf("\tOutput coords   : %s\n", outCoords->legend());
  mprintf("\tSource mask     :");
  mol1Mask.BriefMaskInfo();
  mprintf("\n");
  mprintf("\tTarget mask     :");
  mol0Mask.BriefMaskInfo();
  mprintf("\n");
  if (!bond0ArgStrings.empty()) {
    mprintf("\tBonds will be created between:\n");
    for (unsigned int idx = 0; idx != bond0ArgStrings.size(); idx++)
      mprintf("\t\tAtoms %s and %s\n", bond0ArgStrings[idx].c_str(), bond1ArgStrings[idx].c_str());
  }
  if (hasSrcCharge && (mol1Mask.Nselected() != mol1crd->Top().Natom()))
    mprintf("\tAdjusting source charge to %g\n", srccharge);
  if (hasTgtCharge && (mol0Mask.Nselected() != mol0crd->Top().Natom()))
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
    srcFitFrame.SetupFrameFromMask(srcFitMask, mol1crd->Top().Atoms());
    srcFitFrame.SetCoordinates(mol1frm, srcFitMask);
    Frame tgtFitFrame;
    tgtFitFrame.SetupFrameFromMask(tgtFitMask, mol0crd->Top().Atoms());
    tgtFitFrame.SetCoordinates(mol0frm, tgtFitMask);
    Vec3 refTrans = tgtFitFrame.CenterOnOrigin(false);
    Matrix_3x3 Rot;
    Vec3 Trans;
    srcFitFrame.RMSD_CenteredRef( tgtFitFrame, Rot, Trans, false );
    mol1frm.Trans_Rot_Trans( Trans, Rot, refTrans );
  }

  // Modify the target (mol0) topology
  if (newMol0Top_ != 0) delete newMol0Top_;
  newMol0Top_ = 0;
  Topology* mol0Top = 0;
  if (mol0Mask.Nselected() == mol0crd->Top().Natom()) {
    mol0Top = mol0crd->TopPtr();
  } else {
    mol0Top = modify_top(mol0crd->Top(), mol0Mask, mol0frm);
    newMol0Top_ = mol0Top;
    if (mol0Top == 0) return CpptrajState::ERR;
    // Modify charges if needed
    if (hasTgtCharge) {
      if (mol0Top->RedistributeCharge(tgtcharge)) {
        mprinterr("Error: Redistribute tgt charge failed.\n");
        return CpptrajState::ERR;
      }
    }
  }
  // Modify the source (mol1) topology
  if (newMol1Top_ != 0) delete newMol1Top_;
  newMol1Top_ = 0;
  Topology* mol1Top = 0;
  if (mol1Mask.Nselected() == mol1crd->Top().Natom()) {
    mol1Top = mol1crd->TopPtr();
  } else {
    mol1Top = modify_top(mol1crd->Top(), mol1Mask, mol1frm);
    newMol1Top_ = mol1Top;
    if (mol1Top == 0) return CpptrajState::ERR;
    // Modify charges if needed
    if (hasSrcCharge) {
      if (mol1Top->RedistributeCharge(srccharge)) {
        mprinterr("Error: Redistribute src charge failed.\n");
        return CpptrajState::ERR;
      }
    }
  }

  if (use_ic) {
    if (graft_ic( outCoords, *mol0Top, mol0frm, *mol1Top, mol1frm, bond0ArgStrings, bond1ArgStrings))
      return CpptrajState::ERR;
  } else {
    if (graft_rms( outCoords, *mol0Top, mol0frm, *mol1Top, mol1frm, bond0ArgStrings, bond1ArgStrings))
      return CpptrajState::ERR;
  }
  return CpptrajState::OK;
}


/** Graft using internal coordinates to build the final structure. */
int Exec_Graft::graft_ic(DataSet_Coords* outCoords,
                         Topology const& mol0Top, Frame const& mol0frm,
                         Topology const& mol1Top, Frame const& mol1frm,
                         Sarray const& bond0Atoms, Sarray const& bond1Atoms)
const
{
  // Get bonding atom masks
  if (bond0Atoms.size() != 1 || bond1Atoms.size() != 1) {
    mprinterr("Error: Graft with internal coordinates only works with 1 bond.\n");
    return 1;
  }
  std::string const& tgtbondmask = bond0Atoms[0];
  std::string const& srcbondmask = bond1Atoms[0];

  // Select bond atom indices
  int bondat0 = select_bond_idx(tgtbondmask, mol0Top);
  if (bondat0 < 0) {
    mprinterr("Error: Could not select target bond atom '%s'\n", tgtbondmask.c_str());
    return 1;
  }
  int bondat1 = select_bond_idx(srcbondmask, mol1Top);
  if (bondat1 < 0) {
    mprinterr("Error: Could not select source bond atom '%s'\n", srcbondmask.c_str());
    return 1;
  }

  // Combine topologies.
  Topology combinedTop;
  combinedTop.SetDebug( debug_ );
  combinedTop.SetParmName( outCoords->Meta().Name(), FileName() );
  combinedTop.AppendTop( mol0Top );
  combinedTop.SetParmBox( mol0frm.BoxCrd() );
  if (debug_ > 0)
    combinedTop.Brief("Grafted parm:");

  Frame CombinedFrame = mol0frm;
  Builder builder;
  builder.SetDebug( debug_ );
  if (builder.Combine( combinedTop, CombinedFrame, mol1Top, mol1frm, bondat0, bondat1 )) {
    mprinterr("Error: Fragment combine failed.\n");
    return 1;
  }

  // Add topology to output COORDS set
  if (outCoords->CoordsSetup(combinedTop, CombinedFrame.CoordsInfo())) return 1;

  // Add frame to the output data set
  outCoords->AddFrame( CombinedFrame );

  return 0;
}

/** Graft with RMS-fitting. */
int Exec_Graft::graft_rms(DataSet_Coords* outCoords,
                                    Topology const& mol0Top, Frame const& mol0frm,
                                    Topology const& mol1Top, Frame const& mol1frm,
                                    Sarray const& bond0Atoms, Sarray const& bond1Atoms)
const
{
  // Combine topologies. Use target box info.
  Topology combinedTop;
  combinedTop.SetDebug( debug_ );
  combinedTop.SetParmName( outCoords->Meta().Name(), FileName() );
  combinedTop.AppendTop( mol0Top );
  combinedTop.AppendTop( mol1Top );

  // Add any bonds
  for (unsigned int ii = 0; ii != bond0Atoms.size(); ii++) {
    int bondat0 = select_bond_idx(bond0Atoms[ii], mol0Top);
    if (bondat0 < 0) {
      mprinterr("Error: Could not select target bond atom '%s'\n", bond0Atoms[ii].c_str());
      return 1;
    }
    int bondat1 = select_bond_idx(bond1Atoms[ii], mol1Top);
    if (bondat1 < 0) {
      mprinterr("Error: Could not select source bond atom '%s'\n", bond1Atoms[ii].c_str());
      return 1;
    }
    combinedTop.AddBond( bondat0, bondat1 + mol0Top.Natom() ); // TODO pseudo-parameter?
  }
  // Regenerate the molecule info FIXME should Topology just do this?
  if (combinedTop.DetermineMolecules()) return 1;
  combinedTop.SetParmBox( mol0frm.BoxCrd() );
  combinedTop.Brief("Grafted parm:");
  // Only coords+box for now.
  CoordinateInfo outInfo(mol0frm.BoxCrd(), false, false, false);
  if (outCoords->CoordsSetup(combinedTop, outInfo)) return 1;

  // Combine coords.
  Frame CombinedFrame = outCoords->AllocateFrame();
  std::copy(mol0frm.xAddress(), mol0frm.xAddress()+mol0frm.size(), CombinedFrame.xAddress());
  std::copy(mol1frm.xAddress(), mol1frm.xAddress()+mol1frm.size(), CombinedFrame.xAddress()+mol0frm.size());
  CombinedFrame.SetBox( mol0frm.BoxCrd() );

  // Add to topology
  outCoords->AddFrame( CombinedFrame );

  return 0;
}
