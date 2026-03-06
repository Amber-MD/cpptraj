#include "Exec_Graft.h"
#include "AssociatedData_Connect.h"
#include "CharMask.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords.h"
#include "Structure/Builder.h"
#include "Structure/Zmatrix.h"
#include <algorithm> // std::copy
#include <utility> // std::pair

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Exec_Graft::Exec_Graft() :
  Exec(COORDS),
  debug_(0),
  verbose_(0),
  newMol0Top_(0),
  newMol1Top_(0),
  hasOrient0_(false),
  hasOrient1_(false),
//  orient0_(0),
//  orient1_(0),
  chi0_(0),
  chi1_(0)
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

/** Get connect atoms from associated data */
Exec_Graft::Iarray Exec_Graft::getConnectAtoms( AssociatedData* ad )
{
  if (ad == 0) return Iarray();
  AssociatedData_Connect const& CONN = static_cast<AssociatedData_Connect const&>( *ad );
  return CONN.Connect();
}

/** Print connect atoms to stdout */
void Exec_Graft::print_connect(const char* desc, Iarray const& connect, Topology const& top)
{
  if (connect.empty()) return;
  mprintf("\t  %s :", desc);
  for (Iarray::const_iterator it = connect.begin(); it != connect.end(); ++it)
    if (*it > -1)
      mprintf(" %s", top.AtomMaskName(*it).c_str());
  mprintf("\n");
}

// Exec_Graft::Help()
void Exec_Graft::Help() const
{
  mprintf("\tsrc <source COORDS> [srcframe <#>] [srcmask <mask> [srccharge <charge>]]\n"
          "\ttgt <target COORDS> [tgtframe <#>] [tgtmask <mask> [tgtcharge <charge>]]\n"
          "\t{ic | [srcfitmask <mask>] [tgtfitmask <mask>]} [verbose]\n"
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
  verbose_ = argIn.getKeyInt("verbose", 0);

  // Source (1, fragment)
  Frame mol1frm;
  DataSet_Coords* mol1crd = get_crd(argIn, State.DSL(), "src", "Source COORDS", mol1frm, "srcframe");
  if (mol1crd == 0) return CpptrajState::ERR;
  AssociatedData* ad1 = mol1crd->GetAssociatedData(AssociatedData::CONNECT);
  Iarray connect1 = getConnectAtoms( ad1 );
  // Target (0, base)
  Frame mol0frm;
  DataSet_Coords* mol0crd = get_crd(argIn, State.DSL(), "tgt", "Target COORDS", mol0frm, "tgtframe");
  if (mol0crd == 0) return CpptrajState::ERR;
  AssociatedData* ad0 = mol0crd->GetAssociatedData(AssociatedData::CONNECT);
  Iarray connect0 = getConnectAtoms( ad0 );

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
  print_connect("Source connect atoms", connect1, mol1crd->Top());
  mprintf("\tTarget coords   : %s\n", mol0crd->legend());
  print_connect("Target connect atoms", connect0, mol0crd->Top());
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

  if (use_ic) {
    if (bond0ArgStrings.empty()) {
      // Do we have connect atoms?
      if (!connect0.empty() && !connect1.empty()) {
        // We have connect atoms.
        int tail0 = -1;
        int tail1 = -1;
        int head0 = connect0[0];
        int head1 = connect1[0];
        if (connect0.size() > 1) tail0 = connect0[1];
        if (connect1.size() > 1) tail1 = connect1[1];
        // If we only have one of each that is easy
        if (head0 != -1 && tail0 == -1 && head1 == -1 && tail1 != -1) {
          // Head0, tail1
          bond0ArgStrings.push_back( mol0crd->Top().AtomMaskName(head0) );
          bond1ArgStrings.push_back( mol1crd->Top().AtomMaskName(tail1) );
        } else if (head0 == -1 && tail0 != -1 && head1 != -1 && tail1 == -1) {
          // Head1, tail0
          bond0ArgStrings.push_back( mol0crd->Top().AtomMaskName(tail0) );
          bond1ArgStrings.push_back( mol1crd->Top().AtomMaskName(head1) );
        } else if (head0 != -1 && tail0 != -1 && head1 != -1 && tail1 != -1) {
          // Default head of source (1) to tail of target (0)
          bond0ArgStrings.push_back( mol0crd->Top().AtomMaskName(tail0) );
          bond1ArgStrings.push_back( mol1crd->Top().AtomMaskName(head1) );
        } else {
          mprinterr("Error: No bonds specified and not enough CONNECT atoms:\n");
          mprinterr("Error: head0 %i tail0 %i head1 %i tail1 %i\n", head0+1, tail0+1, head1+1, tail1+1);
          return CpptrajState::ERR;
        }
      }
    }
    // Get internals for both topologies before they are modified.
    get_original_orientations(mol0crd->Top(), mol0frm, mol1crd->Top(), mol1frm,
                              mol0Mask, mol1Mask, bond0ArgStrings, bond1ArgStrings);
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

/** Get internals for molecules before modification. */
int Exec_Graft::get_original_orientations(Topology const& mol0Top, Frame const& mol0frm,
                                          Topology const& mol1Top, Frame const& mol1frm,
                                          AtomMask const& mol0mask, AtomMask const& mol1mask,
                                          Sarray const& bond0Atoms, Sarray const& bond1Atoms)
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

  // Will an atom bonded to this atom disappear?
//  orient0_ = 0.0;
  Atom const& At0 = mol0Top[bondat0];
  CharMask cmask0( mol0mask.ConvertToCharMask(), mol0mask.Nselected() );
  int vanish0idx = -1;
  for (Atom::bond_iterator bat = At0.bondbegin(); bat != At0.bondend(); ++bat) {
    if (!cmask0.AtomInCharMask( *bat )) { // TODO check multiple disappearing atoms
      if (debug_ > 0)
        mprintf("DEBUG: Atom0 %s will vanish.\n", mol0Top.AtomMaskName(*bat).c_str());
      if (vanish0idx == -1) vanish0idx = *bat;
    }
  }
  if (vanish0idx != -1) {
    hasOrient0_ = true;
    chi0_    = Builder::DetermineChiralityAroundAtom(bondat0, mol0frm, mol0Top);
    if (debug_ > 0)
      mprintf("DEBUG: Chirality around %s is %f\n", mol0Top.LeapName(bondat0).c_str(), chi0_);
//    orient0_ = Builder::CalculateOrientationAroundAtom(bondat0, vanish0idx, mol0frm, mol0Top);
  }

//  orient1_ = 0.0;
  Atom const& At1 = mol1Top[bondat1];
  CharMask cmask1( mol1mask.ConvertToCharMask(), mol1mask.Nselected() );
  int vanish1idx = -1;
  for (Atom::bond_iterator bat = At1.bondbegin(); bat != At1.bondend(); ++bat) {
    if (!cmask1.AtomInCharMask( *bat )) { // TODO check multiple disappearing atoms
      if (debug_ > 0)
        mprintf("DEBUG: Atom1 %s will vanish.\n", mol1Top.AtomMaskName(*bat).c_str());
      if (vanish1idx == -1) vanish1idx = *bat;
    }
  }
  if (vanish1idx != -1) {
    hasOrient1_ = true;
    chi1_    = Builder::DetermineChiralityAroundAtom(bondat1, mol1frm, mol1Top);
    if (debug_ > 0)
      mprintf("DEBUG: Chirality around %s is %f\n", mol1Top.LeapName(bondat1).c_str(), chi1_);
//    orient1_ = Builder::CalculateOrientationAroundAtom(bondat1, vanish1idx, mol1frm, mol1Top);
  }

  return 0;
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

  if (debug_ > 0) {
    for (int at = 0; at != mol0Top.Natom(); at++)
      mprintf("\t%6i %s %s\n", at+1, mol0Top.AtomMaskName(at).c_str(), *(mol0Top.Res(mol0Top[at].ResNum()).Name()));
    mprintf("BOND ATOM 0 %i\n", bondat0+1);
    for (int at = 0; at != mol1Top.Natom(); at++)
      mprintf("\t%6i %s %s\n", at+1, mol1Top.AtomMaskName(at).c_str(), *(mol1Top.Res(mol1Top[at].ResNum()).Name()));
    mprintf("BOND ATOM 1 %i\n", bondat1+1);
  }

  // Combine topologies.
  Topology topOut;
  Frame frameOut;
  int total_natom = mol0Top.Natom() + mol1Top.Natom();
  frameOut.SetupFrame( total_natom );
  // Clear frame so that AddXYZ can be used
  frameOut.ClearAtoms();
  // hasPosition - for each atom in topOut, status on whether atom in frameOut needs building
  Cpptraj::Structure::Zmatrix::Barray hasPosition;
  hasPosition.reserve( total_natom );
  // For saving intra-res bonds
  typedef std::pair<int,int> Ipair;
  typedef std::vector<Ipair> IParray;
  IParray intraResBonds;
  // Add mol0 atoms to topology
  topOut.SetDebug( debug_ );
  topOut.SetParmName( outCoords->Meta().Name(), FileName() );
  topOut.SetParmBox( mol0frm.BoxCrd() );
  for (int at = 0; at < mol0Top.Natom(); at++) {
    Atom sourceAtom = mol0Top[at];
    // Save the intra-residue bonds.
    for (Atom::bond_iterator bat = sourceAtom.bondbegin(); bat != sourceAtom.bondend(); ++bat) {
      if (*bat > at) {
        if (debug_ > 1)
          mprintf("Will add bond between %i and %i\n", at+1, *bat+1);
        intraResBonds.push_back( Ipair(at, *bat) );
      }
    }
    sourceAtom.ClearBonds(); // FIXME AddTopAtom should clear bonds
    topOut.AddTopAtom( sourceAtom, mol0Top.Res(mol0Top[at].ResNum()) );
    frameOut.AddVec3( Vec3(mol0frm.XYZ(at)) );
    hasPosition.push_back( true );
  }
  // Add mol1 atoms
  int atomOffset = mol0Top.Natom();
  int resOffset = topOut.Nres();
  if (debug_ > 0)
    mprintf("DEBUG: Atom offset is %i\n", atomOffset);
  for (int itgt = 0; itgt < mol1Top.Natom(); itgt++) {
    Atom sourceAtom = mol1Top[itgt];
    Residue currentRes = mol1Top.Res(sourceAtom.ResNum());
    currentRes.SetOriginalNum( currentRes.OriginalResNum() + resOffset );
    // Save the intra-residue bonds.
    int at0 = itgt + atomOffset;
    for (Atom::bond_iterator bat = sourceAtom.bondbegin(); bat != sourceAtom.bondend(); ++bat) {
      int at1 = *bat + atomOffset;
      if (at1 > at0) {
        if (debug_ > 1)
          mprintf("Will add bond between %i and %i (original %i and %i)\n", at0+1, at1+1, itgt+1, *bat + 1);
        intraResBonds.push_back( Ipair(at0, at1) );
      }
    }
    sourceAtom.ClearBonds(); // FIXME AddTopAtom should clear bonds
    topOut.AddTopAtom( sourceAtom, currentRes );
    frameOut.AddVec3( Vec3(mol1frm.XYZ(itgt)) );
    hasPosition.push_back( false );
  }
  //Add intra-residue bonds
  for (IParray::const_iterator it = intraResBonds.begin(); it != intraResBonds.end(); ++it)
  {
    //mprintf("DEBUG: Intra-res bond: Res %s atom %s to res %s atom %s\n",
    //        topOut.TruncResNameOnumId(topOut[it->first].ResNum()).c_str(), *(topOut[it->first].Name()),
    //        topOut.TruncResNameOnumId(topOut[it->second].ResNum()).c_str(), *(topOut[it->second].Name()));
    topOut.AddBond(it->first, it->second);
  }
  // DEBUG
  //for (int at = 0; at != topOut.Natom(); at++)
  //  mprintf("\t%6i %s %s\n", at+1, topOut.AtomMaskName(at).c_str(), *(topOut.Res(topOut[at].ResNum()).Name()));
  // Build
  Cpptraj::Structure::Builder structureBuilder;
  structureBuilder.SetDebug( debug_ );
  if (structureBuilder.GenerateInternals( mol1frm, mol1Top,
                                          std::vector<bool>(mol1Top.Natom(), true) ))
  {
    mprinterr("Error: Generate internals for %s failed.\n", mol1Top.c_str());
    return 1;
  }
  structureBuilder.UpdateIndicesWithOffset( atomOffset );
  // Connect unit. Head atom of second unit comes first to match LEaP.
  mprintf("\tConnect %s (%i) and %s (%i, original %s and %s)\n",
          topOut.AtomMaskName(bondat1+atomOffset).c_str(), bondat1+atomOffset+1,
          topOut.AtomMaskName(bondat0).c_str(), bondat0+1,
          mol1Top.AtomMaskName(bondat1).c_str(),
          mol0Top.AtomMaskName(bondat0).c_str());
  topOut.AddBond( bondat1 + atomOffset, bondat0 );
  // Set any saved orientations
  if (hasOrient0_)
    structureBuilder.SetAtomChirality(bondat0, chi0_);
  if (hasOrient1_)
    structureBuilder.SetAtomChirality(bondat1 + atomOffset, chi1_);
  // Generate internals around the link
  if (structureBuilder.GenerateInternalsAroundLink(bondat1 + atomOffset, bondat0,
                                                   frameOut, topOut, hasPosition,
                                                   Cpptraj::Structure::Builder::SEQUENCE) )
  {
    mprinterr("Error: Assign torsions around graft bond atoms %s - %s failed.\n",
              topOut.AtomMaskName(bondat1 + atomOffset).c_str(),
              topOut.AtomMaskName(bondat0).c_str());
    return 1;
  }
  // Adjust torsions around link so that longest 'path' is trans
  if (structureBuilder.AdjustIcAroundLink(bondat0, bondat1 + atomOffset, frameOut, topOut))
  {
    mprinterr("Error: Failed to adjust internal coords around the link.\n");
    return 1;
  }
  // Update internal coords from known positions
  // NOTE: By definition, there are no known positions
  //if (structureBuilder.UpdateICsFromFrame( frameOut, topOut, hasPosition )) {
  //  mprinterr("Error: Failed to update Zmatrix with values from existing positions.\n");
  //  return 1;
  //}
  // Convert to Zmatrix and assign missing atom positions
  if (structureBuilder.BuildSequenceFromInternals(frameOut, topOut, hasPosition,
                                                  bondat1 + atomOffset,
                                                  bondat0))
  {
    mprinterr("Error: Grafting %s with %s build from internals failed.\n",
              mol0Top.c_str(), mol1Top.c_str());
    return 1;
  }

  // Finalize topology - determine molecules, dont renumber residues
  topOut.CommonSetup(true, false);
  topOut.Summary();
  // Add to output data set
  if (outCoords->CoordsSetup(topOut, frameOut.CoordsInfo())) return 1;
  outCoords->AddFrame( frameOut );

  return 0;
}

/** Graft with RMS-fitting. */
int Exec_Graft::graft_rms(DataSet_Coords* outCoords,
                                    Topology const& mol0Top, Frame const& mol0frm,
                                    Topology const& mol1Top, Frame const& mol1frm,
                                    Sarray const& bond0Atoms, Sarray const& bond1Atoms)
const
{
  // Combine topologies. Use target box info. Do not merge bond/angle params TODO should they be reduced?
  Topology combinedTop;
  combinedTop.SetDebug( debug_ );
  combinedTop.SetParmName( outCoords->Meta().Name(), FileName() );
  combinedTop.AppendTop( mol0Top, verbose_, false, false );
  combinedTop.AppendTop( mol1Top, verbose_, false, false );

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
