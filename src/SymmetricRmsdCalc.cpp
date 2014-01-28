#include "SymmetricRmsdCalc.h"
#include "DistRoutines.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
SymmetricRmsdCalc::SymmetricRmsdCalc() : debug_(0), fit_(true), useMass_(false) {}

// CONSTRUCTOR - For use when only RMSD is wanted.
SymmetricRmsdCalc::SymmetricRmsdCalc(AtomMask const& maskIn, bool fitIn, 
                                     bool useMassIn, Topology const& topIn) :
  tgtMask_(maskIn), fit_(fitIn), useMass_(useMassIn)
{
  Topology* stripTop = topIn.partialModifyStateByMask( tgtMask_ );
  stripTop->Brief("SymmRMSD"); // DEBUG
  // Since input frames will already be stripped, make target mask have all atoms
  tgtMask_.SetMaskString(0);
  SetupSymmRMSD( *stripTop );
  delete stripTop;
}

// SymmetricRmsdCalc::InitSymmRMSD()
int SymmetricRmsdCalc::InitSymmRMSD(std::string const& tMaskExpr, bool fitIn,
                                    bool useMassIn, int debugIn)
{
  if (tgtMask_.SetMaskString( tMaskExpr )) return 1;
  debug_ = debugIn;
  fit_ = fitIn;
  useMass_ = useMassIn;
  return 0;
}
#ifdef DEBUGSYMMRMSD
static int recursionLevel_;
#endif
/** Recursive function to search for symmetric atoms. */
void SymmetricRmsdCalc::FindSymmetricAtoms(int at, AtomMap const& resmap,
                                           std::string const& Unique,
                                           Iarray& Selected,
                                           Iarray& symmatoms) const
{
  // If this atom has already been selected, leave
  if (Selected[at]) return;
  Selected[at] = 1;
# ifdef DEBUGSYMMRMSD
  ++recursionLevel_;
  for (int i = 0; i < recursionLevel_; i++)
    mprintf("..");
  mprintf("Atom %i(%s)", at+1, resmap[at].c_str());
# endif
  // Does this atom match the unique ID we are looking for?
  if (resmap[at].Unique() == Unique) {
    symmatoms.push_back( at ); // NOTE: This index is relative to the residue 
#   ifdef DEBUGSYMMRMSD
    mprintf(" SYMM");
#   endif
  }
  // Recursively search through all atoms bonded to this atom unless they 
  // are a chiral center
# ifdef DEBUGSYMMRMSD
  mprintf("\n");
# endif
  for (Atom::bond_iterator bndatm = resmap[at].bondbegin();
                           bndatm != resmap[at].bondend();
                           ++bndatm)
  {
    if (!resmap[*bndatm].IsChiral())
      FindSymmetricAtoms( *bndatm, resmap, Unique, Selected, symmatoms );
  }
}

// SymmetricRmsdCalc::SetupSymmRMSD()
/** Find potential symmetric atoms. All residues up to the last selected
  * residue are considered, but all atoms within those residues (even
  * unselected ones) because when symmetric atoms are re-mapped, atoms
  * bonded to the symmetric atoms (which are themselves symmetric) need
  * to be re-mapped as well.
  */
int SymmetricRmsdCalc::SetupSymmRMSD(Topology const& topIn) {
  // Setup target mask.
  if (topIn.SetupIntegerMask( tgtMask_ )) return 1;
  tgtMask_.MaskInfo();
  if (tgtMask_.None()) {
    mprintf("Warning: No atoms selected by mask '%s'\n", tgtMask_.MaskString());
    return 1;
  }
  // Allocate space for selected atoms in target frame. This will also
  // put the correct masses in based on the mask.
  selectedTgt_.SetupFrameFromMask(tgtMask_, topIn.Atoms());
  // Allocate space for remapped frame; same # atoms as original frame
  remapFrame_.SetupFrameV( topIn.Atoms(), topIn.HasVelInfo(), topIn.NrepDim() );
  // Create initial 1 to 1 atom map for all atoms; indices in 
  // SymmetricAtomIndices will correspond to positions in AMap.
  AMap_.resize( topIn.Natom() );
  // Determine last selected residue.
  int last_res = topIn[tgtMask_.back()].ResNum() + 1;
  mprintf("\tResidues up to %s will be considered for symmetry correction.\n",
          topIn.TruncResNameNum(last_res-1).c_str());
  // In each residue, determine which atoms are symmetric.
  SymmetricAtomIndices_.clear();
  AtomMap resmap;
  if (debug_ > 1) resmap.SetDebug(1);
  for (int res = 0; res < last_res; ++res) {
    int res_first_atom = topIn.Res(res).FirstAtom();
    if (debug_>0) mprintf("DEBUG: Residue %s\n", topIn.TruncResNameNum(res).c_str());
    if (resmap.SetupResidue(topIn, res) != 0) return 1;
    if (resmap.CheckBonds() != 0) return 1;
    resmap.DetermineAtomIDs();
    Iarray symmatoms;
    Iarray AtomStatus( resmap.Natom(), UNSELECTED );
    // Loop over all atoms in the residue
    for (int at = 0; at < resmap.Natom(); at++) {
      // If atom is unique in residue, mark non-symmetric 
      if (resmap[at].IsUnique())
        AtomStatus[at] = NONSYMM;
      else if (AtomStatus[at] != SYMM) {
        Iarray Selected( resmap.Natom(), 0 );
        symmatoms.clear();
        // Recursively search for other potentially symmetric atoms in residue.
        // The Selected array is used to keep track of which atoms have been
        // visited in this pass; this is used instead of AtomStatus so that
        // we can travel through atoms already marked as symmetric.
#ifdef  DEBUGSYMMRMSD
        recursionLevel_ = 0;
        mprintf("Starting recursive call for %i(%s)\n", at+1, resmap[at].c_str());
#       endif
        FindSymmetricAtoms(at, resmap, resmap[at].Unique(), Selected, symmatoms);
        if (symmatoms.size() == 1) {
          // Only 1 atom, not symmetric. Reset atom status
          AtomStatus[symmatoms.front()] = NONSYMM;
        } else if (symmatoms.size() > 1) {
          // Shift residue atom #s so they correspond with topology.
          for (Iarray::iterator it = symmatoms.begin(); it != symmatoms.end(); ++it) {
            AtomStatus[*it] = SYMM;
            *it += res_first_atom;
          }
          SymmetricAtomIndices_.push_back( symmatoms );
        }
      }
    }
    if (debug_ > 0) {
      mprintf("DEBUG:\tResidue Atom Status:\n");
      for (int at = 0; at < resmap.Natom(); at++) {
        mprintf("\t%s", topIn.AtomMaskName(at + res_first_atom).c_str());
        switch (AtomStatus[at]) {
          case NONSYMM: mprintf(" Non-symmetric\n"); break;
          case SYMM   : mprintf(" Symmetric\n"); break;
          case UNSELECTED: mprintf(" Unselected\n");
        }
      }
    }
  }
  if (debug_ > 0) {
    mprintf("DEBUG: Potential Symmetric Atom Groups:\n");
    for (AtomIndexArray::const_iterator symmatoms = SymmetricAtomIndices_.begin();
                                        symmatoms != SymmetricAtomIndices_.end();
                                        ++symmatoms)
    {
      mprintf("\t%8u) ", symmatoms - SymmetricAtomIndices_.begin());
      for (Iarray::const_iterator atom = symmatoms->begin();
                                  atom != symmatoms->end(); ++atom)
        mprintf(" %s", topIn.AtomMaskName(*atom).c_str());
      mprintf("\n");
    } 
  }
  return 0;
}

/** It is expected that TGT and REF already correspond to tgtMask. */
double SymmetricRmsdCalc::SymmRMSD(Frame const& TGT, Frame& REF) {
  Vec3 refTrans = REF.CenterOnOrigin( useMass_ );
  return SymmRMSD_CenteredRef(TGT, REF, REF, refTrans);
}

/** REF was already centered at the origin and TGT and REF already
  * correspond to tgtMask.
  */
double SymmetricRmsdCalc::SymmRMSD_CenteredRef( Frame const& TGT, Frame const& REF)
{
  return SymmRMSD_CenteredRef( TGT, REF, REF, Vec3(0.0) );
}

// SymmetricRmsdCalc::SymmRMSD()
double SymmetricRmsdCalc::SymmRMSD_CenteredRef(Frame const& TGT,
                                   Frame const& REF, Frame const& centeredREF,
                                   Vec3 const& refTrans)
{
  // Create initial 1 to 1 atom map for all atoms; indices in 
  // SymmetricAtomIndices will correspond to positions in AMap.
  for (int atom = 0; atom < (int)AMap_.size(); atom++)
    AMap_[atom] = atom;
  // Calculate initial best fit RMSD if necessary
  remapFrame_.SetCoordinates( TGT );
  if (fit_) {
    selectedTgt_.SetCoordinates(TGT, tgtMask_);
    selectedTgt_.RMSD_CenteredRef(centeredREF, rotMatrix_, tgtTrans_, useMass_);
    remapFrame_.Trans_Rot_Trans(tgtTrans_, rotMatrix_, refTrans);
  }
  // Correct RMSD for symmetry
  for (AtomIndexArray::const_iterator symmatoms = SymmetricAtomIndices_.begin();
                                      symmatoms != SymmetricAtomIndices_.end(); ++symmatoms)
  {
    // For each array of symmetric atoms, determine the lowest distance score
#   ifdef DEBUGSYMMRMSD
    mprintf("    Symmetric atoms group %u starting with atom %i\n", 
            symmatoms - SymmetricAtomIndices_.begin(), symmatoms->front() + 1);
#   endif
    cost_matrix_.Initialize( symmatoms->size() );
    for (Iarray::const_iterator tgtatom = symmatoms->begin();
                                tgtatom != symmatoms->end(); ++tgtatom)
    {
      for (Iarray::const_iterator refatom = symmatoms->begin();
                                  refatom != symmatoms->end(); ++refatom)
      {
        double dist2 = DIST2_NoImage( REF.XYZ(*refatom), remapFrame_.XYZ(*tgtatom) );
#       ifdef DEBUGSYMMRMSD
        mprintf("\t\t%i to %i: %f\n", *tgtatom + 1, *refatom + 1, dist2);
#       endif
        cost_matrix_.AddElement( dist2 );
      }
    }
    Iarray resMap = cost_matrix_.Optimize();
#   ifdef DEBUGSYMMRMSD
    mprintf("\tMapping from Hungarian Algorithm:\n");
    for (Iarray::const_iterator ha = resMap.begin(); ha != resMap.end(); ++ha)
      mprintf("\t\tMap col=%u row=%i\n", ha - resMap.begin(), *ha);
#   endif
    // Fill in overall map
    Iarray::const_iterator rmap = resMap.begin();
    for (Iarray::const_iterator atmidx = symmatoms->begin();
                                atmidx != symmatoms->end(); ++atmidx, ++rmap)
    {
      AMap_[*atmidx] = (*symmatoms)[*rmap]; // FIXME: Check indices
#     ifdef DEBUGSYMMRMSD
      mprintf("\tAssigned atom %i to atom %i\n", *atmidx + 1, (*symmatoms)[*rmap] + 1);
#     endif
    }
  }
# ifdef DEBUGSYMMRMSD
  mprintf("    Final Atom Mapping:\n");
  for (unsigned int ref = 0; ref < AMap_.size(); ++ref)
    mprintf("\t%u -> %i\n", ref + 1, AMap_[ref] + 1);
  mprintf("----------------------------------------\n");
# endif
  // Remap the original target frame, then calculate RMSD
  // TODO: Does the topology need to be remapped as well?
  remapFrame_.SetCoordinatesByMap(TGT, AMap_);
  selectedTgt_.SetCoordinates(remapFrame_, tgtMask_);
  double rmsdval;
  if (fit_)
    rmsdval = selectedTgt_.RMSD_CenteredRef( centeredREF, rotMatrix_, tgtTrans_, useMass_ );
  else
    rmsdval = selectedTgt_.RMSD_NoFit( centeredREF, useMass_ );
  return rmsdval;
}
