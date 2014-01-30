#include "SymmetricRmsdCalc.h"
#include "DistRoutines.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
SymmetricRmsdCalc::SymmetricRmsdCalc() : debug_(0), fit_(true),
  useMass_(false), remap_(false) {}

// CONSTRUCTOR - For use when only RMSD is wanted.
SymmetricRmsdCalc::SymmetricRmsdCalc(AtomMask const& maskIn, bool fitIn, 
                                     bool useMassIn, Topology const& topIn, int debugIn) :
  debug_(debugIn), tgtMask_(maskIn), fit_(fitIn), useMass_(useMassIn), remap_(false)
{
  SetupSymmRMSD( topIn );
}

// SymmetricRmsdCalc::InitSymmRMSD()
int SymmetricRmsdCalc::InitSymmRMSD(std::string const& tMaskExpr, bool fitIn,
                                    bool useMassIn, bool remapIn, int debugIn)
{
  if (tgtMask_.SetMaskString( tMaskExpr )) return 1;
  debug_ = debugIn;
  fit_ = fitIn;
  useMass_ = useMassIn;
  remap_ = remapIn;
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
  * residue are considered.
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
  // Allocate space for remapped selected target atoms
  tgtRemap_ = selectedTgt_;
  if (remap_) {
    // Allocate space for remapped frame; same # atoms as original frame
    remapFrame_.SetupFrameV( topIn.Atoms(), topIn.HasVelInfo(), topIn.NrepDim() );
    targetMap_.resize( topIn.Natom() );
  }
  // Create map of original atom numbers to selected indices
  Iarray SelectedIdx( topIn.Natom(), -1 );
  int tgtIdx = 0;
  for (int originalAtom = 0; originalAtom != topIn.Natom(); ++originalAtom)
    if ( originalAtom == tgtMask_[tgtIdx] )
      SelectedIdx[originalAtom] = tgtIdx++;
  if (debug_ > 0) {
    mprintf("DEBUG: Original atom -> Selected Index mapping:\n");
    for (int originalAtom = 0; originalAtom != topIn.Natom(); ++originalAtom)
      mprintf("\t%8i -> %8i\n", originalAtom + 1, SelectedIdx[originalAtom] + 1);
  }
  // Create initial 1 to 1 atom map for all selected atoms; indices in 
  // SymmetricAtomIndices will correspond to positions in AMap.
  AMap_.resize( selectedTgt_.Natom() );
  // Determine last selected residue.
  int last_res = topIn[tgtMask_.back()].ResNum() + 1;
  mprintf("\tResidues up to %s will be considered for symmetry correction.\n",
          topIn.TruncResNameNum(last_res-1).c_str());
  // In each residue, determine which selected atoms are symmetric.
  SymmetricAtomIndices_.clear();
  AtomMap resmap;
  if (debug_ > 1) resmap.SetDebug(1);
  for (int res = 0; res < last_res; ++res) {
    int res_first_atom = topIn.Res(res).FirstAtom();
    // Are any of the residue atoms selected?
    bool atomsAreSelected = false;
    for (int ratom = res_first_atom; ratom != topIn.Res(res).LastAtom(); ++ratom)
      if (SelectedIdx[ratom] != -1) {
        atomsAreSelected = true;
        break;
      }
    if (!atomsAreSelected) continue;
    if (debug_>0) mprintf("DEBUG: Residue %s\n", topIn.TruncResNameNum(res).c_str());
    // Create AtomMap of this residue to determine chiral centers, unique atom IDs etc
    if (resmap.SetupResidue(topIn, res) != 0) return 1;
    if (resmap.CheckBonds() != 0) return 1;
    resmap.DetermineAtomIDs();
    Iarray symmAtoms; // Symmetric atoms, indices relative to resmap
    Iarray selectedSymmAtoms; // Selected Symmetric atoms, indices will be relative to tgtMask
    Iarray AtomStatus( resmap.Natom(), UNSELECTED );
    // Loop over all atoms in the residue
    for (int at = 0; at < resmap.Natom(); at++) {
      // If atom is unique in residue, mark non-symmetric 
      if (resmap[at].IsUnique())
        AtomStatus[at] = NONSYMM;
      else if (AtomStatus[at] != SYMM) {
        Iarray Selected( resmap.Natom(), 0 );
        symmAtoms.clear();
        // Recursively search for other potentially symmetric atoms in residue.
        // The Selected array is used to keep track of which atoms have been
        // visited in this pass; this is used instead of AtomStatus so that
        // we can travel through atoms already marked as symmetric.
#       ifdef DEBUGSYMMRMSD
        recursionLevel_ = 0;
        mprintf("Starting recursive call for %i(%s)\n", at+1, resmap[at].c_str());
#       endif
        FindSymmetricAtoms(at, resmap, resmap[at].Unique(), Selected, symmAtoms);
#       ifdef DEBUGSYMMRMSD
        mprintf("Potentially symmetric:\n");
        for (Iarray::const_iterator sa = symmAtoms.begin(); sa != symmAtoms.end(); ++sa)
          mprintf("\t%8i %4s %8i\n", *sa + res_first_atom + 1, 
                  topIn[*sa + res_first_atom].c_str(),
                  SelectedIdx[ *sa + res_first_atom ] + 1);
#       endif
        // Which of the symmetric atoms are selected?
        selectedSymmAtoms.clear();
        for (Iarray::const_iterator sa = symmAtoms.begin(); sa != symmAtoms.end(); ++sa)
          if (SelectedIdx[ *sa + res_first_atom ] != -1)
            selectedSymmAtoms.push_back( *sa );
        if (selectedSymmAtoms.size() == 1) {
          // Only 1 atom, not symmetric. Reset atom status
          AtomStatus[selectedSymmAtoms.front()] = NONSYMM;
        } else if (selectedSymmAtoms.size() > 1) {
          // Shift residue atom #s so they correspond with tgtMask.
          for (Iarray::iterator it = selectedSymmAtoms.begin(); 
                                it != selectedSymmAtoms.end(); ++it)
          {
            AtomStatus[*it] = SYMM;
            *it = SelectedIdx[ *it + res_first_atom ];
          }
          SymmetricAtomIndices_.push_back( selectedSymmAtoms );
        }
      }
    }
    // If remapping and not all atoms in a residue are selected, warn user.
    if (remap_) {
      for (int at = 0; at < resmap.Natom(); at++) {
        if (AtomStatus[at] == UNSELECTED) {
          mprintf("Warning: Not all atoms selected in residue '%s'. Re-mapped\n"
                  "Warning:   structures may appear distorted.\n", 
                  topIn.TruncResNameNum(res).c_str());
          break;
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
        mprintf(" %s(%i)", topIn.AtomMaskName(tgtMask_[*atom]).c_str(), tgtMask_[*atom] + 1);
      mprintf("\n");
    } 
  }
  return 0;
}

/** It is expected that TGT and REF already correspond to tgtMask. */
double SymmetricRmsdCalc::SymmRMSD(Frame const& TGT, Frame& REF) {
  REF.CenterOnOrigin( useMass_ );
  return SymmRMSD_CenteredRef(TGT, REF);
}

/** Calculate symmetric RMSD, potentially with coordinate remapping.*/
double SymmetricRmsdCalc::SymmRMSD_TGT( Frame const& TGT, Frame const& centeredREF)
{
  selectedTgt_.SetCoordinates( TGT, tgtMask_ );
  double rmsdval = SymmRMSD_CenteredRef( selectedTgt_, centeredREF);
  if (remap_) {
    // Now re-map the target frame
    for (int atom = 0; atom < (int)targetMap_.size(); atom++)
      targetMap_[atom] = atom;
    for (unsigned int ref = 0; ref < AMap_.size(); ++ref)
      targetMap_[ tgtMask_[ref] ] = tgtMask_[AMap_[ref]];
    remapFrame_.SetCoordinatesByMap( TGT, targetMap_ );
  }
  return rmsdval;
}

// SymmetricRmsdCalc::SymmRMSD()
/** selectedTgt and centeredREF must correspond to each other. */
double SymmetricRmsdCalc::SymmRMSD_CenteredRef(Frame const& selectedTgt, Frame const& centeredREF)
{
  // Create initial 1 to 1 atom map for all atoms; indices in 
  // SymmetricAtomIndices will correspond to positions in AMap.
  for (int atom = 0; atom < (int)AMap_.size(); atom++)
    AMap_[atom] = atom;
  tgtRemap_.SetCoordinates(selectedTgt);
  // Calculate initial best fit RMSD if necessary
  if (fit_) {
    tgtRemap_.RMSD_CenteredRef(centeredREF, rotMatrix_, tgtTrans_, useMass_);
    // Since tgtRemap is moved to origin during RMSD calc and centeredREF
    // should already be at the origin, just rotate.
    tgtRemap_.Rotate( rotMatrix_ );
  }
  // Correct RMSD for symmetry
  for (AtomIndexArray::const_iterator symmatoms = SymmetricAtomIndices_.begin();
                                      symmatoms != SymmetricAtomIndices_.end(); ++symmatoms)
  {
    // For each array of symmetric atoms, determine the lowest distance score
#   ifdef DEBUGSYMMRMSD
    mprintf("    Symmetric atoms group %u starting with atom %i\n", 
            symmatoms - SymmetricAtomIndices_.begin(), tgtMask_[symmatoms->front()] + 1);
#   endif
    cost_matrix_.Initialize( symmatoms->size() );
    for (Iarray::const_iterator ta = symmatoms->begin(); ta != symmatoms->end(); ++ta)
    {
      for (Iarray::const_iterator ra = symmatoms->begin(); ra != symmatoms->end(); ++ra)
      { 
        double dist2 = DIST2_NoImage( centeredREF.XYZ(*ra), tgtRemap_.XYZ(*ta) );
#       ifdef DEBUGSYMMRMSD
        mprintf("\t\t%i to %i: %f\n", tgtMask_[*ta] + 1, tgtMask_[*ra] + 1, dist2);
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
      AMap_[*atmidx] = (*symmatoms)[*rmap];
#     ifdef DEBUGSYMMRMSD
      mprintf("\tAssigned atom %i to atom %i\n", tgtMask_[*atmidx] + 1,
              tgtMask_[(*symmatoms)[*rmap]] + 1);
#     endif
    }
  }
# ifdef DEBUGSYMMRMSD
  mprintf("    Final Atom Mapping:\n");
  for (unsigned int ref = 0; ref < AMap_.size(); ++ref)
    mprintf("\t%u -> %i\n", tgtMask_[ref] + 1, tgtMask_[AMap_[ref]] + 1);
  mprintf("----------------------------------------\n");
# endif
  // Remap the target frame for symmetry, then calculate new RMSD.
  // TODO: Does the topology need to be remapped as well?
  double rmsdval;
  tgtRemap_.SetCoordinatesByMap(selectedTgt, AMap_);
  if (fit_)
    rmsdval = tgtRemap_.RMSD_CenteredRef( centeredREF, rotMatrix_, tgtTrans_, useMass_ );
  else
    rmsdval = tgtRemap_.RMSD_NoFit( centeredREF, useMass_ );
  return rmsdval;
}
