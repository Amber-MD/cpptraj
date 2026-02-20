#include "Merge.h"
#include "ParmHolder.h"
#include "../Atom.h"
#include "../CpptrajStdio.h"
#include "../ParameterTypes.h"
#include "../Residue.h"
#include "../StringRoutines.h" // integerToString

using namespace Cpptraj::Parm;

// ----- Bonds -------------------------
#ifdef CPPTRAJ_DEBUG_MERGE
static inline void printIdx(BondType const& bnd1, unsigned int atomOffset) {
  mprintf("DEBUG: Bond from top1 %i - %i will be %i - %i in top0\n",
          bnd1.A1()+1, bnd1.A2()+1, bnd1.A1()+1+atomOffset, bnd1.A2()+1+atomOffset);
}
#endif
static inline TypeNameHolder getTypes(bool& hasH, BondType const& bnd1, std::vector<Atom> const& atoms, BondParmArray const& parms)
{
  Atom const& A1 = atoms[bnd1.A1()];
  Atom const& A2 = atoms[bnd1.A2()];
  if (A1.Element() == Atom::HYDROGEN ||
      A2.Element() == Atom::HYDROGEN)
    hasH = true;
  else
    hasH = false; 
  TypeNameHolder types(2);
  types.AddName( A1.Type() );
  types.AddName( A2.Type() );
  return types;
}

static inline BondType idxWithOffset(BondType const& bnd1, int idx, unsigned int atomOffset) {
  return BondType(bnd1.A1()+atomOffset, bnd1.A2()+atomOffset, idx);
}

// ----- Angles ------------------------
#ifdef CPPTRAJ_DEBUG_MERGE
static inline void printIdx(AngleType const& ang1, unsigned int atomOffset) {
  mprintf("DEBUG: Angle from top1 %i - %i - %i will be %i - %i - %i in top0\n",
          ang1.A1()+1, ang1.A2()+1, ang1.A3()+1, ang1.A1()+1+atomOffset, ang1.A2()+1+atomOffset, ang1.A3()+1+atomOffset);
}
#endif
static inline TypeNameHolder getTypes(bool& hasH, AngleType const& ang1, std::vector<Atom> const& atoms, AngleParmArray const& parms)
{
  Atom const& A1 = atoms[ang1.A1()];
  Atom const& A2 = atoms[ang1.A2()];
  Atom const& A3 = atoms[ang1.A3()];
  if (A1.Element() == Atom::HYDROGEN ||
      A2.Element() == Atom::HYDROGEN ||
      A3.Element() == Atom::HYDROGEN)
    hasH = true;
  else
    hasH = false; 
  TypeNameHolder types(3);
  types.AddName( A1.Type() );
  types.AddName( A2.Type() );
  types.AddName( A3.Type() );
  return types;
}

static inline AngleType idxWithOffset(AngleType const& ang1, int idx, unsigned int atomOffset) {
  return AngleType(ang1.A1()+atomOffset, ang1.A2()+atomOffset, ang1.A3()+atomOffset, idx);
}

// ----- Dihedrals ---------------------
#ifdef CPPTRAJ_DEBUG_MERGE
static inline void printIdx(DihedralType const& dih1, unsigned int atomOffset) {
  //mprintf("DEBUG: Dihedral from top1 %i - %i - %i - %i will be %i - %i - %i - %i in top0\n",
  //        dih1.A1()+1, dih1.A2()+1, dih1.A3()+1, dih1.A4()+1,
  //        dih1.A1()+1+atomOffset, dih1.A2()+1+atomOffset, dih1.A3()+1+atomOffset, dih1.A4()+1+atomOffset);
  const char* desc;
  if (dih1.IsImproper())
    desc = "Improper";
  else
    desc = "Dihedral";
  mprintf("DEBUG: %s %i - %i - %i - %i ( %i %i %i %i )\n", desc,
          dih1.A1()+1, dih1.A2()+1, dih1.A3()+1, dih1.A4()+1,
          dih1.A1()*3, dih1.A2()*3, dih1.A3()*3, dih1.A4()*3);
}
#endif
/** Dihedrals are tricky because they can have multiple parameters for
  * the same 4 types. Have the 5th be a pseudo-type based on the
  * multiplicity (or improper status).
  */
static inline TypeNameHolder getTypes(bool& hasH, DihedralType const& dih1, std::vector<Atom> const& atoms, DihedralParmArray const& parms)
{
  Atom const& A1 = atoms[dih1.A1()];
  Atom const& A2 = atoms[dih1.A2()];
  Atom const& A3 = atoms[dih1.A3()];
  Atom const& A4 = atoms[dih1.A4()];
  if (A1.Element() == Atom::HYDROGEN ||
      A2.Element() == Atom::HYDROGEN ||
      A3.Element() == Atom::HYDROGEN ||
      A4.Element() == Atom::HYDROGEN)
    hasH = true;
  else
    hasH = false; 
  TypeNameHolder types(5);
  types.AddName( A1.Type() );
  types.AddName( A2.Type() );
  types.AddName( A3.Type() );
  types.AddName( A4.Type() );
  if (dih1.IsImproper()) {
    // Sanity check TODO is this needed?
    if (dih1.Idx() != -1 && (parms[dih1.Idx()].Pn() < 1 || parms[dih1.Idx()].Pn() < 1))
      mprinterr("Internal Error: Improper has a multiplicity != 1\n");
    types.AddName( "0" );
  } else {
    if (dih1.Idx() < 0)
      types.AddName( "1" );
    else
      types.AddName( integerToString( (int)parms[dih1.Idx()].Pn() ) );
  }
  return types;
}

static inline DihedralType idxWithOffset(DihedralType const& dih1, int idx, unsigned int atomOffset) {
  return DihedralType(dih1.A1()+atomOffset, dih1.A2()+atomOffset, dih1.A3()+atomOffset, dih1.A4()+atomOffset, dih1.Type(), idx);
}

// -------------------------------------
static inline void noParmWarning(TypeNameHolder const& types) {
  bool all_empty = true;
  for (unsigned int idx = 0; idx < types.Size(); idx++) {
    if (types[idx].len() > 0) { all_empty = false; break; }
  }
  if (all_empty) return;
  if (types.Size() == 2)
    mprintf("Warning: No bond parameters for types %s - %s\n", *(types[0]), *(types[1]));
  else if (types.Size() == 3)
    mprintf("Warning: No angle parameters for types %s - %s - %s\n", *(types[0]), *(types[1]), *(types[2]));
  else if (types.Size() == 5) {
    if (types[4] == "0")
      mprintf("Warning: No improper parameters for types %s - %s - %s - %s\n", *(types[0]), *(types[1]), *(types[2]), *(types[3]));
    else
      mprintf("Warning: No dihedral parameters for types %s - %s - %s - %s multiplicity %s\n", *(types[0]), *(types[1]), *(types[2]), *(types[3]), *(types[4]));
  } else if (types.Size() == 6)
    mprintf("Warning: No CMAP parameters for res %s atoms %s - %s - %s - %s - %s\n", *(types[0]), *(types[1]), *(types[2]), *(types[3]), *(types[4]), *(types[5]));
}

// -------------------------------------
#ifdef CPPTRAJ_DEBUG_MERGE
static inline void printTypes(TypeNameHolder const& types) {
  for (TypeNameHolder::const_iterator it = types.begin(); it != types.end(); ++it)
    mprintf(" %s", *(*it) );
  mprintf("\n");
}
#endif

// -------------------------------------
template <class IdxType, class ParmType, class IdxArray, class ParmArray>
class MergeTopArray
{
  public:
    MergeTopArray() : merge_with_existing_(false) {}
    /// Indicate parameters from top1 already present in top0 should not be added as new
    void SetMergeWithExisting(bool b) { merge_with_existing_ = b; }
  private:
    /// Append parameter for given types
    /** \return Index of newly-added parameter to p0 */
    int append_param(TypeNameHolder const& types,
                     Cpptraj::Parm::ParmHolder<int>& currentTypes0,
                     ParmArray& p0,
                     Cpptraj::Parm::ParmHolder<int> const& currentTypes1,
                     ParmArray const& p1)
    {
#     ifdef CPPTRAJ_DEBUG_MERGE
      printTypes( types ); // DEBUG
#     endif
      // Do we have an existing parameter in top0
      bool found;
      int idx = currentTypes0.FindParam(types, found);
      if (!found) {
        // No parameter yet.
        // Do we have an existing parameter in top1
#       ifdef CPPTRAJ_DEBUG_MERGE
        mprintf("DEBUG: Not found in top0. Looking in top1.\n");
#       endif
        idx = currentTypes1.FindParam(types, found);
        if (!found) {
#         ifdef CPPTRAJ_DEBUG_MERGE
          mprintf("DEBUG: No parameters.\n");
#         endif
          // No parameter in either top.
          idx = -1;
        } else {
#         ifdef CPPTRAJ_DEBUG_MERGE
          mprintf("DEBUG: Found in top1 index %i, adding to top0.\n", idx+1);
#         endif
          // Found a parameter in top1, add it to top0.
          int oldIdx = -1;
          if (merge_with_existing_) {
            // Check if parameter exists
            for (unsigned int i0 = 0; i0 != p0.size(); i0++) {
              if (p1[idx] == p0[i0]) {
                oldIdx = (int)i0;
                break;
              }
            }
          }
          if (oldIdx == -1) {
            // Do not merge with existing or does not yet exist.
            int newIdx = p0.size();
            p0.push_back( p1[idx] );
            idx = newIdx;
          } else {
#           ifdef CPPTRAJ_DEBUG_MERGE
            mprintf("DEBUG: Parm from top1 already present in top0 at position %i\n", oldIdx+1);
#           endif
            idx = oldIdx;
          }
        }
        // Add to existing parameters in top0.
        // Do this even if a parameter was not found so we dont keep looking.
        if (idx < 0) noParmWarning(types);
        currentTypes0.AddParm(types, idx, false);
      }
#     ifdef CPPTRAJ_DEBUG_MERGE
      mprintf("DEBUG: top0 parameter index is %i\n", idx+1);
#     endif
      return idx;
    }
    /// Append term1 to arrayX0/arrayY0 arrays along with parameters
    void append_term(IdxArray& arrayX0,
                     IdxArray& arrayY0,
                     ParmArray& p0,
                     unsigned int atomOffset,
                     IdxType const& term1,
                     Cpptraj::Parm::ParmHolder<int>& currentTypes0,
                     Cpptraj::Parm::ParmHolder<int> const& currentTypes1,
                     ParmArray const& p1,
                     std::vector<Atom> const& atoms1)
    {
#     ifdef CPPTRAJ_DEBUG_MERGE
      printIdx(term1, atomOffset);
#     endif
      bool hasH;
      TypeNameHolder types = getTypes(hasH, term1, atoms1, p1);
#     ifdef CPPTRAJ_DEBUG_MERGE
      mprintf("DEBUG: Looking for types in top0:");
#     endif
      int idx = append_param( types, currentTypes0, p0, currentTypes1, p1 );
      // At this point we have either found a parameter or not.
      if (hasH)
        arrayY0.push_back( idxWithOffset(term1, idx, atomOffset) );
      else
        arrayX0.push_back( idxWithOffset(term1, idx, atomOffset) );
    }
    /// Append term1 to arrayX0 array along with parameter
    void append_term(IdxArray& arrayX0,
                     ParmArray& p0,
                     unsigned int atomOffset,
                     IdxType const& term1,
                     Cpptraj::Parm::ParmHolder<int>& currentTypes0,
                     Cpptraj::Parm::ParmHolder<int> const& currentTypes1,
                     ParmArray const& p1,
                     std::vector<Atom> const& atoms1)
    {
#     ifdef CPPTRAJ_DEBUG_MERGE
      printIdx(term1, atomOffset);
#     endif
      bool hasH;
      TypeNameHolder types = getTypes(hasH, term1, atoms1, p1);
#     ifdef CPPTRAJ_DEBUG_MERGE
      mprintf("DEBUG: Looking for types in top0:");
#     endif
      int idx = append_param( types, currentTypes0, p0, currentTypes1, p1 );
      // At this point we have either found a parameter or not.
      arrayX0.push_back( idxWithOffset(term1, idx, atomOffset) );
    }

    /// Index existing term types in term arrays
    void index_term_types(Cpptraj::Parm::ParmHolder<int>& currentTypes,
                          IdxArray const& terms,
                          std::vector<Atom> const& atoms,
                          ParmArray const& parms)
    {
      bool hasH;
      for (typename IdxArray::const_iterator term = terms.begin(); term != terms.end(); ++term)
      {
        if (term->Idx() > -1) {
          TypeNameHolder types = getTypes(hasH, *term, atoms, parms);
          bool found;
          currentTypes.FindParam(types, found);
          if (!found) {
            currentTypes.AddParm(types, term->Idx(), false);
          }
        }
      }
    }

  public:
    /// Merge two pairs of term/parm arrays
    void MergeTermArrays(IdxArray& terms0,
                         IdxArray& termsh0,
                         ParmArray& p0,
                         std::vector<Atom> const& atoms0,
                         IdxArray const& terms1,
                         IdxArray const& termsh1,
                         ParmArray const& p1,
                         std::vector<Atom> const& atoms1)
    {
      // First index existing parameters
      Cpptraj::Parm::ParmHolder<int> currentTypes0, currentTypes1;
      index_term_types(currentTypes0, terms0, atoms0, p0);
      index_term_types(currentTypes0, termsh0, atoms0, p0);
      index_term_types(currentTypes1, terms1, atoms1, p1);
      index_term_types(currentTypes1, termsh1, atoms1, p1);
      // Loop over separate term arrays from top1 in the correct order
      unsigned int atomOffset = atoms0.size();
      typename IdxArray::const_iterator bx = terms1.begin();
      typename IdxArray::const_iterator by = termsh1.begin();
      while (bx != terms1.end() && by != termsh1.end()) {
        // Which one goes next?
        Atom const& bx0 = atoms1[bx->A1()];
        Atom const& by0 = atoms1[by->A1()];
        if (bx0.ResNum() == by0.ResNum()) {
          if (bx->A1() == by->A1()) {
            // Same residue, same A1. Lower A2 goes first.
            if (*by < *bx) {
              append_term( terms0, termsh0, p0, atomOffset, *by, currentTypes0, currentTypes1, p1, atoms1 );
              ++by;
            } else {
              append_term( terms0, termsh0, p0, atomOffset, *bx, currentTypes0, currentTypes1, p1, atoms1 );
              ++bx;
            }
          } else {
            // Both terms in same residue, different A1.
            // Higher A1 goes first.
            // FIXME fix for scan direction forwards.
            if (by->A1() > bx->A1()) {
              append_term( terms0, termsh0, p0, atomOffset, *by, currentTypes0, currentTypes1, p1, atoms1 );
              ++by;
            } else {
              append_term( terms0, termsh0, p0, atomOffset, *bx, currentTypes0, currentTypes1, p1, atoms1 );
              ++bx;
            }
          }
        } else {
          // Lower residue goes first.
          if (by0.ResNum() < bx0.ResNum()) {
            append_term( terms0, termsh0, p0, atomOffset, *by, currentTypes0, currentTypes1, p1, atoms1 );
            ++by;
          } else {
            append_term( terms0, termsh0, p0, atomOffset, *bx, currentTypes0, currentTypes1, p1, atoms1 );
            ++bx;
          }
        }
      } // END loop over both term arrays from top1
      if (bx != terms1.end()) {
        for (; bx != terms1.end(); ++bx)
          append_term( terms0, termsh0, p0, atomOffset, *bx, currentTypes0, currentTypes1, p1, atoms1 );
      }
      if (by != termsh1.end()) {
        for (; by != termsh1.end(); ++by)
          append_term( terms0, termsh0, p0, atomOffset, *by, currentTypes0, currentTypes1, p1, atoms1 );
      }
    } // END MergeTermArrays
    /// Merge two pairs of improper arrays
    /** Impropers have different ordering based first on the central (3rd) atom */
    void MergeImproperArrays(DihedralArray& terms0,
                             DihedralArray& termsh0,
                             DihedralParmArray& p0,
                             std::vector<Atom> const& atoms0,
                             DihedralArray const& terms1,
                             DihedralArray const& termsh1,
                             DihedralParmArray const& p1,
                             std::vector<Atom> const& atoms1)
    {
      // First index existing parameters
      Cpptraj::Parm::ParmHolder<int> currentTypes0, currentTypes1;
      index_term_types(currentTypes0, terms0, atoms0, p0);
      index_term_types(currentTypes0, termsh0, atoms0, p0);
      index_term_types(currentTypes1, terms1, atoms1, p1);
      index_term_types(currentTypes1, termsh1, atoms1, p1);
      // Loop over separate term arrays from top1 in the correct order
      unsigned int atomOffset = atoms0.size();
      typename IdxArray::const_iterator bx = terms1.begin();
      typename IdxArray::const_iterator by = termsh1.begin();
      while (bx != terms1.end() && by != termsh1.end()) {
        // Which one goes next?
        Atom const& bx0 = atoms1[bx->A1()];
        Atom const& by0 = atoms1[by->A1()];
        if (bx0.ResNum() == by0.ResNum()) {
          if (bx->A3() == by->A3()) {
            // Same residue, same A3. Lower A1 goes first.
            if (*by < *bx) {
              append_term( terms0, termsh0, p0, atomOffset, *by, currentTypes0, currentTypes1, p1, atoms1 );
              ++by;
            } else {
              append_term( terms0, termsh0, p0, atomOffset, *bx, currentTypes0, currentTypes1, p1, atoms1 );
              ++bx;
            }
          } else {
            // Both terms in same residue, different A3.
            // Higher A3 goes first.
            // FIXME fix for scan direction forwards.
            if (by->A3() > bx->A3()) {
              append_term( terms0, termsh0, p0, atomOffset, *by, currentTypes0, currentTypes1, p1, atoms1 );
              ++by;
            } else {
              append_term( terms0, termsh0, p0, atomOffset, *bx, currentTypes0, currentTypes1, p1, atoms1 );
              ++bx;
            }
          }
        } else {
          // Lower residue goes first.
          if (by0.ResNum() < bx0.ResNum()) {
            append_term( terms0, termsh0, p0, atomOffset, *by, currentTypes0, currentTypes1, p1, atoms1 );
            ++by;
          } else {
            append_term( terms0, termsh0, p0, atomOffset, *bx, currentTypes0, currentTypes1, p1, atoms1 );
            ++bx;
          }
        }
      } // END loop over both term arrays from top1
      if (bx != terms1.end()) {
        for (; bx != terms1.end(); ++bx)
          append_term( terms0, termsh0, p0, atomOffset, *bx, currentTypes0, currentTypes1, p1, atoms1 );
      }
      if (by != termsh1.end()) {
        for (; by != termsh1.end(); ++by)
          append_term( terms0, termsh0, p0, atomOffset, *by, currentTypes0, currentTypes1, p1, atoms1 );
      }
    } // END MergeImproperArrays
    /// Merge two term/parm arrays
    void MergeTermArray(IdxArray& terms0,
                        ParmArray& p0,
                        std::vector<Atom> const& atoms0,
                        IdxArray const& terms1,
                        ParmArray const& p1,
                        std::vector<Atom> const& atoms1)
    {
      // First index existing parameters
      Cpptraj::Parm::ParmHolder<int> currentTypes0, currentTypes1;
      index_term_types(currentTypes0, terms0, atoms0, p0);
      index_term_types(currentTypes1, terms1, atoms1, p1);
      unsigned int atomOffset = atoms0.size();
      for (typename IdxArray::const_iterator it = terms1.begin(); it != terms1.end(); ++it) {
        append_term( terms0, p0, atomOffset, *it, currentTypes0, currentTypes1, p1, atoms1 );
      }
    } // END MergeTermArray

  private:
    bool merge_with_existing_;

}; // END MergeTopArray template class

// -----------------------------------------------------------------------------

/** Given bond/bond parameter arrays from top0 and bond/bond parameter
  * arrays from top1, merge the bond arrays and consolidate the
  * parameters.
  */
void Cpptraj::Parm::MergeBondArrays(bool reduce_bond_params, BondArray& bonds0,
                                    BondArray& bondsh0,
                                    BondParmArray& bp0,
                                    AtArray const& atoms0,
                                    BondArray const& bonds1,
                                    BondArray const& bondsh1,
                                    BondParmArray const& bp1,
                                    AtArray const& atoms1)
{
  MergeTopArray<BondType, BondParmType, BondArray, BondParmArray> mergeBonds;
  mergeBonds.SetMergeWithExisting( reduce_bond_params );
  mergeBonds.MergeTermArrays( bonds0, bondsh0, bp0, atoms0,
                              bonds1, bondsh1, bp1, atoms1 );
}

// -----------------------------------------------------------------------------

/** Given angle/angle parameter arrays from top0 and angle/angle parameter
  * arrays from top1, merge the angle arrays and consolidate the
  * parameters.
  */
void Cpptraj::Parm::MergeAngleArrays(bool reduce_angle_params, AngleArray& angles0,
                                    AngleArray& anglesh0,
                                    AngleParmArray& ap0,
                                    AtArray const& atoms0,
                                    AngleArray const& angles1,
                                    AngleArray const& anglesh1,
                                    AngleParmArray const& ap1,
                                    AtArray const& atoms1)
{
  MergeTopArray<AngleType, AngleParmType, AngleArray, AngleParmArray> mergeAngles;
  mergeAngles.SetMergeWithExisting( reduce_angle_params );
  mergeAngles.MergeTermArrays( angles0, anglesh0, ap0, atoms0,
                               angles1, anglesh1, ap1, atoms1 );
}

// -----------------------------------------------------------------------------
static inline void separate_impropers(DihedralArray const& dihIn, DihedralArray& dihOut, DihedralArray& impOut)
{
  for (DihedralArray::const_iterator it = dihIn.begin(); it != dihIn.end(); ++it)
    if (it->IsImproper())
      impOut.push_back( *it );
    else
      dihOut.push_back( *it );
}

static inline void merge_impropers(DihedralArray& dihOut, DihedralArray& dihIn, DihedralArray& impIn) {
  dihOut.reserve( dihIn.size() + impIn.size() );
  dihOut = dihIn;
  dihIn.clear();
  for (DihedralArray::const_iterator it = impIn.begin(); it != impIn.end(); ++it)
    dihOut.push_back( *it  );
  impIn.clear();
}

/** Given dihedral/dihedral parameter arrays from top0 and dihedral/dihedral
  * parameter arrays from top1, merge the dihedral arrays and consolidate the
  * parameters.
  */
void Cpptraj::Parm::MergeDihedralArrays(DihedralArray& dihedrals0,
                                    DihedralArray& dihedralsh0,
                                    DihedralParmArray& dp0,
                                    AtArray const& atoms0,
                                    DihedralArray const& dihedrals1,
                                    DihedralArray const& dihedralsh1,
                                    DihedralParmArray const& dp1,
                                    AtArray const& atoms1)
{
  // Separate any impropers.
  DihedralArray dih0, dihH0, imp0, impH0;
  separate_impropers(dihedrals0, dih0, imp0);
  dihedrals0.clear();
  separate_impropers(dihedralsh0, dihH0, impH0);
  dihedralsh0.clear();

  DihedralArray dih1, dihH1, imp1, impH1;
  separate_impropers(dihedrals1, dih1, imp1);
  separate_impropers(dihedralsh1, dihH1, impH1);

  MergeTopArray<DihedralType, DihedralParmType, DihedralArray, DihedralParmArray> mergeDihedrals;
  mergeDihedrals.SetMergeWithExisting( true );
  mergeDihedrals.MergeTermArrays( dih0, dihH0, dp0, atoms0,
                                  dih1, dihH1, dp1, atoms1 );
  mergeDihedrals.MergeImproperArrays( imp0, impH0, dp0, atoms0,
                                      imp1, impH1, dp1, atoms1 );
  // Re-merge
  merge_impropers(dihedrals0, dih0, imp0);
  merge_impropers(dihedralsh0, dihH0, impH0);
}

// -----------------------------------------------------------------------------
// ----- CMAPs -------------------------
#ifdef CPPTRAJ_DEBUG_MERGE
static inline void printIdx(CmapType const& cmap1, unsigned int atomOffset) {
  mprintf("DEBUG: CMAP from top1 %i - %i - %i - %i - %i will be %i - %i - %i - %i - %i in top0\n",
          cmap1.A1()+1, cmap1.A2()+1, cmap1.A3()+1, cmap1.A4()+1, cmap1.A5()+1,
          cmap1.A1()+1+atomOffset, cmap1.A2()+1+atomOffset, cmap1.A3()+1+atomOffset, cmap1.A4()+1+atomOffset, cmap1.A5()+1+atomOffset);
}
#endif
/** CMAPs are indexed by residue name and 5 atom names, so need 6. 
  */
static inline TypeNameHolder getCmapTypes(CmapType const& cmap1, std::vector<Atom> const& atoms, std::vector<Residue> const& residues)
{
  Atom const& A1 = atoms[cmap1.A1()];
  Atom const& A2 = atoms[cmap1.A2()];
  Atom const& A3 = atoms[cmap1.A3()];
  Atom const& A4 = atoms[cmap1.A4()];
  Atom const& A5 = atoms[cmap1.A5()];

  TypeNameHolder types(6);
  types.AddName( residues[A2.ResNum()].Name() );
  types.AddName( A1.Name() );
  types.AddName( A2.Name() );
  types.AddName( A3.Name() );
  types.AddName( A4.Name() );
  types.AddName( A5.Name() );

  return types;
}

static inline CmapType idxWithOffset(CmapType const& cmap1, int idx, unsigned int atomOffset) {
  return CmapType(cmap1.A1()+atomOffset, cmap1.A2()+atomOffset, cmap1.A3()+atomOffset, cmap1.A4()+atomOffset, cmap1.A5()+atomOffset, idx);
}

/// Index existing CMAP types in term arrays
void index_cmap_types(Cpptraj::Parm::ParmHolder<int>& currentTypes,
                      CmapArray const& terms,
                      std::vector<Atom> const& atoms,
                      std::vector<Residue> const& residues)
{
  for (typename CmapArray::const_iterator term = terms.begin(); term != terms.end(); ++term)
  {
    if (term->Idx() > -1) {
      TypeNameHolder types = getCmapTypes(*term, atoms, residues);
      bool found;
      currentTypes.FindParam(types, found);
      if (!found) {
        currentTypes.AddParm(types, term->Idx(), false);
      }
    }
  }
}

/** Given CMAP arrays from top0 and top1, merge and consolidate parameters. */
void Cpptraj::Parm::MergeCmapArrays(CmapArray& cmap0,
                                    CmapGridArray& cg0,
                                    AtArray const& atoms0,
                                    ResArray const& residues0,
                                    CmapArray const& cmap1,
                                    CmapGridArray const& cg1,
                                    AtArray const& atoms1,
                                    ResArray const& residues1)
{
  unsigned int atomOffset = atoms0.size();
  unsigned int cgOffset = cg0.size();
  // First index existing parameters
  Cpptraj::Parm::ParmHolder<int> currentTypes0, currentTypes1;
  index_cmap_types(currentTypes0, cmap0, atoms0, residues0);
  index_cmap_types(currentTypes1, cmap1, atoms1, residues1);
  // Loop over incoming terms
  for (CmapArray::const_iterator c1 = cmap1.begin(); c1 != cmap1.end(); ++c1)
  {
#   ifdef CPPTRAJ_DEBUG_MERGE
    printIdx(*c1, atomOffset);
#   endif
    TypeNameHolder types = getCmapTypes(*c1, atoms1, residues1);
#   ifdef CPPTRAJ_DEBUG_MERGE
    mprintf("DEBUG: Looking for CMAP types in top0:");
    printTypes( types ); // DEBUG
#   endif
    // Do we have an existing parameter in top0
    bool found;
    int idx = currentTypes0.FindParam(types, found);
    if (!found) {
      // No parameter yet.
      // Do we have an existing parameter in top1
#     ifdef CPPTRAJ_DEBUG_MERGE
      mprintf("DEBUG: CMAP not found in top0. Looking in top1.\n");
#     endif
      idx = currentTypes1.FindParam(types, found);
      if (!found) {
#       ifdef CPPTRAJ_DEBUG_MERGE
        mprintf("DEBUG: No CMAP parameters.\n");
#       endif
        // No parameter in either top.
        idx = -1;
      } else {
#       ifdef CPPTRAJ_DEBUG_MERGE
        mprintf("DEBUG: Found in top1 index %i, adding to top0.\n", idx+1);
#       endif
        // Found a parameter in top1, add it to top0.
        int oldIdx = -1;
        //if (merge_with_existing_) {
          // Check if parameter exists
          for (unsigned int i0 = 0; i0 != cg0.size(); i0++) {
            if (cg1[idx].GridMatches( cg0[i0] )) {
              oldIdx = (int)i0;
              break;
            }
          }
        //}
        if (oldIdx == -1) {
          // Do not merge with existing or does not yet exist.
          // CMAPs are different than other parameters in that their order
          // in the topology matches the order in which they are specified
          // in the original parameter set. Therefore do not just append,
          // use the existing index as an offset.
          //int newIdx = cg0.size();
          //cg0.push_back( cg1[idx] );
          int newIdx = cgOffset + idx;
#         ifdef CPPTRAJ_DEBUG_MERGE
          mprintf("DEBUG: New CMAP index in top0 is %i\n", newIdx);
#         endif
          if ((unsigned int)newIdx >= cg0.size())
            cg0.resize( newIdx + 1 );
          cg0[newIdx] = cg1[idx];
          idx = newIdx;
        } else {
#         ifdef CPPTRAJ_DEBUG_MERGE
          mprintf("DEBUG: Parm from top1 already present in top0 at position %i\n", oldIdx+1);
#         endif
          idx = oldIdx;
        }
      }
      // Add to existing parameters in top0.
      // Do this even if a parameter was not found so we dont keep looking.
      if (idx < 0) noParmWarning(types);
      currentTypes0.AddParm(types, idx, false);
    }
#   ifdef CPPTRAJ_DEBUG_MERGE
    mprintf("DEBUG: top0 parameter index is %i\n", idx+1);
#   endif
    // At this point we have either found a parameter or not.
    cmap0.push_back( idxWithOffset(*c1, idx, atomOffset) );
  } // END loop over cmap terms from top1
}

// -----------------------------------------------------------------------------

/** Given bond/bond parameter array from top0 and bond/bond parameter
  * array from top1, merge the bond arrays and consolidate the
  * parameters.
  */
void Cpptraj::Parm::MergeBondArray(BondArray& bonds0,
                                   BondParmArray& bp0,
                                   AtArray const& atoms0,
                                   BondArray const& bonds1,
                                   BondParmArray const& bp1,
                                   AtArray const& atoms1)
{
# ifdef CPPTRAJ_DEBUG_MERGE
  mprintf("DEBUG: Enter MergeBondArray()\n");
# endif
  MergeTopArray<BondType, BondParmType, BondArray, BondParmArray> mergeBonds;
  mergeBonds.MergeTermArray( bonds0, bp0, atoms0,
                             bonds1, bp1, atoms1 );
}

// -----------------------------------------------------------------------------
/** Given improper/improper parameter array from top0 and improper/improper
  * parameter array from top1, merge the improper arrays and consolidate the
  * parameters.
  */
void Cpptraj::Parm::MergeImproperArray(DihedralArray& impropers0,
                                       DihedralParmArray& ip0,
                                       AtArray const& atoms0,
                                       DihedralArray const& impropers1,
                                       DihedralParmArray const& ip1,
                                       AtArray const& atoms1)
{
  MergeTopArray<DihedralType, DihedralParmType, DihedralArray, DihedralParmArray> mergeImpropers;
  mergeImpropers.SetMergeWithExisting( true ); // FIXME check this
  // NOTE could be MergeTermArrays as well
  mergeImpropers.MergeTermArray( impropers0, ip0, atoms0,
                                 impropers1, ip1, atoms1 );
}
