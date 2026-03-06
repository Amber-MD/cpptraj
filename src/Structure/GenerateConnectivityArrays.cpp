#include "GenerateConnectivityArrays.h"
#include "../ParameterTypes.h"
#include "../Atom.h"
#include "../Residue.h"
#include "../CpptrajStdio.h"
#include <algorithm> // std::swap

/// The direction in which atoms in residues should be scanned
static Cpptraj::Structure::AtomScanDirectionType CpptrajStructureAtomScanDirection_ = Cpptraj::Structure::SCAN_ATOMS_BACKWARDS;

/** Set default atom scan direction. */
void Cpptraj::Structure::SetAtomScanDirection( AtomScanDirectionType direction ) {
  if (direction == SCAN_ATOMS_BACKWARDS)
    mprintf("\tSetting atom scan direction to backwards.\n");
  else
    mprintf("\tSetting atom scan direction to forwards.\n");
  CpptrajStructureAtomScanDirection_ = direction;
}

/// Set start and end atoms along with offset based on atom scan direction
static inline void set_indices(int& start, int& end, int& offset, int firstatom, int lastatom)
{
  if (CpptrajStructureAtomScanDirection_ == Cpptraj::Structure::SCAN_ATOMS_BACKWARDS) {
    start = lastatom - 1;
    end = firstatom - 1;
    offset = -1;
  } else {
    start = firstatom;
    end = lastatom;
    offset = 1;
  }
}

/** Generate atom array in same order as LEaP. */
std::vector<int> Cpptraj::Structure::GenerateAtomArray(std::vector<Residue> const& residues,
                                                       std::vector<Atom> const& atoms)
{
  std::vector<int> out;
  for (std::vector<Residue>::const_iterator res = residues.begin(); res != residues.end(); ++res)
  {
    int start, end, offset;
    set_indices(start, end, offset, res->FirstAtom(), res->LastAtom());
    for (int iat = start; iat != end; iat += offset)
      out.push_back( iat );
  }
  return out;
}

/** From atom connectivity, generate a bond array in the same order as LEaP. */ // TODO use in GenerateBAT
BondArray Cpptraj::Structure::GenerateBondArray(std::vector<Residue> const& residues,
                                                std::vector<Atom> const& atoms)
{
  BondArray out;
  // BONDS
  //int bidx = 0; // DEBUG
  for (std::vector<Residue>::const_iterator res = residues.begin(); res != residues.end(); ++res)
  {
    int start, end, offset;
     // FIXME - this is a hack to get cpptraj to spit out solvent bonds the same way leap does. Leap lib/mol2 files get forwards direction.
    if (res->NameIsSolvent()) {
      start = res->FirstAtom();
      end = res->LastAtom();
      offset = 1;
      for (int iat = start; iat != end; iat += offset)
      {
        Atom const& At = atoms[iat];
        for (Atom::bond_iterator bat = At.bondbegin(); bat != At.bondend(); ++bat)
        {
          if (iat > *bat) {
            out.push_back( BondType(iat, *bat, -1) );
          }
        }
      }
    } else {
      set_indices(start, end, offset, res->FirstAtom(), res->LastAtom());
      for (int iat = start; iat != end; iat += offset)
      //for (int iat = res->LastAtom()-1; iat >= res->FirstAtom(); iat--)
      {
        Atom const& At = atoms[iat];
        for (Atom::bond_iterator bat = At.bondbegin(); bat != At.bondend(); ++bat)
        {
          if (iat < *bat) {
            //mprintf("DEBUG: BOND  i= %i  %i - %i (%i %i)\n",  bidx++, iat+1, *bat+1, iat*3, *bat*3);
            out.push_back( BondType(iat, *bat, -1) );
          }
          //else
          //  mprintf("DEBUG: X    i= %i  %i - %i (%i %i)\n",   bidx++, iat+1, *bat+1, iat*3, *bat*3);
        }
      }
    }
  }
  return out;
}

/** Generate a spanning tree around two atoms in the same manner as LEaP.
  * If at1 is -1, get the full tree. Otherwise, just get the tree in 
  * the direction of at0.
  */
std::vector<int> Cpptraj::Structure::GenerateSpanningTree(int at0, int at1, int targetDepth,
                                                          std::vector<Atom> const& atoms)
{
  std::vector<int> out;

  std::vector<bool> atomSeen(atoms.size(), false);

  std::vector<int> queue;
  std::vector<int> depth;

  out.push_back( at0 );
  queue.push_back( at0 );
  depth.push_back( 0 );
  atomSeen[at0] = true;
  unsigned int idx = 0;
  while (idx < queue.size()) {
    Atom const& currentAt = atoms[queue[idx]];
    for (Atom::bond_iterator bat = currentAt.bondbegin(); bat != currentAt.bondend(); ++bat)
    {
      if (*bat != at1 && !atomSeen[*bat]) {
        out.push_back( *bat );
        atomSeen[*bat] = true;
        if (targetDepth < 1) {
          if (atoms[*bat].Nbonds() > 1) {
            queue.push_back( *bat );
            depth.push_back( depth[idx]+1 );
          }
        } else {
          if (depth[idx] + 1 < targetDepth && atoms[*bat].Nbonds() > 1) {
            queue.push_back( *bat );
            depth.push_back( depth[idx]+1 );
          }
        }
      }
    }
    idx++;
  }

  return out;
}

/** From atom connectiviy, generate an angle array in the same order as LEaP. */ // TODO use in GenerateBAT
AngleArray Cpptraj::Structure::GenerateAngleArray(std::vector<Residue> const& residues,
                                                  std::vector<Atom> const& atoms)
{
  AngleArray out;
  // ANGLES TODO combine above
//  int aidx = 0;
  for (std::vector<Residue>::const_iterator res = residues.begin(); res != residues.end(); ++res)
  {
    int start, end, offset;
    set_indices(start, end, offset, res->FirstAtom(), res->LastAtom());
    for (int iat1 = start; iat1 != end; iat1 += offset)
    //for (int iat1 = res->LastAtom()-1; iat1 >= res->FirstAtom(); iat1--)
    {
      Atom const& At1 = atoms[iat1];
      for (int bidx1 = 0; bidx1 < At1.Nbonds(); bidx1++) {
        int iat2 = At1.Bond(bidx1);
        Atom const& At2 = atoms[iat2];
        for (int bidx2 = 0; bidx2 < At2.Nbonds(); bidx2++) {
          int iat3 = At2.Bond(bidx2);
          if (iat1 < iat3) {
            //mprintf("DEBUG: ANGLE  i= %i  %i - %i - %i (%i %i %i)\n", aidx++, iat1+1, iat2+1, iat3+1, iat1*3, iat2*3, iat3*3);
            out.push_back( AngleType(iat1, iat2, iat3, -1) );
          }
        }
      }
    }
  }
  return out;
}

/** From atom connectivity, generate a dihedral array in the same order as LEaP. */ // TODO use in GenerateBAT
DihedralArray Cpptraj::Structure::GenerateDihedralArray(std::vector<Residue> const& residues,
                                                        std::vector<Atom> const& atoms)
{
  DihedralArray out;
  // TORSIONS TODO combine above
//  int didx = 0;
  for (std::vector<Residue>::const_iterator res = residues.begin(); res != residues.end(); ++res)
  {
    int start, end, offset;
    set_indices(start, end, offset, res->FirstAtom(), res->LastAtom());
    for (int iat1 = start; iat1 != end; iat1 += offset)
    //for (int iat1 = res->LastAtom()-1; iat1 >= res->FirstAtom(); iat1--)
    {
      Atom const& At1 = atoms[iat1];
      for (int bidx1 = 0; bidx1 < At1.Nbonds(); bidx1++) {
        int iat2 = At1.Bond(bidx1);
        Atom const& At2 = atoms[iat2];
        for (int bidx2 = 0; bidx2 < At2.Nbonds(); bidx2++) {
          int iat3 = At2.Bond(bidx2);
          if (iat3 != iat1) {
            Atom const& At3 = atoms[iat3];
            for (int bidx3 = 0; bidx3 < At3.Nbonds(); bidx3++) {
              int iat4 = At3.Bond(bidx3);
              if (iat4 != iat2 && iat1 < iat4) {
                //mprintf("DEBUG: DIHEDRAL  i= %i  %i - %i - %i - %i (%i %i %i %i)\n", didx++, iat1+1, iat2+1, iat3+1, iat4+1, iat1*3, iat2*3, iat3*3, iat4*3);
                out.push_back( DihedralType( iat1, iat2, iat3, iat4, -1 ) );
              }
            }
          }
        }
      }
    }
  }
  return out;
}

/** Try to order an improper the same way that LEaP does.
  * LEaP has wild card names first, followed by atom types
  * in alphabetical order.
  */
static void order_improper_atoms(int* indices, std::vector<Atom> const& atoms)
{
  if (atoms[indices[0]].Type() > atoms[indices[1]].Type()) std::swap( indices[0], indices[1] );
  if (atoms[indices[1]].Type() > atoms[indices[2]].Type()) std::swap( indices[1], indices[2] );
  if (atoms[indices[0]].Type() > atoms[indices[1]].Type()) std::swap( indices[0], indices[1] );
  if (atoms[indices[1]].Type() > atoms[indices[2]].Type()) std::swap( indices[1], indices[2] );
}

// DEBUG
static inline void printName(Atom const& AJ) {
  mprintf(" :%i@%s (%s)", AJ.ResNum()+1, AJ.Name().Truncated().c_str(), AJ.Type().Truncated().c_str());
}

/** From atom connectivity, generate an improper array in the same order as LEaP.
  * No attempt is made to determine if this is an sp2 center; that is done
  * during parameterization.
  */
DihedralArray Cpptraj::Structure::GenerateImproperArray(std::vector<Residue> const& residues,
                                                        std::vector<Atom> const& atoms)
{
  DihedralArray out;
//  int iidx = 0;
  for (std::vector<Residue>::const_iterator res = residues.begin(); res != residues.end(); ++res)
  {
    int start, end, offset;
    set_indices(start, end, offset, res->FirstAtom(), res->LastAtom());
    for (int iat3 = start; iat3 != end; iat3 += offset)
    //for (int iat3 = res->LastAtom()-1; iat3 >= res->FirstAtom(); iat3--)
    {
      Atom const& AJ = atoms[iat3];
      if (AJ.Nbonds() >= 3) {
        for (int bidx0 = 0; bidx0 < AJ.Nbonds(); bidx0++) {
          for (int bidx1 = bidx0 + 1; bidx1 < AJ.Nbonds(); bidx1++) {
            for (int bidx2 = bidx1 + 1; bidx2 < AJ.Nbonds(); bidx2++) {
              int iat1 = AJ.BondIdxArray()[bidx0];
              int iat2 = AJ.BondIdxArray()[bidx1];
              int iat4 = AJ.BondIdxArray()[bidx2];
//              mprintf("DEBUG: IMPROPER  i= %i  %i - %i - %i - %i (%i %i %i %i)\n", iidx++, iat1+1, iat2+1, iat3+1, iat4+1, iat1*3, iat2*3, iat3*3, iat4*3);
              int indices[3];
              indices[0] = iat1;
              indices[1] = iat2;
              indices[2] = iat4;
              order_improper_atoms(indices, atoms);
              out.push_back( DihedralType(indices[0], indices[1], iat3, indices[2], DihedralType::BOTH) );
              // DEBUG
              //mprintf("DEBUG:\tOriginal order :");
              //printName(atoms[iat1]);
              //printName(atoms[iat2]);
              //printName(atoms[iat3]);
              //printName(atoms[iat4]);
              //mprintf("\nDEBUG:\tLeap order     :");
              //printName(atoms[indices[0]]);
              //printName(atoms[indices[1]]);
              //printName(atoms[iat3]);
              //printName(atoms[indices[2]]);
              //mprintf("\n");
            }
          }
        } // END outer loop over bond indices
      }
    } // END loop over residue atoms
  } // END loop over residues
  return out;
}

/// Given 2 bonded atoms at1 and at2, get all angles at1-at2-X
static inline void enumerateAngles(AngleArray& angles, int at1, int at2, std::vector<Atom> const& atoms)
{
  Atom const& A2 = atoms[at2];
  if (A2.Nbonds() > 1) {
    for (Atom::bond_iterator bat = A2.bondbegin(); bat != A2.bondend(); ++bat)
    {
      if (*bat != at1 && at1 < *bat) {
        mprintf("DEBUG: ANGLE  i= %zu  %i - %i - %i (%i %i %i)\n",  angles.size(), at1+1, at2+1, *bat+1, at1*3, at2*3, *bat*3);
        angles.push_back( AngleType(at1, at2, *bat, -1) );
      }
    }
  }
}

/// Given 2 bonded atoms at1 and at2, get all torsions X-at1-at2-Y
static inline void enumerateDihedrals(DihedralArray& dihedrals, int at1, int at2, std::vector<Atom> const& atoms)
{
  Atom const& A1 = atoms[at1];
  Atom const& A2 = atoms[at2];
  if (A1.Nbonds() > 1 && A2.Nbonds() > 1) {
    for (Atom::bond_iterator bat1 = A1.bondbegin(); bat1 != A1.bondend(); ++bat1)
    {
      if (*bat1 != at2) {
        for (Atom::bond_iterator bat2 = A2.bondbegin(); bat2 != A2.bondend(); ++bat2)
        {
          if (*bat2 != at1) {
            // LEaP convention appears to be first atom less than last atom
            if (*bat1 < *bat2)
              dihedrals.push_back( DihedralType(*bat1, at1, at2, *bat2, -1) );
            else
              dihedrals.push_back( DihedralType(*bat2, at2, at1, *bat1, -1) );
          }
        }
      }
    }
  }
}

/** Generate angle and torsion arrays from given bond array. */
void Cpptraj::Structure::GenerateAngleAndTorsionArraysFromBonds(AngleArray& angles,
                                                                DihedralArray& dihedrals,
                                                                std::vector<Atom> const& atoms,
                                                                BondArray const& allBonds)
{
  for (BondArray::const_iterator it = allBonds.begin(); it != allBonds.end(); ++it)
  {
    // Forward direction. A1-A2-X
    enumerateAngles( angles, it->A1(), it->A2(), atoms );
    // Reverse direction. A2-A1-X
    enumerateAngles( angles, it->A2(), it->A1(), atoms );
    // Dihedrals
    enumerateDihedrals( dihedrals, it->A1(), it->A2(), atoms );
  }
}

/** Merge two separate bond arrays. */
void Cpptraj::Structure::MergeBondArrays(BondArray& bondsOut,
                                         BondArray const& bonds,
                                         BondArray const& bondsh,
                                         std::vector<Atom> const& atoms)
{
  bondsOut.clear();
  bondsOut.reserve( bonds.size() + bondsh.size() );

  BondArray::const_iterator bx = bonds.begin();
  BondArray::const_iterator by = bondsh.begin();

  while (bx != bonds.end() && by != bondsh.end()) {
    // Which one goes next?
    Atom const& bx0 = atoms[bx->A1()];
    Atom const& by0 = atoms[by->A1()];
    if (bx0.ResNum() == by0.ResNum()) {
      if (bx->A1() == by->A1()) {
        // Same residue, same A1. Lower A2 goes first.
        if (by->A2() < bx->A2()) {
          bondsOut.push_back( *by );
          ++by;
        } else {
          bondsOut.push_back( *bx );
          ++bx;
        }
      } else {
        // Same residue. Higher A1 goes first. FIXME fix for scan direction forwards
        if (by->A1() > bx->A1()) {
          bondsOut.push_back( *by );
          ++by;
        } else {
          bondsOut.push_back( *bx );
          ++bx;
        }
      }
    } else {
      // Lower residue goes first.
      if (by0.ResNum() < bx0.ResNum()) {
        bondsOut.push_back( *by );
        ++by;
      } else {
        bondsOut.push_back( *bx );
        ++bx;
      }
    }
  } // END loop over both
  if (bx != bonds.end()) {
    for (; bx != bonds.end(); ++bx)
      bondsOut.push_back( *bx );
  }
  if (by != bondsh.end()) {
    for (; by != bondsh.end(); ++by)
      bondsOut.push_back( *by );
  }
}
