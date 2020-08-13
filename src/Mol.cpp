#include "Mol.h"
#include "CpptrajStdio.h"
#include "Topology.h"

/** \return An array of unique molecule types out of those specified by molNums.
  */
Mol::Marray Mol::UniqueCount(Topology const& top, std::vector<int> const& molNums)
{
  Marray mols;
  // Hold begin residue for each molecule type
  Iarray Res0;
  // Hold end residue for each molecule type
  Iarray Res1;
  for (std::vector<int>::const_iterator mnum = molNums.begin();
                                        mnum != molNums.end(); ++mnum)
  {
    Molecule const& CurMol = top.Mol(*mnum);
    // Has this molecule been seen before?
    int natom = CurMol.NumAtoms();
    int res0 = top[ CurMol.BeginAtom() ].ResNum(); // CurMol first residue
    int res1 = top[ CurMol.EndAtom()-1 ].ResNum(); // CurMol last residue
    int nres = res1 - res0 + 1;
    int matchIdx = -1;
    for (int idx = 0; idx != (int)mols.size(); idx++)
    {
      // First check number of atoms, then number residues, then residue names.
      if ( mols[idx].natom_ == natom ) {
        if ( mols[idx].nres_ == nres ) {
          matchIdx = idx;
          int cridx = res0; // current molecule residue index
          for (int rridx = Res0[idx]; rridx <= Res1[idx]; rridx++, cridx++)
          {
            if ( top.Res(rridx).Name() != top.Res(cridx).Name() ) {
              // Residue name mismatch.
              matchIdx = -1;
              break;
            }
          } // END loop over all residues
        }
      }
      if (matchIdx != -1) break;
    } // END loop over found molecule types
    if (matchIdx == -1) {
      // New molecule
      mols.push_back( Type(*mnum, natom, nres, top.Res(res0).Name().Truncated()) );
      Res0.push_back( res0 );
      Res1.push_back( res1 );
    } else {
      // Existing molecule. Update count.
      mols[matchIdx].UpdateCount(*mnum);
    }
  } // END loop over all molecules
  return mols;
}

// Mol::UniqueCount()
Mol::Marray Mol::UniqueCount(Topology const& top) {
  if (top.Nmol() < 1) {
    mprintf("\t'%s' No molecule info.\n", top.c_str());
    return Marray();
  }
  AtomMask mask( "*" );
  if (top.SetupIntegerMask( mask )) return Marray();
  std::vector<int> molNums = top.MolnumsSelectedBy( mask );
  return UniqueCount(top, molNums);
}
