#include "Mol.h"
#include "CpptrajStdio.h"
#include "Topology.h"

/** \return An array of unique molecule types out of those specified by molNums.
  */
Mol::Marray Mol::UniqueCount(Topology const& top, std::vector<int> const& molNums)
{
  Marray mols;
  // Loop over all molecule numbers
  for (std::vector<int>::const_iterator mnum = molNums.begin();
                                        mnum != molNums.end(); ++mnum)
  {
    Molecule const& CurMol = top.Mol(*mnum);
    // Has this molecule been seen before? Check # atoms and residues.
    int natom = CurMol.NumAtoms();
    int nres  = top.NresInMol(*mnum);
    std::string curName = top.Res(top[ CurMol.MolUnit().Front() ].ResNum()).Name().Truncated();
    int matchIdx = -1;
    for (int idx = 0; idx != (int)mols.size(); idx++)
    {
      // First check number of atoms, then number residues.
      if ( mols[idx].natom_ == natom ) {
        if ( mols[idx].nres_ == nres ) {
          Molecule const& PrevMol = top.Mol( mols[idx].idxs_.front() );
          // Check that # segments match
          if ( CurMol.MolUnit().nSegments() == PrevMol.MolUnit().nSegments() ) {
            matchIdx = idx;
            // Now check that residue names in CurMol match those in PrevMol
            for (unsigned int iseg = 0; iseg != CurMol.MolUnit().nSegments(); iseg++)
            {
              int curRes0  = top[ CurMol.MolUnit()[iseg].Begin()  ].ResNum();
              int curRes1  = top[ CurMol.MolUnit()[iseg].End()-1  ].ResNum();
              int prevRes0 = top[ PrevMol.MolUnit()[iseg].Begin() ].ResNum();
              int prevRes1 = top[ PrevMol.MolUnit()[iseg].End()-1 ].ResNum();
              if (curRes1 - curRes0 != prevRes1 - prevRes0) {
                matchIdx = -1;
                break;
              }
              while (curRes0 <= curRes1) {
                if (top.Res(curRes0).Name() != top.Res(prevRes0).Name()) {
                  // Residue name mismatch.
                  matchIdx = -1;
                  break;
                }
                curRes0++;
                prevRes0++;
              }
              if (matchIdx == -1) break;
            } // END loop over segments
          } // END # segments match 
        } // END # residues match
      } // END # atoms match
      if (matchIdx != -1) break;
    } // END loop over found molecule types
    if (matchIdx == -1) {
      // New molecule
      mols.push_back( Type(*mnum, natom, nres, curName) );
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
