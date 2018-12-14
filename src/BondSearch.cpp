#include <algorithm> // std::min,max
#include "BondSearch.h"
#include "DistRoutines.h"
#include "CpptrajStdio.h"
#include "PairList.h"
#include "StringRoutines.h" // ByteString
#ifdef TIMER
# include "Timer.h"
#endif

/** Create a bounding box around atoms in frameIn. */
Box CreateBoundingBox(Frame const& frameIn, Vec3& min)
{
  Box box;
  mprintf("\tCreating bounding box.\n");
  min = Vec3( frameIn.XYZ(0) );
  Vec3 max = min;
  for (int at = 1; at != frameIn.Natom(); at++)
  {
    const double* xyz = frameIn.XYZ(at);
    for (int i = 0; i < 3; i++) {
      min[i] = std::min( min[i], xyz[i] );
      max[i] = std::max( max[i], xyz[i] );
    }
  }
  // Make the offset 4 angstroms in each direction
  static const double boffset = 4.0;
  max += boffset;
  min -= boffset;
  Vec3 len = max - min;
  box.SetBetaLengths(90.0, len[0], len[1], len[2]);
  return box;
}

Box CreateBoundingBox(Frame const& frameIn) {
  Vec3 min;
  return CreateBoundingBox(frameIn, min);
}

/** Add bonds within residues to top using coords in frameIn. */
void BondsWithinResidues(Topology& top, Frame const& frameIn, double offset) {
  // ----- STEP 1: Determine bonds within residues
  for (Topology::res_iterator res = top.ResStart(); res != top.ResEnd(); ++res)
  {
    int stopatom = res->LastAtom();
    // Check for bonds between each atom in the residue.
    for (int atom1 = res->FirstAtom(); atom1 != stopatom; ++atom1) {
      Atom::AtomicElementType a1Elt = top[atom1].Element();
      // If this is a hydrogen and it already has a bond, move on.
      if (a1Elt==Atom::HYDROGEN && top[atom1].Nbonds() > 0 )
        continue;
      for (int atom2 = atom1 + 1; atom2 != stopatom; ++atom2) {
        Atom::AtomicElementType a2Elt = top[atom2].Element();
        double D2 = DIST2_NoImage(frameIn.XYZ(atom1), frameIn.XYZ(atom2) );
        double cutoff2 = Atom::GetBondLength(a1Elt, a2Elt) + offset;
        cutoff2 *= cutoff2;
        if (D2 < cutoff2) {
          top.AddBond(atom1, atom2);
          // Once a bond has been made to hydrogen move on.
          if (a1Elt==Atom::HYDROGEN) break;
        }
      }
    }
  }
}

/** Search for bonds within residues, then bonds between residues using a Grid. */
int BondSearch_Grid(Topology& top, Frame const& frameIn, double offset, int debug) {
  mprintf("\tDetermining bond info from distances using Grid.\n");
  if (frameIn.empty()) {
    mprinterr("Internal Error: No coordinates set; cannot search for bonds.\n");
    return 1;
  }
# ifdef TIMER
  Timer time_total, time_within, time_between, time_box;
  time_total.Start();
  time_box.Start();
# endif
  //Box box = frameIn.BoxCrd();
  Vec3 min(0.0);
  Box box = CreateBoundingBox( frameIn, min );
  box.PrintInfo();
  // Create grid indices.
  static const double spacing = 3.0;
  int nx = (int)(box.BoxX() / spacing);
  int ny = (int)(box.BoxY() / spacing);
  int nz = (int)(box.BoxZ() / spacing);
  int nynz = ny*nz;
  typedef std::vector<int> Iarray;
  typedef std::vector<Iarray> I2array;
  I2array GridIdx(nx * ny * nz);
  unsigned int nentries = 0;
  mprintf("\t%6s %6s %6s %6s\n", "X", "Y", "Z", "IDX");
  for (int x = 0; x < nx; x++) {
    int mx = std::max(0,  x-1);
    int px = std::min(nx, x+2);
    for (int y = 0; y < ny; y++) {
      int my = std::max(0,  y-1);
      int py = std::min(ny, y+2);
      int ynz = y*nz;
      for (int z = 0; z < nz; z++) {
        int idx = (x*(nynz))+(ynz)+z;
        int mz = std::max(0,  z-1);
        int pz = std::min(nz, z+2);
        mprintf("\t%6i %6i %6i %6i {", x, y, z, idx);
        // Who are my neighbors? TODO make faster by only looking forwards
        for (int xx = mx; xx < px; xx++) {
          for (int yy = my; yy < py; yy++) {
            int yynz = yy*nz;
            for (int zz = mz; zz < pz; zz++) {
              // NOTE: SELF NOT EXCLUDED HERE
              int ndx = (xx*(nynz))+(yynz)+zz;
              GridIdx[idx].push_back( ndx );
              nentries++;
              mprintf(" %5i", ndx);
            }
          }
        }
        mprintf(" }\n");
      } 
    }
  }
  // Approximate the memory used by the grid
  size_t totalgmem = (GridIdx.size()*sizeof(Iarray) + nentries*sizeof(int));
  mprintf("\tGrid index memory: %s\n", ByteString(totalgmem, BYTE_DECIMAL).c_str());
# ifdef TIMER
  time_box.Stop();
  time_within.Start();
# endif
  // ----- STEP 1: Determine bonds within residues.
  BondsWithinResidues(top, frameIn, offset);
# ifdef TIMER
  time_within.Stop();
  time_between.Start();
# endif
  // ----- STEP 2: Determine which atoms are eligible to form bonds
  //               with other residues.
  mprintf("DEBUG: Min= %8.3f %8.3f %8.3f\n", min[0], min[1], min[2]);
  I2array GridAtm(nx * ny * nz);
  for (int at = 0; at != top.Natom(); at++) {
    Atom const& atm = top[at];
    if (atm.Element() != Atom::HYDROGEN) {
      // TODO determine sane number of bonds based on element
      const double* XYZ = frameIn.XYZ( at );
      int ix = (int)((XYZ[0]-min[0]) / spacing); // TODO wrap
      int iy = (int)((XYZ[1]-min[1]) / spacing);
      int iz = (int)((XYZ[2]-min[2]) / spacing);
      if (ix < 0 || ix >= nx) mprintf("Warning: Atom %i X out of bounds\n", at+1);
      if (iy < 0 || iy >= ny) mprintf("Warning: Atom %i Y out of bounds\n", at+1);
      if (iz < 0 || iz >= nz) mprintf("Warning: Atom %i Z out of bounds\n", at+1);
      int idx = (ix*(nynz))+(iy*nz)+iz;
      //mprintf("GRID ATOM %6i {%8.3f %8.3f %8.3f} {%6i %6i %6i} %6i\n",
      //        at+1, XYZ[0], XYZ[1], XYZ[2], ix, iy, iz, idx);
      GridAtm[idx].push_back( at );
    }
  }
  mprintf("DEBUG: Populated grid bins:\n");
  for (I2array::const_iterator bin = GridAtm.begin(); bin != GridAtm.end(); ++bin)
  {
    if (!bin->empty()) {
      mprintf("\t%6u", bin-GridAtm.begin());
      for (Iarray::const_iterator atm = bin->begin(); atm != bin->end(); ++atm)
        mprintf(" %6i", *atm);
      mprintf("\n");
    }
  }
  // Loop over grid bins
  for (unsigned int ig = 0; ig != GridAtm.size(); ig++)
  {
    Iarray const& Bin1 = GridAtm[ig];
    if (!Bin1.empty()) {
      // Nbr contains indices for this and neighboring cells. 
      Iarray const& Nbr = GridIdx[ig];
      for (Iarray::const_iterator at1 = Bin1.begin(); at1 != Bin1.end(); ++at1)
      {
        Atom::AtomicElementType a1Elt = top[*at1].Element();
        for (Iarray::const_iterator nidx = Nbr.begin(); nidx != Nbr.end(); ++nidx)
        {
          Iarray const& Bin2 = GridAtm[*nidx];
          for (Iarray::const_iterator at2 = Bin2.begin(); at2 != Bin2.end(); ++at2)
          {
            if (*at1 != *at2) {
              Atom::AtomicElementType a2Elt = top[*at2].Element();
              double D2 = DIST2_NoImage(frameIn.XYZ(*at1), frameIn.XYZ(*at2) );
              double cutoff2 = Atom::GetBondLength(a1Elt, a2Elt) + offset;
              cutoff2 *= cutoff2;
              if (D2 < cutoff2)
                //top.AddBond(*at1, *at2);
                // DEBUG - add low, high
                top.AddBond(std::min(*at1, *at2), std::max(*at1, *at2));
            }
          } // END loop over Bin2 atoms
        } // END inner loop over neighbor+self bins
      } // END loop over Bin1 atoms
    } // END if Bin1 not empty
  } // END loop over all grid points

      
# ifdef TIMER
  time_between.Stop();
  time_total.Stop();
  time_box.WriteTiming(2, "Box creation", time_total.Total());
  time_within.WriteTiming(2, "Distances within residues", time_total.Total());
  time_between.WriteTiming(2, "Distances between residues", time_total.Total());
  time_total.WriteTiming(1, "Total for determining bonds via distances");
# endif

  return 0;
}

/** Search for bonds between atoms in residues and atoms in adjacent residues
  * using distance-based criterion that depends on atomic elements.
  * \param top Topology to add bonds to.
  * \param frameIn Frame containing atomic coordinates.
  * \param offset Offset to add when determining if a bond is present.
  * \param debug If > 0 print extra info.
  */
int BondSearch( Topology& top, Frame const& frameIn, double offset, int debug) {
  mprintf("\tDetermining bond info from distances.\n");
  if (frameIn.empty()) {
    mprinterr("Internal Error: No coordinates set; cannot search for bonds.\n");
    return 1;
  }
# ifdef TIMER
  Timer time_total, time_within, time_between;
  time_total.Start();
  time_within.Start();
# endif
  // ----- STEP 1: Determine bonds within residues
  BondsWithinResidues( top, frameIn, offset );
# ifdef TIMER
  time_within.Stop();
  time_between.Start();
# endif
  // ----- STEP 2: Determine bonds between adjacent residues
  Topology::mol_iterator nextmol = top.MolStart();
  if (top.Nmol() > 0)
    ++nextmol;
  for (Topology::res_iterator res = top.ResStart() + 1; res != top.ResEnd(); ++res)
  {
    // If molecule information is already present, check if first atom of 
    // this residue >= first atom of next molecule, which indicates this
    // residue and the previous residue are in different molecules.
    if ( (nextmol != top.MolEnd()) &&
         (res->FirstAtom() >= nextmol->BeginAtom()) )
    {
      ++nextmol;
      continue;
    }
    // If this residue is recognized as solvent, no need to check previous or
    // next residue
    if ( res->NameIsSolvent() ) {
      ++res;
      if (res == top.ResEnd()) break;
      continue;
    }
    // Get previous residue
    Topology::res_iterator previous_res = res - 1;
    // If previous residue is recognized as solvent, no need to check previous.
    if ( previous_res->NameIsSolvent() ) continue;
    // Get previous residue start atom
    int startatom = previous_res->FirstAtom();
    // Previous residue stop atom, this residue start atom
    int midatom = res->FirstAtom();
    // This residue stop atom
    int stopatom = res->LastAtom();
    // Check for bonds between adjacent residues
    for (int atom1 = startatom; atom1 != midatom; atom1++) {
      Atom::AtomicElementType a1Elt = top[atom1].Element();
      if (a1Elt==Atom::HYDROGEN) continue;
      for (int atom2 = midatom; atom2 != stopatom; atom2++) {
        Atom::AtomicElementType a2Elt = top[atom2].Element();
        if (a2Elt==Atom::HYDROGEN) continue;
        double D2 = DIST2_NoImage(frameIn.XYZ(atom1), frameIn.XYZ(atom2) );
        double cutoff2 = Atom::GetBondLength(a1Elt, a2Elt) + offset;
        cutoff2 *= cutoff2;
        if (D2 < cutoff2)
          top.AddBond(atom1, atom2);
      }
    }
  }
# ifdef TIMER
  time_between.Stop();
  time_total.Stop();
  time_within.WriteTiming(2, "Distances within residues", time_total.Total());
  time_between.WriteTiming(2, "Distances between residues", time_total.Total());
  time_total.WriteTiming(1, "Total for determining bonds via distances");
# endif
  if (debug > 0)
    mprintf("\t%s: %zu bonds to hydrogen, %zu other bonds.\n", top.c_str(),
            top.BondsH().size(), top.Bonds().size());
  return 0;
}

/** Bond search with pair list. */
int BondSearch_PL( Topology& top, Frame const& frameIn, double offset, int debug) {
  mprintf("\tDetermining bond info from distances using pair list.\n");
  if (frameIn.empty()) {
    mprinterr("Internal Error: No coordinates set; cannot search for bonds.\n");
    return 1;
  }
# ifdef TIMER
  Timer time_total;
  time_total.Start();
# endif
  // Pair list setup requires a box. Will need to create one if not present.
  Box box = frameIn.BoxCrd();
  if (box.Type() == Box::NOBOX)
   box = CreateBoundingBox(frameIn); 
  box.PrintInfo();

  // Initialize and set up pair list. TODO determine cutoff from spacing, maybe 3x3x3 ang voxels?
  static const double cutoff = 8.0;
  static const double skinnb = 0.1;
  PairList PL;
  PL.InitPairList( cutoff, skinnb, debug );
  PL.SetupPairList( box );
  Matrix_3x3 ucell, recip;
  box.ToRecip( ucell, recip );
  PL.CreatePairList( frameIn, ucell, recip, AtomMask(0, frameIn.Natom()) );

  // Determine bonds
  int cidx;
  for (cidx = 0; cidx < PL.NGridMax(); cidx++)
  {
    PairList::CellType const& thisCell = PL.Cell( cidx );
    if (thisCell.NatomsInGrid() > 0)
    {
      // cellList contains this cell index and all neighbors.
      PairList::Iarray const& cellList = thisCell.CellList();
      // transList contains index to translation for the neighbor.
      PairList::Iarray const& transList = thisCell.TransList();
      // Loop over all atoms of thisCell.
      for (PairList::CellType::const_iterator it0 = thisCell.begin();
                                              it0 != thisCell.end(); ++it0)
      {
        Vec3 const& xyz0 = it0->ImageCoords();
        Atom::AtomicElementType a0Elt = top[it0->Idx()].Element();
        // If hydrogen and already bonded, skip it.
        if (a0Elt == Atom::HYDROGEN && top[it0->Idx()].Nbonds() > 0) continue;
        // Check bonds to all other atoms in thisCell
        for (PairList::CellType::const_iterator it1 = it0 + 1;
                                                it1 != thisCell.end(); ++it1)
        {
          Vec3 const& xyz1 = it1->ImageCoords();
          Atom::AtomicElementType a1Elt = top[it1->Idx()].Element();
          // If hydrogen and already bonded, skip it.
          if (a1Elt == Atom::HYDROGEN && top[it1->Idx()].Nbonds() > 0) continue;
          Vec3 dxyz = xyz1 - xyz0;
          double rij2 = dxyz.Magnitude2();
          double cutoff2 = Atom::GetBondLength(a0Elt, a1Elt) + offset;
          cutoff2 *= cutoff2;
          if (rij2 < cutoff2) {
            //mprintf("DEBUG: IN CELL BOND: %s - %s\n",
            //  top.TruncResAtomNameNum(it0->Idx()).c_str(),
            //  top.TruncResAtomNameNum(it1->Idx()).c_str());
            //top.AddBond(it0->Idx(), it1->Idx());
            // DEBUG - Add low, high
            top.AddBond(std::min(it0->Idx(), it1->Idx()), std::max(it0->Idx(), it1->Idx()));
          }
        } // END loop over other atoms in thisCell
        // Loop over all neighbor cells
        for (unsigned int nidx = 1; nidx != cellList.size(); nidx++)
        {
          PairList::CellType const& nbrCell = PL.Cell( cellList[nidx] );
          // Translate vector for neighbor cell
          Vec3 const& tVec = PL.TransVec( transList[nidx] );
          //mprintf("\tNEIGHBOR %i (idxs %i - %i)\n", nbrCell, beg1, end1);
          // Loop over every atom in nbrCell
          for (PairList::CellType::const_iterator it1 = nbrCell.begin();
                                                  it1 != nbrCell.end(); ++it1)
          {
            Vec3 const& xyz1 = it1->ImageCoords();
            Atom::AtomicElementType a1Elt = top[it1->Idx()].Element();
            // If hydrogen and already bonded, skip it.
            if (a1Elt == Atom::HYDROGEN && top[it1->Idx()].Nbonds() > 0) continue;
            Vec3 dxyz = xyz1 + tVec - xyz0;
            double rij2 = dxyz.Magnitude2();
            double cutoff2 = Atom::GetBondLength(a0Elt, a1Elt) + offset;
            cutoff2 *= cutoff2;
            if (rij2 < cutoff2) {
              //mprintf("DEBUG: BETWEEN CELL BOND: %s - %s\n",
              //  top.TruncResAtomNameNum(it0->Idx()).c_str(),
              //  top.TruncResAtomNameNum(it1->Idx()).c_str());
              //top.AddBond(it0->Idx(), it1->Idx());
              // DEBUG - Add low, high
              top.AddBond(std::min(it0->Idx(), it1->Idx()), std::max(it0->Idx(), it1->Idx()));
            }
          } // END loop over atoms in neighbor cell
        } // END loop over neighbor cells
      } // END loop over atoms in this cell
    } // END if cell is populated
  } // END loop over cells
# ifdef TIMER
  time_total.Stop();
  time_total.WriteTiming(1, "Total for determining bonds via distances (pair list)");
# endif
  return 0;
}

// BondSearch()
int BondSearch(Topology& top, BondSearchType type, Frame const& frameIn, double offset, int debug)
{
  int err = 0;
  switch (type) {
    case SEARCH_REGULAR  : err = BondSearch(top, frameIn, offset, debug); break;
    case SEARCH_PAIRLIST : err = BondSearch_PL(top, frameIn, offset, debug); break;
    case SEARCH_GRID     : err = BondSearch_Grid(top, frameIn, offset, debug); break;
  }
  return err;
}
