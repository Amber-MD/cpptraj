#ifndef INC_PAIRLISTTEMPLATE_H
#define INC_PAIRLISTTEMPLATE_H
#include "PairList.h"
#include "ExclusionArray.h"
namespace Cpptraj {
/// Template for doing pair list calculations
/** This template is designed to make it simpler to use the PairList class.
  * It takes a PairList, an excluded atom array (ExclusionArray), a cutoff
  * (squared), and an "Engine".
  * The Engine determines what is done when the cutoff is satisfied vs
  * when the interaction is excluded. It should also be a template and have
  * the following functions:
  * void FrameBeginCalc() : What to do at the beginning of each frame/calc.
  * void SetupAtom0(int idx0) : What to do for the outer loop atom.
  * void SetupAtom1(int idx1) : What to do for the inner loop atoms.
  * void CutoffSatisfied(double rij2, int idx0, idx1) : What to do when the
  *      cutoff squared is satisfied for atom indices 0 and 1.
  * void AtomPairExcluded(double rij2, int idx0, idx1) : What to do when the
  *      pair interaction between atom indices 0 and 1 is excluded.
  * 
  */
template <typename T, template <typename> class EngineClass>
void PairListTemplate(PairList const& PL, ExclusionArray const& Excluded, T cut2,
                      EngineClass<T>& engine)
{
  engine.FrameBeginCalc();

  int cidx;
# ifdef _OPENMP
# pragma omp parallel private(cidx) reduction(+: engine)
  {
# pragma omp for
# endif
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
        engine.SetupAtom0( *it0 );
        Vec3 const& xyz0 = it0->ImageCoords();
#       ifdef DEBUG_PAIRLIST
        mprintf("DBG: Cell %6i (%6i atoms):\n", cidx, thisCell.NatomsInGrid());
#       endif
        // Exclusion list for this atom
        ExclusionArray::ExListType const& excluded = Excluded[it0->Idx()];
        // Calc interaction of atom to all other atoms in thisCell.
        for (PairList::CellType::const_iterator it1 = it0 + 1;
                                                it1 != thisCell.end(); ++it1)
        {
          engine.SetupAtom1( *it1 );
          Vec3 const& xyz1 = it1->ImageCoords();
          Vec3 dxyz = xyz1 - xyz0;
          T rij2 = dxyz.Magnitude2();
#         ifdef DEBUG_PAIRLIST
          mprintf("\tAtom %6i to atom %6i (%f)\n", it0->Idx()+1, it1->Idx()+1, sqrt(rij2));
#         endif
          // If atom excluded, calc adjustment, otherwise calc elec. energy.
          if (excluded.find( it1->Idx() ) == excluded.end())
          {
            if ( rij2 < cut2 ) {
#             ifdef NBDBG
              if (it0->Idx() < it1->Idx())
                mprintf("NBDBG %6i%6i\n", it0->Idx()+1, it1->Idx()+1);
              else
                mprintf("NBDBG %6i%6i\n", it1->Idx()+1, it0->Idx()+1);
#             endif
              engine.CutoffSatisfied(rij2, *it0, *it1);
            }
          } else {
            engine.AtomPairExcluded(rij2, *it0, *it1);
          }
        } // END loop over other atoms in thisCell
        // Loop over all neighbor cells
        for (unsigned int nidx = 1; nidx != cellList.size(); nidx++)
        {
          PairList::CellType const& nbrCell = PL.Cell( cellList[nidx] );
#         ifdef DEBUG_PAIRLIST
          if (nbrCell.NatomsInGrid()>0) mprintf("\tto neighbor cell %6i\n", cellList[nidx]+1);
#         endif
          // Translate vector for neighbor cell
          Vec3 const& tVec = PL.TransVec( transList[nidx] );
#         ifdef DEBUG_PAIRLIST
          if (nbrCell.NatomsInGrid()>0) mprintf("DBG:\tto neighbor cell %6i (%6i atoms) tVec= %f %f %f\n", cellList[nidx], nbrCell.NatomsInGrid(), tVec[0], tVec[1], tVec[2]);
#         endif
          //mprintf("\tNEIGHBOR %i (idxs %i - %i)\n", nbrCell, beg1, end1);
          // Loop over every atom in nbrCell
          for (PairList::CellType::const_iterator it1 = nbrCell.begin();
                                                  it1 != nbrCell.end(); ++it1)
          {
            engine.SetupAtom1( *it1 );
            Vec3 const& xyz1 = it1->ImageCoords();
            Vec3 dxyz = xyz1 + tVec - xyz0;
            T rij2 = dxyz.Magnitude2();
            //mprintf("\t\tAtom %6i {%f %f %f} to atom %6i {%f %f %f} = %f Ang\n", it0->Idx()+1, xyz0[0], xyz0[1], xyz0[2], it1->Idx()+1, xyz1[0], xyz1[1], xyz1[2], sqrt(rij2));
#           ifdef DEBUG_PAIRLIST
            mprintf("\t\tAtom %6i to atom %6i (%f)\n", it0->Idx()+1, it1->Idx()+1, sqrt(rij2));
#           endif
            //mprintf("\t\tNbrAtom %06i\n",atnum1);
            // If atom excluded, calc adjustment, otherwise calc elec. energy.
            // TODO Is there better way of checking this?
            if (excluded.find( it1->Idx() ) == excluded.end())
            {

              //mprintf("\t\t\tdist= %f\n", sqrt(rij2));
              if ( rij2 < cut2 ) {
#               ifdef NBDBG
                if (it0->Idx() < it1->Idx())
                  mprintf("NBDBG %6i%6i\n", it0->Idx()+1, it1->Idx()+1);
                else
                  mprintf("NBDBG %6i%6i\n", it1->Idx()+1, it0->Idx()+1);
#               endif
                engine.CutoffSatisfied(rij2, *it0, *it1);
              }
            } else {
              engine.AtomPairExcluded(rij2, *it0, *it1);
            }
          } // END loop over neighbor cell atoms
        } // END Loop over neighbor cells
      } // Loop over thisCell atoms
    } // END if thisCell is not empty
  } // Loop over cells
# ifdef _OPENMP
  } // END pragma omp parallel
# endif
} // END PairListTemplate
} // END namespace Cpptraj
#endif
