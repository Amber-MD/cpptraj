#include "EnergyDecomp_Ewald.h"
#include "../ExclusionArray.h"
#include "../ParameterTypes.h"
#include <cmath>

using namespace Cpptraj::Energy;

EnergyDecomp_Ewald::EnergyDecomp_Ewald() {}

// -----------------------------------------------------------------------------

double EnergyDecomp_Ewald::ERFC(double rIn) {
# ifdef _OPENMP
  return EW_.ErfcEW( rIn );
# else
  t_erfc_.Start();
  double erfcval = EW_.ErfcEW( rIn );
  t_erfc_.Stop();
  return erfcval;
# endif
}

double EnergyDecomp_Ewald::adjust(double q0, double q1, double rij) {
  t_adjust_.Start();
  //double erfc = erfc_func(ew_coeff_ * rij);
  double erfc = ERFC(rij);
  double d0 = (erfc - 1.0) / rij;
  t_adjust_.Stop();
  return (q0 * q1 * d0);
}

void EnergyDecomp_Ewald::calcAdjust(double& e_adjust, double& Eljpme_correction_excl,
                                    PairList::AtmType const& atom0,
                                    PairList::AtmType const& atom1,
                                    double q0, double q1, double rij2)
{
  e_adjust += adjust(q0, q1, sqrt(rij2));
  if (EW_.LW_Coeff() > 0.0) {
    // LJ PME direct space exclusion correction
    // NOTE: Assuming excluded pair is within cutoff
    double kr2 = EW_.LW_Coeff() * EW_.LW_Coeff() * rij2;
    double kr4 = kr2 * kr2;
    //double kr6 = kr2 * kr4;
    double expterm = exp(-kr2);
    double r4 = rij2 * rij2;
    double r6 = rij2 * r4;
    double Cij = Cparam_[atom0.Idx()] * Cparam_[atom1.Idx()];
    Eljpme_correction_excl += (1.0 - (1.0 +  kr2 + kr4/2.0)*expterm) / r6 * Cij;
  }
}

/** Nonbonded energy */
void EnergyDecomp_Ewald::ene_nb(double& Eelec, double& Evdw, double& Eljpme_correction,
                                double rij2, double q0, double q1,
                                PairList::AtmType const& atom0,
                                PairList::AtmType const& atom1)
{
  double rij = sqrt( rij2 );
  double qiqj = q0 * q1;
  //double erfc = erfc_func(ew_coeff_ * rij);
  double erfc = ERFC(rij);
  double e_elec = qiqj * erfc / rij;
  Eelec += e_elec;
  //mprintf("EELEC %4i%4i%12.5f%12.5f%12.5f%3.0f%3.0f%3.0f\n",
  //int ta0, ta1;
  //if (it0->Idx() < it1->Idx()) {
  //  ta0=it0->Idx(); ta1=it1->Idx();
  //} else {
  //  ta1=it0->Idx(); ta0=it1->Idx();
  //}
  //mprintf("PELEC %6i%6i%12.5f%12.5f%12.5f\n", ta0, ta1, rij, erfc, e_elec);
  int nbindex = NB_->GetLJindex(TypeIndices_[atom0.Idx()],
                                TypeIndices_[atom1.Idx()]);
  if (nbindex > -1) {
    double vswitch = EW_.Switch_Fn(rij2);
    NonbondType const& LJ = NB_->NBarray()[ nbindex ];
    double r2    = 1.0 / rij2;
    double r6    = r2 * r2 * r2;
    double r12   = r6 * r6;
    double f12   = LJ.A() * r12;  // A/r^12
    double f6    = LJ.B() * r6;   // B/r^6
    double e_vdw = f12 - f6;      // (A/r^12)-(B/r^6)
    Evdw += (e_vdw * vswitch);
    //mprintf("PVDW %8i%8i%20.6f%20.6f\n", ta0+1, ta1+1, e_vdw, r2);
    if (EW_.LW_Coeff() > 0.0) {
      // LJ PME direct space correction
      double kr2 = EW_.LW_Coeff() * EW_.LW_Coeff() * rij2;
      double kr4 = kr2 * kr2;
      //double kr6 = kr2 * kr4;
      double expterm = exp(-kr2);
      double Cij = Cparam_[atom0.Idx()] * Cparam_[atom1.Idx()];
      Eljpme_correction += (1.0 - (1.0 +  kr2 + kr4/2.0)*expterm) * r6 * vswitch * Cij;
    }
  }
}

/** Direct space nonbonded energy calculation for Ewald */
void EnergyDecomp_Ewald::ene_ewald_direct(double& Eelec, double& evdw_out, double& eadjust_out,
                                          Frame const& frameIn, PairList const& PL,
                                          ExclusionArray const& Excluded)
{
  t_direct_.Start();
  Eelec = 0.0;
  double e_adjust = 0.0;
  double Evdw = 0.0;
  double Eljpme_correction = 0.0;
  double Eljpme_correction_excl = 0.0;
  int cidx;
# ifdef _OPENMP
# pragma omp parallel private(cidx) reduction(+: Eelec, Evdw, e_adjust, Eljpme_correction,Eljpme_correction_excl)
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
        Vec3 const& xyz0 = it0->ImageCoords();
        double q0 = Charge_[it0->Idx()];
#       ifdef DEBUG_PAIRLIST
        mprintf("DBG: Cell %6i (%6i atoms):\n", cidx, thisCell.NatomsInGrid());
#       endif
        // Exclusion list for this atom
        ExclusionArray::ExListType const& excluded = Excluded[it0->Idx()];
        // Calc interaction of atom to all other atoms in thisCell.
        for (PairList::CellType::const_iterator it1 = it0 + 1;
                                                it1 != thisCell.end(); ++it1)
        {
          Vec3 const& xyz1 = it1->ImageCoords();
          double q1 = Charge_[it1->Idx()];
          Vec3 dxyz = xyz1 - xyz0;
          double rij2 = dxyz.Magnitude2();
#         ifdef DEBUG_PAIRLIST
          mprintf("\tAtom %6i to atom %6i (%f)\n", it0->Idx()+1, it1->Idx()+1, sqrt(rij2));
#         endif
          // If atom excluded, calc adjustment, otherwise calc elec. energy.
          if (excluded.find( it1->Idx() ) == excluded.end())
          {

            if ( rij2 < EW_.Cut2() ) {
#             ifdef NBDBG
              if (it0->Idx() < it1->Idx())
                mprintf("NBDBG %6i%6i\n", it0->Idx()+1, it1->Idx()+1);
              else
                mprintf("NBDBG %6i%6i\n", it1->Idx()+1, it0->Idx()+1);
#             endif
              ene_nb(Eelec, Evdw, Eljpme_correction, rij2, q0, q1, *it0, *it1); // FIXME
            }
          } else {
            calcAdjust(e_adjust, Eljpme_correction_excl, *it0, *it1, q0, q1, rij2); //FIXME
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
            Vec3 const& xyz1 = it1->ImageCoords();
            double q1 = Charge_[it1->Idx()];
            Vec3 dxyz = xyz1 + tVec - xyz0;
            double rij2 = dxyz.Magnitude2();
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
              if ( rij2 < EW_.Cut2() ) {
#               ifdef NBDBG
                if (it0->Idx() < it1->Idx())
                  mprintf("NBDBG %6i%6i\n", it0->Idx()+1, it1->Idx()+1);
                else
                  mprintf("NBDBG %6i%6i\n", it1->Idx()+1, it0->Idx()+1);
#               endif
                ene_nb(Eelec, Evdw, Eljpme_correction, rij2, q0, q1, *it0, *it1); // FIXME
              }
            } else {
              calcAdjust(e_adjust, Eljpme_correction_excl, *it0, *it1, q0, q1, rij2); //FIXME
            }
          } // END loop over neighbor cell atoms
        } // END Loop over neighbor cells
      } // Loop over thisCell atoms
    } // END if thisCell is not empty
  } // Loop over cells
  t_direct_.Stop();
# ifdef DEBUG_PAIRLIST
  mprintf("DEBUG: Elec                             = %16.8f\n", Eelec);
  mprintf("DEBUG: Eadjust                          = %16.8f\n", e_adjust);
  mprintf("DEBUG: LJ vdw                           = %16.8f\n", Evdw);
  mprintf("DEBUG: LJ vdw PME correction            = %16.8f\n", Eljpme_correction);
  mprintf("DEBUG: LJ vdw PME correction (excluded) = %16.8f\n", Eljpme_correction_excl);
# endif
  evdw_out = Evdw + Eljpme_correction + Eljpme_correction_excl;
  eadjust_out = e_adjust;
} // END nb ewald


