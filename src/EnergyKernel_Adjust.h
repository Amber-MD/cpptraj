              e_adjust += Adjust(q0, q1, sqrt(rij2));
/*
              // Electrostatic exclusion adjustment
#             ifndef _OPENMP
              t_adjust_.Start();
              t_erfc_.Start();
#             endif
              double rij = sqrt(rij2);
              //double erfc = erfc_func(ew_coeff_ * rij);
              double erfc = ERFC(ew_coeff_ * rij);
#             ifndef _OPENMP
              t_erfc_.Stop();
#             endif
              double d0 = (erfc - 1.0) / rij;
#             ifndef _OPENMP
              t_adjust_.Stop();
#             endif
              e_adjust += q0 * q1 * d0;
*/
#             ifdef CPPTRAJ_EKERNEL_LJPME
              // LJ PME direct space correction
              // NOTE: Assuming excluded pair is within cutoff
              double kr2 = lw_coeff_ * lw_coeff_ * rij2;
              double kr4 = kr2 * kr2;
              //double kr6 = kr2 * kr4;
              double expterm = exp(-kr2);
              double r4 = rij2 * rij2;
              double r6 = rij2 * r4;
              double Cij = Cparam_[it0->Idx()] * Cparam_[it1->Idx()];
              Eljpme_correction_excl += (1.0 - (1.0 +  kr2 + kr4/2.0)*expterm) / r6 * Cij;
#             endif
