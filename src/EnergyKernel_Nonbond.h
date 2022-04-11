                double rij = sqrt( rij2 );
                double qiqj = q0 * q1;
#               ifndef _OPENMP
                t_erfc_.Start();
#               endif
                //double erfc = erfc_func(ew_coeff_ * rij);
                double erfc = ERFC(ew_coeff_ * rij);
#               ifndef _OPENMP
                t_erfc_.Stop();
#               endif
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
                int nbindex = NB_->GetLJindex(TypeIndices_[it0->Idx()],
                                              TypeIndices_[it1->Idx()]);
                if (nbindex > -1) {
                  double vswitch = switch_fn(rij2, cut2_0_, cut2_);
                  NonbondType const& LJ = NB_->NBarray()[ nbindex ];
                  double r2    = 1.0 / rij2;
                  double r6    = r2 * r2 * r2;
                  double r12   = r6 * r6;
                  double f12   = LJ.A() * r12;  // A/r^12
                  double f6    = LJ.B() * r6;   // B/r^6
                  double e_vdw = f12 - f6;      // (A/r^12)-(B/r^6)
                  Evdw += (e_vdw * vswitch);
                  //mprintf("PVDW %8i%8i%20.6f%20.6f\n", ta0+1, ta1+1, e_vdw, r2);
#                 ifdef CPPTRAJ_EKERNEL_LJPME
                  // LJ PME direct space correction
                  double kr2 = lw_coeff_ * lw_coeff_ * rij2;
                  double kr4 = kr2 * kr2;
                  //double kr6 = kr2 * kr4;
                  double expterm = exp(-kr2);
                  double Cij = Cparam_[it0->Idx()] * Cparam_[it1->Idx()];
                  Eljpme_correction += (1.0 - (1.0 +  kr2 + kr4/2.0)*expterm) * r6 * vswitch * Cij;
#                 endif
                }
