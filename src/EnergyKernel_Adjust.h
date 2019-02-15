              e_adjust += Adjust(q0, q1, sqrt(rij2));
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
