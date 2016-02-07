#ifndef INC_ESANDER_H
#define INC_ESANDER_H
#ifdef USE_SANDERLIB
#include "Topology.h"
#include "sander.h"
class Energy_Sander {
  public:
    Energy_Sander() : top_pindex_(-1) {}
    ~Energy_Sander();
    int Initialize(Topology const&, Frame&);
    int CalcEnergy(Frame&);
    double Ebond()     const { return energy_.bond;     }
    double Eangle()    const { return energy_.angle;    }
    double Edihedral() const { return energy_.dihedral; }
    double Evdw14()    const { return energy_.vdw_14;   }
    double Eelec14()   const { return energy_.elec_14;  }
    double Evdw()      const { return energy_.vdw;      }
    double Eelec()     const { return energy_.elec;     }
    const double *EbondPtr()  const { return &(energy_.bond);     }
  private:
    sander_input input_;
    pot_ene energy_;
    FileName top_filename_;
    std::vector<double> forces_;
    int top_pindex_;
};
#endif /* USE_SANDERLIB */
#endif
