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
    double Etotal()    const {
      return (energy_.bond + energy_.angle + energy_.dihedral +
              energy_.vdw_14 + energy_.elec_14 + energy_.vdw +
              energy_.elec);
    }
    const double* EbondPtr()     const { return &(energy_.bond);     }
    const double* EanglePtr()    const { return &(energy_.angle);    }
    const double* EdihedralPtr() const { return &(energy_.dihedral); }
    const double* Evdw14Ptr()    const { return &(energy_.vdw_14);   }
    const double* Eelec14Ptr()   const { return &(energy_.elec_14);  }
    const double* EvdwPtr()      const { return &(energy_.vdw);      }
    const double* EelecPtr()     const { return &(energy_.elec);     }
  private:
    sander_input input_;
    pot_ene energy_;
    FileName top_filename_;
    std::vector<double> forces_;
    int top_pindex_;
};
#endif /* USE_SANDERLIB */
#endif
