#ifndef INC_PARM_LJ1264_PARAMS_H
#define INC_PARM_LJ1264_PARAMS_H
#include "ParmEnum.h"
#include <map>
#include <string>
namespace Cpptraj {
namespace Parm {
/// Used to assign LJ 12-6-4 C coefficients
class LJ1264_Params {
  public:
    LJ1264_Params();
    /// Init - mask, c4 file, water model, polarization file, tuning factor
    int Init_LJ1264(std::string const&, std::string const&, WaterModelType, std::string const&, double);
  private:
    typedef std::pair<std::string, double> NameMapPair;
    typedef std::map<std::string, double> NameMapType;

    /// Set C4 parameters for TIP3P
    void set_tip3p_params();
    void set_tip4pew_params();
    void set_spce_params();
    void set_opc3_params();
    void set_opc_params();
    void set_fb3_params();
    void set_fb4_params();
    /// Set C4 params for given water model
    void setupForWaterModel(WaterModelType);

    int read_2col_file(std::string const&, NameMapType&);
    int read_pol(std::string const&);
    int read_c4(std::string const&);

    WaterModelType waterModel_; ///< Water model that C3 parameters have been set for
    NameMapType c4params_; ///< C4 parameters
    NameMapType pol_; ///< Polarizabilities
    double tunfactor_; ///< Tuning factor
    std::string mask_; ///< Mask for selecting metal centers
};
}
}
#endif
