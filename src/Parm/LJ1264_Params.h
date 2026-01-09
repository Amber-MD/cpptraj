#ifndef INC_PARM_LJ1264_PARAMS_H
#define INC_PARM_LJ1264_PARAMS_H
#include "../AtomMask.h"
#include "ParmEnum.h"
#include <map>
#include <string>
class Topology;
namespace Cpptraj {
namespace Parm {
/// Used to assign LJ 12-6-4 C coefficients
/** Adds the LENNARD_JONES_CCOEF term for the ion 12-6-4 Lennard-Jones potential
  * term. If provided, the mask will allow you to specify which ions. If not
  * provided, :ZN will be used for the mask by default. The C4 parameter between
  * the ion and water can either be taken from the references for the requested
  * [watermodel] (TIP3P, TIP4PEW, SPCE, OPC3, OPC, FB3, and FB4) or provided in
  * the file specified by the c4file keyword. The polarizabilities must be
  * present in the the polfile file. The chemical symbol of the element will be
  * used to determine the atom type.
  * Parameters are expected in a file with 2 columns:
  *     <atom type>   <parameter>

  * All defaults come from Ref. [1], [2], [3], [4], [5] and [6]

  * [1] Pengfei Li and Kenneth M. Merz, J. Chem. Theory Comput., 2014, 10,
  *     289-297.
  * [2] Pengfei Li, Lin F. Song and Kenneth M. Merz, J. Phys. Chem. B, 2015,
  *     119, 883-895.
  * [3] Pengfei Li, Lin F. Song and Kenneth M. Merz, J. Chem. Theory Comput.,
  *     2015, 11, 1645-1657.
  * [4] Zhen Li, Lin Frank Song, Pengfei Li, and Kenneth M. Merz Jr. J. Chem.
  *     Theory Comput., 2020, 16, 4429-4442.
  * [5] Arkajyoti Sengupta, Zhen Li, Lin Frank Song, Pengfei Li, and Kenneth M.
  *     Merz Jr., J. Chem. Inf. Model., 2021, 61, 869-880.
  * [6] Zhen Li, Lin Frank Song, Pengfei Li, and Kenneth M. Merz Jr. J. Chem.
  *     Theory Comput., 2021, 17, 2342-2354.
  */
class LJ1264_Params {
  public:
    LJ1264_Params();
    /// Init - mask, c4 file, water model, polarization file, tuning factor, debug
    int Init_LJ1264(std::string const&, std::string const&, WaterModelType, std::string const&, double, int);
    /// Assign Topology with LJ 12-6-4 params
    int AssignLJ1264(Topology&);

    /// \return true if C4 parameters have been set up
    bool HasC4Params() const { return !c4params_.empty(); }
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

    static const double DEFAULT_WATER_POL_; ///< Default Polarizability of water

    WaterModelType waterModel_; ///< Water model that C4 parameters have been set for
    NameMapType c4params_; ///< C4 parameters
    NameMapType pol_; ///< Polarizabilities
    double tunfactor_; ///< Tuning factor
    double WATER_POL_; ///< Polarizability of water
    AtomMask mask_; ///< Mask for selecting metal centers
    int debug_;
};
}
}
#endif
