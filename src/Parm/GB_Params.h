#ifndef INC_PARM_GBPARAMS_H
#define INC_PARM_GBPARAMS_H
#include <string>
class ArgList;
class Topology;
namespace Cpptraj {
namespace Parm {
/// Known radii types. KEEP SYNCED WITH GB_RadiiTypeStr[] in GB_Params.cpp
enum GB_RadiiType { 
  BONDI=0,          ///< 0, bondi, Bondi radii
  BONDI_AMBER6,     ///< 1, amber6, Amber6 modified Bondi radii
  MBONDI,           ///< 2, mbondi, Modified bondi radii
  PB_AMBER,         ///< 3, pbamber, PB radii from Huo and Kollman, currently unused
  MBONDI2,          ///< 6, mbondi2, H(N)-modified Bondi radii
  PARSE,            ///< 7, parse, PARSE radii
  MBONDI3,          ///< 8, ArgH and AspGluO modified Bondi2 radii
  UNKNOWN_GB
};
/// \return Keyword corresponding to Amber radii flag
const char* GbTypeKey(GB_RadiiType);
/// \return string corresponding to Amber Radii Flag
std::string GbAmberFlag(GB_RadiiType);
/// \return string corresponding to gb radii set
std::string GbTypeStr(GB_RadiiType);
/// \return GB radii type corresponding to string
GB_RadiiType GbTypeFromKey(std::string const&);

/// Used to assign GB radii to a Topology
class GB_Params {
  public:
    /// CONSTRUCTOR
    GB_Params();
    /// \return recognized GB radii keywords.
    static std::string HelpText();
    /// Initialize 
    int Init_GB_Radii(ArgList&, GB_RadiiType);
    /// Initialize  no args
    int Init_GB_Radii(GB_RadiiType);
    // Print current setup
    void GB_Info() const;
    /// Assign GB radii
    int Assign_GB_Radii(Topology&) const;
  private:
    static void assign_gb_radii(Topology&, GB_RadiiType);
    static void assign_gb_screen(Topology&, GB_RadiiType);
    GB_RadiiType gbradii_; ///< Type of GB radii to assign
};
}
}
#endif 
