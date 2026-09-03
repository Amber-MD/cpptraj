#ifndef INC_AMBERETERM_H
#define INC_AMBERETERM_H
#include <vector>
#include <string>
#include <map>
namespace Cpptraj {
class AmberEterm {
  public:
    typedef std::vector<double> Darray;

    enum FieldType { ETOT= 0, EPTOT, GMAX, BOND,
                     ANGLE, DIHED, VDWAALS, EEL, EGB, EPB, ECAVITY, EDISPER,
                     VDW14, EEL14, RESTRAINT, EAMBER, DENSITY,
                     RMS, EKTOT, ESURF, EAMD_BOOST, VOLUME, TEMP,
                     PRESS, DVDL, CMAP, N_FIELDTYPES };

    AmberEterm();
    Darray AllocEnergyArray();
    std::vector<bool> AllocExistsArray();
    int GetAmberEterms(const char*, Darray&, std::vector<bool>&) const;
  private:
    static const char* Enames_[];

    typedef std::map<std::string, unsigned int> NameIdxMap;
    typedef std::pair<std::string, unsigned int> NameIdxPair;

    FieldType getTermIdx(std::string const&) const;

    /// Map field names to indices into energy sets.
    NameIdxMap termIdxMap_;
};
}
#endif
