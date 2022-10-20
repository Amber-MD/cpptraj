#ifndef INC_STRUCTURE_SUGARTOKEN_H
#define INC_STRUCTURE_SUGARTOKEN_H
#include <string>
// Forward declares
class ArgList;

namespace Cpptraj {
namespace Structure {
/// Hold general info for a specific sugar
class SugarToken {
  public:
    /// Base ring type
    enum RingTypeEnum { PYRANOSE = 0,  ///< Ring is 5 carbons, 1 oxygen
                        FURANOSE,      ///< Ring is 4 carbons, 1 oxygen
                        UNKNOWN_RING   ///< Some unknown ring type
                      };
    /// Sugar form
    enum FormTypeEnum { ALPHA = 0, BETA, UNKNOWN_FORM };
    /// Sugar chirality
    enum ChirTypeEnum { IS_D, IS_L, UNKNOWN_CHIR };

    /// CONSTRUCTOR
    SugarToken();
    /// CONSTRUCTOR - name, glycam code, form, chirality, ring type
    SugarToken(std::string const&, std::string const&, FormTypeEnum, ChirTypeEnum, RingTypeEnum);
    /// CONSTRUCTOR - ring type
    SugarToken(RingTypeEnum);
    /// /return <res>, set up from line: '<res> <code> <form> <chir> <ring> <name>'
    std::string SetFromLine(ArgList const&);
    /// \return string containing token info
    std::string InfoStr() const;

    std::string const& FullName()   const { return name_; }
    std::string const& GlycamCode() const { return glycamCode_; }
    FormTypeEnum Form()             const { return form_; }
    ChirTypeEnum Chirality()        const { return chir_; }
    RingTypeEnum RingType()         const { return ring_; }
    const char* formStr()           const { return formstr_[form_]; }
    const char* chirStr()           const { return chirstr_[chir_]; }
    const char* ringStr()           const { return ringstr_[ring_]; }
    static const char* formStr(FormTypeEnum f) { return formstr_[f]; }
    static const char* chirStr(ChirTypeEnum c) { return chirstr_[c]; }
    static const char* ringStr(RingTypeEnum r) { return ringstr_[r]; }

    void SetChirality(ChirTypeEnum c) { chir_ = c; }
    void SetForm(FormTypeEnum f)      { form_ = f; }
    void SetRingType(RingTypeEnum r)  { ring_ = r; }
  private:
    /// Keep synced with RingTypeEnum
    static const char* ringstr_[];
    /// Keep synced with FormTypeEnum
    static const char* formstr_[];
    /// Keep synced with ChirTypeEnum
    static const char* chirstr_[];


    std::string name_;       ///< Full sugar name
    //std::string resname_;    ///< PDB residue name
    std::string glycamCode_; ///< Glycam residue code
    FormTypeEnum form_;      ///< Sugar form
    ChirTypeEnum chir_;      ///< Sugar chirality
    RingTypeEnum ring_;      ///< Sugar ring type
};

}
}
#endif
