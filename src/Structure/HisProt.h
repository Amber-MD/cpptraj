#ifndef INC_STRUCTURE_HISPROT_H
#define INC_STRUCTURE_HISPROT_H
#include <string>
class ArgList;
class Topology;
class NameType;
namespace Cpptraj {
namespace Structure {
/// For determining Histidine protonation
class HisProt {
  public:
    /// CONSTRUCTOR
    HisProt();
    /// \return keywords recognized by init
    static const char* keywords_;
    /// Initialize
    int InitHisProt(ArgList&, int);
    /// Print info to stdout
    void HisProtInfo() const;
    /// \return Epsilon-protonated name
    std::string const& EpsilonProtHisName() const { return hiename_; }
    /// \return Delta-protonated name
    std::string const& DeltaProtHisName() const { return hidname_; }
    /// \return Doubly-protonated name
    std::string const& DoubleProtHisName() const { return hipname_; }
    /// Determine histidine protonation from hydrogen atom placement
    int DetermineHisProt(Topology&) const;
  private:
    /// Try to determine protonation state of histidines from any hydrogens present.
    static int determineHisProt(Topology&,
                                NameType const&, NameType const&,
                                NameType const&, NameType const&,
                                NameType const&, NameType const&,
                                std::string const&);

    std::string nd1name_; ///< Delta nitrogen atom name
    std::string ne2name_; ///< Epsilon nitrogen atom name
    std::string hisname_; ///< Histidine original residue name
    std::string hiename_; ///< Epsilon-protonated residue name
    std::string hidname_; ///< Delta-protonated residue name
    std::string hipname_; ///< Doubly-protonated residue name
    std::string default_; ///< Default name if nothing detected
};
}
}
#endif
