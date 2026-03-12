#ifndef INC_STRUCTURE_MODXNA_INFO_H
#define INC_STRUCTURE_MODXNA_INFO_H
#include <string>
//namespace Cpptraj {
//namespace Structure {
/// Hold info used by the modXNA FF for joining fragments
class ModXNA_Info {
  public:
    /// CONSTRUCTOR
    ModXNA_Info();

    enum StatType { OK = 0, NOT_MODXNA, ERR };

    /// Parse the given modXNA title string
    StatType ParseModxnaStr(std::string const&);
    /// \return true if there is modXNA info
    bool HasModxna() const { return (!fragmentName_.empty()); }
    /// Print to stdout
    void PrintModxna() const;
    /// Generate ModXNA metadata string for e.g. writing ModXNA mol2 files
    std::string GenMetadataString() const;

    std::string const& FragmentName() const { return fragmentName_; }
    std::string const& Head()         const { return head_; }
    std::string const& HeadStrip()    const { return headStrip_; }
    std::string const& Tail()         const { return tail_; }
    std::string const& TailStrip()    const { return tailStrip_; }
    std::string const& Anchor()       const { return anchor_; }
    std::string const& AnchorStrip()  const { return anchorStrip_; }
  private:
    std::string fragmentName_;
    std::string head_;
    std::string headStrip_;
    std::string tail_;
    std::string tailStrip_;
    std::string anchor_;
    std::string anchorStrip_;
};
//}
//}
#endif
