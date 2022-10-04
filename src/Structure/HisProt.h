#ifndef INC_STRUCTURE_HISPROT_H
#define INC_STRUCTURE_HISPROT_H
class Topology;
class NameType;
namespace Cpptraj {
namespace Structure {
/// Try to determine protonation state of histidines from any hydrogens present.
int DetermineHisProt(Topology&,
                     NameType const&, NameType const&,
                     NameType const&, NameType const&, NameType const&, NameType const&);
}
}
#endif
