#ifndef INC_STRUCTURE_DISULFIDE_H
#define INC_STRUCTURE_DISULFIDE_H
#include <string>
class CpptrajFile;
class Frame;
class Topology;
namespace Cpptraj {
namespace Structure {
class ResStatArray;
/// Search for disulfide bonds
int SearchForDisulfides(ResStatArray&,
                        double, std::string const&, std::string const&, bool,
                        Topology&, Frame const&, std::string const&, CpptrajFile*);
}
}
#endif
