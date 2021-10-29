#ifndef INC_CHIRALITY_H
#define INC_CHIRALITY_H
class Topology;
class Frame;
namespace Cpptraj {
namespace Chirality {

enum ChiralType { ERR = 0, IS_S, IS_R };

/// \return Chirality at specified atom
ChiralType DetermineChirality(double&, int, Topology const&, Frame const&);

}
}
#endif
