#ifndef INC_ENERGY_RECIPCOORDS_H
#define INC_ENERGY_RECIPCOORDS_H
#include <vector>
class AtomMask;
class Frame;
namespace Cpptraj {
namespace Energy {
/// Hold selected XYZ coords for PME recip calcs with helPME
class RecipCoords {
  public:
    typedef std::vector<double> Darray;

    RecipCoords();

    int ReserveForSelected(AtomMask const&);
    void FillRecipCoords(Frame const&, AtomMask const&);

  private:
    Darray coordsD_;
};
}
}
#endif
