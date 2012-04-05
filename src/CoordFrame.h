#ifndef INC_COORDFRAME_H
#define INC_COORDFRAME_H
#include <vector>
#include "Vec3.h"
#include "AtomMask.h"
/// Test Frame replacement
class CoordFrame {
  public:
    CoordFrame();
    CoordFrame(int,const double*);
    CoordFrame(const CoordFrame&);
    CoordFrame &operator=(const CoordFrame&);

    double DIST2_NoImage(int, int);
    void Translate( Vec3 & );
    Vec3 GeometricCenter(AtomMask &Mask);

    inline bool empty() {
      return (coords_.empty());
    }
  private:
    std::vector<Vec3> coords_;
};
#endif
