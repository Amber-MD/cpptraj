#ifndef INC_BOX_H
#define INC_BOX_H
#include <vector>
class Box {
  public:
    enum BoxType { NOBOX=0, ORTHO, NONORTHO }; 

    Box();
    Box(const Box&);
    Box &operator=(const Box&);
    void PrintBoxType();
    void SetTruncOct();
    void SetNoBox();
    void SetBetaLengths(std::vector<double> &);
    void SetAngles(double*);
    int AmberIfbox();
    double ToRecip(double*,double*);
    inline BoxType Type() {
      return btype_;
    }
    inline void ToDouble(double *boxOut) {
      boxOut[0] = box_[0];
      boxOut[1] = box_[1];
      boxOut[2] = box_[2];
      boxOut[3] = box_[3];
      boxOut[4] = box_[4];
      boxOut[5] = box_[5];
    }

  private:
    static const double TRUNCOCTBETA;

    int debug_;
    double box_[6];
    BoxType btype_;

    void SetBoxType();
};
#endif
