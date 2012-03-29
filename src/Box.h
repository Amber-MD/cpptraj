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
    void SetBetaLengths(std::vector<double> &);
    int AmberIfbox();
    double ToRecip(double*,double*);

  private:
    static const double TRUNCOCTBETA;

    int debug_;
    double box_[6];
    BoxType btype_;

    void SetBoxType();
};
#endif
