#ifndef INC_BOX_H
#define INC_BOX_H
/// Hold box information.
class Box {
  public:
    enum BoxType { NOBOX=0, ORTHO, TRUNCOCT, RHOMBIC, NONORTHO }; 

    Box();
    Box(const Box&);
    Box& operator=(const Box&);

    const char* TypeName(); 

    void SetBetaLengths(double,double,double,double);
    void SetAngles(const double*);
    void SetTruncOct();
    void SetNoBox();
    void SetMissingInfo(const Box&);

    double ToRecip(double*,double*);

    void SetX(double xin)     { box_[0] = xin; }
    void SetY(double yin)     { box_[1] = yin; }
    void SetZ(double zin)     { box_[2] = zin; }
    void SetAlpha(double ain) { box_[3] = ain; }
    void SetBeta(double bin)  { box_[4] = bin; }
    void SetGamma(double gin) { box_[5] = gin; }

    BoxType Type() const { return btype_;  }
    double BoxX()  { return box_[0]; }
    double BoxY()  { return box_[1]; }
    double BoxZ()  { return box_[2]; }
    double Alpha() { return box_[3]; }
    double Beta()  { return box_[4]; }
    double Gamma() { return box_[5]; }

  private:
    static const double TRUNCOCTBETA;
    static const char BoxNames[][15];
    int debug_;
    BoxType btype_;
    double box_[6];

    void SetBoxType();
};
#endif
