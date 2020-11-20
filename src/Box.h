#ifndef INC_BOX_H
#define INC_BOX_H
#include "Matrix_3x3.h"
#include "Vec3.h"
#ifdef MPI
# include "Parallel.h"
#endif
/// Hold box information; 3xlengths, 3xangles.
class Box {
  public:
    enum BoxType { NOBOX=0, ORTHO, TRUNCOCT, RHOMBIC, NONORTHO }; 

    Box();
    //Box(const double*);
    //Box(const float*);
    //Box(Matrix_3x3 const&);
    Box(const Box&);
    Box& operator=(const Box&);
    void swap(Box&);
#   ifdef MPI
    int SyncBox(Parallel::Comm const&);
#   endif

    void SetupFromShapeMatrix(const double*);

    void SetupFromXyzAbg(const double*);

    //void SetBetaLengths(double,double,double,double);
    //void SetBox(const double*);
    //void SetBox(const float*);
    //void SetBox(Matrix_3x3 const&);
    //void SetBox(float,float,float,float,float,float);
    //void SetTruncOct();
    //void SetNoBox();
    //void SetMissingInfo(const Box&);
    /// Calculate Frac->Cart and Cart->Frac matrices.
    //double ToRecip(Matrix_3x3&, Matrix_3x3&) const;
    /// Calculate unit cell matrix, optionally scaling lengths.
    //Matrix_3x3 UnitCell(double) const;
    /// Print Box info to STDOUT
    void PrintInfo() const;

    //void SetX(double xin)     { box_[0] = xin; }
    //void SetY(double yin)     { box_[1] = yin; }
    //void SetZ(double zin)     { box_[2] = zin; }
    //void SetAlpha(double ain) { box_[3] = ain; }
    //void SetBeta(double bin)  { box_[4] = bin; }
    //void SetGamma(double gin) { box_[5] = gin; }

    const char* TypeName() const { return BoxNames_[btype_]; }
    BoxType Type() const { return btype_;  }
    double BoxX()  const { return box_[0]; }
    double BoxY()  const { return box_[1]; }
    double BoxZ()  const { return box_[2]; }
    double Alpha() const { return box_[3]; }
    double Beta()  const { return box_[4]; }
    double Gamma() const { return box_[5]; }
    bool HasBox()  const { return (btype_ != NOBOX); }
    Vec3 Center()  const { return Vec3(box_[0]/2.0, box_[1]/2.0, box_[2]/2.0); }
    Vec3 Lengths() const { return Vec3(box_[0], box_[1], box_[2]);             }

    // For interfacing with file IO
    //double* boxPtr()             { return box_; }
    //const double* boxPtr() const { return box_; }

    //double const& operator[](int idx) const { return box_[idx]; }
    //double&       operator[](int idx)       { return box_[idx]; }
  private:
    static const double TRUNCOCTBETA_;
    static const double TruncOctDelta_;
    static const double TruncOctMin_;
    static const double TruncOctMax_;
    static const double TruncOctEps_;
    static const char* BoxNames_[];

    static inline bool IsTruncOct(double);
    static inline bool BadTruncOctAngle(double);
    static inline bool IsAngle(double,double);

    void SetBoxType();

    /// Calculate fractional matrix from unit cell matrix.
    static inline double CalcFracFromUcell(Matrix_3x3&, Matrix_3x3 const&);
    /// Calculate unit cell matrix from XYZ ABG array.
    static inline void CalcUcellFromXyzAbg(Matrix_3x3&, const double*);
    /// Calculate unit cell matrix from XYZ ABG array based on box type, optionally scaling lengths.
    static void CalcUcellFromXyzAbg(Matrix_3x3&, BoxType, const double*, double);
    /// Calculate XYZ ABG array from unit cell matrix
    static void CalcXyzAbgFromUcell(double*, Matrix_3x3 const&);
    /// Calculate XYZ ABG array from symmetric shape matrix
    static void CalcXyzAbgFromShape(double*, const double*);

    //static Vec3 RecipLengths(Matrix_3x3 const&);
    //static void XyzAbgToUcell(Matrix_3x3&, const double*);


    //int debug_; // TODO: Replace with ifdefs or just comment out?
    BoxType btype_;       ///< Box type
    double box_[6];       ///< Box X Y Z alpha beta gamma
    Matrix_3x3 unitCell_; ///< Unit cell (Cartesian) matrix
    Matrix_3x3 fracCell_; ///< Fractional coord. cell matrix
    double cellVolume_;   ///< Unit cell volume
};
#endif
