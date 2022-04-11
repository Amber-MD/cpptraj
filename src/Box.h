#ifndef INC_BOX_H
#define INC_BOX_H
#include "Matrix_3x3.h"
#include "Vec3.h"
#ifdef MPI
# include "Parallel.h"
#endif
/// Hold box information; unit and fractional cell vectors, 3xlengths, 3xangles.
class Box {
  public:
    /// Various box parameters; corresponds to XYZ ABG array.
    enum ParamType { X=0, Y, Z, ALPHA, BETA, GAMMA };
    /// Cell shape types
    enum CellShapeType { NO_SHAPE = 0,
                         CUBIC,        TETRAGONAL, ORTHORHOMBIC,
                         MONOCLINIC,   TRICLINIC,  HEXAGONAL,
                         RHOMBOHEDRAL, OCTAHEDRAL, RHOMBIC_DODECAHEDRON };
    /// CONSTRUCTOR
    Box();
    /// COPY CONSTRUCTOR
    Box(const Box&);
    /// ASSIGNMENT
    Box& operator=(const Box&);
    /// SWAP
    void swap(Box&);
#   ifdef MPI
    int BroadcastBox(Parallel::Comm const&);
    int SendBox(int, Parallel::Comm const&) const;
    int RecvBox(int, Parallel::Comm const&);
#   endif
    /// Remove all box information
    void SetNoBox();
    /// Rotate unit/frac cell vectors
    void RotateUcell(Matrix_3x3 const&);
    /// Inverse rotate unit/frac cell vectors
    void InverseRotateUcell(Matrix_3x3 const&);
    /// \return Current cell shape based on lengths and angles
    CellShapeType CellShape() const;
    // -------------------------------------------
    // Setup routines; will set the box type and perform checks.
    int SetupFromShapeMatrix(const double*);

    int SetupFromUcell(const double*);

    int SetupFromUcell(Matrix_3x3 const& ucell) { return SetupFromUcell(ucell.Dptr()); }

    int SetupFromXyzAbg(double,double,double,double,double,double);

    int SetupFromXyzAbg(const double*);

    // -------------------------------------------
    // Assign routines; just assign box lengths, not type; no checks.
    void AssignFromUcell(const double*);

    void AssignFromXyzAbg(double,double,double,double,double,double);

    void AssignFromXyzAbg(const double*);

    void AssignFromShapeMatrix(const double*);

    // -------------------------------------------
    /// Set incoming array of length 6 with shape matrix; suitable for Charmm
    void GetSymmetricShapeMatrix(double*) const;

    /// \return vector containing reciprocal lengths from fractional cell matrix
    Vec3 RecipLengths() const;

    /// Print Box info to STDOUT
    void PrintInfo() const;
    /// Print box debug info to stdout
    void PrintDebug(const char*) const;
    /// \return Name of current cell shape
    const char* CellShapeName() const { return CellShapeStr_[CellShape()]; }
    /// \return Specified XYZ ABG parameter
    double Param(ParamType p) const { return box_[p]; }
    /// \return True if box info present
    bool HasBox()             const { return (box_[0] > 0); }
    //bool HasBox()             const { return (btype_ != NOBOX); }
    /// \return Vector containing box center.
    Vec3 Center()             const { return Vec3(box_[0]/2.0, box_[1]/2.0, box_[2]/2.0); }
    /// \return Vector containing box center.
    // TODO Pretty sure this is the right way to calculate box center, but it breaks
    //      tests so more testing is needed.
    //Vec3 Center()             const { return (UnitCell().TransposeMult(Vec3(0.5))); }
    /// \return Vector containing box lengths
    Vec3 Lengths()            const { return Vec3(box_[0], box_[1], box_[2]);             }
    /// \return the unit cell matrix.
    /** The rows of this matrix represent the cell axes A, B C.
      * Fractional coords F can be converted to Cartesian coords C via
      * C = UnitCell().TransposeMult( F )
      */
    Matrix_3x3 const& UnitCell() const { return unitCell_; }
    /// \return the fractional cell matrix.
    /** Cartesian coords C can be converted to fractional coords F via
      * F = FracCell() * C
      */
    Matrix_3x3 const& FracCell() const { return fracCell_; }
    /// \return the cell volume
    double CellVolume()          const { return cellVolume_; }
    /// \return True if A is aligned X axis (xyz abg reference).
    bool Is_X_Aligned() const;
    /// \return True if X-aligned and cell is orthogonal.
    bool Is_X_Aligned_Ortho() const;
    /// \return True if unit cell off-diagonal elements are symmetric
    bool Is_Symmetric() const;

    /// \return true if the given angle is suitable for a truncated octahedron
    static bool IsTruncOct(double);
    // \return Truncated oct angle in degrees TODO should this be in Constants?
    static double TruncatedOctAngle() { return TRUNCOCTBETA_; }
 
   // For interfacing with file IO
    /// Double pointer starting at lengths (XYZ)
    const double* XyzPtr() const { return box_; }
    /// Double pointer starting at angles (ABG)
    const double* AbgPtr() const { return box_+3; }
  private:
    static const double TRUNCOCTBETA_;
    static const double TruncOctDelta_;
    static const double TruncOctMin_;
    static const double TruncOctMax_;
    static const double TruncOctEps_;
    static const double EqEps_;
    /// Names corresponding to ParamType
    static const char* ParamStr_[];
    /// Names corresponding to CellShapeType
    static const char* CellShapeStr_[];

    static inline bool BadTruncOctAngle(double);
    static inline bool IsEq(double,double);

    void printBoxStatus(const char*) const;
    /// \return 1 if box has issues, 0 otherwise.
    int CheckBox() const;

    /// Calculate fractional matrix from unit cell matrix.
    static inline double CalcFracFromUcell(Matrix_3x3&, Matrix_3x3 const&);
    /// Calculate unit cell matrix from XYZ ABG array.
    static inline void CalcUcellFromXyzAbg(Matrix_3x3&, const double*);
    /// Calculate XYZ ABG array from unit cell matrix
    static void CalcXyzAbgFromUcell(double*, Matrix_3x3 const&);
    /// Calculate XYZ ABG array from symmetric shape matrix
    static void CalcXyzAbgFromShape(double*, const double*);
    /// Calculate symmetric shape matrix from XYZ ABG array
    static void CalcShapeFromXyzAbg(double*, const double*);

    double box_[6];       ///< Box X Y Z alpha beta gamma
    Matrix_3x3 unitCell_; ///< Unit cell (Cartesian) matrix
    Matrix_3x3 fracCell_; ///< Fractional coord. cell matrix
    double cellVolume_;   ///< Unit cell volume
};
#endif
