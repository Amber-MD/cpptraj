#ifndef INC_FRAME_H
#define INC_FRAME_H
#include "Atom.h"
#include "AtomMask.h"
// Class: Frame
/// Hold coordinates, perform various operations/transformations on them.
/** Intended to hold coordinates e.g. from a trajectory or reference frame,
  * along with box coordinates (used in imaging calculations) and optionally 
  * with mass information and/or velocity information. Mass is stored since 
  * several functions (like COM, RADGYR etc) have the option to factor in 
  * the mass of the atoms involved, and this avoids having to pass a mass 
  * pointer in, which takes the burden of keeping track of mass away from 
  * actions etc. Mass is stored when the frame is initially created, and is 
  * modified if necessary by SetFrame (which is the case when e.g. calculating
  * per-residue RMSD).
  */
class Frame {
    //friend class TrajectoryFile;
  public:
    /// Potential imaging types 
    enum ImageType { NOIMAGE=0, ORTHO, NONORTHO };
    // Construction/Destruction/Assignment
    Frame();
    virtual ~Frame(); // Destructor is virtual since this class can be inherited
    Frame(int);
    Frame(std::vector<Atom> const&);
    Frame(Frame &, AtomMask &);
    Frame(const Frame&);
    Frame& operator=(Frame);
    // Convert to/from arrays
    Frame &operator=(const std::vector<float>&);
    std::vector<float> ConvertToFloat(AtomMask &);
    // Access internal data
    void printAtomCoord(int);
    void Info(const char*);
    void GetAtomXYZ(double*, int);
    void AddXYZ(const double *);
    int Natom();
    int size() { return (int)X_.size(); }
    double Temperature() { return T_; }
    bool empty();
    double MaxImagedDistance();
    const double* XYZ(int atnum) { return &(X_[0]) + (atnum*3); } // TODO: Replace?
    double& operator[](int idx)  { return X_[idx];              }
    // Box routines
    void BoxXYZ(double* XYZ) { XYZ[0]=box_[0]; XYZ[1]=box_[1]; XYZ[2]=box_[2]; }
    double BoxX() { return box_[0]; }
    double BoxY() { return box_[1]; }
    double BoxZ() { return box_[2]; }
    const double* Box() { return box_; }
    // Routines for accessing internal data pointers
    inline double* xAddress() { return xaddress_; }
    inline double* vAddress() { return vaddress_; }
    inline double* bAddress() { return box_;      }
    inline double* tAddress() { return &T_;       }
    // Frame memory allocation/reallocation
    int SetupFrame(int);
    int SetupFrameM(std::vector<Atom> const&);
    int SetupFrameV(std::vector<Atom> const&,bool);
    int SetupFrameFromMask(AtomMask &, std::vector<Atom> const&);
    // Frame setup of coords (no memory realloc)
    void SetCoordinatesByMask(const double*, AtomMask const&); 
    void SetCoordinates(Frame&,AtomMask&);
    void SetCoordinates(Frame&);
    void SetFrame(Frame&,AtomMask&);
    Frame* FrameCopy(); // TODO: Obsolete
    // Frame Setup with Atom Mapping
    void SetCoordinatesByMap(Frame &, std::vector<int>const&);
    void SetReferenceByMap(Frame&,std::vector<int>const&);
    void SetTargetByMap(Frame&,std::vector<int>const&);
    // Basic Arithmetic
    void ZeroCoords();
    Frame & operator+=(const Frame&);
    Frame & operator-=(const Frame&);
    Frame & operator*=(const Frame&);
    const Frame operator*(const Frame&) const;
    int Divide(Frame&, double); 
    void Divide(double);
    void AddByMask(Frame &, AtomMask &); 
    // Center of mass / Geometric Center
    double CenterOfMass(double*, AtomMask&);
    double GeometricCenter(double*, AtomMask&);
    double CenterOfMass(double*,int,int);
    double GeometricCenter(double*,int,int);
    // Coordinate manipulation
    void SCALE(AtomMask&, double, double, double);
    void Translate(double *);
    void Translate(double *, int,int);
    void Translate(double *, int);
    void Trans_Rot_Trans(double *, double *);
    void Rotate(double *);
    void InverseRotate(double *);
    void Center(AtomMask &, bool,bool);
    void CenterReference(double *, bool);
    void ShiftToGeometricCenter();
    // Imaging
    void SetupImageTruncoct(double*, AtomMask*,bool,bool);
    void ImageNonortho(bool, double*, double*, double*, bool, bool, bool, std::vector<int> const&);
    void ImageNonortho(double*, double*, bool, bool, double*, double*, double*);
    void SetupImageOrtho(double*, double*, bool);
    void ImageOrtho(double*,double*, bool, bool, std::vector<int> const&);
    void ImageOrtho(double*, double*, double*, double*);
    void UnwrapNonortho( Frame&, AtomMask& );
    void UnwrapOrtho( Frame&, AtomMask& );
    // Coordinate calculation
    double BoxToRecip(double *, double *);
    double DIST2(AtomMask&, AtomMask&, bool, ImageType, double *, double *);
    double DIST2(double *, double *, ImageType, double *, double *);
    double DIST2(int, int, ImageType, double *, double *);
    double DIST2(double*, int, ImageType, double *, double *);
    double DIST(int, int);
    double DIST2(int, int);
    double COORDDIST(int, int);
    double COORDDIST2(int, int);
    void COORDVECTOR(double*, int, int);
    double ANGLE(AtomMask&, AtomMask&, AtomMask&,bool);
    double ANGLE(int, int, int);
    double DIHEDRAL(AtomMask&, AtomMask&, AtomMask&, AtomMask&,bool);
    double DIHEDRAL(int,int,int,int);
    double PUCKER(AtomMask&,AtomMask&,AtomMask&,AtomMask&,AtomMask&,int,bool,bool);
    double RADGYR(AtomMask &, bool, double *);
    double RMSD(Frame&, double*, double*,bool);
    double RMSD_CenteredRef( Frame &, double[9], double[6], bool);
    double RMSD(Frame&,bool);
    double DISTRMSD( Frame& );

    void SetAxisOfRotation(double *, int, int);
    void RotateAroundAxis(double *, double, AtomMask &);
    void CalculateInertia(AtomMask&, double*, double*);

  private:
    typedef std::vector<double> Darray;

    int natom_;     ///< Number of atoms stored in frame.
    Darray X_;      ///< Coord array, X0 Y0 Z0 X1 Y1 Z1 ...
    Darray::iterator lastx_;
    Darray V_;      ///< Velocities
    Darray Mass_;   ///< Masses
    double box_[6]; ///< Box coords, 3xlengths, 3xangles
    double T_;      ///< Temperature
    double* xaddress_; ///< For direct input into frame coordinate array.
    double* vaddress_; ///< For direct input into frame velocity array.

    void swap(Frame&, Frame&);
};
#endif
