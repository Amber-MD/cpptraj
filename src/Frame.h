#ifndef INC_FRAME_H
#define INC_FRAME_H
#include "Atom.h"
#include "AtomMask.h"
#include "Vec3.h"
// Class: Frame
/// Hold coordinates, perform various operations/transformations on them.
/** Intended to hold coordinates e.g. from a trajectory or reference frame,
  * along with box coordinates (used in imaging calculations) and optionally 
  * with mass information and/or velocity information. Frame can be set up
  * coords only, coords and masses, or coords/masses/velocities. Mass is stored
  * since several functions (like COM, RADGYR etc) have the option to factor in 
  * the mass of the atoms involved, and this avoids having to pass a mass 
  * pointer in, which takes the burden of keeping track of mass away from 
  * actions etc. Mass is stored when the frame is initially created, and is 
  * modified if necessary by SetFrame (which is the case when e.g. calculating
  * per-residue RMSD).
  *
  * - Implementation Details:
  *
  * In addition to the constructors, there are two classes of routine that
  * can be used to set up Frames. The SetupX routines do any memory allocation,
  * and assign masses, and the SetX routines assign coordinates. The SetX 
  * routines will dynamically adjust the size of the frame up to maxnatom, but
  * no reallocation will occur so the frame should be set up for the largest
  * possible # of atoms it will hold. This avoids expensive reallocations.
  * The representation of coordinates (X) and velocities (V) are double*
  * instead of STL vectors so as to easily interface with the FileIO routines
  * which are much faster than iostream ops. 
  */
class Frame {
  public:
    // Construction/Destruction/Assignment
    Frame();
    virtual ~Frame(); // Destructor is virtual since this class can be inherited
    Frame(int);
    Frame(std::vector<Atom> const&);
    Frame(Frame const&, AtomMask const&);
    Frame(const Frame&);
    Frame& operator=(Frame);
    // Convert to/from arrays
    Frame &operator=(std::vector<float> const&);
    std::vector<float> ConvertToFloat(AtomMask const&);
    // Access internal data
    void printAtomCoord(int);
    void Info(const char*);
    void GetAtomXYZ(double*, int);
    void AddXYZ(const double *);
    bool empty()                 { return (natom_ == 0);        }
    bool HasVelocity()           { return (V_ != NULL);         }
    int Natom()                  { return natom_;               }
    int size()                   { return ncoord_;              }
    double Temperature()         { return T_;                   }
    const double* XYZ(int atnum) { return X_ + (atnum*3);       } // TODO: Replace?
    const double* CRD(int idx)   { return X_ + idx;             } // TODO: Replace?
    double& operator[](int idx)  { return X_[idx];              } // TODO: Make const?
    // Box routines
    void BoxXYZ(double* XYZ) { XYZ[0]=box_[0]; XYZ[1]=box_[1]; XYZ[2]=box_[2]; }
    double BoxX() { return box_[0]; }
    double BoxY() { return box_[1]; }
    double BoxZ() { return box_[2]; }
    const double* Box() { return box_; }
    // Routines for accessing internal data pointers
    inline double* xAddress() { return X_;   }
    inline double* vAddress() { return V_;   }
    inline double* bAddress() { return box_; }
    inline double* tAddress() { return &T_;  }
    // Frame memory allocation/reallocation
    int SetupFrame(int);
    int SetupFrameM(std::vector<Atom> const&);
    int SetupFrameV(std::vector<Atom> const&,bool);
    int SetupFrameFromMask(AtomMask const&, std::vector<Atom> const&);
    // Frame setup of coords (no memory realloc)
    void SetCoordinatesByMask(const double*, AtomMask const&); 
    void SetCoordinates(Frame const&, AtomMask const&);
    void SetCoordinates(Frame const&);
    void SetFrame(Frame const&, AtomMask const&);
    Frame* FrameCopy(); // TODO: Obsolete
    // Frame Setup with Atom Mapping
    void SetCoordinatesByMap(Frame const&, std::vector<int>const&);
    void SetReferenceByMap(Frame const&, std::vector<int>const&);
    void SetTargetByMap(Frame const&, std::vector<int>const&);
    // Basic Arithmetic
    void ZeroCoords();
    Frame & operator+=(const Frame&);
    Frame & operator-=(const Frame&);
    Frame & operator*=(const Frame&);
    const Frame operator*(const Frame&) const;
    int Divide(Frame const&, double); 
    void Divide(double);
    void AddByMask(Frame const&, AtomMask const&); 
    // Center of mass / Geometric Center
    Vec3 VCenterOfMass(AtomMask const&);
    Vec3 VGeometricCenter(AtomMask const&);
    double CenterOfMass(double*, AtomMask const&);
    double GeometricCenter(double*, AtomMask const&);
    double CenterOfMass(double*,int,int);
    double GeometricCenter(double*,int,int);
    // Coordinate manipulation
    void SCALE(AtomMask const&, double, double, double);
    void Translate(const double *);
    void Translate(const double *, int,int);
    void Translate(const double *, int);
    void Trans_Rot_Trans(const double *, const double *);
    void Rotate(const double *);
    void InverseRotate(const double *);
    void Center(AtomMask const&, bool,bool);
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
    double RADGYR(AtomMask &, bool, double *);
    double RMSD(Frame &, double*, double*,bool);
    double RMSD_CenteredRef( Frame const&, double[9], double[6], bool);
    double RMSD(Frame const&,bool);
    double DISTRMSD( Frame& );

    void SetAxisOfRotation(double *, int, int);
    void RotateAroundAxis(double *, AtomMask &);
    void CalculateInertia(AtomMask&, double*, double*);

  private:
    typedef std::vector<double> Darray;
    static const size_t COORDSIZE_;
    static const size_t BOXSIZE_;

    int natom_;     ///< Number of atoms stored in frame.
    int maxnatom_;  ///< Maximum number of atoms this frame can store.
    int ncoord_;    ///< Number of coordinates stored in frame (natom * 3).
    double box_[6]; ///< Box coords, 3xlengths, 3xangles
    double T_;      ///< Temperature
    double* X_;     ///< Coord array, X0 Y0 Z0 X1 Y1 Z1 ...
    double* V_;     ///< Velocities (same arrangement as Coords).
    Darray Mass_;   ///< Masses.

    void swap(Frame&, Frame&);
};
#endif
