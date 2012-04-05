#ifndef INC_FRAME_H
#define INC_FRAME_H
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
    friend class TrajectoryFile;
  public:
    // Construction/Destruction/Assignment
    Frame();
    Frame(int);
    Frame(double*, int);
    Frame(int, double *);
    Frame(AtomMask &, double *);
    Frame(Frame &, AtomMask &);
    Frame(const Frame&);
    void Info(const char*);
    virtual ~Frame(); // Destructor is virtual since this class can be inherited
    Frame & operator=(const Frame&);
    // Convert to/from arrays
    Frame &operator=(const std::vector<float>&);
    std::vector<float> ConvertToFloat(AtomMask &);
    std::vector<double> ConvertToDouble();
    double *DoubleArray();
    // Access internal data
    const double *CoordPtr();
    void ConvertToPtrajXYZ(double *, double *, double *, double*);
    void SetFromPtrajXYZ(double *, double *, double *);
    void GetAtomXYZ(double*, int);
    int Natom();
    bool empty();
    double MaxImagedDistance();
    // Frame memory allocation/reallocation
    int SetupFrame(int);
    int SetupFrame(int, double*);
    int SetupFrameV(int,double*,bool);
    int SetupFrameFromMask(AtomMask&, double *);
    // Frame setup of coords (mass/velo)
    void SetCoordinates(Frame&,AtomMask&);
    void SetCoordinates(Frame &, std::vector<int>&);
    void SetReferenceByMap(Frame&,std::vector<int>&);
    void SetTargetByMap(Frame&,std::vector<int>&);
    void SetCoordinates(Frame&);
    void SetFrame(Frame&,AtomMask&);

    Frame *FrameCopy();

    // Basic Arithmetic
    void ZeroCoords();
    Frame & operator+=(const Frame&);
    Frame & operator-=(const Frame&);
    Frame & operator*=(const Frame&);
    const Frame operator*(const Frame&) const;
    int Divide(Frame&, double); 
    void Divide(double);
    void AddByMask(Frame &, AtomMask &); 
    // Coordinate manipulation
    void Translate(double *);
    void Translate(double *, int,int);
    void Trans_Rot_Trans(double *, double *);
    void Rotate(double *);
    void InverseRotate(double *);
    void Center(AtomMask *, bool,bool);
    void CenterReference(double *, bool);
    void ShiftToGeometricCenter();
    void ImageNonortho(bool, AtomMask *, bool, bool, bool, std::vector<int> &);
    void ImageOrtho(bool, bool, bool, std::vector<int> &);

    void printAtomCoord(int);

    // Center of mass
    double CenterOfMass(AtomMask*, double *);
    double GeometricCenter(AtomMask*, double *);
    double CenterOfMass(double*,int,int);
    double GeometricCenter(double*,int,int);
    // Coordinate calculation
    double BoxToRecip(double *, double *);
    double DIST2(AtomMask*, AtomMask*, bool, int, double *, double *);
    double DIST2(int, int, int, double *, double *);
    double DIST2(double*, int, int, double *, double *);
    double DIST(int, int);
    double DIST2(int, int);
    double COORDDIST(int, int);
    double COORDDIST2(int, int);
    void COORDVECTOR(double*, int, int);
    double ANGLE(AtomMask*, AtomMask*, AtomMask*,bool);
    double ANGLE(int, int, int);
    double DIHEDRAL(AtomMask *, AtomMask *, AtomMask *, AtomMask *,bool);
    double DIHEDRAL(int,int,int,int);
    double PUCKER(AtomMask*,AtomMask*,AtomMask*,AtomMask*,AtomMask*,int,bool,bool);
    double RADGYR(AtomMask *, bool, double *);
    double RMSD(Frame*, double*, double*,bool);
    double RMSD_CenteredRef( Frame &, double[9], double[6], bool);
    double RMSD(Frame*,bool);
    double DISTRMSD( Frame * );

    void SetAxisOfRotation(double *, int, int);
    void RotateAroundAxis(double *, double, AtomMask &);

    // TEST: Iterator to coordinates
    class Iterator : public std::iterator<std::forward_iterator_tag, double*> {
      public:
        Iterator() : xptr_(0) {}
        Iterator(double *p) : xptr_(p) {}
        ~Iterator() {}
        // Assignment
        Iterator &operator=(const Iterator &rhs) {
          xptr_ = rhs.xptr_;
          return *this;
        }
        // Relations
        bool operator==(const Iterator &rhs) {
          return (xptr_ == rhs.xptr_);
        }
        bool operator!=(const Iterator &rhs) {
          return (xptr_ != rhs.xptr_);
        }
        // Pre-increment
        Iterator &operator++() {
          xptr_ += 3;
          return *this;
        }
        // Post-increment
        Iterator operator++(int) { // Return copy, not reference
          Iterator tmp(*this);
          ++(*this);
          return tmp;
        }
        // Value - Make const?
        double* operator*() {
          return xptr_;
        }
        double* operator->() {
          return xptr_;
        }
      private:
        double *xptr_;
    }; // END class Iterator
    Iterator begin() {
      return( Iterator(X_) );
    }
    Iterator end() {
      return( Iterator(X_ + Ncoord_) );
    }
    
  protected:
    static const size_t COORDSIZE_;
    static const size_t BOXSIZE_;

    int natom_;     ///< Number of atoms
    int maxnatom_;  ///< Number of atoms for which space has been allocated
    int Ncoord_;    ///< Number of coords, natom*3
    double *X_;     ///< Coord array, X0 Y0 Z0 X1 Y1 Z1 ...
    double *V_;     ///< Velocities
    double *Mass_;  ///< Mass
    double box_[6]; ///< Box coords, 3xlengths, 3xangles
    double T_;      ///< Temperature

    //void swap(Frame &, Frame &);
    int ReallocateFrame(int, bool, bool);
};
#endif
