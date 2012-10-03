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
#ifdef NASTRUCTDEBUG
    Frame(int,const double*);
#endif
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
    double BoxX() { return box_[0]; }
    double BoxY() { return box_[1]; }
    double BoxZ() { return box_[2]; }
    Vec3 BoxLengths() { return Vec3( box_[0], box_[1], box_[2] ); }
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
    // DEBUG -------------------------------------------------------------------
    // NOTE: Placing these functions in the header since most modern compilers
    //       will actually try to inline them which results in a decent
    //       speedup for most routines (e.g. when imaging).
    Vec3 VCenterOfMass( AtomMask const& Mask ) {
      double Coord0 = 0.0;
      double Coord1 = 0.0;
      double Coord2 = 0.0;
      double sumMass = 0.0;
      for (AtomMask::const_iterator atom = Mask.begin(); atom != Mask.end(); ++atom)
      {
        unsigned int xidx = (*atom) * 3;
        double mass = Mass_[*atom];
        sumMass += mass;
        Coord0 += ( X_[xidx  ] * mass );
        Coord1 += ( X_[xidx+1] * mass );
        Coord2 += ( X_[xidx+2] * mass );
      }
      if (sumMass == 0.0) return Vec3();
      return Vec3( Coord0 / sumMass, Coord1 / sumMass, Coord2 / sumMass );
    }
    Vec3 VGeometricCenter( AtomMask const& Mask ) {
      double Coord0 = 0.0;
      double Coord1 = 0.0;
      double Coord2 = 0.0;
      for (AtomMask::const_iterator atom = Mask.begin(); atom != Mask.end(); ++atom)
      {
        unsigned int xidx = (*atom) * 3;
        Coord0 += X_[xidx  ];
        Coord1 += X_[xidx+1];
        Coord2 += X_[xidx+2];
      }
      double sumMass = (double)Mask.Nselected();
      if (sumMass == 0) return Vec3();
      return Vec3( Coord0 / sumMass, Coord1 / sumMass, Coord2 / sumMass );
    }
    Vec3 VCenterOfMass(int startAtom, int stopAtom) {
      double Coord0 = 0.0;
      double Coord1 = 0.0;
      double Coord2 = 0.0;
      double sumMass = 0.0;
      Darray::iterator mass = Mass_.begin() + startAtom;
      int startAtom3 = startAtom * 3;
      int stopAtom3 = stopAtom * 3;
      for (int i = startAtom3; i < stopAtom3; i += 3) {
        sumMass += (*mass);
        Coord0 += ( X_[i  ] * (*mass) );
        Coord1 += ( X_[i+1] * (*mass) );
        Coord2 += ( X_[i+2] * (*mass) );
        ++mass;
      }
      if (sumMass == 0.0) return Vec3();
      return Vec3( Coord0 / sumMass, Coord1 / sumMass, Coord2 / sumMass );
    }
    Vec3 VGeometricCenter(int startAtom, int stopAtom) {
      double Coord0 = 0.0;
      double Coord1 = 0.0;
      double Coord2 = 0.0;
      int startAtom3 = startAtom * 3;
      int stopAtom3 = stopAtom * 3;
      for (int i = startAtom3; i < stopAtom3; i += 3) {
        Coord0 += X_[i  ];
        Coord1 += X_[i+1];
        Coord2 += X_[i+2];
      }
      double sumMass = (double)(stopAtom - startAtom);
      if (sumMass == 0) return Vec3();
      return Vec3( Coord0 / sumMass, Coord1 / sumMass, Coord2 / sumMass );
    }
    // END DEBUG ---------------------------------------------------------------
    double CenterOfMass(double*, AtomMask const&);
    double GeometricCenter(double*, AtomMask const&);
    double CenterOfMass(double*,int,int);
    double GeometricCenter(double*,int,int);
    // Coordinate manipulation
    void SCALE(AtomMask const&, double, double, double);
    void Translate(const double *);
    void Translate(const double *, int,int);
    void Translate(const double *, int);
    void Translate(const Vec3& Vec, int firstAtom, int lastAtom) {
      double Vec0 = Vec[0];
      double Vec1 = Vec[1];
      double Vec2 = Vec[2];
      int startatom3 = firstAtom * 3;
      int stopatom3 = lastAtom * 3;
      for (int i = startatom3; i < stopatom3; i += 3) {
        X_[i  ] += Vec0;
        X_[i+1] += Vec1;
        X_[i+2] += Vec2;
      }
    }
    void Trans_Rot_Trans(const double *, const double *);
    void Rotate(const double *);
    void InverseRotate(const double *);
    void Center(AtomMask const&, bool,bool);
    void CenterReference(double *, bool);
    void ShiftToGeometricCenter();
    // Coordinate calculation
    double BoxToRecip(double *, double *);
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
