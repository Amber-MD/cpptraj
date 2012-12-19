#ifndef INC_FRAME_H
#define INC_FRAME_H
#include "Atom.h"
#include "AtomMask.h"
#include "Box.h"
// Class: Frame
/// Hold coordinates, perform various operations/transformations on them.
/** Intended to hold coordinates e.g. from a trajectory or reference frame,
  * along with box coordinates (used in imaging calculations) and optionally 
  * with mass information and/or velocity information. Frame can be set up
  * coords only, coords and masses, or coords/masses/velocities. Mass is stored
  * since several functions (like COM, RMSD, Inertia etc) have the option to 
  * factor in the mass of the atoms involved, and this avoids having to pass a 
  * mass pointer in, which takes the burden of keeping track of mass away from 
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
    typedef std::vector<float> CRDtype;
    // Construction/Destruction/Assignment
    Frame();
    ~Frame();
    /// Set up empty frame for given # of atoms.
    Frame(int);
    /// Set up to be the size of given atom array (including masses).
    Frame(std::vector<Atom> const&);
    /// Copy input frame according to input mask.
    Frame(Frame const&, AtomMask const&);
    Frame(const Frame&);
    Frame& operator=(Frame);
    // Convert to/from CRDtype arrays
    void SetFromCRD(CRDtype const&, int);
    void SetFromCRD(CRDtype const&, int, AtomMask const&);
    CRDtype ConvertToCRD(int) const;
    // Access internal data
    void printAtomCoord(int);
    void Info(const char*);
    void AddXYZ(const double *);
    void AddVec3(Vec3 const&);
    double& operator[](int idx)        { return X_[idx];        } // TODO: Make const?
    const double& operator[](int idx) const { return X_[idx];   }
    bool empty()                 const { return (natom_ == 0);  }
    bool HasVelocity()           const { return (V_ != NULL);   }
    int Natom()                  const { return natom_;         }
    int size()                   const { return ncoord_;        }
    double Temperature()         const { return T_;             }
    const double* XYZ(int atnum) const { return X_ + (atnum*3); } 
    const double* CRD(int idx)   const { return X_ + idx;       } 
    double Mass(int atnum)       const { return Mass_[atnum];   }
    const Box& BoxCrd()          const { return box_;           }
    // Routines for accessing internal data pointers
    inline double* xAddress() { return X_;          }
    inline double* vAddress() { return V_;          }
    inline double* bAddress() { return box_.Dptr(); }
    inline double* tAddress() { return &T_;         }
    // Frame memory allocation/reallocation
    int SetupFrame(int);
    int SetupFrameM(std::vector<Atom> const&);
    int SetupFrameV(std::vector<Atom> const&,bool);
    int SetupFrameFromMask(AtomMask const&, std::vector<Atom> const&);
    // Frame setup of coords (no memory realloc)
    void SetCoordinates(Frame const&, AtomMask const&);
    void SetCoordinates(Frame const&);
    void SetFrame(Frame const&, AtomMask const&);
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
    // -------------------------------------------------------------------------
    // NOTE: These functions are placed in the header since most modern 
    //       compilers will try to inline them which results in a decent
    //       speedup for most routines (e.g. when imaging).
    /// \return true if 1st two coord sets are 0; indicates possible corruption.
    bool CheckCoordsInvalid() {
      if (natom_ > 1) {
        return (X_[0] == 0.0 && X_[1] == 0.0 && X_[2] == 0.0 &&
                X_[3] == 0.0 && X_[4] == 0.0 && X_[5] == 0.0   );
      }
      return false;
    }
    /// \return Center of mass of atoms in mask.
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
    /// \return Geometric center of atoms in mask.
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
    /// \return Center of mass of atoms in range.
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
    /// \return Geometric center of atoms in range.
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
    /// Scale coordinates of atoms in mask by given X|Y|Z constants
    void Scale(AtomMask const&, double, double, double);
    /// Translate atoms in range by Vec
    void Translate(const Vec3& Vec, int firstAtom, int lastAtom) {
      int startatom3 = firstAtom * 3;
      int stopatom3 = lastAtom * 3;
      for (int i = startatom3; i < stopatom3; i += 3) {
        X_[i  ] += Vec[0];
        X_[i+1] += Vec[1];
        X_[i+2] += Vec[2];
      }
    }
    /// Translate atom by Vec
    void Translate(const Vec3& Vec, int atom) {
      int icrd = atom * 3;
      X_[icrd  ] += Vec[0];
      X_[icrd+1] += Vec[1];
      X_[icrd+2] += Vec[2];
    }
    /// Translate all atoms by Vec
    void Translate(const Vec3& Vec) {
      for (int i = 0; i < ncoord_; i += 3) {
        X_[i  ] += Vec[0];
        X_[i+1] += Vec[1];
        X_[i+2] += Vec[2];
      }
    }
    /// Translate all atoms by negative Vec
    void NegTranslate(const Vec3& Vec) {
      for (int i = 0; i < ncoord_; i += 3) {
        X_[i  ] -= Vec[0];
        X_[i+1] -= Vec[1];
        X_[i+2] -= Vec[2];
      }
    }
    /// Rotate all coords by matrix
    void Rotate(Matrix_3x3 const& T) {
      for (int i = 0; i < ncoord_; i += 3) {
        double x = X_[i  ];
        double y = X_[i+1];
        double z = X_[i+2];
        X_[i  ] = (x*T[0]) + (y*T[1]) + (z*T[2]);
        X_[i+1] = (x*T[3]) + (y*T[4]) + (z*T[5]);
        X_[i+2] = (x*T[6]) + (y*T[7]) + (z*T[8]);
      }
    }
    /// Rotate all atoms in mask by matrix
    void Rotate(Matrix_3x3 const& RotMatrix, AtomMask& mask) {
      for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom) {
        double* XYZ = X_ + (*atom * 3) ;
        double x = XYZ[0];
        double y = XYZ[1];
        double z = XYZ[2];
        XYZ[0] = (x*RotMatrix[0]) + (y*RotMatrix[1]) + (z*RotMatrix[2]);
        XYZ[1] = (x*RotMatrix[3]) + (y*RotMatrix[4]) + (z*RotMatrix[5]);
        XYZ[2] = (x*RotMatrix[6]) + (y*RotMatrix[7]) + (z*RotMatrix[8]);
      }
    }
    /// Apply translation followed by rotation followed by second translation
    void Trans_Rot_Trans(Vec3 const& t1, Matrix_3x3 const& R, Vec3 const& t2) {
      for (int i = 0; i < ncoord_; i+=3) {
        double x = X_[i  ] + t1[0];
        double y = X_[i+1] + t1[1];
        double z = X_[i+2] + t1[2];
        X_[i  ] = x*R[0] + y*R[1] + z*R[2] + t2[0];
        X_[i+1] = x*R[3] + y*R[4] + z*R[5] + t2[1];
        X_[i+2] = x*R[6] + y*R[7] + z*R[8] + t2[2];
      }
    }
    // -------------------------------------------------------------------------
    void Center(AtomMask const&, bool,bool);
    Vec3 CenterOnOrigin(bool);
    // Coordinate calculation
    double RMSD(Frame &, bool );
    double RMSD(Frame &, Matrix_3x3&, Vec3&, Vec3&, bool);
    double RMSD_CenteredRef( Frame const&, bool);
    double RMSD_CenteredRef( Frame const&, Matrix_3x3&, Vec3&, bool);
    double RMSD_NoFit(Frame const&,bool);
    double DISTRMSD( Frame const& );

    Vec3 SetAxisOfRotation(int, int);
    Vec3 CalculateInertia(AtomMask const&, Matrix_3x3&);

  private:
    typedef std::vector<double> Darray;
    static const size_t COORDSIZE_;
    static const size_t BOXSIZE_;

    int natom_;     ///< Number of atoms stored in frame.
    int maxnatom_;  ///< Maximum number of atoms this frame can store.
    int ncoord_;    ///< Number of coordinates stored in frame (natom * 3).
    Box box_;       ///< Box coords, 3xlengths, 3xangles
    double T_;      ///< Temperature
    double* X_;     ///< Coord array, X0 Y0 Z0 X1 Y1 Z1 ...
    double* V_;     ///< Velocities (same arrangement as Coords).
    Darray Mass_;   ///< Masses.

    void swap(Frame&, Frame&);
    void ReallocateX();
};
#endif
