#ifndef INC_FRAME_H
#define INC_FRAME_H
#include <vector>
#include "AtomMask.h"
#include "Box.h"
#include "CoordinateInfo.h"
// Forward declarations
class Atom;
/// Hold coordinates, perform various operations/transformations on them.
/** Intended to hold coordinates e.g. from a trajectory or reference frame,
  * along with box coordinates (used in imaging calculations), mass information,
  * and optionally velocity information. Frame can be set up coords only (all 
  * masses set to 1.0), coords and masses, or coords/masses/velocities. Mass is 
  * stored since several functions (like COM, RMSD, Inertia etc) have the option
  * to factor in the mass of the atoms involved, and this avoids having to pass
  * a mass pointer in, which takes the burden of keeping track of mass away from 
  * actions etc. Mass is stored when the frame is initially created, and is 
  * modified if necessary by SetFrame (which is the case when e.g. calculating
  * per-residue RMSD).
  *
  * - Implementation Details:
  *
  * In addition to the constructors, there are two classes of routine that
  * can be used to set up Frames. The SetupX routines do any memory allocation,
  * and assign masses, and the SetX routines assign coordinates/velocities. The
  * SetX routines will dynamically adjust the size of the frame up to maxnatom,
  * but no reallocation will occur so the frame should be set up for the largest
  * possible # of atoms it will hold. This avoids expensive reallocations.
  * The representation of coordinates (X) and velocities (V) are double*
  * instead of STL vectors so as to easily interface with the FileIO routines
  * which tend to be much faster than iostream ops. 
  */
class Frame {
  public:
    typedef std::vector<double> Darray;
    // Construction/Destruction/Assignment
    Frame();
    ~Frame();
    /// Set up empty frame for given # of atoms.
    Frame(int);
    /// Set up to be the size of given atom array (including masses).
    Frame(std::vector<Atom> const&);
    /// Copy input frame according to input mask.
    Frame(Frame const&, AtomMask const&);
    /// Set up frame for # of atoms with external coordinates data.
    Frame(int, double*);
    Frame(const Frame&);
    Frame& operator=(Frame);
    /// For dealing with replica indices TODO put in ReplicaInfo
    typedef std::vector<int> RemdIdxType;
    /// For reading replica values
    typedef std::vector<double> RemdValType;
    /// \return Size of Frame in memory
    size_t DataSize() const;
    /// Print XYZ coordinates for given atom.
    void printAtomCoord(int) const;
    /// Print information about the frame.
    void Info(const char*) const;
    /// Set atom/coord count to zero but do not clear memory.
    void ClearAtoms();
    /// Add the given XYZ coord to this frame.
    void AddXYZ(const double *);
    /// Add the given Vec3 to this frame.
    void AddVec3(Vec3 const&);
    /// Swap the coordinates and velocities of two atoms
    void SwapAtoms(int, int);
    // Access internal data
    double& operator[](int idx)             { return X_[idx];        }
    const double& operator[](int idx) const { return X_[idx];        }
    bool empty()                      const { return (natom_ == 0);  }
    bool HasVelocity()                const { return (V_ != 0);      }
    bool HasForce()                   const { return (F_ != 0);      }
    int Natom()                       const { return natom_;         }
    int size()                        const { return ncoord_;        }
    int NrepDims()                    const { return (int)remd_indices_.size(); } // TODO: deprecate
    double Temperature()              const { return T_;             }
    int Step()                        const { return step_;          }
    double pH()                       const { return pH_;            }
    double RedOx()                    const { return redox_;         }
    double Time()                     const { return time_;          }
    /// \return CoordinateInfo that describes the Frame
    CoordinateInfo CoordsInfo()       const;
    /// \return pointer to start of XYZ coords for given atom.
    const double* XYZ(int atnum)      const { return X_ + (atnum*3); } 
    /// \return pointer to specified coordinate.
    const double* CRD(int idx)        const { return X_ + idx;       }
    /// \return pointer to start of velocity XYZ for given atom.
    const double* VelXYZ(int atnum)     const { return V_ + (atnum*3); }
    /// \return pointer to start of force XYZ for given atom.
    const double* FrcXYZ(int atnum)     const { return F_ + (atnum*3); }
    /// \return mass of specified atom.
    double Mass(int atnum)            const { return Mass_[atnum];   }
    /// \return Box information
    const Box& BoxCrd()               const { return box_;           }
    /// \return replica indices
    RemdIdxType const& RemdIndices()  const { return remd_indices_;  }
    /// \return overall replica index
    int RepIdx()                      const { return repidx_;        }
    /// \return overall coordinate index
    int CrdIdx()                      const { return crdidx_;        }
    /// Set box from another box
    void SetBox( Box const& b ) { box_ = b; }
    /// Modify box in place
    Box& ModifyBox() { return box_; }
    /// Set temperature
    void SetTemperature(double tIn) { T_ = tIn;     }
    /// Set step
    void SetStep( int sIn )         { step_ = sIn;  }
    /// Set pH
    void Set_pH(double phIn)        { pH_ = phIn;   }
    /// Set RedOx potential
    void SetRedOx(double rIn)       { redox_ = rIn; }
    /// Set time
    void SetTime(double tIn)        { time_ = tIn;  }
    /// Set replica index
    void SetRepIdx(int rIn)         { repidx_ = rIn; }
    /// Set coordinate index
    void SetCrdIdx(int cIn)         { crdidx_ = cIn; }
    /// Set masses
    void SetMass(std::vector<Atom> const&);
    /// Copy atoms from input frame to here
    void CopyFrom(Frame const&, int, int);
    /// Copy unit from input frame to here
    void CopyFrom(Frame const&, Unit const&);
    // ----- Access to internal data pointers ----
    inline double* xAddress() { return X_;                }
    inline double* vAddress() { return V_;                }
    inline double* fAddress() { return F_;                }
    inline double* tAddress() { return &T_;               }
    inline double* mAddress() { return &time_;            }
    inline int* iAddress()    { return &remd_indices_[0]; }
    inline int* repidxPtr()   { return &repidx_;          }
    inline int* crdidxPtr()   { return &crdidx_;          }
    inline const double* xAddress() const { return X_;                }
    inline const double* vAddress() const { return V_;                }
    inline const double* fAddress() const { return F_;                }
    inline const double* tAddress() const { return &T_;               }
    inline const double* mAddress() const { return &time_;            }
    inline const int* iAddress()    const { return &remd_indices_[0]; }
    inline const int* repidxPtr()   const { return &repidx_;          }
    inline const int* crdidxPtr()   const { return &crdidx_;          }
    // ----- Frame memory allocation routines ----
    /// Allocate frame for given # atoms, no mass or velocity.
    int SetupFrame(int);
    /// Allocate frame for given # atoms with mass, no velocity. 
    int SetupFrameM(std::vector<Atom> const&);
    /// Allocate frame with given XYZ coords and masses, no velocity.
    int SetupFrameXM(Darray const&, Darray const&);
    /// Allocate frame based on given frame.
    int SetupFrame(Frame const&);
    /// Allocate frame for given # atoms with mass and opt. velocity/indices.
    int SetupFrameV(std::vector<Atom> const&, CoordinateInfo const&);
    /// Allocate frame for selected # atoms, coords/mass only.
    int SetupFrameFromMask(AtomMask const&, std::vector<Atom> const&);
    // ----- Add/remove components from Frame ----
    int AddVelocities(Darray const&);
    void RemoveVelocities();
    int AddForces(Darray const&);
    void RemoveForces();
    int AddMasses(Darray const&);
    void RemoveMasses();
    // ----- Frame coords set routines -----------
    /// Copy coordinates, box, and temp. from input frame according to mask. 
    void SetCoordinates(Frame const&, AtomMask const&);
    /// Copy only coordinates from input frame to this frame.
    void SetCoordinates(Frame const&);
    /// Copy only coordinates and box info from input frame to this frame.
    void SetCoordAndBox(Frame const&);
    /// Set coordinates from external memory.
    int SetCoordinates(int, double*);
    /// Copy entire input frame according to mask.
    void SetFrame(Frame const&, AtomMask const&);
    /// Copy entire input frame
    void SetFrame(Frame const&);
    /// Zero the force array
    void ZeroForces();
    /// Zero the velocity array
    void ZeroVelocities();
    // ----- Frame coordinate remapping ----------
    /// Copy entire input frame, reorder according to input map. 
    void SetCoordinatesByMap(Frame const&, std::vector<int>const&);
    /// Modify this frame to include only mapped atoms from input frame. TODO use mask?
    void StripUnmappedAtoms(Frame const&, std::vector<int>const&);
    /// Copy only input coordinates, reorder according to input map.
    void ModifyByMap(Frame const&, std::vector<int>const&);
    // ----- Basic Arithmetic --------------------
    void ZeroCoords();
    Frame& operator+=(const Frame&);
    Frame& operator-=(const Frame&);
    Frame& operator*=(const Frame&);
    const Frame operator*(const Frame&) const;
    const Frame operator-(const Frame&) const;
    int Divide(Frame const&, double); 
    void Divide(double);
    void Multiply(double);
    /// Increment atoms in this frame by the selected atoms in given frame.
    int AddByMask(Frame const&, AtomMask const&);
    /// \return Total momentum vector for atoms in mask
    Vec3 VMomentum(AtomMask const&) const;
    /// \return Total momentum vector for atoms in mask; also set total mass.
    Vec3 VMomentum(AtomMask const&, double&) const;
    // -------------------------------------------------------------------------
    // NOTE: These functions are placed in the header since most modern 
    //       compilers will try to inline them which results in a decent
    //       speedup for most routines (e.g. when imaging).
    /// \return true if 1st two coord sets are 0; indicates possible corruption.
    inline bool CheckCoordsInvalid() const;
    /// \return Center of mass of atoms in mask.
    inline Vec3 VCenterOfMass( AtomMask const& ) const;
    /// \return Geometric center of atoms in mask.
    inline Vec3 VGeometricCenter( AtomMask const& ) const;
    /// \return Center of mass of atoms in range.
    inline Vec3 VCenterOfMass(int, int) const;
    /// \return Geometric center of atoms in range.
    inline Vec3 VGeometricCenter(int, int) const;
    /// \return Center of mass of atoms in Unit
    inline Vec3 VCenterOfMass(Unit const&) const;
    /// \return Geometric center of atoms in Unit
    inline Vec3 VGeometricCenter(Unit const&) const;
    /// Translate atoms in mask by Vec
    inline void Translate(Vec3 const&, AtomMask const&);
    /// Translate atoms in range by Vec
    inline void Translate(Vec3 const&, int, int);
    /// Translate atoms in Unit by Vec
    inline void Translate(Vec3 const&, Unit const&);
    /// Translate atom by Vec
    inline void Translate(Vec3 const&, int);
    /// Translate all atoms by Vec
    inline void Translate(Vec3 const&);
    /// Translate all atoms by negative Vec
    inline void NegTranslate(Vec3 const&);
    /// Rotate all coords by matrix
    inline void Rotate(Matrix_3x3 const&);
    /// Rotate all atoms in range by matrix
    inline void Rotate(Matrix_3x3 const&, int, int);
    /// Rotate all atoms in unit by matrix
    inline void Rotate(Matrix_3x3 const&, Unit const&);
    /// Rotate all atoms in mask by matrix
    inline void Rotate(Matrix_3x3 const&, AtomMask const&);
    /// Apply inverse of rotation defined by matrix to all atoms in mask
    inline void InverseRotate(Matrix_3x3 const&, AtomMask const&);
    /// Apply translation followed by rotation followed by second translation
    inline void Trans_Rot_Trans(Vec3 const&, Matrix_3x3 const&, Vec3 const&);
    /// Apply translation, rotation, 2nd translation for atoms in mask
    inline void Trans_Rot_Trans(AtomMask const&, Vec3 const&, Matrix_3x3 const&, Vec3 const&);
    // -------------------------------------------------------------------------
    /// Scale coordinates of atoms in mask by given X|Y|Z constants
    void Scale(AtomMask const&, double, double, double);
    /// Translate atoms to origin.
    Vec3 CenterOnOrigin(bool);
    // Align on reference
    void Align(Frame const&, AtomMask const&);
    // Coordinate calculation
    double RMSD(Frame &, bool );
    double RMSD(Frame &, Matrix_3x3&, Vec3&, Vec3&, bool);
    double RMSD_CenteredRef( Frame const&, bool);
    double RMSD_CenteredRef( Frame const&, Matrix_3x3&, Vec3&, bool);
    double RMSD_NoFit(Frame const&,bool) const;
    inline double RMSD_FitToRef(Frame const&, Vec3 const&);
    double DISTRMSD( Frame const& ) const;
    /// Set axis of rotation to be around line connecting given atoms.
    Vec3 SetAxisOfRotation(int, int);
    /// Set axis of rotation to be around line connecting given points.
    Vec3 SetAxisOfRotation(Vec3 const&, Vec3 const&);
    /// Calculate inertia matrix.
    Vec3 CalculateInertia(AtomMask const&, Matrix_3x3&) const;
    /// Calculate temperature of atoms in mask.
    double CalcTemperature(AtomMask const&,int) const;
    /// Set an orthogonal bounding box
    void SetOrthoBoundingBox(std::vector<double> const& Radii, double);
#   ifdef MPI
    // ----- Parallel Routines -------------------
    int SendFrame(int, Parallel::Comm const&) const;
    int RecvFrame(int, Parallel::Comm const&);
    int SumToMaster(Parallel::Comm const&);
#   endif
  private:
    static const size_t COORDSIZE_;
    static const size_t BOXSIZE_;

    int natom_;     ///< Number of atoms stored in frame.
    int maxnatom_;  ///< Maximum number of atoms this frame can store.
    int ncoord_;    ///< Number of coordinates stored in frame (natom * 3).
    int step_;      ///< Simulation step
    Box box_;       ///< Box coords, 3xlengths, 3xangles
    double T_;      ///< Temperature
    double pH_;     ///< pH
    double redox_;  ///< RedOx potential
    double time_;   ///< Time
    double* X_;     ///< Coord array, X0 Y0 Z0 X1 Y1 Z1 ...
    double* V_;     ///< Velocities (same arrangement as Coords).
    double* F_;     ///< Frame (same arrangement as Coords).
    RemdIdxType remd_indices_; ///< replica indices.
    int repidx_;    ///< overall replica index
    int crdidx_;    ///< overall coordinate index.
    Darray Mass_;   ///< Masses.
    bool memIsExternal_; ///< True if Frame is not responsible for freeing memory.

    void swap(Frame&, Frame&);
    void IncreaseX();
    inline bool ReallocateX(int);
    /// Allocate coords/velo/force based on given num atoms and coordinate info.
    bool setupFrame(unsigned int, CoordinateInfo const&);
};
// ---------- INLINE FUNCTION DEFINITIONS --------------------------------------

bool Frame::CheckCoordsInvalid() const {
  if (natom_ > 1) {
    return (X_[0] == 0.0 && X_[1] == 0.0 && X_[2] == 0.0 &&
            X_[3] == 0.0 && X_[4] == 0.0 && X_[5] == 0.0   );
  }
  return false;
}

Vec3 Frame::VCenterOfMass( AtomMask const& Mask ) const {
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
  if (sumMass == 0.0) return Vec3(0,0,0);
  return Vec3( Coord0 / sumMass, Coord1 / sumMass, Coord2 / sumMass );
}

Vec3 Frame::VGeometricCenter( AtomMask const& Mask ) const {
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
  if (sumMass == 0) return Vec3(0,0,0);
  return Vec3( Coord0 / sumMass, Coord1 / sumMass, Coord2 / sumMass );
}

Vec3 Frame::VCenterOfMass(int startAtom, int stopAtom) const {
  double Coord0 = 0.0;
  double Coord1 = 0.0;
  double Coord2 = 0.0;
  double sumMass = 0.0;
  Darray::const_iterator mass = Mass_.begin() + startAtom;
  int startAtom3 = startAtom * 3;
  int stopAtom3 = stopAtom * 3;
  for (int i = startAtom3; i < stopAtom3; i += 3) {
    sumMass += (*mass);
    Coord0 += ( X_[i  ] * (*mass) );
    Coord1 += ( X_[i+1] * (*mass) );
    Coord2 += ( X_[i+2] * (*mass) );
    ++mass;
  }
  if (sumMass == 0.0) return Vec3(0,0,0);
  return Vec3( Coord0 / sumMass, Coord1 / sumMass, Coord2 / sumMass );
}

Vec3 Frame::VCenterOfMass(Unit const& unit) const {
  Vec3 out(0.0);
  for (Unit::const_iterator seg = unit.segBegin(); seg != unit.segEnd(); ++seg)
    out += VCenterOfMass(seg->Begin(), seg->End());
  return out;
}

Vec3 Frame::VGeometricCenter(int startAtom, int stopAtom) const {
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
  if (sumMass == 0) return Vec3(0,0,0);
  return Vec3( Coord0 / sumMass, Coord1 / sumMass, Coord2 / sumMass );
}

Vec3 Frame::VGeometricCenter(Unit const& unit) const {
  Vec3 out(0.0);
  for (Unit::const_iterator seg = unit.segBegin(); seg != unit.segEnd(); ++seg)
    out += VGeometricCenter(seg->Begin(), seg->End());
  return out;
}

void Frame::Translate(Vec3 const& Vec, AtomMask const& maskIn) {
  for (AtomMask::const_iterator atm = maskIn.begin(); atm != maskIn.end(); ++atm)
  {
    int idx = *atm * 3;
    X_[idx  ] += Vec[0];
    X_[idx+1] += Vec[1];
    X_[idx+2] += Vec[2];
  }
}

void Frame::Translate(Vec3 const& Vec, int firstAtom, int lastAtom) {
  int startatom3 = firstAtom * 3;
  int stopatom3 = lastAtom * 3;
  for (int i = startatom3; i < stopatom3; i += 3) {
    X_[i  ] += Vec[0];
    X_[i+1] += Vec[1];
    X_[i+2] += Vec[2];
  }
}

void Frame::Translate(Vec3 const& Vec, Unit const& unit) {
  for (Unit::const_iterator seg = unit.segBegin(); seg != unit.segEnd(); ++seg)
    Translate(Vec, seg->Begin(), seg->End());
}

void Frame::Translate(Vec3 const& Vec, int atom) {
  int icrd = atom * 3;
  X_[icrd  ] += Vec[0];
  X_[icrd+1] += Vec[1];
  X_[icrd+2] += Vec[2];
}

void Frame::Translate(Vec3 const& Vec) {
  for (int i = 0; i < ncoord_; i += 3) {
    X_[i  ] += Vec[0];
    X_[i+1] += Vec[1];
    X_[i+2] += Vec[2];
  }
}

void Frame::NegTranslate(Vec3 const& Vec) {
  for (int i = 0; i < ncoord_; i += 3) {
    X_[i  ] -= Vec[0];
    X_[i+1] -= Vec[1];
    X_[i+2] -= Vec[2];
  }
}

void Frame::Rotate(Matrix_3x3 const& T) {
  for (int i = 0; i < ncoord_; i += 3) {
    double x = X_[i  ];
    double y = X_[i+1];
    double z = X_[i+2];
    X_[i  ] = (x*T[0]) + (y*T[1]) + (z*T[2]);
    X_[i+1] = (x*T[3]) + (y*T[4]) + (z*T[5]);
    X_[i+2] = (x*T[6]) + (y*T[7]) + (z*T[8]);
  }
}

void Frame::Rotate(Matrix_3x3 const& T, int firstAtom, int lastAtom) {
  int startatom3 = firstAtom * 3;
  int stopatom3 = lastAtom * 3;
  for (int i = startatom3; i < stopatom3; i += 3) {
    double x = X_[i  ];
    double y = X_[i+1];
    double z = X_[i+2];
    X_[i  ] = (x*T[0]) + (y*T[1]) + (z*T[2]);
    X_[i+1] = (x*T[3]) + (y*T[4]) + (z*T[5]);
    X_[i+2] = (x*T[6]) + (y*T[7]) + (z*T[8]);
  }
}

void Frame::Rotate(Matrix_3x3 const& T, Unit const& unit) {
  for (Unit::const_iterator seg = unit.segBegin(); seg != unit.segEnd(); ++seg)
    Rotate(T, seg->Begin(), seg->End());
}

void Frame::Rotate(Matrix_3x3 const& RotMatrix, AtomMask const& mask) {
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

void Frame::InverseRotate(Matrix_3x3 const& RotMatrix, AtomMask const& mask) {
  for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom) {
    double* XYZ = X_ + (*atom * 3) ;
    double x = XYZ[0];
    double y = XYZ[1];
    double z = XYZ[2];
    XYZ[0] = (x*RotMatrix[0]) + (y*RotMatrix[3]) + (z*RotMatrix[6]);
    XYZ[1] = (x*RotMatrix[1]) + (y*RotMatrix[4]) + (z*RotMatrix[7]);
    XYZ[2] = (x*RotMatrix[2]) + (y*RotMatrix[5]) + (z*RotMatrix[8]);
  }
}

void Frame::Trans_Rot_Trans(Vec3 const& t1, Matrix_3x3 const& R, Vec3 const& t2) {
  for (int i = 0; i < ncoord_; i+=3) {
    double x = X_[i  ] + t1[0];
    double y = X_[i+1] + t1[1];
    double z = X_[i+2] + t1[2];
    X_[i  ] = x*R[0] + y*R[1] + z*R[2] + t2[0];
    X_[i+1] = x*R[3] + y*R[4] + z*R[5] + t2[1];
    X_[i+2] = x*R[6] + y*R[7] + z*R[8] + t2[2];
  }
}

void Frame::Trans_Rot_Trans(AtomMask const& mask, Vec3 const& t1, Matrix_3x3 const& R,
                            Vec3 const& t2)
{
  for (AtomMask::const_iterator at = mask.begin(); at != mask.end(); ++at) {
    int i = *at * 3;
    double x = X_[i  ] + t1[0];
    double y = X_[i+1] + t1[1];
    double z = X_[i+2] + t1[2];
    X_[i  ] = x*R[0] + y*R[1] + z*R[2] + t2[0];
    X_[i+1] = x*R[3] + y*R[4] + z*R[5] + t2[1];
    X_[i+2] = x*R[6] + y*R[7] + z*R[8] + t2[2];
  }
}

double Frame::RMSD_FitToRef(Frame const& Ref, Vec3 const& reftrans) {
  Matrix_3x3 U;
  Vec3 tgttrans;
  double rval = RMSD_CenteredRef( Ref, U, tgttrans, false );
  // Frame will have already been translated to origin.
  // Rotate and translate back to reference.
  for (int i = 0; i < ncoord_; i += 3) {
    double x = X_[i  ];
    double y = X_[i+1];
    double z = X_[i+2];
    X_[i  ] = x*U[0] + y*U[1] + z*U[2] + reftrans[0];
    X_[i+1] = x*U[3] + y*U[4] + z*U[5] + reftrans[1];
    X_[i+2] = x*U[6] + y*U[7] + z*U[8] + reftrans[2];
  }
  return rval;
}
#endif
