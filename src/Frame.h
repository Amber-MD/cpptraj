#ifndef INC_FRAME_H
#define INC_FRAME_H
/// Class: Frame
/// Hold the coordinates of a trajectory frame or reference frame, along 
/// with box coordinates (used in imaging calculations) and optionally with 
/// mass information and/or velocity information. Mass is stored since several 
/// functions (like COM, RADGYR etc) have the option to factor in the mass of 
/// the atoms involved, and this avoids having to pass a mass pointer in, 
/// which takes the burden of keeping track of mass away from actions etc.
/// Mass is stored when the frame is initially created, and is modified if 
/// necessary by SetFrameFromMask (which is the case when e.g. calculating 
/// RMSD).
#include "AtomMask.h"
class Frame {
  public:
    double *X;     // Coord array, X0 Y0 Z0 X1 Y1 Z1 ...
    int natom;     // Number of atoms
    int maxnatom;  // Number of atoms for which space has been allocated
    int N;         // Number of coords, natom*3
    double box[6]; // Box coords, 3xlengths, 3xangles
    double T;      // Temperature
    double *V;     // Velocities
    double *Mass;  // Mass

    Frame();
    Frame(int, double*);
    Frame(int,double*,bool);
    Frame(AtomMask *, double *);
    virtual ~Frame();             // Destructor is virtual since this class can be inherited
    Frame *Copy();
    int Resize(int,bool,bool);
    // Coordinate manipulation
    void ZeroCoords();
    void AddCoord(Frame*);
    void Divide(double);
    void Translate(double *);
    void Translate(double *, int);
    void Rotate(double *);
    void InverseRotate(double *);
    void Center(AtomMask *, double *,bool);
    void ShiftToCenter( Frame * );
    // Coordinate assignment/extraction
    void printAtomCoord(int);
    void GetCoord(double *, int);
    void SetCoord(int, double *);
    double *Coord(int);
    void SetFrameFromMask(Frame*, AtomMask *);
    int SetFrameCoordsFromMask(double *, AtomMask *);
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
    double COORDDIST(int, int);
    double ANGLE(AtomMask*, AtomMask*, AtomMask*,bool);
    double ANGLE(int, int, int);
    double DIHEDRAL(AtomMask *, AtomMask *, AtomMask *, AtomMask *,bool);
    double PUCKER(AtomMask*,AtomMask*,AtomMask*,AtomMask*,AtomMask*,int,bool,bool);
    double RADGYR(AtomMask *, bool, double *);
    double RMSD(Frame*, double*, double*,bool);
    double RMSD(Frame*,bool);
    double DISTRMSD( Frame * );
};
#endif
