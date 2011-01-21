#ifndef INC_FRAME_H
#define INC_FRAME_H
// Frame.h
// Hold the coordinates of a trajectory frame or reference frame
#include "AtomMask.h"

class Frame {
  public:
    double *X;     // Coord array, X0 Y0 Z0 X1 Y1 Z1 ...
    int natom;     // Number of atoms
    int N;         // Number of coords in X, natom*3
    double box[6]; // Box coords, 3xlengths, 3xangles
    double T;      // Temperature
    Frame *V;      // Velocities
    double *Mass;  // Hold mass. Several functions (like COM, RADGYR etc) have 
                   // the option to factor in the mass of the atoms involved. 
                   // Mass is passed in when the frame is created, and is 
                   // modified if necessary by SetFrameFromMask (which is the
                   // case when e.g. calculating RMSD).
    Frame(int, double*);
    Frame(AtomMask *, double *);
    ~Frame();
    void printAtomCoord(int);
    void SetCoord(double *, int);
    Frame *Copy();

    char *BufferToBox(char *, int, int);
    char *BufferToFrame(char *, int);
    char *FrameToBuffer(char *, const char*, int,int);
    char *BoxToBuffer(char *, int, const char *, int);
    void floatToFrame(float *);
    void frameToFloat(float *);

    double COM(AtomMask*, double *, bool);
    double COM(double*,bool,int,int);
    double COM(double*,bool);

    void BoxToRecip(double *, double *);
    double DIST2(AtomMask*, AtomMask*, bool, int, double *, double *);
    double DIST2(int, int, int, double *, double *);
    double DIST(int, int);

    double ANGLE(AtomMask*, AtomMask*, AtomMask*);
    double ANGLE(int, int, int);
    double DIHEDRAL(AtomMask *, AtomMask *, AtomMask *, AtomMask *);
    double RADGYR(AtomMask *, bool, double *);

    void SetFrameFromMask(Frame*, AtomMask *);
    int SetFrameCoordsFromMask(double *, AtomMask *);

    void Translate(double *);
    void Translate(double *, int, int);
    void Rotate(double *);
    void Center(AtomMask *, double *,bool);
    void ShiftToCenter( Frame * );

    double RMSD(Frame*, double*, double*,bool);
    double RMSD(Frame*,bool);
};
#endif
