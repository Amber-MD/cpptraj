#ifndef INC_FRAME_H
#define INC_FRAME_H
// Frame.h
// Hold the coordinates of a trajectory frame or reference frame
#include "AtomMask.h"

class Frame {
  public:
    double *X;     // Coord array, X0 Y0 Z0 X1 Y1 Z1 ...
    int natom;     // Number of atoms
    int N;         // Number of entries in X, natom*3
    double box[6]; // Box coords, 3xlengths, 3xangles
    double T;      // Temperature
    Frame *V;      // Velocities
    double *Mass;  // Hold mass. For most functions right now mass is passed in,
                   // however for functions like RMSD where the coordinates were
                   // probably modified the frame will need its own copy of the
                   // modified mass array. Could potentially make it so that
                   // information such as mass is passed in during frame creation.

    Frame(int, double*);
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
    void ClosestImage(double *, double *, int *);
    double MinImageNonOrtho(double *, double *, double *, double *, bool, int *);
    double DIST_ImageNonOrtho(AtomMask *, AtomMask *, bool, double *, double *);
    double DIST_ImageNonOrtho(double *, double *, double, double *, int *);
    double DIST_ImageOrtho(AtomMask *, AtomMask *, bool);
    double DIST(AtomMask*, AtomMask*,bool);
    double DIST(int, int);
    double ANGLE(AtomMask*, AtomMask*, AtomMask*);
    double ANGLE(int, int, int);
    double DIHEDRAL(AtomMask *, AtomMask *, AtomMask *, AtomMask *);
    void SetFrameFromMask(Frame*, AtomMask *);
    void Translate(double *);
    void Translate(double *, int, int);
    void Rotate(double *);
    void Center(AtomMask *, double *,bool);
    void ShiftToCenter( Frame * );
    double RMSD(Frame*, double*, double*,bool);
    double RMSD(Frame*,bool);
};
#endif
