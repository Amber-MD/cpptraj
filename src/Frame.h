#ifndef INC_FRAME_H
#define INC_FRAME_H
// Frame.h
// Hold the coordinates of a trajectory frame or reference frame
#include "AtomMask.h"

class Frame {
    void ClosestImage(double *, double *, int *);
    double DIST2_ImageNonOrtho(double *, double *);
    double DIST2_ImageNonOrtho(double *, double *, double, int *);
    double DIST2_ImageOrtho(double *, double *);
    double DIST2_NoImage(double *, double *);

  public:
    double *X;     // Coord array, X0 Y0 Z0 X1 Y1 Z1 ...
    int natom;     // Number of atoms
    int N;         // Number of entries in X, natom*3
    double box[6]; // Box coords, 3xlengths, 3xangles
    double ucell[9];
    double recip[9];
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
    void BoxToRecip();
    double MinImageNonOrtho2(double *, double *, bool, int *);
    double DIST2(AtomMask*, AtomMask*, bool, int);
    double DIST2(int, int, int);
    double DIST(int, int);
    double ANGLE(AtomMask*, AtomMask*, AtomMask*);
    double ANGLE(int, int, int);
    double DIHEDRAL(AtomMask *, AtomMask *, AtomMask *, AtomMask *);
    double RADGYR(AtomMask *, bool, double *);
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
