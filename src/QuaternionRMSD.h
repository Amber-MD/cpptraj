#ifndef INC_QUATERNIONRMSD_H
#define INC_QUATERNIONRMSD_H
#include "Frame.h"
/// Print quaternion RMSD citation to stdout
void PrintQrmsdCitation();
/// \return Quaternion RMSD, set rotation matrix (in cols) and target translation vector.
double QuaternionRMSD_CenteredRef(Frame const&, Frame&, Matrix_3x3&, Vec3&, bool);
/// \return Quaternion RMSD
double QuaternionRMSD_CenteredRef(Frame const&, Frame&, bool);
#endif
