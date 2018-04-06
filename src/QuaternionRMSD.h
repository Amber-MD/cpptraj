#ifndef INC_QUATERNIONRMSD_H
#define INC_QUATERNIONRMSD_H
#include "Frame.h"
double QuaternionRMSD_CenteredRef(Frame const&, Frame&, Matrix_3x3&, Vec3&, bool, double);
double QuaternionRMSD_CenteredRef(Frame const&, Frame&, Matrix_3x3&, Vec3&, bool);
double QuaternionRMSD_CenteredRef(Frame const&, Frame&, bool);
#endif
