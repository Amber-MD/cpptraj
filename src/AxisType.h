#ifndef INC_AXISTYPE_H
#define INC_AXISTYPE_H
#include "Frame.h"
#include "AmberParm.h"
#ifdef NASTRUCTDEBUG
#  include "CpptrajFile.h"
#endif
/*! \file AxisType.h
    \brief Hold classes and functions used for NA structure analysis.
 */

// Class: AxisType
/// Frame for NA bases, intended for use with NAstruct Action.
/** AxisType is a special kind of Frame. It will be used in 2 cases: 
  * - To hold the coordinates of reference bases for RMS fitting onto input
  *   coordinates in order to obtain reference frames.
  * - For holding the coordinates of the origin and XYZ axes of the 
  *   reference frames themselves.
  */
class AxisType : public Frame {
    // The following enumerated types are used to index the standard
    // NA reference frame coordinates for each base.
    //enum ADEatoms { C1p, N9, C8, N7, C5, C6, N6, N1, C2, N3, C4 };
    //enum CYTatoms { C1p, N1, C2, O2, N3, C4, N4, C5, C6 };
    //enum GUAatoms { C1p, N9, C8, N7, C5, C6, O6, N1, C2, N2, N3, C4 };
    //enum THYatoms { C1p, N1, C2, O2, N3, C4, O4, C5, C7, C6 };
    //enum URAatoms { C1p, N1, C2, O2, N3, C4, O4, C5, C6 };
    // The following arrays contain the corresponding atom names.
    static const NAME ADEnames[];
    static const double ADEcoords[][3];
    static const int ADEhbonds[];
    static const NAME CYTnames[];
    static const double CYTcoords[][3];
    static const int CYThbonds[];
    static const NAME GUAnames[];
    static const double GUAcoords[][3];
    static const int GUAhbonds[];
    static const NAME THYnames[];
    static const double THYcoords[][3];
    static const int THYhbonds[];
    static const NAME URAnames[];
    static const double URAcoords[][3];
    static const int URAhbonds[];
  public:
    /// Type for each standard NA base.
    enum NAbaseType { UNKNOWN_BASE, ADE, CYT, GUA, THY, URA };
  private:
    /// Identify NA base from residue name
    NAbaseType ID_base(char*);
    static const char NAbaseName[][4];
    NAME *Name;         ///< Atom/Axis names
    int residue_number; ///< Original residue number

    int AllocAxis(int);
  public:
    NAbaseType ID;
    double R[9];
    double *HbondCoord[3];
    int HbondAtom[3];

    AxisType();
    AxisType(const AxisType&);
    AxisType & operator=(const AxisType&);
    ~AxisType();
    void RX(double Vin[3]);
    void RY(double Vin[3]);
    void RZ(double Vin[3]);
    double *Origin();
    char *BaseName();
    int BaseNum() { return residue_number; }
    bool AtomNameIs(int, char *);
    char *AtomName(int);
    void PrintAtomNames();
    void SetFromFrame(AxisType *);
    void StoreRotMatrix(double*);
    void SetPrincipalAxes();
    enum RefReturn { NA_OK, NA_UNKNOWN, NA_ERROR };
    RefReturn SetRefCoord(AmberParm *, int, AtomMask &);
    void FlipYZ();
    void FlipXY();
#ifdef NASTRUCTDEBUG
    void WritePDB(CpptrajFile *, int, char *, int *); // DEBUG
#endif
};
#endif  
