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
  * - For holding the coordinates and rotation matrix/ origin of input 
  *   frames.
  */
class AxisType : public Frame {
  public:
    /// Type for each standard NA base.
    enum NAbaseType { UNKNOWN_BASE, ADE, CYT, GUA, THY, URA };
  private:
    /// Identify NA base from residue name
    NAbaseType ID_base(char*);
    /// Strings corresponding to NAbaseType
    static const char NAbaseName[][4];
    /// Atom Names
    NAME *Name;
    /// Original residue number
    int residue_number;
    /// Second base number if this is a base pair
    int second_resnum;
    /// Origin coordinates
    double origin[3];
    /// DEBUG - Storage for writing out BaseName + residue_number 
    char basename_num[32];

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
    void OXYZ(double Vin[3]);
    double *Origin();

    char *BaseName();
    char *ResName();
    int ResNum() { return residue_number; }
    int ResNum2() { return second_resnum; }
    bool AtomNameIs(int, char *);
    char *AtomName(int);
    void PrintAtomNames();
    void PrintAxisInfo(const char *);

    void SetAxisFromMask(AxisType &, AtomMask &);
    void StoreRotMatrix(double*,double*);
    void StoreBPresnums(int,int);

    enum RefReturn { NA_OK, NA_UNKNOWN, NA_ERROR };
    RefReturn SetRefCoord(AmberParm *, int, AtomMask &,AtomMask&,NAbaseType);
    void FlipYZ();
    void FlipXY();
#ifdef NASTRUCTDEBUG
    void WritePDB(CpptrajFile &, int, char *, int *); // DEBUG
    void WriteAxesPDB(CpptrajFile &, int, char *, int *);
#endif
};
#endif  
