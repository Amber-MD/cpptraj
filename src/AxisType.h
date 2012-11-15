#ifndef INC_AXISTYPE_H
#define INC_AXISTYPE_H
#include "Topology.h"
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
class AxisType {
  public:
    /// Type for each standard NA base.
    enum NAbaseType { UNKNOWN_BASE, ADE, CYT, GUA, THY, URA };

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

    const char* ResName();
    int ResNum() { return residue_number; }
    int ResNum2() { return second_resnum; }
    const double* xAddress() { return X_; }
    bool AtomNameIs(int, char *);
    const char* AtomName(int);
    void PrintAtomNames();
    void PrintAxisInfo(const char *);

    void SetCoordsFromFrame( Frame& ); // TODO: Make const
    void StoreRotMatrix(double*,double*);
    void StoreBPresnums(int,int);

    enum RefReturn { NA_OK, NA_UNKNOWN, NA_ERROR };
    RefReturn SetRefCoord(Topology *, int, AtomMask &,AtomMask&,NAbaseType);
    void FlipYZ();
    void FlipXY();
    // P/O4' atom routines
    bool HasPatom() { return patomidx_ > -1; }
    bool HasO4atom() { return o4atomidx_ > -1; }
    int Pidx() { return patomidx_; }
    int O4idx() { return o4atomidx_; }
    void SetPcrd( const double* Xin ) { 
      atomcrd_[0] = Xin[0];
      atomcrd_[1] = Xin[1];
      atomcrd_[2] = Xin[2];
    }
    void SetO4crd( const double* Xin ) {
      atomcrd_[3] = Xin[0];
      atomcrd_[4] = Xin[1];
      atomcrd_[5] = Xin[2];
    }
    const double* Pcrd() { return atomcrd_; }
    const double* O4crd() { return atomcrd_+3; }
#ifdef NASTRUCTDEBUG
    const char* BaseName();
    int Natom() { return natom_; }
    double operator[](int idx) { return X_[idx]; }
#endif
  private:
    static const int ADENATOM;
    static const char ADEnames[][5];
    static const int ADEhbonds[];
    static const double ADEcoords[][3];
    static const int CYTNATOM;
    static const char CYTnames[][5];
    static const int CYThbonds[];
    static const double CYTcoords[][3];
    static const int GUANATOM;
    static const char GUAnames[][5];
    static const int GUAhbonds[];
    static const double GUAcoords[][3];
    static const int THYNATOM;
    static const char THYnames[][5];
    static const int THYhbonds[];
    static const double THYcoords[][3];
    static const int URANATOM;
    static const char URAnames[][5];
    static const int URAhbonds[];
    static const double URAcoords[][3];
    /// Current number of atoms
    int natom_;
    /// Maximum number of atoms
    int maxnatom_;
    /// Number of coordinates
    int Ncoord_;
    /// Hold coordinates
    double* X_;
    /// Strings corresponding to NAbaseType
    static const char NAbaseName[][4];
    /// Atom Names
    std::vector<NameType> Name;
    /// Original residue number
    int residue_number;
    /// Second base number if this is a base pair
    int second_resnum;
    /// Origin coordinates
    double origin[3];
    /// Index of phosphorus atom if present
    int patomidx_;
    /// Index of O4' atom if present
    int o4atomidx_;
    /// Hold phosphorus/O4' atom coordinates
    double atomcrd_[6];
#ifdef NASTRUCTDEBUG
    /// DEBUG - Storage for writing out BaseName + residue_number 
    std::string basename_num_;
#endif
    /// Identify NA base from residue name
    NAbaseType ID_base(NameType const&);
    /// Allocate memory
    int AllocAxis(int);
};
#endif  
