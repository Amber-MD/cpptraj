#ifndef INC_ACTION_NASTRUCT_H
#define INC_ACTION_NASTRUCT_H
// Action_NAstruct
#include <vector>
#include <map>
#include "Action.h"
#include "AxisType.h"
#include "Range.h"
#ifdef NASTRUCTDEBUG
#include "PDBfile.h"
#endif
/// Basic Nucleic acid structure analysis. 
/** Calculate nucleic acid base/base pair structural parameters.
  * Algorithms for calculation of base/base pair structural parameters
  * adapted from:
  *   Babcock MS, Pednault EPD, Olson WK, "Nucleic Acid Structure Analysis: 
  *   Mathematics for Local Cartesian and Helical Structure Parameters That
  *   Are Truly Comparable Between Structures", J. Mol. Biol. (1994) 237,
  *   125-156.
  * NA base reference frame coordinates taken from:
  *   Olson WK, Bansal M, Burley SK, Dickerson RE, Gerstein M, Harvey SC,
  *   Heinemann U, Lu XJ, Neidle S, Shekked Z, Sklenar H, Suzuki M, Tung CS,
  *   Westhof E, Wolberger C, Berman H, "A Standard Reference Frame for the 
  *   Description of Nucleic Acid Base-pair Geometry", J. Mol. Biol. (2001)
  *   313, 229-237.
  */
class Action_NAstruct: public Action {
  public:
    Action_NAstruct();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_NAstruct(); }
    static void Help();

    ~Action_NAstruct();

    void Print();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    // Variables
    std::vector<AxisType> RefCoords;    ///< Hold reference frame coordinates for each base
    std::vector<AxisType> BaseAxes;     ///< Hold axis coordinates for each base
    std::vector<AxisType> BasePairAxes; ///< Hold axis coordinates for each base pair
    std::vector<AtomMask> ExpMasks;     ///< Hold atom masks of each base for input coordinates
    std::vector<AtomMask> FitMasks;     ///< Hold atom masks used to RMS-fit each base
    std::vector<int> BasePair;          ///< [Base1-0,Base2-0,IsAnti-0], [Base1-1...
    std::vector<int> NumberOfHbonds_;   ///< # hbonds between each base pair
    double HBdistCut2_;                 ///< distance Cutoff^2 for determining hydrogen bonds
    //double HBangleCut2_;                ///< Angle Cutoff^2 for determining if bases can h bond
    double originCut2_;                 ///< Cutoff^2 for determining base-pairing vi origins
    int maxResSize_;                    ///< Max residue size, used to set up frames for RMS fit.
    int debug_;
    Range resRange;                     ///< Range to search for NA residues.
    bool printheader_;                  ///< If true, print header to naout files.
    bool useReference_;                 ///< If true, use reference to determine base pairing.
    std::string outputsuffix_;          ///< Output file suffix (BP.<suffix> etc)
    std::string dataname_;              ///< NA DataSet name (default NA).

    typedef std::map<std::string,AxisType::NAbaseType> ResMapType;
    ResMapType CustomMap;
    // Datasets - 1 entry per BP/BPstep
    typedef std::vector<DataSet*> Darray;
    Darray SHEAR_;
    Darray STRETCH_;
    Darray STAGGER_;
    Darray BUCKLE_;
    Darray PROPELLER_;
    Darray OPENING_;
    Darray BPHBONDS_;
    Darray MAJOR_;
    Darray MINOR_;
    Darray SHIFT_;
    Darray SLIDE_;
    Darray RISE_;
    Darray TILT_;
    Darray ROLL_;
    Darray TWIST_;
    Darray XDISP_;
    Darray YDISP_;
    Darray HRISE_;
    Darray INCL_;
    Darray TIP_;
    Darray HTWIST_;
    // TODO: Replace these with new DataSet type
    DataSetList* masterDSL_;
    // Functions
    void ClearLists();
    int GCpair(AxisType *, AxisType *);
    int ATpair(AxisType *, AxisType *);
    int basesArePaired(AxisType *, AxisType *);
    int determineBasePairing();
    int setupBaseAxes(Frame *);
    int calculateParameters(AxisType &, AxisType &, AxisType*, double*);
    int helicalParameters(AxisType &, AxisType &, double *);
    int determineBaseParameters(int);
    int determineBasepairParameters(int);
#   ifdef NASTRUCTDEBUG
    // DEBUG - Class to hold PDB output
    class AxisPDBwriter : protected PDBfile {
      public:
        AxisPDBwriter() : pdbatom_(1) { }
        ~AxisPDBwriter() { CloseFile(); }
        void Open(const char *fname) { OpenWrite( fname ); }
        void Write(AxisType &axis, const double* XYZin, int res, const char *resname) {
          const double* XYZ = XYZin;
          for (int i = 0; i < axis.Natom(); ++i) {
            WriteRec( PDBfile::ATOM, pdbatom_++, axis.AtomName(i), resname, 'X', res+1,
                      XYZ[0], XYZ[1], XYZ[2], 1.0, 0.0, "", false);
            XYZ += 3;
          }
        }
        void WriteAxes(AxisType &axis, int res, const char *resname) {
          double OXYZ[3], V[3];
          axis.OXYZ( OXYZ );
          // Origin
          WriteRec( PDBfile::ATOM, pdbatom_++, "Orig", resname, 'X', res+1,
                    OXYZ[0], OXYZ[1], OXYZ[2], 1.0, 0.0, "", false);
          // X vector
          axis.RX( V );
          WriteRec( PDBfile::ATOM, pdbatom_++, "X", resname, 'X', res+1,
                    OXYZ[0]+V[0], OXYZ[1]+V[1], OXYZ[2]+V[2], 1.0, 0.0, "", false);
          // Y vector
          axis.RY( V );
          WriteRec( PDBfile::ATOM, pdbatom_++, "Y", resname, 'X', res+1,
                    OXYZ[0]+V[0], OXYZ[1]+V[1], OXYZ[2]+V[2], 1.0, 0.0, "", false);
          // Z vector
          axis.RZ( V );
          WriteRec( PDBfile::ATOM, pdbatom_++, "Z", resname, 'X', res+1,
                    OXYZ[0]+V[0], OXYZ[1]+V[1], OXYZ[2]+V[2], 1.0, 0.0, "", false);
        }
      private:
        int pdbatom_;
    };
    // DEBUG - used to trigger AxisPDBwriter for first call of calculateParameters
    bool calcparam;
#   endif
};
#endif
