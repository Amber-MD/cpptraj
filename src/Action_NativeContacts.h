#ifndef INC_ACTION_NATIVECONTACTS_H
#define INC_ACTION_NATIVECONTACTS_H
#include <set>
#include "Action.h"
#include "ImagedAction.h"
#include "DataSet_MatrixDbl.h"
/// Calculate the number of native/non-native contacts based on distance
/** Intended to combine and replace contacts, mindist, and maxdist actions.
  * Go through AtomMask arrays and calculate 
  */
class Action_NativeContacts : public Action {
  public:
    Action_NativeContacts();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_NativeContacts(); }
    static void Help();
  private:
    typedef std::vector<AtomMask> Marray;
    typedef std::vector<int> Iarray;
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    Marray SetupList(Topology const&, AtomMask const&, Iarray&) const;
    int SetupContactLists(Marray&, Marray&, Iarray&, Iarray&, 
                          Topology const&, Frame const&);
    int DetermineNativeContacts(Topology const&, Frame const&);
    inline double GetMinMax(AtomMask const&, AtomMask const&, Frame const&,
                            double&) const;

    double distance_;     ///< Cutoff distance
    int debug_;           ///< Action debug level.
    int ensembleNum_;
    int matrix_min_;      ///< Used for map output
    unsigned int nframes_;///< Number of frames, for normalizing map
    bool first_;          ///< If true use first frame as reference
    bool byResidue_;      ///< If true calculate distances by residue
    bool includeSolvent_; ///< If true include solvent residues
    ImagedAction image_;  ///< Hold imaging-related info/routines.
    AtomMask Mask1_;      ///< First mask in which to search
    AtomMask Mask2_;      ///< Second mask in which to search
    Marray contactList1_; ///< First list of potential contacts to search
    Marray contactList2_; ///< Second list of potential contacts to search
    Iarray contactIdx1_;  ///< Hold atom/residue indices corresponding to contactList1
    Iarray contactIdx2_;  ///< Hold atom/residue indices corresponding to contactList2
    std::string cfile_;   ///< File to write native contact list to.
    DataSet* numnative_;  ///< Hold # of native contacts
    DataSet* nonnative_;  ///< Hold # of non-native contacts
    DataSet* mindist_;    ///< Hold minimum observed distance among contacts
    DataSet* maxdist_;    ///< Hold maximum observed distance among contacts
    DataSet_MatrixDbl* map_; ///< Hold contacts map
    Topology* CurrentParm_;
    Matrix_3x3 ucell_, recip_;
    /// Define contact, either atom or residue number pair.
    typedef std::pair<int,int> contactType;
    /// Define list of contacts, stored as a set for lookup efficiency.
    typedef std::set<contactType> contactListType;
    contactListType nativeContacts_; ///< List of native contacts.
};
#endif
