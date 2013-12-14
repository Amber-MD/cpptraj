#ifndef INC_ACTION_NATIVECONTACTS_H
#define INC_ACTION_NATIVECONTACTS_H
#include <set>
#include "Action.h"
#include "ImagedAction.h"
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
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    Marray SetupList(Topology const&, AtomMask const&) const;
    int SetupContactLists(Topology const&, Frame const&);
    double GetMinMax(AtomMask const&, AtomMask const&, Frame const&,
                     double&, double&) const;

    double distance_;    ///< Cutoff distance
    int debug_;          ///< Action debug level.
    bool first_;         ///< If true use first frame as reference
    bool byResidue_;     ///< If true calculate distances by residue
    ImagedAction image_; ///< Hold imaging-related info/routines.
    AtomMask Mask1_;
    AtomMask Mask2_;
    Marray contactList1_;
    Marray contactList2_;
    DataSet* numnative_;
    DataSet* nonnative_;
    DataSet* mindist_;
    DataSet* maxdist_;
    Topology* CurrentParm_;
    Matrix_3x3 ucell_, recip_;
    /// Define contact, either atom or residue number pair.
    typedef std::pair<int,int> contactType;
    /// Define list of contacts, stored as a multimap for efficiency.
    typedef std::set<contactType> contactListType;
    contactListType nativeContacts_; ///< List of native contacts.
    
};
#endif
