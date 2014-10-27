#ifndef INC_ACTION_NMRRST_H
#define INC_ACTION_NMRRST_H
#include <map>
#include "Action.h"
#include "ImagedAction.h"
#include "BufferedLine.h"
// Class: Action_NMRrst
class Action_NMRrst: public Action {
  public:
    Action_NMRrst();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_NMRrst(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    int ReadNmrRestraints( std::string const& );
    int ReadXplor( BufferedLine& );
    int ReadAmber( BufferedLine& );
    /// Type to hold NOE data.
    struct noeDataType {
      int resNum1_;     ///< First residue number.
      int resNum2_;     ///< Second residue number.
      std::string aName1_; ///< First atom name.
      std::string aName2_; ///< Second atom name.
      AtomMask dMask1_; ///< First mask for distance pair.
      AtomMask dMask2_; ///< Second mask for distance pair.
      double bound_;    ///< Lower bound.
      double boundh_;   ///< Upper bound.
      double rexp_;     ///< Expected distance
      DataSet* dist_;   ///< Distance DataSet.
      bool active_;     ///< True if NOE was properly set up.
    };
    typedef std::vector<noeDataType> noeArray;
    noeArray NOEs_;

    typedef std::vector<int> Iarray;

    /// Potential NOE site.
    class Site {
      public:
        Site() : resNum_(-1) {}
        Site(int r, Iarray const& i) : resNum_(r), indices_(i), shortestCount_(i.size(), 0) {}
        int ResNum() const { return resNum_; }
        typedef Iarray::const_iterator Idx_it;
        Idx_it IdxBegin() const { return indices_.begin(); }
        Idx_it IdxEnd() const { return indices_.end(); }
      private:
        int resNum_; ///< Site residue number.
        Iarray indices_; ///< Site atom indices.
        Iarray shortestCount_; ///< # times atom was part of shortest distance.
    };
    typedef std::vector<Site> SiteArray;
    SiteArray potentialSites_;

    /// Used to map NOEs to unique values, res1 always < res2. 
    typedef std::pair<int,int> Ptype;

    /// NOE between two sites.
    class NOEtype {
      public:
        NOEtype() : dist_(0) {}
      private:
        Site site1_; ///< First site, lower resNum
        Site site2_; ///< Second site, higher resNum
        DataSet* dist_; ///< Distance data.
    };
    typedef std::map<Ptype, NOEtype> NOEmap;
    NOEmap FoundNOEs_;
    
    ImagedAction Image_;
    bool useMass_;
    bool findNOEs_;
    int resOffset_;
    int debug_;
};
#endif
