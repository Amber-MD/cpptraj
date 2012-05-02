#ifndef INC_ACTIONREFERENCE_H
#define INC_ACTIONREFERENCE_H
#include "TopologyList.h"
#include "Frame.h"
#include "TrajectoryFile.h"
#include "FrameList.h"
/// Hold reference frame/trajectory information for an action
class ActionReference {
  public:
   enum RefModeType {
      NOREF=0, FIRST, REF, REFTRAJ, FIRST_DONE
    }; 

    ActionReference();

    RefModeType RefMode() { return refMode_; }

    void ForbidRefMode( RefModeType );
    void SetFirst(bool,char*,bool);
    void RefInfo();
    int RefNselected();
    Topology *GetRefParm();

    int RefInit(bool, bool, char*, ArgList&, FrameList*, TopologyList*, double*);
    int RefSetup(Topology*);
    void RefAction(Frame *,double*);
  protected:
    Frame SelectedRef_;
    Frame RefFrame_;
  private:
    std::vector<bool> allowed_;
    RefModeType refMode_;
    AtomMask RefMask_;
    TrajectoryFile RefTraj_;
    Topology *RefParm_;
    bool fitReference_;
    bool useRefMass_;
    std::string refname_;
    bool refParmSetup_;

    int SetRefMask();
    void SetRefStructure(Frame&,double*,const char*);
    int SetRefTraj(char*);
};
#endif
