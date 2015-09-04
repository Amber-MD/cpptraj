#ifndef INC_TRAJIN_H
#define INC_TRAJIN_H
#include "InputTrajCommon.h"
/// Read in 1 frame at a time.
class Trajin {
  public:
    Trajin() : debug_(0) {}
    virtual ~Trajin() {}
    virtual int SetupTrajRead(std::string const&, ArgList&, Topology*) = 0;
    virtual int ReadTrajFrame(int, Frame&) = 0;
    virtual int BeginTraj() = 0;
    virtual void EndTraj() = 0;
    virtual void PrintInfo(int) const = 0;
    virtual CoordinateInfo const& TrajCoordInfo() const = 0;

    inline int GetNextFrame(Frame&);

    InputTrajCommon const& Traj() const { return traj_; }

    void SetDebug(int d)                { debug_ = d;         }
  protected:
    InputTrajCommon& SetTraj() { return traj_; }

    int debug_;
  private:
    InputTrajCommon traj_;
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
int Trajin::GetNextFrame(Frame& frameIn) {
  if (traj_.Counter().CheckFinished()) return 0;
  if (ReadTrajFrame( traj_.Counter().Current(), frameIn )) return 0;
  traj_.Counter().UpdateCounters();
  return 1;
}
#endif 
