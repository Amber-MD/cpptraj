#ifndef INC_ACTION_DIFFUSION_H
#define INC_ACTION_DIFFUSION_H
#include "Action.h"
class Action_Diffusion : public Action {
  public:
    Action_Diffusion();
  private:
    int init();
    int setup();
    int action();

    Frame initial_;
    std::vector<double> previousx_;
    std::vector<double> previousy_;
    std::vector<double> previousz_;
    bool printIndividual_;
    double time_;
    bool hasBox_;
    std::vector<double> distancex_;
    std::vector<double> distancey_;
    std::vector<double> distancez_;
    std::vector<double> distance_;
    std::vector<double> deltax_;
    std::vector<double> deltay_;
    std::vector<double> deltaz_;
    AtomMask mask_;
    CpptrajFile outputx_;
    CpptrajFile outputy_;
    CpptrajFile outputz_;
    CpptrajFile outputr_;
    CpptrajFile outputa_;
    double boxcrd_[3];    ///< Hold box coordinates each frame
    double boxcenter_[3]; ///< Hold center of box each frame
};
#endif
