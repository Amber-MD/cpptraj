#ifndef INC_ACTION_STFC_DIFFUSION_H
#define INC_ACTION_STFC_DIFFUSION_H
#include "Action.h"
/** \author Hannes H. Loeffler
  * \author C++ adaptation by Daniel R. Roe
  */
class Action_STFC_Diffusion : public Action {
  public:
    Action_STFC_Diffusion();
    ~Action_STFC_Diffusion();
  private:
    int init();
    int setup();
    int action();

    void calculateMSD(double*,double*);

    bool printDistances_; // iarg1
    enum CalcType { DEFAULT = 0, COM, DIST };
    CalcType calcType_; // iarg2
    enum DirectionType { DX = 0, DY, DZ, DXY, DXZ, DYZ, DXYZ };
    DirectionType direction_; // iarg3
    AtomMask mask_;
    AtomMask mask2_;
    CpptrajFile output_;
    CpptrajFile outputad_;
    CpptrajFile outputnw_;
    std::string outputAverDist_;
    std::string outputNumWat_;
    double time_;
    double lowerCutoff_;
    double upperCutoff_;

    bool hasBox_;
    int Ninitial_atm_;
    int Ninitial_crd_;
    double* initialxyz_;
    //std::vector<double> initialx_;
    //std::vector<double> initialy_;
    //std::vector<double> initialz_;
    double* distancexyz_;
    //std::vector<double> distancex_;
    //std::vector<double> distancey_;
    //std::vector<double> distancez_;
    std::vector<double> distance_;
    double* deltaxyz_;
    //std::vector<double> deltax_;
    //std::vector<double> deltay_;
    //std::vector<double> deltaz_;
    double* previousxyz_;
    //std::vector<double> previousx_;
    //std::vector<double> previousy_;
    //std::vector<double> previousz_;

    std::vector<double> dSum1_;
    std::vector<double> dSum2_;
    int elapsedFrames_;
};
#endif    
