#ifndef INC_ACTION_DIFFUSION_H
#define INC_ACTION_DIFFUSION_H
#include "Action.h"
#include "ImagedAction.h"
class Action_Diffusion : public Action {
  public:
    Action_Diffusion();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Diffusion(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();

    typedef DataSetList::DataListType Dlist;
    typedef std::vector<double> Darray;

    inline void LoadInitial(Frame const&);
    void CalcDiffForSet(unsigned int&, Dlist const&, int, std::string const&) const;
    void CalcDiffusionConst(unsigned int&, DataSet*, int, std::string const&) const;

    ImagedAction image_; ///< Imaging routines
    Frame initial_;   ///< Initial frame (all atoms)
    Darray previous_; ///< Previous coordinates (selected atoms)
    DataSet* avg_x_;  ///< Hold average diffusion in X direction each frame
    DataSet* avg_y_;  ///< Hold average diffusion in Y direction each frame
    DataSet* avg_z_;  ///< Hold average diffusion in Z direction each frame
    DataSet* avg_r_;  ///< Hold average MSD each frame
    DataSet* avg_a_;  ///< Hold average distance each frame
    Dlist atom_x_; ///< Hold individual distances in X direction each frame
    Dlist atom_y_; ///< Hold individual distances in Y direction each frame
    Dlist atom_z_; ///< Hold individual distances in Z direction each frame
    Dlist atom_r_; ///< Hold individual MSDs each frame
    Dlist atom_a_; ///< Hold individual distances each frame
    bool printIndividual_; ///< If true print diffusion for individual atoms
    bool calcDiffConst_;   ///< If true calculate diffusion const from MSD vs time
    double time_;          ///< Time step between frames
    DataSet* diffConst_;   ///< Hold diffusion constants.
    DataSet* diffLabel_;   ///< Hold diffusion constant labels.
    DataSet* diffSlope_;   ///< Hold MSD vs time line slopes.
    DataSet* diffInter_;   ///< Hold MSD vs time line intercepts.
    DataSet* diffCorrl_;   ///< Hold MSD vs time line correlation.
    int debug_;
    Darray delta_;         ///< Hold current distances from initial frame for selected atoms
    AtomMask mask_;        ///< Selected atoms
    DataFile* outputx_;
    DataFile* outputy_;
    DataFile* outputz_;
    DataFile* outputr_;
    DataFile* outputa_;
    DataFile* diffout_;
    Vec3 boxcenter_; ///< Hold center of box each frame
    DataSetList* masterDSL_;
    std::string dsname_;
    Dimension Xdim_;
#   ifdef MPI
    Parallel::Comm trajComm_;
#   endif
};
#endif
