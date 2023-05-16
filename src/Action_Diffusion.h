#ifndef INC_ACTION_DIFFUSION_H
#define INC_ACTION_DIFFUSION_H
#include "Action.h"
#include "ImageOption.h"
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
#   ifdef MPI
    int SyncAction();

    void average_multiple_time_origins(DataSet*, std::vector<int> const&, int) const;
#   endif

    typedef DataSetList::DataListType Dlist;
    typedef std::vector<Vec3> Varray;

    inline void LoadInitial(Frame const&);
    void CalcDiffForSet(unsigned int&, Dlist const&, int, std::string const&) const;
    void CalcDiffusionConst(unsigned int&, DataSet*, int, std::string const&) const;

    ImageOption imageOpt_; ///< Used to determine if imaging should be used.
    Frame initial_;   ///< Initial frame (all atoms)
    Varray previousFrac_; ///< Previous fractional coordinates for imaging (selected atoms)
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
    AtomMask mask_;        ///< Selected atoms
    DataFile* outputx_;
    DataFile* outputy_;
    DataFile* outputz_;
    DataFile* outputr_;
    DataFile* outputa_;
    DataFile* diffout_;
    DataSetList* masterDSL_;
    std::string dsname_;
    Dimension Xdim_;
    DataSet* avgucell_;   ///< Hold average unit cell parameters for removing box fluctuations
    Box avgbox_;          ///< Hold average box for removing box fluctuations
#   ifdef MPI
    Parallel::Comm trajComm_;
    bool multipleTimeOrigins_; ///< If true, parallel diffusion calc with imaging, multiple time originsi
#   endif
};
#endif
