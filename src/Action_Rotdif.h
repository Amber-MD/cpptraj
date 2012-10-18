#ifndef INC_ACTION_ROTDIF_H
#define INC_ACTION_ROTDIF_H
#include "Action.h"
#include "Random.h"
// Class: Action_Rotdif
/// Estimate rotational diffusion tensors from MD simulations
/** To estimate rotational diffusion tensors from MD simulations along the
  * lines described by Wong & Case, (Evaluating rotational diffusion from
  * protein MD simulations, J. Phys. Chem. B 112:6013, 2008) using the 
  * following procedure:
  * - Create random vectors
  * - Calculate rotation matrices via RMS fitting to a reference
  * - Apply the rotation matrices to the random vectors
  * - computes a best-fit rotational diffusion tensor, plus various statistics 
  *   on this. There are two optimizations: one giving the optimal diffusion
  *   tensor in the small-anisotropy limit, and the second optimizing the 
  *   diffusion tensor with the full expression.
  * \author Original code: Vance Wong, George Giambasu
  * \author Adapted by: DRR
  */
class Action_Rotdif: public Action {
  public:
    Action_Rotdif();
    ~Action_Rotdif();

  private:
    int init();
    int setup();
    int action();
    void print();

    int rseed_;       ///< Random seed
    int nvecs_;       ///< Number of random vectors to generate
    double tfac_;     ///< time step
    double ti_;       ///< initial time (ns)
    double tf_;       ///< final time (ns), should be less than ncorr * tfac
    int NmeshPoints_; ///< # cubic spline mesh points 
    int itmax_;       ///< Max number of iterations
    double delmin_;   ///< Convergence criterion
    double d0_;       ///< Initial guess for iso diffusion tensor
    int olegendre_;   ///< order of Legendre polynomial in the correlation function
    int ncorr_;       ///< Max length to compute time correlation fns (# frames)
    double delqfrac_; ///< how to scale simplexes
    double amoeba_ftol_;
    int amoeba_itmax_;
    bool do_gridsearch_;
    bool useMass_;
    bool usefft_;

    // Workspace for LAPACK functions
    double *work_;
    int lwork_;
    double D_tensor_[9];
    double D_XYZ_[3];

    std::string randvecOut_;
    std::string randvecIn_;
    std::string rmOut_;
    std::string deffOut_;
    std::string corrOut_;

    Frame SelectedRef_;
    AtomMask TargetMask_;
    Frame SelectedTgt_;
    CpptrajFile outfile_;

    // Variables used by the random number generator
    Random_Number RNgen_;

    std::vector<double*> Rmatrices_; ///< Store rotation matrices
    double *random_vectors_;         ///< Hold nvecs random vectors
    double *D_eff_;                  ///< Hold calculated effective D values for each vector
    std::vector<double> tau1_;       ///< Hold tau for l=1, full anisotropy
    std::vector<double> tau2_;       ///< Hold tau for l=2, full anisotropy
    std::vector<double> *Tau_;       ///> Hold tau being compared based on olegendre
    std::vector<double> sumc2_;      

    double *randvec();
    int compute_corr(double *, int, int, double *, double *);
    int fft_compute_corr(double*, int, int, double*, int);
    double calcEffectiveDiffusionConst(double );

    int calc_Asymmetric(double *, double *);
    double chi_squared(double *);
    double Amotry(double[][6], double *, double *, int, double); 
    int Amoeba(double[][6], double *);
    int Simplex_min(double*);
    int Grid_search(double *, int);
    int Tensor_Fit(double*);
    int DetermineDeffs();
};
#endif  
