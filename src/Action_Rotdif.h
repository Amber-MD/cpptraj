#ifndef INC_ACTION_ROTDIF_H
#define INC_ACTION_ROTDIF_H
#include "Action.h"
#include "Random.h"
// Class: Rotdif
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
  * \author Original code: George Giambasu
  * \author Adapted by: DRR
  */
class Rotdif: public Action {
    int rseed;     ///< Random seed
    int nvecs;     ///< Number of random vectors to generate
    double tfac;   ///< time step
    double ti;     ///< initial time (ns)
    double tf;     ///< final time (ns), should be less than ncorr * tfac
    int itmax;     ///< Max number of iterations
    double delmin; ///< Convergence criterion
    double d0;     ///< Initial guess for iso diffusion tensor
    int olegendre; ///< order of Legendre polynomial in the correlation function
    int ncorr;     ///< Max length to compute time correlation fns (# frames)
    double delqfrac; ///< how to scale simplexes
    double amoeba_ftol;
    int amoeba_itmax;
    bool do_gridsearch;

    // Workspace for LAPACK functions
    double *work;
    int lwork;
    double D_tensor[9];
    double D_XYZ[3];

    char *randvecOut;
    char *randvecIn;
    char *rmOut;
    char *deffOut;

    Frame RefFrame;
    AmberParm *RefParm;
    Frame SelectedRef;
    AtomMask RefMask;
    AtomMask TargetMask;
    Frame SelectedTarget;
    CpptrajFile outfile;

    // Variables used by the random number generator
    Random_Number RNgen;

    std::vector<double*> Rmatrices; ///< Store rotation matrices
    double *random_vectors;         ///< Hold nvecs random vectors
    double *D_eff;                  ///< Hold calculated effective D values for each vector
    std::vector<double> tau1;       ///< Hold tau for l=1, full anisotropy
    std::vector<double> tau2;       ///< Hold tau for l=2, full anisotropy
    std::vector<double> *Tau;       ///> Hold tau being compared based on olegendre
    std::vector<double> sumc2;      

    double *randvec();
    int compute_corr(double *, int, int, double *, double *);
    double calcEffectiveDiffusionConst(double );

    int calc_Asymmetric(double *, double *);
    double chi_squared(double *);
    double Amotry(double[][6], double *, double *, int, double); 
    int Amoeba(double[][6], double *);
    int Simplex_min(double*);
    int Grid_search(double *, int);
    int Tensor_Fit(double*);
    int DetermineDeffs();
    
  public:
    Rotdif();
    ~Rotdif();

    int init();
    int setup();
    int action();
    void print();
};
#endif  
