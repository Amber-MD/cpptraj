#ifndef INC_ACTION_ROTDIF_H
#define INC_ACTION_ROTDIF_H
#include "Action.h"
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

    Frame RefFrame;
    AmberParm *RefParm;
    Frame SelectedRef;
    AtomMask RefMask;
    AtomMask TargetMask;
    Frame SelectedTarget;

    // Variables used by the random number generator
    //int iff;
    //int inext;
    //int inextp;
    //float ma[55];

    //double random_number();
    double *randvec();

    std::vector<double*> Rmatrices; ///< Store rotation matrices
    
  public:
    Rotdif();
    ~Rotdif();

    int init();
    int setup();
    int action();
    void print();
};
#endif  
