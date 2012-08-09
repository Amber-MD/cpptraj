#ifndef INC_ACTION_RADIAL_H
#define INC_ACTION_RADIAL_H
#include "Action.h"
#include "Histogram.h"
// Class: Action_Radial
/// Calculate the radial distribution (pair correlation) function.
class Action_Radial: public Action {
  public:
    Action_Radial();
    ~Action_Radial();

    void print();
  private:
    int init();
    int setup();
    int action();

    int* RDF_;                ///< Hold bin counts.
    int **rdf_thread_;        ///< Hold bin count on each thread.
    AtomMask Mask1_;          ///< Atoms to calculate RDF for.
    AtomMask Mask2_;          ///< Optional mask to calc RDF to atoms in Mask1.
    AtomMask OuterMask_;      ///< Mask with the most atoms.
    AtomMask InnerMask_;      ///< Mask with the fewest atoms.
    bool center1_;            ///< If true calculate RDF of atoms in Mask2 to COM of Mask1.
    bool useVolume_;          ///< If true normalize based on input volume.
    double volume_;           ///< Hold sum of volume for averaging.
    char *outfilename_;       ///< Output file name.
    double maximum2_;         ///< Largest distance squared that can be binned.
    double spacing_;          ///< Bin spacing.
    double one_over_spacing_; ///< 1/spacing, used to avoid man division ops.
    int numBins_;             ///< The number of bins.
    int numthreads_;          ///< Number of threads.
    int numFrames_;           ///< Number of frames for which RDF is calcd.
    int numDistances_;        ///< Number of distances binned, only informational.
    double density_;          ///< Particle density (molecules/Ang^3).
};
#endif  
