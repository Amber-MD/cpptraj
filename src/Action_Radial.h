#ifndef INC_ACTION_RADIAL_H
#define INC_ACTION_RADIAL_H
#include "Action.h"
#include "Histogram.h"
// Class: Radial
/// Calculate the radial distribution (pair correlation) function.
class Radial: public Action {
    int *rdf;    ///< Hold bin counts
    AtomMask Mask1;
    AtomMask Mask2;
    AtomMask OuterMask;
    AtomMask InnerMask;
    bool center1;
    bool useVolume;
    double volume;
    char *outfilename;
    double maximum;
    double maximum2;
    double spacing;
    double one_over_spacing;
    int numBins;
    int numthreads;
    int **rdf_thread;
    //DataSetList rdfdata;     
    int numFrames;
    int numDistances;
    double density;
  public:
    Radial();
    ~Radial();

    int init();
    int setup();
    int action();
    void print();
};
#endif  
