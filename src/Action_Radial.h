#ifndef INC_ACTION_RADIAL_H
#define INC_ACTION_RADIAL_H
#include "Action.h"
#include "Histogram.h"
// Class: Radial
/// Calculate the radial distribution (pair correlation) function.
class Radial: public Action {
    Histogram rdf;
    bool noimage;
    int imageType; 
    AtomMask Mask1, Mask2;
    AtomMask OuterMask, InnerMask;
    bool center1;
    bool useVolume;
    double volume;
    char *outfilename;
    double maximum;
    double maximum2;
    double spacing;
    DataSetList rdfdata;
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
