#ifndef INC_ACTION_RADIAL_H
#define INC_ACTION_RADIAL_H
/// Class: Radial
/// Used to calculate the radial distribution function (pair correlation)
/// of atom(s) in mask0 to atoms in mask2.
#include "Action.h"
#include "Histogram.h"
class Radial: public Action {
    Histogram rdf;
    bool noimage;
    bool hasSecondMask;
    int imageType; 
    AtomMask Mask1, Mask2;
    char *outfilename;
    double maximum;
    double maximum2;
    double spacing;
    DataSetList rdfdata;
    int numFrames;
    int numDistances;
  public:
    Radial();
    ~Radial();

    int init();
    int setup();
    int action();
    void print();
};
#endif  
