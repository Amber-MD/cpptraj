#ifndef INC_ACTION_VOLMAP_H
#define INC_ACTION_VOLMAP_H
#include "Action.h"
#include "Grid.h"
class Action_Volmap : public Action {
  public:
    Action_Volmap();
   ~Action_Volmap();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Volmap(); }
    static void Help();
    static void RawHelp();

    void Print();
  private:
    /** \brief grid resolutions */
    double dx_, dy_, dz_;
    /** \brief X, Y, and Z-coordinates of the grid center */
    double xcenter_, ycenter_, zcenter_;
    /** \brief minimum values in the x-, y-, and z-dimensions */
    double xmin_, ymin_, zmin_;
    /** \brief number of frames we analyzed so we can average at the end */
    int Nframes_;
    /** \brief mask to center the grid on */
    AtomMask centermask_;
    /** \brief mask of atoms to grid */
    AtomMask densitymask_;
    /** \brief the grid we are using */
    Grid grid_;
    /** \brief the grid with the peak locations */
    Grid peakgrid_;
    /** \brief file name with the peak locations as Carbons in XYZ file format
     */
    std::string peakfilename_;
    /// \brief The value below which to ignore all peaks
    double peakcut_;
    /** \brief file name for the output density */
    std::string filename_;
    /** \brief the atomic radii of each atom in the gridded selection */
    float *halfradii_;
    /** \brief the clearance between the edges of our grid and centermask_ */
    double buffer_;
    /** \brief the scaling factor to divide all radii by */
    double radscale_;
    /** \brief size of the grid in x,y,z dims */
    double xsize_, ysize_, zsize_;
    /** \brief declares if we are outputting a DX-formatted density file */
    bool dxform_;
    /** \brief gets the LJ radius for a given atom from a topology */
    double GetRadius_(Topology const&, int);

    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
};
#endif
