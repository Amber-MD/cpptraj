#ifndef INC_ACTION_AUTOIMAGE_H
#define INC_ACTION_AUTOIMAGE_H
#include "Action.h"
namespace Image {
  class List_Unit;
}
/// Perform imaging, attempting to keep solute in 1 configuration
class Action_AutoImage : public Action {
  public:
    Action_AutoImage();
    ~Action_AutoImage();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_AutoImage(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    enum Mode { BY_DISTANCE=0, BY_VECTOR, UNSPECIFIED };

    typedef std::vector<Vec3> Varray;

    /// \return vector that can be used to move an imaged fixed molecule back to its reference location
    static inline Vec3 calc_frac_image_vec(Vec3 const&, bool&);
    /// \return vector needed to image a molecule back into the cell
    static inline Vec3 wrap_frac(Vec3 const&, Vec3 const&, Vec3 const&, bool&);
    /// Translate given unit
    static inline void translate_unit(Varray&, Vec3 const&, Unit const&);
    /// \return center of unit according to given array of coords
    static inline Vec3 unit_center(Varray const&, Unit const&);
    /// Move anchor molecule to the center
    Vec3 center_anchor_molecule(ActionFrame& frm, bool, bool) const;
    /// Autoimage using vectors 
    Action::RetType autoimage_by_vector(int, ActionFrame&);
    /// Autoimage using distances
    Action::RetType autoimage_by_distance(int, ActionFrame&);

    Varray RefVecs_;      ///< Fixed molecule reference vectors for BY_VECTOR
    AtomMask anchorMask_; ///< Used to center anchor region.
    std::string anchor_;  ///< Mask expression for anchor region.
    std::string fixed_;   ///< Mask expression for fixed region.
    std::string mobile_;  ///< Mask expression for mobile region.
    int debug_;
    Mode mode_;           ///< Autoimage mode
    bool origin_;         ///< If true imaging occurs w.r.t. coordinate origin.
    bool usecom_;         ///< If true imaging of mobile region uses molecule center.
    bool truncoct_;       ///< If true image into truncated octahedron shape.
    bool useMass_;        ///< If true use center of mass
    bool movingAnchor_;   ///< If true anchor position set to previous fixed molecule
    enum TriclinicArg {OFF, FORCE, FAMILIAR};
    TriclinicArg triclinic_; ///< Determine whether triclinic code should be used.

    Image::List_Unit* fixedList_;  ///< Contain atom indices for fixed elements
    Image::List_Unit* mobileList_; ///< Contain atom indices for mobile elements.
};
#endif
