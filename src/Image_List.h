#ifndef INC_IMAGE_LIST_H
#define INC_IMAGE_LIST_H
#include <string>
class Topology;
class Frame;
class Vec3;
class Unit;
namespace Image {
/// Abstract base class for holding entities to be imaged
class List {
  public:
    List() {}
    virtual ~List() {}
    /// \return Number of entities to be imaged
    virtual unsigned int nEntities() const = 0;
    /// Set up the list according to Topology and mask expression
    virtual int SetupList(Topology const&, std::string const&) = 0;
    /// \return Coordinates of specified entity
    virtual Vec3 GetCoord(unsigned int, Frame const&) const = 0;
    /// Do translation of specified entity according to given translation vector.
    virtual void DoTranslation(Frame&, unsigned int, Vec3 const&) const = 0;
    /// \return Unit containing all entities in the list
    virtual Unit AllEntities() const = 0;
    /// Print entities to STDOUT
    virtual void PrintEntities() const = 0;
};
}
#endif
