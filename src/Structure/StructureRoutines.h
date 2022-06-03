#ifndef INC_STRUCTURE_STRUCTUREROUTINES_H
#define INC_STRUCTURE_STRUCTUREROUTINES_H
class Atom;
class NameType;
class Residue;
namespace Cpptraj {
namespace Structure {

void SetStructureDebugLevel(int); // TODO use debug level everywhere throughout Structure
void ChangeResName(Residue&, NameType const&);
void ChangeAtomName(Atom&, NameType const&);

}
}
#endif
