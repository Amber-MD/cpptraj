#ifndef INC_STRUCTURE_STRUCTUREENUM_H
#define INC_STRUCTURE_STRUCTUREENUM_H
namespace Cpptraj {
namespace Structure {

enum ChiralType { CHIRALITY_ERR = 0, IS_S, IS_R, IS_UNKNOWN_CHIRALITY };
/// \return String corresponding to ChiralType
const char* chiralStr(ChiralType);

}
}
#endif
