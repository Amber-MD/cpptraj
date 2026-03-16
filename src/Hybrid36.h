#ifndef INC_HYBRID36_H
#define INC_HYBRID36_H
namespace Cpptraj {
/// Implement the PDB hybrid36 extensions
namespace Hybrid36 {
/// Convert integer value to string result
const char* hy36encode(unsigned, int, char*);
/// Convert string to integer result
const char* hy36decode(unsigned, const char*, unsigned, int*);
}
}
#endif
