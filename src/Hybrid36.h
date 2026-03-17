#ifndef INC_HYBRID36_H
#define INC_HYBRID36_H
namespace Cpptraj {
/// Implement the PDB hybrid36 extensions
class Hybrid36 {
  public:
    Hybrid36() {}
    /// Encode to H36: width (4/5), value to convert, resulting string
    static int Encode(unsigned, int, char*);
    /// Decode from H36: width (4/5), string to convert (size must be width), resulting int
    static int Decode(unsigned, const char*, int&);
    /// Convert integer value to string result
  private:
    static const char* hy36encode(unsigned, int, char*);
    /// Convert string to integer result
    static const char* hy36decode(unsigned, const char*, unsigned, int*);
};
}
#endif
