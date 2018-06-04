#ifndef INC_CHARMMPARAM_H
#define INC_CHARMMPARAM_H
/// Used to read in CHARMM parameters from CHARMM parameter file.
class CharmmParam {
  public:
    CharmmParam() {}
    int ReadParams(std::string const&);
  private:
};
#endif
