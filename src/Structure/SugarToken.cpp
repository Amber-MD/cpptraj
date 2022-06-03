#include "SugarToken.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"

using namespace Cpptraj::Structure;

const char* SugarToken::ringstr_[] = {"pyranose", "furanose", "?"};
const char* SugarToken::formstr_[] = {"alpha", "beta", "?"};
const char* SugarToken::chirstr_[] = {"D", "L", "?"};

/** CONSTRUCTOR */
SugarToken::SugarToken() :
  form_(UNKNOWN_FORM),
  chir_(UNKNOWN_CHIR),
  ring_(UNKNOWN_RING)
{}

/** CONSTRUCTOR - name, glycam code, form, chirality, ring type */
SugarToken::SugarToken(std::string const& fn, std::string const& gc,
                                            FormTypeEnum ft, ChirTypeEnum ct, RingTypeEnum rt) :
  name_(fn),
  glycamCode_(gc),
  form_(ft),
  chir_(ct),
  ring_(rt)
{}

/** CONSTRUCTOR - ring type */
SugarToken::SugarToken(RingTypeEnum rt) :
  form_(UNKNOWN_FORM),
  chir_(UNKNOWN_CHIR),
  ring_(rt)
{}

/** \return String containing name, glycam code, and form-chirality-ring type. */
std::string SugarToken::InfoStr() const {
  return std::string("\"" + name_ + "\" " + glycamCode_ + " " +
                     std::string(formstr_[form_]) + "-" +
                     std::string(chirstr_[chir_]) + "-" +
                     std::string(ringstr_[ring_]));
}

/** Set up from line: <res> <code> <form> <chir> <ring> <name>
  *                   0     1      2      3      4      5
  * \return Residue name <res>
  */
std::string SugarToken::SetFromLine(ArgList const& line) {
  const char* lineIn = line.ArgLine();
  if (line.Nargs() != 6) {
    mprinterr("Error: Expected 6 columns, got %i\n"
              "Error: %s\n", line.Nargs(), lineIn);
    return std::string("");
  }
  name_ = line[5];
  glycamCode_ = line[1];
  if (line[2] == "A")
    form_ = ALPHA;
  else if (line[2] == "B")
    form_ = BETA;
  else {
    mprinterr("Error: Unrecognized anomer type: %s\n"
              "Error: Line: %s\n", line[2].c_str(), lineIn);
    return std::string("");
  }
  if (line[3] == "D")
    chir_ = IS_D;
  else if (line[3] == "L")
    chir_ = IS_L;
  else {
    mprinterr("Error: Unrecognized configuration: %s\n"
              "Error: Line: %s\n", line[3].c_str(), lineIn);
    return std::string("");
  }
  if (line[4] == "P")
    ring_ = PYRANOSE;
  else if (line[4] == "F")
    ring_ = FURANOSE;
  else {
    mprinterr("Error: Unrecognized ring: %s\n"
              "Error: Line: %s\n", line[4].c_str(), lineIn);
    return std::string("");
  }

  return line[0];
}

