#ifndef INC_NC_ROUTINES_H
#define INC_NC_ROUTINES_H
#include <string>
// Routines that can exist without NetCDF
namespace NC {
  /// Recognized NetCDF Conventions types (previously NetcdfFile::NCTYPE). MUST MATCH NC_ConventionsStr_ in NC_Routines.cpp
  enum ConventionsType {
    NC_AMBERTRAJ = 0,  ///< Amber NetCDF trajectory.
    NC_AMBERRESTART,   ///< Amber NetCDF restart.
    NC_AMBERENSEMBLE,  ///< Amber NetCDF ensemble trajectory.
    NC_CPPTRAJCMATRIX, ///< Cpptraj NetCDF clustering pairwise distance matrix.
    NC_CPPTRAJDATA,    ///< Cpptraj NetCDF data file.
    NC_UNKNOWN
  };
  /// \return NetCDF conventions of file with given name.
  ConventionsType GetConventions(std::string const&);
}
#ifdef BINTRAJ
// Routines that require NetCDF
namespace NC {
  /// \return true if given code is error and print message, false otherwise.
  bool CheckErr(int);
  /// \return Text for attribute of given variable ID.
  std::string GetAttrText(int, int, const char*);
  /// \return Text for given global attribute.
  std::string GetAttrText(int, const char*);
  /// \return dimension ID of given attribute and set dimension length.
  int GetDimInfo(int, const char*, unsigned int&);
  // FIXME This version here for backwards compatibility.
  int GetDimInfo(int, const char*, int&);
  /// Print debug info to STDOUT for given NetCDF id.
  void Debug(int);
  /// \return NetCDF conventions for given NetCDF id, silent if not present/recognized.
  ConventionsType GetConventions(int);
  /// \return NetCDF conventions for given NetCDF id, potentially verbose.
  ConventionsType GetConventions(int, bool);
  /// \return Character string corresponding to given ConventionsType
  const char* conventionsStr(ConventionsType);
  /// Put specified conventions string to given NetCDF id.
  int PutConventions(int, ConventionsType);
}
#endif /* BINTRAJ */
#endif
