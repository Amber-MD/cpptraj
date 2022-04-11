
// Header to control exporting of symbols under Windows.
// Currently, for pytraj to work, only global variables need explicit export specifications.
// (see https://github.com/Amber-MD/cmake-buildscripts/issues/29)

// NOTE: if you are on Windows and linking to a dynamically compiled cpptraj, then
// you must define CPPTRAJ_USE_DLL for it to link properly.

#ifndef INC_SYMBOL_EXPORTING_H
#define INC_SYMBOL_EXPORTING_H

#ifdef _WIN32
# ifdef libcpptraj_EXPORTS // this definition is auto-added by CMake
#  define CPPTRAJ_EXPORT __declspec(dllexport)
# elif defined(CPPTRAJ_USE_DLL)
#  define CPPTRAJ_EXPORT __declspec(dllimport)
# else
#  define CPPTRAJ_EXPORT
# endif
#else
# define CPPTRAJ_EXPORT
#endif

#endif
