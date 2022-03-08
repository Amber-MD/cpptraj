# Parts of the build system that happen after enable_language()

include(TargetArch)
include(CheckFunctionExists)
include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)
include(CheckLinkerFlag)
include(TryLinkLibrary)

include(LibraryTracking)
include(BuildReport)
include(CompilationOptions)
include(CopyTarget)
include(LibraryUtils)
