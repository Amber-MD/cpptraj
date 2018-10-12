#ifndef CPPTRAJ_VERSION_STRING
/* The external versioning for CPPTRAJ currently goes like this:
 *   V<AmberTools version>.<AmberTools patch level>
 */
#define CPPTRAJ_VERSION_STRING "V18.01"
/// This is used in NetcdfFile
#define CPPTRAJ_VERSION_STRLEN 6
#endif
#ifndef CPPTRAJ_INTERNAL_VERSION
/* The internal versioning for CPPTRAJ is supposed to go like this:
 *   V<major>.<minor>.<revision>
 * <major>    - Incremented whenever there is a major API change, e.g.
 *              changing the Action base class, etc.
 * <minor>    - Incremented whenever there are changes to behavior, e.g.
 *              syntax, functionality, or output. The same input should
 *              produce exactly the same output if the major and minor
 *              version numbers are the same.
 * <revision> - All other changes (so at least each PR).
 *
 * Whenever a number that precedes <revision> is incremented, all subsequent
 * numbers should be reset to 0.
 */
#define CPPTRAJ_INTERNAL_VERSION "V4.9.2"
#endif
