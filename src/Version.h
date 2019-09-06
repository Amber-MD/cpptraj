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
#define CPPTRAJ_INTERNAL_VERSION "V4.19.0"
/// PYTRAJ relies on this
#define CPPTRAJ_VERSION_STRING CPPTRAJ_INTERNAL_VERSION
#endif
