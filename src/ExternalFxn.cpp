#include "Random.h"

/** Allow external programs to change the default random number generator. 
  * Wrapping it here avoids e.g. pytraj having to use namespace Cpptraj
  * which will fail since I made both a Cpptraj namespace and class.
  * Its primary function right now is to allow pytraj tests to set
  * the default RNG back to Marsaglia.
  * It also avoids having to import the Random_Number::RngType enum.
  * TODO this can be deprecated when the Cpptraj class is renamed.
  */
void EXT_SetDefaultRng(int rtypeIn) {
  Random_Number::RngType rt = (Random_Number::RngType)rtypeIn;
  Random_Number::SetDefaultRng( rt );
}
