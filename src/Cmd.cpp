#include "Cmd.h"

void Cmd::Clear() { if (object_ != 0) delete object_; }

bool Cmd::KeyMatches(const char* keyIn) const {
  for (key_iterator key = keywords_.begin(); key != keywords_.end(); ++key)
    if ( key->compare( keyIn ) == 0 ) return true;
  return false;
}
