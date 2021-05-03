#include "Pucker_PuckerToken.h"

using namespace Cpptraj;

Pucker::PuckerToken::PuckerToken() {}

Pucker::PuckerToken::PuckerToken(NameArray const& namesIn) :
  atomNames_(namesIn)
{}
