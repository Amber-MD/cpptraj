#include "CpptrajState.h"

/** Select lists from ArgList */
std::vector<bool> CpptrajState::ListsFromArg( ArgList& argIn, bool allowEnableAll ) {
  std::vector<bool> enabled( (int)N_LISTS );
  enabled[L_ACTION]   = argIn.hasKey("actions");
  enabled[L_TRAJIN]   = argIn.hasKey("trajin");
  enabled[L_REF]      = argIn.hasKey("ref");
  enabled[L_TRAJOUT]  = argIn.hasKey("trajout");
  enabled[L_PARM]     = argIn.hasKey("parm");
  enabled[L_ANALYSIS] = argIn.hasKey("analysis");
  enabled[L_DATAFILE] = argIn.hasKey("datafile");
  enabled[L_DATASET]  = argIn.hasKey("dataset");
  if (!allowEnableAll) return enabled;
  // If nothing is enabled, set all enabled
  bool nothing_enabled = true;
  for (std::vector<bool>::iterator en = enabled.begin(); en != enabled.end(); ++en)
    if (*en) {
      nothing_enabled = false;
      break;
    }
  if (nothing_enabled) enabled.assign( (int)N_LISTS, true );
  return enabled;
}

/** List all members of specified lists */
int CpptrajState::ListAll( ArgList& argIn ) {
  std::vector<bool> enabled = ListsFromArg( argIn, true );
  if ( enabled[L_ACTION]   ) actionList_.List();
  if ( enabled[L_TRAJIN]   ) trajinList_.List();
  if ( enabled[L_REF]      ) refFrames_.List();//{mprintf("\nREFERENCE COORDS:\n");refFrames_.List();}
  if ( enabled[L_TRAJOUT]  ) trajoutList_.List();//{mprintf("\nOUTPUT TRAJECTORIES:\n");trajoutList_.List();}
  if ( enabled[L_PARM]     ) parmFileList_.List();
  if ( enabled[L_ANALYSIS] ) analysisList_.List();
  if ( enabled[L_DATAFILE] ) DFL_.List();
  if ( enabled[L_DATASET]  ) DSL_.List();//{mprintf("\nDATASETS:\n");DSL_.List();}
  return 0;
}

/** Set debug level of specified lists. */
int CpptrajState::SetListDebug( ArgList& argIn ) {
  std::vector<bool> enabled = ListsFromArg( argIn, true );
  debug_ = argIn.getNextInteger(0);
  if ( enabled[L_ACTION]   ) actionList_.SetDebug( debug_ );
  if ( enabled[L_TRAJIN]   ) trajinList_.SetDebug( debug_ );
  if ( enabled[L_REF]      ) refFrames_.SetDebug( debug_ );
  if ( enabled[L_TRAJOUT]  ) trajoutList_.SetDebug( debug_ );
  if ( enabled[L_PARM]     ) parmFileList_.SetDebug( debug_ );
  if ( enabled[L_ANALYSIS] ) analysisList_.SetDebug( debug_ );
  if ( enabled[L_DATAFILE] ) DFL_.SetDebug( debug_ );
  if ( enabled[L_DATASET]  ) DSL_.SetDebug( debug_ );
  return 0;
}

/** Clear specified lists */
int CpptrajState::ClearList( ArgList& argIn ) {
  std::vector<bool> enabled = ListsFromArg( argIn, argIn.hasKey("all") );
  if ( enabled[L_ACTION]   ) actionList_.Clear();
  if ( enabled[L_TRAJIN]   ) trajinList_.Clear();
  if ( enabled[L_REF]      ) refFrames_.Clear();
  if ( enabled[L_TRAJOUT]  ) trajoutList_.Clear();
  if ( enabled[L_PARM]     ) parmFileList_.Clear();
  if ( enabled[L_ANALYSIS] ) analysisList_.Clear();
  if ( enabled[L_DATAFILE] ) DFL_.Clear();
  if ( enabled[L_DATASET]  ) DSL_.Clear();
  return 0;
}
