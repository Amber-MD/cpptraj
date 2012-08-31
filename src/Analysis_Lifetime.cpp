#include "Analysis_Lifetime.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Analysis_Lifetime::Analysis_Lifetime() :
  windowSize_(0),
  cut_(0.5),
  averageonly_(false),
  cumulative_(false)
{}

/** Usage: lifetime [out <filename>] <dsetarg0> [ <dsetarg1> ... ]
  *                 [window <windowsize> [name <setname>]] [averageonly]
  *                 [cumulative] [cut <cutoff>]
  */
int Analysis_Lifetime::Setup( DataSetList* datasetlist ) {
  // Get Keywords
  outfilename_ = analyzeArgs_.GetStringKey("out");
  std::string setname_ = analyzeArgs_.GetStringKey("name");
  windowSize_ = analyzeArgs_.getKeyInt("window", -1);
  averageonly_ = analyzeArgs_.hasKey("averageonly");
  cumulative_ = analyzeArgs_.hasKey("cumulative");
  cut_ = analyzeArgs_.getKeyDouble("cut", 0.5);
  // Select datasets
  while (analyzeArgs_.ArgsRemain())
    inputDsets_ += datasetlist->GetMultipleSets( analyzeArgs_.GetStringNext() );
  if (inputDsets_.empty()) {
    mprinterr("Error: lifetime: No data sets selected.\n");
    return 1;
  }

  // Create output datasets
  if ( windowSize_ != -1) {
    if (setname_.empty()) 
      setname_ = datasetlist->GenerateDefaultName( "lifetime" );
    int didx = 0;
    for (DataSetList::const_iterator set = inputDsets_.begin(); set != inputDsets_.end(); ++set)
    {
      DataSet* outSet = datasetlist->AddSetIdx( DataSet::FLOAT, setname_, didx );
      if (outSet==NULL) {
        mprinterr("Error: lifetime: Could not allocate output set for %s\n", 
                  (*set)->Legend().c_str());
        return 1;
      }
      outSet->SetLegend( (*set)->Legend() );
      outputDsets_.push_back( outSet );
      if (!averageonly_) {
        // MAX
        // FIXME: CHeck for NULLS
        outSet = datasetlist->AddSetIdxAspect( DataSet::INT, setname_, didx, "max" );
        outSet->SetLegend( (*set)->Legend() );
        maxDsets_.push_back( outSet );
        // AVG
        outSet = datasetlist->AddSetIdxAspect( DataSet::FLOAT, setname_, didx, "avg" );
        outSet->SetLegend( (*set)->Legend() );
        avgDsets_.push_back( outSet );
      }
      ++didx;
    }
  }

  if (!averageonly_)
    mprintf("    LIFETIME: Calculating average lifetime using a cutoff of %f", cut_);
  else
    mprintf("    LIFETIME: Calculating only averages");
  mprintf(" of data in %i sets:\n", inputDsets_.size());
  inputDsets_.Info();
  if (windowSize_ != -1) {
    mprintf("\tAverage of data over windows will be saved to sets named %s\n",
            setname_.c_str());
    mprintf("\tWindow size for averaging: %i\n", windowSize_);
    if (cumulative_)
      mprintf("\tCumulative averages will be saved.\n");
    if (!outfilename_.empty()) {
      mprintf("\tOutfile: %s", outfilename_.c_str());
      if (!averageonly_)
        mprintf(", max.%s, avg.%s", outfilename_.c_str(), outfilename_.c_str());
      mprintf("\n");
    }
  }

  return 0;
}

// Analysis_Lifetime::Analyze()
int Analysis_Lifetime::Analyze() {
  float favg;
  std::vector<DataSet*>::iterator outSet = outputDsets_.begin();
  std::vector<DataSet*>::iterator maxSet = maxDsets_.begin();
  std::vector<DataSet*>::iterator avgSet = avgDsets_.begin();
  for (DataSetList::const_iterator inSet = inputDsets_.begin(); 
                                   inSet != inputDsets_.end(); ++inSet)
  {
    mprintf("\t\tCalculating lifetimes for set %s\n", (*inSet)->Legend().c_str());
    // Loop over all values in set.
    int setSize = (*inSet)->Size();
    double sum = 0;
    int windowcount = 0; // Used to trigger averaging
    int Ncount = 0;      // Used in averaging; if !cumulative, == windowcount
    int frame = 0;       // Frame to add data at.
    int currentLifetimeCount = 0; // # of frames value has been present this lifetime
    int maximumLifetimeCount = 0; // Max observed lifetime
    int Nlifetimes = 0;           // # of separate lifetimes observed
    int sumLifetimes = 0;         // sum of lifetimeCount for each lifetime observed
    for (int i = 0; i < setSize; ++i) {
      double dval = (*inSet)->Dval(i);
      //mprintf("\t\t\tValue[%i]= %.2f", i,dval);
      if (averageonly_) 
        // Average only
        sum += dval;
      else {
        // Lifetime calculation
        if ( dval > cut_ ) {
          // Value is present at time i
          ++sum;
          ++currentLifetimeCount;
          //mprintf(" present; sum=%i LC=%i\n", sum, currentLifetimeCount);
        } else {
          //mprintf(" not present");
          // Value is not present at time i
          if (currentLifetimeCount > 0) {
            if (currentLifetimeCount > maximumLifetimeCount)
              maximumLifetimeCount = currentLifetimeCount;
            sumLifetimes += currentLifetimeCount;
            ++Nlifetimes;
            //mprintf("; LC=%i maxLC=%i NL=%i\n", currentLifetimeCount, maximumLifetimeCount, Nlifetimes);
            currentLifetimeCount = 0;
          }
          //mprintf("\n");
        }
      }
      //sum += (*inSet)->Dval(i);
      ++Ncount;
      ++windowcount;
      if (windowcount == windowSize_) {
        double windowavg = sum / (double)Ncount;
        float fval = (float)windowavg;
        (*outSet)->Add( frame, &fval );
        if (!averageonly_) {
          // Store lifetime information for this window
          // Update current lifetime total
          if (currentLifetimeCount > 0) {
            if (currentLifetimeCount > maximumLifetimeCount)
              maximumLifetimeCount = currentLifetimeCount;
            sumLifetimes += currentLifetimeCount;
            ++Nlifetimes;
          }
          // If Nlifetimes is 0 then value was never present. 
          if (Nlifetimes == 0) 
            favg = 0.0;
          else
            favg = (float)sumLifetimes / (float)Nlifetimes;
          //mprintf("\t\t\t[%i]Max lifetime observed: %i frames\n", frame,maximumLifetimeCount);
          //mprintf("\t\t\t[%i]Avg lifetime: %f frames\n", frame, favg);
          (*maxSet)->Add( frame, &maximumLifetimeCount );
          (*avgSet)->Add( frame, &favg );
        }
        frame += windowcount;
        // Window counter is always reset
        windowcount = 0;
        if (!cumulative_) {
          // Reset average counters
          sum = 0;
          Ncount = 0;
          // Reset lifetime counters
          currentLifetimeCount = 0;
          maximumLifetimeCount = 0;
          Nlifetimes = 0;
          sumLifetimes = 0;
        }
      }
    }
    // Print lifetime information if no window
    if ( !averageonly_ && windowSize_ == -1 ) {
      // Update current lifetime total
      if (currentLifetimeCount > 0) {
        if (currentLifetimeCount > maximumLifetimeCount)
          maximumLifetimeCount = currentLifetimeCount;
        sumLifetimes += currentLifetimeCount;
        ++Nlifetimes;
      }
      // If Nlifetimes is 0 then value was never present. 
      if (Nlifetimes == 0) 
        favg = 0.0;
      else
        favg = (float)sumLifetimes / (float)Nlifetimes;
      mprintf("\t\t\tMax lifetime observed: %i frames\n", maximumLifetimeCount);
      //mprintf("\t\tSumLifeTimes=%i  Nlifetimes=%i\n",sumLifetimes,Nlifetimes);
      mprintf("\t\t\tAvg lifetime: %f frames\n", favg);
    }
    ++outSet;
    ++maxSet;
    ++avgSet;
  }
  return 0;
}

void Analysis_Lifetime::PrintListToFile(DataFileList *dfl, std::vector<DataSet*>& list,
                                        std::string const& outname)
{
  DataFile *outfile = NULL;
  if (outname.empty()) return;
  for (std::vector<DataSet*>::iterator set = list.begin(); set != list.end(); ++set)
  {
    outfile = dfl->Add( outname.c_str(), *set );
    if (outfile == NULL) {
      mprinterr("Error adding set %s to file %s\n", (*set)->Legend().c_str(), outname.c_str());
      return;
    }
  }
  outfile->ProcessArgs("noemptyframes");
}

void Analysis_Lifetime::Print(DataFileList* datafilelist) {
  if (!outfilename_.empty()) {
    PrintListToFile(datafilelist, outputDsets_, outfilename_);
    if (!averageonly_) {
      PrintListToFile(datafilelist, maxDsets_, "max." + outfilename_);
      PrintListToFile(datafilelist, avgDsets_, "avg." + outfilename_);
    }
  }
}

