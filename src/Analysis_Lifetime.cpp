#include "Analysis_Lifetime.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"

// CONSTRUCTOR
Analysis_Lifetime::Analysis_Lifetime() :
  windowSize_(0),
  cut_(0.5),
  averageonly_(false),
  cumulative_(false),
  deltaAvg_(false)
{}

void Analysis_Lifetime::Help() {
  mprintf("\t[out <filename>] <dsetarg0> [ <dsetarg1> ... ]\n");
  mprintf("\t[window <windowsize> [name <setname>]] [averageonly]\n");
  mprintf("\t[cumulative] [cut <cutoff>]\n");
  mprintf("\tCalculate lifetimes for specified data set(s) (time data is > <cutoff>,\n");
  mprintf("\tdefault 0.5) over windows of given size.\n");
}

Analysis::RetType Analysis_Lifetime::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  // Get Keywords
  DataFile* outfile = DFLin->AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  if (outfile != 0) outfile->ProcessArgs("noemptyframes");
  DataFile* maxfile = 0;
  DataFile* avgfile = 0;
  std::string setname = analyzeArgs.GetStringKey("name");
  windowSize_ = analyzeArgs.getKeyInt("window", -1);
  averageonly_ = analyzeArgs.hasKey("averageonly");
  if (!averageonly_ && outfile != 0) {
    maxfile = DFLin->AddDataFile("max." + outfile->DataFilename().Full(), analyzeArgs);
    maxfile->ProcessArgs("noemptyframes");
    avgfile = DFLin->AddDataFile("avg." + outfile->DataFilename().Full(), analyzeArgs);
    avgfile->ProcessArgs("noemptyframes");
  }
  cumulative_ = analyzeArgs.hasKey("cumulative");
  deltaAvg_ = analyzeArgs.hasKey("delta");
  cut_ = analyzeArgs.getKeyDouble("cut", 0.5);
  // Select datasets from remaining args
  ArgList dsetArgs = analyzeArgs.RemainingArgs();
  for (ArgList::const_iterator dsa = dsetArgs.begin(); dsa != dsetArgs.end(); ++dsa)
    inputDsets_ += datasetlist->GetMultipleSets( *dsa );
  if (inputDsets_.empty()) {
    mprinterr("Error: lifetime: No data sets selected.\n");
    return Analysis::ERR;
  }
  // Sort input datasets
  inputDsets_.sort();

  // Create output datasets
  if ( windowSize_ != -1) {
    if (setname.empty()) 
      setname = datasetlist->GenerateDefaultName( "lifetime" );
    int didx = 0;
    for (DataSetList::const_iterator set = inputDsets_.begin(); set != inputDsets_.end(); ++set)
    {
      DataSet* outSet = datasetlist->AddSetIdx( DataSet::FLOAT, setname, didx );
      if (outSet==0) {
        mprinterr("Error: lifetime: Could not allocate output set for %s\n", 
                  (*set)->Legend().c_str());
        return Analysis::ERR;
      }
      outSet->SetLegend( (*set)->Legend() );
      outputDsets_.push_back( outSet );
      if (outfile != 0) outfile->AddSet( outSet );
      if (!averageonly_) {
        // MAX
        // FIXME: CHeck for nullS
        outSet = datasetlist->AddSetIdxAspect( DataSet::INT, setname, didx, "max" );
        outSet->SetLegend( (*set)->Legend() );
        maxDsets_.push_back( outSet );
        if (maxfile != 0) maxfile->AddSet( outSet );
        // AVG
        outSet = datasetlist->AddSetIdxAspect( DataSet::FLOAT, setname, didx, "avg" );
        outSet->SetLegend( (*set)->Legend() );
        avgDsets_.push_back( outSet );
        if (avgfile != 0) avgfile->AddSet( outSet );
      }
      ++didx;
    }
  } else if (outfile != 0) {
    mprinterr("Error: Output file name specified but no window size given ('window <N>')\n");
    return Analysis::ERR;
  }

  if (!averageonly_)
    mprintf("    LIFETIME: Calculating average lifetime using a cutoff of %f", cut_);
  else
    mprintf("    LIFETIME: Calculating only averages");
  mprintf(" of data in %i sets\n", inputDsets_.size());
  if (debugIn > 0)
    inputDsets_.List();
  if (windowSize_ != -1) {
    mprintf("\tAverage of data over windows will be saved to sets named %s\n",
            setname.c_str());
    mprintf("\tWindow size for averaging: %i\n", windowSize_);
    if (cumulative_)
      mprintf("\tCumulative averages will be saved.\n");
    if (deltaAvg_)
      mprintf("\tChange of average from previous average will be saved.\n");
    if (outfile != 0) {
      mprintf("\tOutfile: %s", outfile->DataFilename().base());
      if (!averageonly_)
        mprintf(", %s, %s", maxfile->DataFilename().base(), avgfile->DataFilename().base());
      mprintf("\n");
    }
  }

  return Analysis::OK;
}

// Analysis_Lifetime::Analyze()
Analysis::RetType Analysis_Lifetime::Analyze() {
  float favg;
  int current = 0;
  ProgressBar progress( inputDsets_.size() );
  std::vector<DataSet*>::iterator outSet = outputDsets_.begin();
  std::vector<DataSet*>::iterator maxSet = maxDsets_.begin();
  std::vector<DataSet*>::iterator avgSet = avgDsets_.begin();
  for (DataSetList::const_iterator inSet = inputDsets_.begin(); 
                                   inSet != inputDsets_.end(); ++inSet)
  {
    //mprintf("\t\tCalculating lifetimes for set %s\n", (*inSet)->Legend().c_str());
    progress.Update( current++ );
    // Loop over all values in set.
    int setSize = (*inSet)->Size();
    double sum = 0.0;
    double previous_windowavg = 0.0;
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
        float fval = (float)(windowavg - previous_windowavg);
        if (deltaAvg_) previous_windowavg = windowavg;
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
  return Analysis::OK;
}
