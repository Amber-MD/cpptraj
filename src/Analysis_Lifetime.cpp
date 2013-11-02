#include "Analysis_Lifetime.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"
#include "StringRoutines.h" // integerToString

// CONSTRUCTOR
Analysis_Lifetime::Analysis_Lifetime() :
  windowSize_(0),
  cut_(0.5),
  averageonly_(false),
  cumulative_(false),
  deltaAvg_(false),
  Compare_(Compare_GreaterThan)
{}

void Analysis_Lifetime::Help() {
  mprintf("\t[out <filename>] <dsetarg0> [ <dsetarg1> ... ]\n");
  mprintf("\t[window <windowsize> [name <setname>]] [averageonly]\n");
  mprintf("\t[cumulative] [cut <cutoff>] [greater | less]\n");
  mprintf("\tCalculate lifetimes for specified data set(s) (time data is > <cutoff>,\n");
  mprintf("\tdefault 0.5) over windows of given size.\n");
}

Analysis::RetType Analysis_Lifetime::Setup(Array1D const& dsArray) {
  if (dsArray.empty()) return Analysis::ERR;
  outfileName_.clear();
  inputDsets_ = dsArray;
  windowSize_ = -1;
  averageonly_ = false;
  cumulative_ = false;
  deltaAvg_ = false;
  cut_ = 0.5;
  Compare_ = Compare_GreaterThan;
  return Analysis::OK;
}

Analysis::RetType Analysis_Lifetime::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  // Get Keywords
  outfileName_ = analyzeArgs.GetStringKey("out");
  std::string setname = analyzeArgs.GetStringKey("name");
  windowSize_ = analyzeArgs.getKeyInt("window", -1);
  averageonly_ = analyzeArgs.hasKey("averageonly");
  cumulative_ = analyzeArgs.hasKey("cumulative");
  deltaAvg_ = analyzeArgs.hasKey("delta");
  cut_ = analyzeArgs.getKeyDouble("cut", 0.5);
  if (analyzeArgs.hasKey("greater"))
    Compare_ = Compare_GreaterThan;
  else if (analyzeArgs.hasKey("less"))
    Compare_ = Compare_LessThan;
  else
    Compare_ = Compare_GreaterThan;
  // Select datasets from remaining args
  if (inputDsets_.AddSetsFromArgs( analyzeArgs.RemainingArgs(), *datasetlist )) {
    mprinterr("Error: lifetime: Could not add data sets.\n");
    return Analysis::ERR;
  }

  // Create output datasets
  DataFile* outfile = 0;
  DataFile* maxfile = 0;
  DataFile* avgfile = 0;
  if ( windowSize_ != -1) {
    outfile = DFLin->AddDataFile(outfileName_, analyzeArgs);
    if (!averageonly_ && outfile != 0) {
      maxfile = DFLin->AddDataFile("max." + outfile->DataFilename().Full(), analyzeArgs);
      avgfile = DFLin->AddDataFile("avg." + outfile->DataFilename().Full(), analyzeArgs);
    }
    if (setname.empty()) 
      setname = datasetlist->GenerateDefaultName( "lifetime" );
    int didx = 0;
    for (Array1D::const_iterator set = inputDsets_.begin(); set != inputDsets_.end(); ++set)
    {
      DataSet_1D* outSet = (DataSet_1D*)datasetlist->AddSetIdx( DataSet::FLOAT, setname, didx );
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
        outSet = (DataSet_1D*)datasetlist->AddSetIdxAspect( DataSet::INTEGER, setname, didx, "max" );
        outSet->SetLegend( (*set)->Legend() );
        maxDsets_.push_back( outSet );
        if (maxfile != 0) maxfile->AddSet( outSet );
        // AVG
        outSet = (DataSet_1D*)datasetlist->AddSetIdxAspect( DataSet::FLOAT, setname, didx, "avg" );
        outSet->SetLegend( (*set)->Legend() );
        avgDsets_.push_back( outSet );
        if (avgfile != 0) avgfile->AddSet( outSet );
      }
      ++didx;
    }
    // Set step to window size.
    std::string fileArgs = "xstep " + integerToString( windowSize_ ); 
    if (outfile != 0) outfile->ProcessArgs( fileArgs );
    if (maxfile != 0) maxfile->ProcessArgs( fileArgs );
    if (avgfile != 0) avgfile->ProcessArgs( fileArgs );
  }

  if (!averageonly_)
    mprintf("    LIFETIME: Calculating average lifetime using a cutoff of %f", cut_);
  else
    mprintf("    LIFETIME: Calculating only averages");
  mprintf(" of data in %i sets\n", inputDsets_.size());
  if (debugIn > 0)
    for (Array1D::const_iterator set = inputDsets_.begin(); set != inputDsets_.end(); ++set)
      mprintf("\t%s\n", (*set)->Legend().c_str());
  if (Compare_ == Compare_GreaterThan) 
    mprintf("\tValues greater than %f are considered present.\n", cut_);
  else
    mprintf("\tValues less than %f are considered present.\n", cut_);
  if (windowSize_ != -1) {
    mprintf("\tAverage of data over windows will be saved to sets named %s\n",
            setname.c_str());
    mprintf("\tWindow size for averaging: %i\n", windowSize_);
    if (cumulative_)
      mprintf("\tCumulative averages will be saved.\n");
    if (deltaAvg_)
      mprintf("\tChange of average from previous average will be saved.\n");
  }
  if (!outfileName_.empty()) {
    mprintf("\tOutfile: %s", outfileName_.c_str());
    if (!averageonly_ && outfile != 0)
      mprintf(", %s, %s", maxfile->DataFilename().base(), avgfile->DataFilename().base());
    mprintf("\n");
  }


  return Analysis::OK;
}

// Analysis_Lifetime::Analyze()
Analysis::RetType Analysis_Lifetime::Analyze() {
  float favg;
  int current = 0;
  CpptrajFile standalone_out;
  bool standalone = (!averageonly_ && windowSize_ == -1);
  if (standalone) {
    if (standalone_out.OpenWrite( outfileName_ ))
      return Analysis::ERR;
    standalone_out.Printf("%-10s %10s %10s %10s %10s %s\n","#Set","Nlifetimes",
                          "MaxLT","AvgLT","TotFrames","SetName");
  }
  ProgressBar progress( inputDsets_.size() );
  for (unsigned int setIdx = 0; setIdx < inputDsets_.size(); setIdx++) {
    DataSet_1D const& DS = static_cast<DataSet_1D const&>( *inputDsets_[setIdx] );
    if (standalone)
      mprintf("\t\tCalculating lifetimes for set %s\n", DS.Legend().c_str());
    else
      progress.Update( current++ );
    // Loop over all values in set.
    int setSize = (int)DS.Size();
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
      double dval = DS.Dval(i);
      //mprintf("\t\t\tValue[%i]= %.2f", i,dval);
      if (averageonly_) 
        // Average only
        sum += dval;
      else {
        // Lifetime calculation
        if ( Compare_(dval, cut_) ) {
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
      //sum += inputDsets_[setIdx]->Dval(i);
      ++Ncount;
      ++windowcount;
      if (windowcount == windowSize_) {
        double windowavg = sum / (double)Ncount;
        float fval = (float)(windowavg - previous_windowavg);
        if (deltaAvg_) previous_windowavg = windowavg;
        outputDsets_[setIdx]->Add( frame, &fval );
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
          maxDsets_[setIdx]->Add( frame, &maximumLifetimeCount );
          avgDsets_[setIdx]->Add( frame, &favg );
        }
        //frame += windowcount;
        frame++;
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
      standalone_out.Printf("%10u %10i %10i %10.4f %10.0f %s\n",setIdx,
                            Nlifetimes, maximumLifetimeCount, favg, sum,
                            DS.Legend().c_str());
    }
  }
  return Analysis::OK;
}
