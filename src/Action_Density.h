// -*- mode: c++; -*-

#ifndef INC_ACTION_DENSITY_H
#define INC_ACTION_DENSITY_H

#include <vector>
#include <map>

#include "Action.h"
#include "OnlineVarT.h"
#include "ImagedAction.h"

#define ROUTINE_VERSION_STRING "1.0.2"

/** Calculate density along a coordinate.
  * \author Hannes H. Loeffler.
  */
class Action_Density : public Action {
public:
  Action_Density();

    DispatchObject* Alloc() const { return (DispatchObject*) new Action_Density(); }
    
    void Help() const;

private:
  Action::RetType Init(ArgList&, ActionInit&, int);
  Action::RetType Setup(ActionSetup&);
  Action::RetType DoAction(int, ActionFrame&);
  void Print();

  Action::RetType HistSetup(ActionSetup&);
  Action::RetType DensitySetup(ActionSetup&);
  Action::RetType HistAction(int, ActionFrame&);
  Action::RetType DensityAction(int, ActionFrame&);
  void PrintHist();
  void PrintDensity();

  /// Used to accumulate histogram bins.
  typedef Stats<double> BinType;
  /// Used to map histogram indices to histogram bins.
  typedef std::map<long,BinType> HistType;
  /// Used to hold all histograms
  typedef std::vector<HistType> HistArray;
  /// Used to hold output DataSets
  typedef std::vector<DataSet*> DSarray;
  /// Array of double values
  typedef std::vector<double> Darray;
  /// Used to hold properties for each mask.
  typedef std::vector<Darray> PropArray;

  static const std::string emptystring;
  static const char* PropertyStr_[];
  static const char* AxisStr_[];
  static const double AMU_ANG_TO_G_CM3;

  enum DirectionType {DX = 0, DY, DZ};
  enum PropertyType {NUMBER = 0, MASS, CHARGE, ELECTRON};

  DirectionType axis_;          ///< Which axis to bin along.
  DirectionType area_coord_[2]; ///< Hold which two axes used to calc. area
  PropertyType property_;       ///< Property being binned.

  double delta_;                ///< Histogram spacing
  Stats<double> area_;          ///< Used to accumulate average area

  std::vector<AtomMask> masks_; ///< Hold masks of things to bin.

  DSarray AvSets_;              ///< Hold normalized histogram bin average data sets for each mask
  DSarray SdSets_;              ///< Hold histogram bin SD data sets for each mask
  HistArray histograms_;        ///< Hold raw histograms for each mask

  PropArray properties_;        ///< Hold properties for each mask.
  
  DataSet* density_;            ///< Hold total system density (if not binning)
  ImagedAction image_;          ///< Used to calculate system volume for total density.
};
#endif    
