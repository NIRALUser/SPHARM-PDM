#ifndef _MRMLColorableHelper_h
#define _MRMLColorableHelper_h

#include <vtkMRMLColorTableNode.h>
#include "MRMLNodeHelper.h"

//Copied from vtkMRMLColorTableNode.h
  /// 
  /// The list of valid table types
 
  /// Grey - greyscale ramp
  /// Iron - neutral
  /// Rainbow - red-orange-yellow-blue-indigo-violet
  /// Ocean - bluish ramp
  /// Desert - orange ramp
  /// InvGrey - inverted greyscale ramp
  /// ReverseRainbow - inverted Rainbow
  /// FMRI - fMRI map
  /// FMRIPA - fMRI Positive Activation map
  /// Labels - the Slicer2 FullRainbow editor labels
  /// Random - 255 random colors
  /// User - user defined in the GUI
  /// File - read in from file
  /// Red - red ramp (like greyscale but with red, meant for layering with cyan)
  /// Green - green ramp (like greyscale but with green, layering with magenta)
  /// Blue - blue ramp (like greyscale but with blue, layering with yellow)
  /// Yellow - yellow ramp (complementary ramp to blue, layering yeilds gray)
  /// Cyan - cyan ramp (complementary ramp to red, layering yeilds gray)
  /// Magenta - magenta ramp (complementary ramp to green, layering yeilds gray)
  /// Warm# - ramps of warm colors that are complimentary to Cool#
  /// WarmShade# - ramps of warm colors with variation in value that are
  ///       complimentary to CoolShade# 
  /// WarmTint# - ramps of warm colors with variation in saturation that are
  ///       complimentary to CoolTint# 
  /// Cool# - ramps of cool colors that are complimentary to Warm#
  /// CoolShade# - ramps of cool colors with variation in value that are
  ///       complimentary to WarmShade# 
  /// CoolTint# - ramps of cool colors with variation in saturation that are
  ///       complimentary to WarmSTint# 

/*
      FullRainbow = 0,
      Grey = 1,
      Iron = 2,
      Rainbow = 3,
      Ocean = 4,
      Desert = 5,
      InvGrey = 6,
      ReverseRainbow = 7,
      FMRI = 8,
      FMRIPA = 9,
      Labels = 10,
      Random = 12,
      User = 13,
      File = 14,
      Red = 15,
      Green = 16,
      Blue = 17,
      Yellow = 18,
      Cyan = 19,
      Magenta = 20,
      Warm1 = 21,
      Warm2 = 22,
      Warm3 = 23,
      Cool1 = 24,
      Cool2 = 25,
      Cool3 = 26,
      WarmShade1 = 27,
      WarmShade2 = 28,
      WarmShade3 = 29,
      CoolShade1 = 30,
      CoolShade2 = 31,
      CoolShade3 = 32,
      WarmTint1 = 33,
      WarmTint2 = 34,
      WarmTint3 = 35,
      CoolTint1 = 36,
      CoolTint2 = 37,
      CoolTint3 = 38
*/


class MRMLColorableHelper : public MRMLNodeHelper
{
   public:
      MRMLColorableHelper() ;
      ~MRMLColorableHelper() ;
      int SetOpacity( double value ) ;
      float GetOpacity() ;
      int SetRGB( double R , double G , double B ) ;
      double GetR() ;
      double GetG() ;
      double GetB() ;
      void Print() ;
      std::string GetColor() ;
      void SetColor( int color ) ;
      void SetColorString( std::string color ) ;
      const char* GetColorString() ;
      int GetFirstColor() ;
      int GetLastColor() ;
      void PrintColors() ;
      void SetActiveScalarName( std::string activeScalar ) ;
      std::string GetActiveScalarName( ) ;
   private:
      vtkMRMLColorTableNode *colorNode ;
      std::string m_ColorString ;
      std::string m_ActiveScalar ;
      double m_Opacity ;
      double m_R ;
      double m_G ;
      double m_B ;
};

#endif
