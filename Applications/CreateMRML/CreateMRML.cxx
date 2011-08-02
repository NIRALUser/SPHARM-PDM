#include "CreateMRMLSceneHelper.h"
#include "MRMLModelHelper.h"
#include "MRMLTransformHelper.h"
#include "MRMLVolumeHelper.h"
#include "MRMLFiducialHelper.h"

#include <vector>
#include <string.h>
/*
//files
-t transform
-v volume
-m model
-h help
-p parent
-n Node name
-f File name
-cc color code
-dc displayed color
-op opacity
-y volume type
format:
t: transform
v: volume
m: model
q: fiducial
-t or -v or -m or -h
if -t or -v or -m
    -f compulsory (file)
    -p optional (default: no parent) ( parent)
    -n optional (default: name of the file stripped from extension and directory) (node name)
if -v
    -c optional (default: grey)
    -y type

If error: prints exec_name -h to print help
//visualization
****
*/

void PrintHelp( const char* arg0, bool extended )
{
  std::cout << "usage: " << std::endl;
  std::cout << arg0 << " SceneFileName [-t/-v/-m/-q/-s][-f/-p/-n][-dc/-op][-cc][-y][-h]" << std::endl;
  if( !extended )
    {
    return;
    }
  std::cout << std::endl;
  std::cout << "-t: Adds a transform" << std::endl;
  std::cout << "-v: Adds a volume" << std::endl;
  std::cout << "-m: Adds a model" << std::endl;
  std::cout << "-q: Adds a fiducial" << std::endl;
  std::cout << "-s filename: Adds an input scene" << std::endl;
  std::cout << "-h: Prints help" << std::endl;
  std::cout << std::endl;
  std::cout << "For each object added to the scene:" << std::endl;
  std::cout << "-f FileName (compulsory![except for fiducials])[relative path from object file to scene file]"
            << std::endl;
  std::cout << "-p ParentNodeName (default: no parent)" << std::endl;
  std::cout << "-n NodeName (default: name of the file stripped from extension and directory)" << std::endl;
  std::cout << std::endl;
  std::cout << "If node is a volume, a fiducial or a model:" << std::endl;
  std::cout << "-op Opacity (default: 1)" << std::endl;
  std::cout << "-dc r,g,b (default: 0.5,0.5,0.5)" << std::endl;
  std::cout << "-as ActiveScalar (default: \"\")" << std::endl;
  std::cout << std::endl;
  std::cout << "If node is a volume or a model:" << std::endl;
  std::cout << "-cc ColorCode (default: 1 [grey])" << std::endl;
  std::cout << std::endl;
  std::cout << "If node is a volume:" << std::endl;
  std::cout << "-y VolumeType (default: scalar)" << std::endl;
  std::cout << "-l: Label map (only for scalar volumes)" << std::endl;
  std::cout << std::endl;
  std::cout << "If node is a transform:" << std::endl;
  std::cout << "-l: matrix (12 values of the 3x3 transform matrix and vector, ie: 1,0,0,0,1,0,0,0,1,2,3,4)"
            << std::endl;
  std::cout << std::endl;
  std::cout << "If node is a fiducial:" << std::endl;
  std::cout << "-id fiducialID (compulsory!)" << std::endl;
  std::cout << "-lbl label" << std::endl;
  std::cout << "-pos x,y,z (position)" << std::endl;
  std::cout << "-o w,x,y,z (orientation)" << std::endl;
  std::cout << "-ts textScale" << std::endl;
  std::cout << "-sc r,g,b (SelectedColor)" << std::endl;
  std::cout << std::endl;
  std::cout << "Color Table:" << std::endl;
  MRMLVolumeHelper volume;
  volume.PrintColors();
  std::cout << std::endl;
  std::cout << "Volume types:" << std::endl;
  volume.PrintVolumeTypes();
}

int ReadCommonSubArguments( std::string arg, std::string value, MRMLNodeHelper* ptr )
{
  if( !arg.compare("-f") )
    {
    if( ptr->GetFileName().compare( "" ) )
      {
      std::cerr << "Error: Multiple filenames for one object" << std::endl;
      return 1;
      }
    if( value[0] == '/' || ( value.size() >= 3  && value[1] == ':'  && value[2] == '\\' ) )
      {
      std::cerr << "Error: Path to file has to be relative" << std::endl;
      }
    ptr->SetFileName( value );
    return -1;
    }
  else if( !arg.compare( "-p" ) )
    {
    if( ptr->GetParentName().compare( "" ) )
      {
      std::cerr << "Error: Multiple parents for one object" << std::endl;
      return 1;
      }
    ptr->SetParentName( value );
    return -1;
    }
  else if( !arg.compare( "-n" ) )
    {
    if( ptr->GetNodeName().compare( "" ) )
      {
      std::cerr << "Error: Multiple node names for one object" << std::endl;
      return 1;
      }
    ptr->SetNodeName( value );
    return -1;
    }
  else if( !arg.compare( "-t" )
           || !arg.compare( "-v" )
           || !arg.compare( "-m" )
           || !arg.compare( "-q" )
           || !arg.compare( "-h" )
           || !arg.compare( "-s" )
           )
    {
    return 2;
    }
  return 0;
}

int ReadColorCodeArgument( bool & colorCodeSet,
                           const char *argv,
                           MRMLColorableHelper* ptr
                           )
{
  if( colorCodeSet )
    {
    std::cerr << "Error: Multiple color codes given for one object" << std::endl;
    return 1;
    }
  colorCodeSet = true;
  std::istringstream stream;
  stream.str( argv );
  int nb;
  stream >> nb;
  if( !stream.fail() )
    {
    if( nb < ptr->GetFirstColor() || nb > ptr->GetLastColor() )
      {
      std::cerr << "Color Table: " << std::endl;
      ptr->PrintColors();
      return 1;
      }
    ptr->SetColor( nb );
    }
  else
    {
    ptr->SetColorString( argv );
    }
  return 0;
}

int ReadValue( std::string argv, double & val, const char* text )
{
  std::istringstream stream( argv );

  stream >> val;
  if( stream.fail() )
    {
    std::cerr << "Error: Problem reading " << text << " values" << std::endl;
    return 1;
    }
  return 0;
}

int ReadVectors( std::string argv, std::vector<double> & vec, unsigned int size, const char* text)
{
//   double RGB[ 3 ] ;
  vec.resize( size );
  std::vector<std::string> RGBstring( size, "" );
  std::size_t              pos;
  bool                     error = false;
  for( unsigned int i = 0; i < size - 1; i++ )
    {
    pos = argv.find( "," );
    if( std::string::npos != pos )
      {
      RGBstring[i] = argv.substr( 0, pos );
      argv = argv.substr( pos + 1, argv.size() - pos  - 1 );
      }
    else
      {
      error = 1;
      break;
      }
    }
  pos = argv.find( "," );
  if( pos != std::string::npos || error )
    {
    std::cerr << "Error: Their should be only " << size << " values given for " << text << std::endl;
    return 1;
    }
  RGBstring[size - 1] = argv;
  for( unsigned int i = 0; i < size; i++ )
    {
    if( ReadValue( RGBstring[i], vec[i], text ) )
      {
      return 1;
      }
    }
  return 0;
}

int ReadPositionArgument( bool & positionSet, std::string argv, MRMLFiducialHelper* ptr )
{
  if( positionSet )
    {
    std::cerr << "Error: Multiple positions given for one object" << std::endl;
    return 1;
    }
  std::vector<double> position;
  if( ReadVectors( argv, position, 3, "position" ) )
    {
    return 1;
    }
  ptr->SetPosition( position[0], position[1], position[2] );
  positionSet = true;
  return 0;
}

int ReadSelectedColorArgument( bool & colorSet, std::string argv, MRMLFiducialHelper* ptr )
{
  if( colorSet )
    {
    std::cerr << "Error: Multiple colors given for one object" << std::endl;
    return 1;
    }
  std::vector<double> color;
  if( ReadVectors( argv, color, 3, "selected color" ) )
    {
    return 1;
    }
  ptr->SetSelectedColor( color[0], color[1], color[2] );
  colorSet = true;
  return 0;
}

int ReadOrientationArgument( bool & orientationSet, std::string argv, MRMLFiducialHelper* ptr )
{
  if( orientationSet )
    {
    std::cerr << "Error: Multiple orientations given for one object" << std::endl;
    return 1;
    }
  std::vector<double> orientation;
  if( ReadVectors( argv, orientation, 4, "orientation" ) )
    {
    return 1;
    }
  ptr->SetOrientation( orientation[0], orientation[1], orientation[2], orientation[3] );
  orientationSet = true;
  return 0;
}

int ReadColorArgument( bool & colorSet,
                       std::string argv,
                       MRMLColorableHelper* ptr
                       )
{
  if( colorSet )
    {
    std::cerr << "Error: Multiple colors given for one object" << std::endl;
    return 1;
    }
  std::vector<double> RGB;
  if( ReadVectors( argv, RGB, 3, "color" ) )
    {
    return 1;
    }
  if( ptr->SetRGB( (float)RGB[0], (float)RGB[1], (float)RGB[2] ) )
    {
    std::cerr << "Error: Color values must be between 0.0 and 1.0" << std::endl;
    return 1;
    }
  colorSet = true;
  return 0;
}

int ReadOpacityArgument( bool & opacitySet, const char* str, MRMLColorableHelper *object )
{
  if( opacitySet )
    {
    std::cerr << "Multiple opacities given for one object" << std::endl;
    return 1;
    }
  opacitySet = true;
  std::istringstream stream;
  stream.str( str );
  double val;
  stream >> val;
  if( stream.fail() || object->SetOpacity( val ) )
    {
    std::cerr << "Error: Opacity value must be between 0.0 and 1.0" << std::endl;
    return 1;
    }
  return 0;
}

int ReadModelSubArguments( int argc,
                           const char *argv[],
                           int & pos,
                           std::vector<MRMLNodeHelper *> & arguments
                           )
{
  MRMLModelHelper* ptr = new MRMLModelHelper;
  int              exit = 0;
  bool             colorCodeSet = false;
  bool             colorSet = false;
  bool             opacitySet = false;
  bool             activeModelSet = false;

  while( pos < argc - 2 )
    {
    pos++;
    exit = ReadCommonSubArguments( argv[pos], argv[pos + 1], ptr );
    if( exit > 0 )
      {
      break;
      }
    else if( exit < 0 )
      {
      pos++;
      continue;
      }
    if( !strcmp( argv[pos], "-cc" ) )
      {
      if( ReadColorCodeArgument( colorCodeSet, argv[pos + 1], ptr ) )
        {
        exit = 1;
        break;
        }
      pos++;
      }
    else if( !strcmp( argv[pos], "-as" ) )
      {
      if( activeModelSet )
        {
        std::cerr << "Error: Multiple active model given for one model" << std::endl;
        exit = 1;
        break;
        }
      ptr->SetActiveScalarName( argv[pos + 1] );
      activeModelSet = true;
      pos++;
      }
    else if( !strcmp( argv[pos], "-op" ) )
      {
      if( ReadOpacityArgument( opacitySet, argv[pos + 1], ptr ) )
        {
        exit = 1;
        break;
        }
      pos++;
      }
    else if( !strcmp( argv[pos], "-dc" ) )
      {
      if( ReadColorArgument( colorSet, argv[pos + 1], ptr ) )
        {
        exit = 1;
        break;
        }
      pos++;
      }
    else
      {
      std::cerr << "Error: No attributes found for one object" << std::endl;
      exit = 1;
      break;
      }
    }

  if( exit == 1 || !ptr->GetFileName().compare( "" ) )
    {
    delete ptr;
    return 1;
    }
  if( exit == 2 )
    {
    pos--;
    }
  arguments.push_back( ptr );
  return 0;
}

int ReadFiducialSubArguments( int argc,
                              const char *argv[],
                              int & pos,
                              std::vector<MRMLNodeHelper *> & arguments
                              )
{
  MRMLFiducialHelper* fiducial = new MRMLFiducialHelper;
  int                 exit = 0;
  bool                idSet = false;
  bool                labelSet = false;
  bool                positionSet = false;
  bool                orientationSet = false;
  bool                colorSet = false;
  bool                opacitySet = false;
  bool                scaleSet = false;
  bool                selectedColorSet = false;
  bool                activeModelSet = false;

  while( pos < argc - 2 )
    {
    pos++;
    exit = ReadCommonSubArguments( argv[pos], argv[pos + 1], fiducial );
    if( exit > 0 )
      {
      break;
      }
    else if( exit < 0 )
      {
      pos++;
      continue;
      }
    if( !strcmp( argv[pos], "-id" ) )
      {
      if( idSet )
        {
        std::cerr << "Error: Multiple Ids given for one fiducial" << std::endl;
        exit = 1;
        break;
        }
      fiducial->SetId( argv[pos + 1] );
      idSet = true;
      pos++;
      }
    else if( !strcmp( argv[pos], "-as" ) )
      {
      if( activeModelSet )
        {
        std::cerr << "Error: Multiple active model given for one fiducial" << std::endl;
        exit = 1;
        break;
        }
      fiducial->SetActiveScalarName( argv[pos + 1] );
      activeModelSet = true;
      pos++;
      }
    else if( !strcmp( argv[pos], "-lbl" ) )
      {
      if( labelSet )
        {
        std::cerr << "Error: Multiple labels given for one fiducial" << std::endl;
        exit = 1;
        break;
        }
      fiducial->SetLabelText( argv[pos + 1] );
      labelSet = true;
      pos++;
      }
    else if( !strcmp( argv[pos], "-pos" ) )
      {
      if( ReadPositionArgument( positionSet, argv[pos + 1], fiducial ) )
        {
        exit = 1;
        break;
        }
      pos++;
      }
    else if( !strcmp( argv[pos], "-sc" ) )
      {
      if( ReadSelectedColorArgument( selectedColorSet, argv[pos + 1], fiducial ) )
        {
        exit = 1;
        break;
        }
      pos++;
      }
    else if( !strcmp( argv[pos], "-op" ) )
      {
      if( ReadOpacityArgument( opacitySet, argv[pos + 1], fiducial ) )
        {
        exit = 1;
        break;
        }
      pos++;
      }
    else if( !strcmp( argv[pos], "-dc" ) )
      {
      if( ReadColorArgument( colorSet, argv[pos + 1], fiducial ) )
        {
        exit = 1;
        break;
        }
      pos++;
      }
    else if( !strcmp( argv[pos], "-o" ) )
      {
      if( ReadOrientationArgument( orientationSet, argv[pos + 1], fiducial ) )
        {
        exit = 1;
        break;
        }
      pos++;
      }
    else if( !strcmp( argv[pos], "-ts" ) )
      {
      if( scaleSet )
        {
        std::cerr << "Error: Multiple scales given for one fiducial" << std::endl;
        exit = 1;
        break;
        }
      double scale;
      if( ReadValue( argv[pos + 1], scale, "scale value" ) )
        {
        exit = 1;
        break;
        }
      fiducial->SetTextScale( scale );
      scaleSet = true;
      pos++;
      }
    else
      {
      std::cerr << "Error: No attributes found for one object" << std::endl;
      exit = 1;
      break;
      }
    }

  if( exit == 1 )
    {
    delete fiducial;
    return 1;
    }
  if( fiducial->GetFileName().compare( "" ) )
    {
    std::cerr << "Error: No file name should be given to a fiducial" << std::endl;
    delete fiducial;
    return 1;
    }
  if( !idSet )
    {
    std::cerr << "Error: An id should be given to a fiducial" << std::endl;
    delete fiducial;
    return 1;

    }
  if( exit == 2 )
    {
    pos--;
    }
  arguments.push_back( fiducial );
  return 0;
}

int ReadTransformArgument( bool & transformSet,
                           std::string argv,
                           MRMLTransformHelper* ptr
                           )
{
  if( transformSet )
    {
    std::cerr << "Error: Multiple transforms given for one object" << std::endl;
    return 1;
    }
  std::vector<double> matrix;
  if( ReadVectors( argv, matrix, 12, "transform" ) )
    {
    return 1;
    }
  if( ptr->SetTransform( matrix ) )
    {
    return 1;
    }
  transformSet = true;
  return 0;
}

int ReadTransformSubArguments( int argc,
                               const char *argv[],
                               int & pos,
                               std::vector<MRMLNodeHelper *> & arguments
                               )
{
  MRMLTransformHelper* ptr = new MRMLTransformHelper;
  int                  exit = 0;
  bool                 transformSet = false;

  while( pos < argc - 2 )
    {
    pos++;
    exit = ReadCommonSubArguments( argv[pos], argv[pos + 1], ptr );
    if( exit > 0 )
      {
      break;
      }
    else if( exit < 0 )
      {
      pos++;
      continue;
      }
    if( !strcmp( argv[pos], "-l" ) )
      {
      if( ReadTransformArgument( transformSet, argv[pos + 1], ptr ) )
        {
        exit = 1;
        break;
        }
      pos++;
      }
    else
      {
      std::cerr << "Error: No attributes found for one object" << std::endl;
      exit = 1;
      break;
      }
    }

  if( exit == 1 || !ptr->GetFileName().compare( "" ) )
    {
    delete ptr;
    return 1;
    }
  if( exit == 2 )
    {
    pos--;
    }
  arguments.push_back( ptr );
  return 0;
}

int ReadVolumeSubArguments( int argc,
                            const char *argv[],
                            int & pos,
                            std::vector<MRMLNodeHelper *> & arguments
                            )
{
  MRMLVolumeHelper* ptr = new MRMLVolumeHelper;
  bool              colorCodeSet = false;
  bool              colorSet = false;
  bool              typeSet = false;
  bool              labelSet = false;
  bool              opacitySet = false;
  bool              activeModelSet = false;
  int               exit = 0;

  while( pos < argc - 1 )
    {
    pos++;
    if( !strcmp( argv[pos], "-l" ) )
      {
      if( labelSet )
        {
        std::cerr << "Error: label map flag given mutliple times for one object" << std::endl;
        exit = 1;
        break;
        }
      labelSet = true;
      ptr->LabelMap( true );
      continue;
      }
    if( pos >= argc - 1 )
      {
      break;
      }
    exit = ReadCommonSubArguments( argv[pos], argv[pos + 1], ptr );
    if( exit > 0 )
      {
      break;
      }
    else if( exit < 0 )
      {
      pos++;
      continue;
      }
    if( !strcmp( argv[pos], "-cc" ) )
      {
      if( ReadColorCodeArgument( colorCodeSet, argv[pos + 1], ptr ) )
        {
        exit = 1;
        break;
        }
      pos++;
      }
    else if( !strcmp( argv[pos], "-as" ) )
      {
      if( activeModelSet )
        {
        std::cerr << "Error: Multiple active model given for one volume" << std::endl;
        exit = 1;
        break;
        }
      ptr->SetActiveScalarName( argv[pos + 1] );
      activeModelSet = true;
      pos++;
      }
    else if( !strcmp( argv[pos], "-op" ) )
      {
      if( ReadOpacityArgument( opacitySet, argv[pos + 1], ptr ) )
        {
        exit = 1;
        break;
        }
      pos++;
      }
    else if( !strcmp( argv[pos], "-dc" ) )
      {
      if( ReadColorArgument( colorSet, argv[pos + 1], ptr ) )
        {
        exit = 1;
        break;
        }
      pos++;
      }
    else if( !strcmp( argv[pos], "-y" ) )
      {
      if( typeSet )
        {
        std::cerr << "Error: Multiple types given for one object" << std::endl;
        exit = 1;
        break;
        }
      exit = ptr->SetVolumeType( argv[pos + 1] );
      if( exit )
        {
        ptr->PrintVolumeTypes();
        break;
        }
      typeSet = true;
      pos++;
      }
    else
      {
      std::cerr << "Error: No attributes found for one object" << std::endl;
      exit = 1;
      break;
      }
    }

  if( exit == 1 || !ptr->GetFileName().compare( "" ) )
    {
    delete ptr;
    return 1;
    }
  if( exit == 2 )
    {
    pos--;
    }
  arguments.push_back( ptr );
  return 0;
}

int ReadArguments( int argc,
                   const char *argv[],
                   std::vector<MRMLNodeHelper *> & arguments,
                   std::string & input
                   )
{
  if( argc < 2 )
    {
    PrintHelp( argv[0], 0 );
    return 1;
    }
  if( argc == 2 && strcmp( argv[1], "-h" ) )
    {
    PrintHelp( argv[0], 0 );
    return 1;
    }
  int i = 1;
  if( strcmp( argv[i], "-h" ) )
    {
    i = 2;
    }
  if( !strcmp( argv[1], "-v" )
      || !strcmp( argv[1], "-t" )
      || !strcmp( argv[1], "-m" )
      || !strcmp( argv[1], "-q" )
      || !strcmp( argv[1], "-s" )
      )
    {
    std::cerr << "First argument must be the scene file name [or -h]" << std::endl;
    PrintHelp( argv[0], 0 );
    return 1;
    }
  for( ; i < argc; i++ )
    {
    if( !strcmp( argv[i], "-v" ) )
      {
      if( ReadVolumeSubArguments( argc, argv, i, arguments ) )
        {
        PrintHelp( argv[0], 0 );
        return 1;
        }
      }
    else if( !strcmp( argv[i], "-t" ) )
      {
      if( ReadTransformSubArguments( argc, argv, i, arguments ) )
        {
        PrintHelp( argv[0], 0 );
        return 1;
        }
      }
    else if( !strcmp( argv[i], "-m" ) )
      {
      if( ReadModelSubArguments( argc, argv, i, arguments ) )
        {
        PrintHelp( argv[0], 0 );
        return 1;
        }
      }
    else if( !strcmp( argv[i], "-q" ) )
      {
      if( ReadFiducialSubArguments( argc, argv, i, arguments ) )
        {
        PrintHelp( argv[0], 0 );
        return 1;
        }
      }
    else if( !strcmp( argv[i], "-h" ) )
      {
      PrintHelp( argv[0], 1 );
      return -1;
      }
    else if( !strcmp( argv[i], "-s" ) )
      {
      if( i + 1 >= argc )
        {
        PrintHelp( argv[0], 1 );
        return -1;
        }
      input = argv[i + 1];
      i++;
      }
    else
      {
      PrintHelp( argv[0], 0 );
      return 1;
      }
    }
  return 0;
}

void CheckNodeName( std::vector<MRMLNodeHelper *> & arguments )
{
  for( unsigned int i = 0; i < arguments.size(); i++ )
    {
    if( !arguments[i]->GetNodeName().compare( "" ) )
      {
      std::string name = arguments[i]->GetFileName();
      // remove folder
      std::size_t pos = name.find_last_of( "/\\" );
      if( pos != std::string::npos )
        {
        name = name.substr( pos + 1 );
        }
      // remove extension
      pos = name.find_last_of( "." );
      if( pos != std::string::npos )
        {
        name = name.substr( 0, pos );
        }
      arguments[i]->SetNodeName( name );
      }
    }
}

void DeleteArguments( std::vector<MRMLNodeHelper *> arguments )
{
  for( unsigned int i = 0; i < arguments.size(); i++ )
    {
    delete arguments[i];
    }
}

int main( int argc, const char* argv[] )
{
  int         output;
  std::string input;

  std::vector<MRMLNodeHelper *> arguments;
  output = ReadArguments( argc, argv, arguments, input );
  if( output > 0 )
    {
    DeleteArguments( arguments );
    return EXIT_FAILURE;
    }
  else if( output < 0 )
    {
    DeleteArguments( arguments );
    return EXIT_SUCCESS;
    }
  output = EXIT_SUCCESS;
  CheckNodeName( arguments );
  CreateMRMLSceneHelper MRMLCreator;
  MRMLCreator.SetOutputSceneName( argv[1] );
  MRMLCreator.SetInputSceneName( input );
  MRMLCreator.SetInputs( arguments );
  output = MRMLCreator.Write();
  DeleteArguments( arguments );
  return output;
}
