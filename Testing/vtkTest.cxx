#include "vtkTestMain.h"
#include "itkImageIOFactoryRegisterManager.h"
#include "itkTransformIOFactoryRegisterManager.h"

#ifdef WIN32
#ifdef BUILD_SHARED_LIBS
#define MODULE_IMPORT __declspec(dllimport)
#else
#define MODULE_IMPORT
#endif
#else
#define MODULE_IMPORT
#endif

extern "C" MODULE_IMPORT int ModuleEntryPoint(int, char *[]);

void RegisterTests()
{
  StringToTestFunctionMap["ModuleEntryPoint"] = ModuleEntryPoint;
}
