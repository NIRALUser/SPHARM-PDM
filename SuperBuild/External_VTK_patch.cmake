set(ftsystem "${VTKSource}/ThirdParty/freetype/vtkfreetype/builds/unix/ftsystem.c")
file(READ ${ftsystem} code)
string(REPLACE
"#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif"
"#include <fcntl.h>" code "${code}")

file(WRITE ${ftsystem} "${code}")
