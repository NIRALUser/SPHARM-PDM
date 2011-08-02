#include <QApplication>
#include "shapeAnalysisMANCOVA_Wizard.h"
#include "shapeAnalysisMANCOVA_WizardCLP.h"

using namespace std;
int main(int argc, char *argv[])
{

  PARSE_ARGS;
  QApplication app(argc, argv);

  std::string file;
  file = infile;
  shapeAnalysisMANCOVA_Wizard window(file);

  window.show();

  return app.exec();
}
