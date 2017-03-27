#include "ParticleModuleCLP.h"
#include "ParticleModuleComputation.h"
#include "ParticleModuleParameters.h"

using namespace std;
int main(int argc, char * argv [])
{
  PARSE_ARGS;
  std::cout << " " << std::endl;
  cout << "------ ParticleModule start -----" << endl;

  ParticleModuleComputation m_computation;
  m_computation.SetColumnMeshFile(columMeshFile);
  m_computation.SetTemplate(TemplateFile);
  m_computation.SetParameterFile(GroupeProjectInputFile.c_str() );
  m_computation.SetOutputDirectory(GroupeProjectOutputDirectory.c_str() );
  m_computation.SetSmooting(SmoothingValue);
  m_computation.SetEnforcedSpace(sx, sy, sz);
  m_computation.SetEndRegularization(EndingRegularization);
  m_computation.SetStartRegularization(StartingRegularization);
  m_computation.SetOptimizationIteration(OptimizationIterations);
  m_computation.SetCheckpointingInterval(CheckpointingInterval);
  m_computation.SetProcrustesScalingState(ProcrustesScalingState);
  m_computation.SetProcrustesIntervalState(ProcrustesIntervalState);
  m_computation.SetAdaptivityStrengthState(AdaptivityStrengthState);
  m_computation.SetRelativeWeightingState(RelativeWeightingState);
  m_computation.SetProcrustesScaling(ProcrustesScaling);
  m_computation.SetProcrustesInterval(ProcrustesInterval);
  m_computation.SetAdaptivityStrength(AdaptivityStrength);
  m_computation.SetRelativeWeighting(RelativeWeighting);

  int InfileProb;
  InfileProb = m_computation.ReadFile(m_computation.GetParameterFile() );
  if( !InfileProb )
    {
    m_computation.OrganizeOutputDirectory();
    m_computation.CreateLptsFiles();
    m_computation.CreatePreprocessFiles();
    m_computation.RunParticleCorrespondencePreprocessing();
    m_computation.CreateCorrespondenceFiles();
    m_computation.RunShapeWorksRun();
    m_computation.CreateVTKFiles();
    m_computation.CreateMrml();
    }

  return 0;

}
