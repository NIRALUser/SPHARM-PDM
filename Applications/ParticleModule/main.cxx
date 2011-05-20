#include "ParticleModuleParameters.h"



using namespace std;
int main()
{
	//TODO PARSE_ARGS; 
	cout<<"------ Particle Module start -----"<<endl;
	
	ParticleModuleParameters m_ParticleModuleParameters;

	m_ParticleModuleParameters.SetOutputDirectory();
	
	m_ParticleModuleParameters.CreateLptsFiles();

	m_ParticleModuleParameters.CreatePreprocessFiles();

	//m_ParticleModuleParameters.RunParticleCorrespondencePreprocessing();

	m_ParticleModuleParameters.CreateMhaFiles();
	m_ParticleModuleParameters.CreateCorrespondanceFiles();

	/*std::string Correspondance_file;
	CreateCorrespondanceFiles(mhaFile,Correspondance_file);*/
	
	//m_ParticleModuleParameters.RunShapeWorksShop();

	//LptsToVTK();
}
