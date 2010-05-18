#include "ShapeAnalysisModuleCLP.h"
#include "ShapeAnalysisModuleComputation.h"
#include "Parameters.h"

using namespace std;
// int main(int argc, char * argv [])
// {
//    	PARSE_ARGS; 
// 	cout<<"------ Shape Analysis Module start -----"<<endl;
// 	
// 	ShapeAnalysisModuleComputation m_computation;
// 
// 	m_computation.SetModulePath(argv[0]);
// 	m_computation.SetChooseInputColumnState(inputColumn);
// 	if(inputColumn) m_computation.SetColumnVolumeFile(columVolumeFile);
// 	
// 	m_computation.SetParameterFile(GroupeProjectInputFile.c_str());
// 	m_computation.SetOutputDirectory(GroupeProjectOutputDirectory.c_str());
// 	m_computation.ReadFile(m_computation.GetParameterFile());
// 	
// 	m_computation.SetOverwriteSegPostProcess(OverwriteSegPostProcess);
// 	m_computation.SetOverwriteGenParaMesh(OverwriteGenParaMesh);
// 	m_computation.SetOverwriteParaToSPHARMMesh(OverwriteParaToSPHARMMesh);
//  	m_computation.SetEnforcedSpace(sx,sy,sz);
// 	m_computation.SetLabelState(LabelState);
// 	m_computation.SetLabel(ValueLabelNumber);
// 	m_computation.SetNumberOfIterations(NumberofIterations);
// 	m_computation.SetSubdivLevel(SubdivLevelValue);
// 	m_computation.SetSPHARMDegree(SPHARMDegreeValue);
// 	
// 	m_computation.SetGaussianFilteringState(GaussianFiltering);
// 	if(GaussianFiltering)
// 		m_computation.SetVarianceBox(VarianceX,VarianceY,VarianceZ);
// 
//  	m_computation.SetTemplateState(Template);
// 	if(Template)
// 	{
// 		m_computation.SetFlipTemplate(flipTemplate.c_str());
//    		m_computation.SetRegTemplate(regTemplate.c_str());
// 	}
// 
// 	m_computation.SetTemplateMState(MTemplate);
// 	m_computation.SetParaOut1State(ParaOut1Template);
// 
// 	if(NoFLip) m_computation.SetFinalFlip(1,0,0,0,0,0,0,0);
// 	else if(FlipALongAxisX) m_computation.SetFinalFlip(0,1,0,0,0,0,0,0);
// 	else if (FlipALongAxisY) m_computation.SetFinalFlip(0,0,1,0,0,0,0,0);
// 	else if (FlipALongAxisZ) m_computation.SetFinalFlip(0,0,0,1,0,0,0,0);
// 	else if (FlipALongAxisXandY) m_computation.SetFinalFlip(0,0,0,0,1,0,0,0);
// 	else if (FlipALongAxisYandZ) m_computation.SetFinalFlip(0,0,0,0,0,1,0,0);
// 	else if (FlipALongAxisXandZ) m_computation.SetFinalFlip(0,0,0,0,0,0,1,0);
// 	else if (FlipALongAxisXandYandZ) m_computation.SetFinalFlip(0,0,0,0,0,0,0,1);
// 	
// 	m_computation.Computation();
// 
// 	return EXIT_SUCCESS ;
// }

int main(int argc, char * argv [])
{

	cout<<"------ Shape Analysis Module start -----"<<endl;
	ShapeAnalysisModuleComputation m_computation;
	m_computation.SetModulePath("/home/hamelc/ShapeAnalysisTool/build64/ShapeAnalysisModule");
	m_computation.SetChooseInputColumnState(false);
	m_computation.SetColumnVolumeFile(5);
	//m_computation.SetParameterFile("../ShapeAnalysisTest/ShapeLocal/Output/ShapeAnalysisModule_OutputFile.csv");
	m_computation.SetParameterFile("/home/hamelc/ShapeAnalysisTool/InputFile.csv");
	m_computation.SetOutputDirectory("/home/hamelc/ShapeAnalysisTool/ShapeAnalysisTest/ShapeLocal2");
	m_computation.ReadFile(m_computation.GetParameterFile());
	
	m_computation.SetOverwriteSegPostProcess(true);
 	m_computation.SetOverwriteGenParaMesh(true);
 	m_computation.SetOverwriteParaToSPHARMMesh(true);
	
 	m_computation.SetEnforcedSpace(0.5,0.5,0.5);

	m_computation.SetLabelState(false);
	m_computation.SetLabel(1);


	m_computation.SetNumberOfIterations(1);
	m_computation.SetSubdivLevel(11);
	m_computation.SetSPHARMDegree(15);
	
	m_computation.SetGaussianFilteringState(true);
	m_computation.SetVarianceBox(12,12,12);

	m_computation.SetTemplateState(false);
	m_computation.SetFlipTemplate("/home/hamelc/ShapeAnalysisModule/ShapeAnalysisTest/test_surfSPHARM.coef");
	m_computation.SetRegTemplate("/home/hamelc/ShapeAnalysisModule/ShapeAnalysisTest/test_surfSPHARM.vtk");
	m_computation.SetTemplateMState(false);

	m_computation.SetParaOut1State(true);

	m_computation.SetFinalFlip(1,0,0,0,0,0,0,0);
	
	m_computation.Computation();

	return EXIT_SUCCESS ;
}

