	#include "miscMANCOVA.h"
	#include <sstream>
	#include <boost/math/distributions/students_t.hpp>
	
	using namespace std;
	vnl_vector<double> TiedRankCalc(vnl_vector<double>& data)
	{
	// computes ranking of values in data vector, but resolves ties by
	// returning the average rank
	
	vnl_vector<double> output(data.size());
	
	typedef double T;
	std::vector<std::pair<T, int> > pairVector;
	for (unsigned int i=0; i<data.size();++i)
	pairVector.push_back(std::pair<T, int>(data(i), i));
	
	std::sort(pairVector.begin(), pairVector.end());
	
	// need to go through values and average if there are ties
	
	unsigned int lastIndex = 0;
	double lastValue = (pairVector[lastIndex]).first;
	unsigned int currentIndex = 0;
	
	if ( data.size()<=1 ) // there is no data or only one data point, so just return
	{
	for ( unsigned int i=0; i<data.size(); ++i )
	output(i) = i;
	return output;
	}
	
	double nextValue = (pairVector[currentIndex+1]).first;
	
	bool keepGoing = true;
	
	do 
	{
	
	while (lastValue == nextValue) // check for ties
	{
	currentIndex++;
	if ( currentIndex+1>=data.size() ) // we are at the end
		{
		break;
		}
	else
		{
		nextValue = (pairVector[currentIndex+1]).first; // check the
								// next value
		}
	}
	
	unsigned int nrOfTies = currentIndex-lastIndex + 1;
	double averageRankSum = 0;
	
	for (unsigned int i=lastIndex; i<=currentIndex;++i )
	{
	//averageRankSum += (pairVector[i]).second;
	averageRankSum += i;
	}  
	
	averageRankSum /= nrOfTies;
	
	for (unsigned int i=lastIndex; i<=currentIndex;++i )
	{
	output((pairVector[i]).second) = averageRankSum;
	}  
	
	if ( currentIndex>=data.size()-1 ) // we are at the end, so stop processing
	{
	keepGoing = false;
	}
	else
	{
	currentIndex++; // let's look at the next values
	lastIndex = currentIndex;
	lastValue = (pairVector[lastIndex]).first;
	
	if ( currentIndex+1>=data.size() ) // we are now at the end, so
						// simply return the result
		{
		output( currentIndex ) = (pairVector[currentIndex]).second;
		keepGoing = false;
		}
	else
		{ // there is a next value and it is the following ...
	nextValue = (pairVector[currentIndex+1]).first;
		}
	}
	
	} while ( keepGoing ); // that's all folks ...
	
	return output;
	
	}
	
	double computePearsonCorrelation( vnl_vector<double>& x, vnl_vector<double>& y )
	{
	unsigned int numSubjects = x.size();
	
	// first compute the means
	
	double dXM = 0;
	double dYM = 0;
	
	for ( unsigned int sub=0; sub<numSubjects; sub++ )
	{
	dXM += x[sub];
	dYM += y[sub];
	}
	
	dXM/=numSubjects;
	dYM/=numSubjects;
	
	// now compute the variances
	
	double dXV = 0;
	double dYV = 0;
	
	for ( unsigned int sub=0; sub<numSubjects; sub++ )
	{
	dXV+=(x[sub]-dXM)*(x[sub]-dXM);
	dYV+=(y[sub]-dYM)*(y[sub]-dYM);
	}
	
	dXV/=(numSubjects-1);
	dYV/=(numSubjects-1);
	
	// which gives the standard deviations
	
	double dXS = sqrt( dXV );
	double dYS = sqrt( dYV );
	
	// from which we can compute the correlation coefficient
	
	double dR = 0;
	
	for ( unsigned int sub=0; sub<numSubjects; sub++ )
	{
	dR+= (x[sub]-dXM)*(y[sub]-dYM);
	}
	dR/=(numSubjects-1)*dXS*dYS;
	
	return dR;
	
	}
	
	double computePearsonCorrelationWithP( vnl_vector<double>& x, vnl_vector<double>& y, double &dPDiff, double &dPGreater, double &dPSmaller )
	{
	
	// Equations from the book
	// Handbook of Parametric and Nonparametric Statistical Procedures
	// David J. Sheskin, Chapman & Hall/CRC
	//
	// as well as from the BOOST documentation
	
	unsigned int numSubjects = x.size();
	
	double dR = computePearsonCorrelation( x, y );
	
	// Degrees of freedom:
	unsigned v = numSubjects - 1;
	// t-statistic:
	double t_stat = dR*sqrt(numSubjects-2)/sqrt(1-dR*dR);
	
	boost::math::students_t dist(v);
	double q = boost::math::cdf(boost::math::complement(dist, fabs(t_stat)));
	
	dPDiff = 2*q; // need to multiply by 2, because this is a two-sided test
	// just want to reject the possibility that there is no correlation
	
	dPGreater = boost::math::cdf(boost::math::complement(dist, t_stat));
	
	dPSmaller = boost::math::cdf(dist, t_stat);
	
	return dR;
	
	}
	
	void output_vector(vnl_vector<double> data, std::string outbase, std::string toAppend)
	{
	// Create and write out difference vector file:
	std::string outFile=outbase+toAppend;
	
	std::ofstream outputFile;
	outputFile.open(outFile.c_str(), std::ios::out); 
	
	// Set up header:
	outputFile<<"NUMBER_OF_POINTS = "<<data.size()<<std::endl;
	outputFile<<"DIMENSION = 1"<<std::endl;
	outputFile<<"TYPE = Scalar"<<std::endl;
	
	for (unsigned int feat=0;feat<data.size();++feat)
	{
	outputFile<<data(feat)<<std::endl;
	}
	outputFile.close();
	}
	
	void output_matrix(vnl_matrix<double> data, std::string outbase, std::string toAppend)
	{
	// Create and write out difference vector file:
	std::string outFile=outbase+toAppend;
	
	std::ofstream outputFile;
	outputFile.open(outFile.c_str(), std::ios::out); 
	
	// Set up header:
	outputFile<<"NUMBER_OF_POINTS = "<<data.rows()<<std::endl;
	outputFile<<"DIMENSION = "<<data.columns()<<std::endl;
	outputFile<<"TYPE = Vector"<<std::endl;
	
	for (unsigned int feat=0;feat<data.rows();++feat)
	{
	for(unsigned int dim=0;dim<data.columns();++dim)
		outputFile<<data(feat, dim)<<" ";
	outputFile<<std::endl;
	}
	outputFile.close();
	}
	
	void write_SubjectPoints( std::string fileName, unsigned int sub, 
				vnl_matrix<double>  * &featureValue, 
				unsigned int numFeatures,
				MeshType::Pointer & surfaceMesh, 
				MeshSpatialObjectType::Pointer & SOMesh )
	{
	// get the points
	PointsContainerPointer currentPoints;
	currentPoints = surfaceMesh->GetPoints();
	
	PointType currentPoint;
	
	for (unsigned int iI=0; iI<numFeatures; ++iI )
	{
	for ( unsigned int iJ=0; iJ<dimension; iJ++ )
	{
	currentPoint[iJ] = (*featureValue)[sub][iI*dimension+iJ];
	}
	currentPoints->SetElement(iI,currentPoint);
	}
	
	// Write out total mean:
	surfaceMesh->SetPoints(currentPoints); 
	SOMesh->SetMesh(surfaceMesh);
	MeshWriterType::Pointer writer = MeshWriterType::New();
	writer->SetInput(SOMesh);
	writer->SetFileName(fileName.c_str());
	writer->Update();
	
	}
	
	
void
load_MeshList_file( std::string filename, int surfaceColumn, unsigned int numIndependent, 
			unsigned int numGroupTypes,std::vector<int>  groupTypeColumns,std::vector<int> independentColumns, bool scaleOn, int scaleColumn, unsigned int &numSubjects, unsigned int &numA, unsigned int &numB,
			vnl_matrix<int>* &groupLabel, vnl_matrix<double>* &scaleFactor, 
			std::string* &meshFileNames,  vnl_matrix<double> * &indValue, bool computeScaleFactorFromVolumes )
{
	
	
	cout<<"filename: "<<filename<<" numIndependent: "<<numIndependent<<endl;
	//bool interactionTest =false; // Replace
	//unsigned int testColumn=1;
	const int MAXLINE  = 5000; 
	static char line [ MAXLINE ];
	
	char *extension=const_cast<char*>(strrchr(filename.c_str(),'.'));
	std::ifstream datafile(filename.c_str(),std::ios::in); 
	
	if (!datafile.is_open()) 
	{
	std::cerr << "ERROR: Mesh list file does not exist" << filename << std::endl;
	exit(-1);
	}
	
	numSubjects = 0;
	
	datafile.clear(); 
	datafile.getline(line,MAXLINE);
	while (!datafile.eof())
	{
	if (line[0] != '#') // skip over comments
	{
		numSubjects++;
	}
	datafile.getline(line,MAXLINE);  
	
	}
	
	// if the input file is a csv file
	if(!strcmp(extension,".csv")) 
	{
		//because of the header
		numSubjects=numSubjects-1;
		
		if (debug) std::cout << "Num Subjects: " << numSubjects << std::endl;
		
		scaleFactor =  new vnl_matrix<double>; scaleFactor->set_size(numSubjects,1);
		meshFileNames = new std::string [numSubjects];
		
		// Create the data matrix:
		indValue = new vnl_matrix<double>;
		indValue->set_size(numSubjects, numIndependent); // Each subjects data is flattened, this data will go into X
		
		groupLabel = new vnl_matrix<int>;
		groupLabel->set_size(numSubjects, numGroupTypes);
		
		// read the list
		unsigned int curLine = 0; 	
		//char * cur_token;
		string Line;
		
		
		datafile.clear();
		datafile.seekg(ios::beg); // go back to the beginning of the
						// file to read the actual data
		
		getline(datafile,Line); //skip the csv file header

		int temp_int=0;
		double temp_double=0;
		int indexColumn=0;
		
		while (!datafile.eof())
		{
			getline(datafile,Line);
			
			if (Line[0] != '#' && !(Line.empty()) ) 
			{
				char filename[MAXLINE];
				//int retval;
			
				istringstream buffer(Line);
				string cur_token;
			
		
			while(getline(buffer,cur_token,','))
			{
			if ( !cur_token.empty() )
			{
			for (unsigned int gr=0; gr<numGroupTypes; ++gr) // Read all of the group types
			{
			if(indexColumn==groupTypeColumns[gr]) // if the curent column contains a group type
			{	
				if ( cur_token.empty() )
				{
					std::cerr << "ERROR: Could not read group types properly." << std::endl;
					exit(-1);
				}
				temp_int=atoi(cur_token.c_str());
				(*groupLabel)[curLine][gr]=temp_int;
			}
			}
			if(scaleOn)
			{
				if(indexColumn==scaleColumn)//if the current column contains a scale factor
				{
					if ( cur_token.empty() )
					{
						std::cerr << "ERROR: Could not read scaling factor properly." << std::endl;
						exit(-1);
					}
					temp_double=atof(cur_token.c_str());
					(*scaleFactor)[curLine][0]=temp_double; 
				}
			}
			
			if(indexColumn==surfaceColumn) // if the current column contains the surface file
			{		
				if ( cur_token.empty() )
				{
					std::cerr << "ERROR: Could not read filename properly." << std::endl;
					exit(-1);
				}
				strcpy(filename,cur_token.c_str());
				meshFileNames[curLine] = std::string(filename);
			}
		
			for (unsigned int ind_vars=0; ind_vars<numIndependent; ++ind_vars) // Read all of the independent variables
			{
			
			if(indexColumn==independentColumns[ind_vars]) // if the current column contain a independent variable
			{ 
				if ( cur_token.empty() )
				{
					std::cerr << "ERROR: Could not read independent variable properly." << std::endl;
					exit(-1);
				}
				temp_double=atof(cur_token.c_str());
				(*indValue)[curLine][ind_vars]=temp_double;
			}
			}  
			indexColumn++;
			}
			
			else 
			{
			std::cerr << "WARNING: Could not read group token. Skipping line" << std::endl;
			}
		
			}
			curLine++;
			indexColumn=0;
		
			} 
			
		}
		datafile.close();	
	}
	
	else
	{
	if (debug) std::cout << "Num Subjects: " << numSubjects << std::endl;
	
	scaleFactor =  new vnl_matrix<double>; scaleFactor->set_size(numSubjects,1);
	meshFileNames = new std::string [numSubjects];
	
	// Create the data matrix:
	indValue = new vnl_matrix<double>;
	indValue->set_size(numSubjects, numIndependent); // Each subjects data is flattened, this data will go into X
	
	groupLabel = new vnl_matrix<int>;
	groupLabel->set_size(numSubjects, numGroupTypes);
	
	// read the list
	unsigned int curLine = 0; 	
	char * cur_token;
	
	datafile.clear();
	datafile.seekg(0, std::ios::beg); // go back to the beginning of the
					// file to read the actual data
	datafile.getline(line,MAXLINE);
	
	
	int temp_int=0;
	double temp_double=0;
	
	while (!datafile.eof())
	{
	if (line[0] != '#' && !(line[0]==0) ) 
		{
		char filename[MAXLINE];
		int retval;
		
		cur_token = strtok (line," \t");
		
		if ( cur_token!=NULL )
		{
		for (unsigned int gr=0; gr<numGroupTypes; ++gr) // Read all of the group types
		{
		if ( cur_token==NULL )
			{
			std::cerr << "ERROR: Could not read group types properly." << std::endl;
			exit(-1);
			}
		retval = sscanf(cur_token,"%d",&temp_int);
		(*groupLabel)[curLine][gr]=temp_int;
		cur_token = strtok (NULL," \t");
		}
	
		if ( cur_token==NULL )
		{
		std::cerr << "ERROR: Could not read scaling factor properly." << std::endl;
		exit(-1);
		}
		retval = sscanf(cur_token,"%lf",&temp_double); cur_token = strtok (NULL," \t");
		(*scaleFactor)[curLine][0]=temp_double; 
		
		if ( cur_token==NULL )
		{
		std::cerr << "ERROR: Could not read filename properly." << std::endl;
		exit(-1);
		}
		retval = sscanf(cur_token,"%s",filename); cur_token = strtok (NULL," \t");
	
		meshFileNames[curLine] = std::string(filename);
	
		for (unsigned int ind_vars=0; ind_vars<numIndependent; ++ind_vars) // Read all of the independent variables
		{
		if ( cur_token==NULL )
		{
		std::cerr << "ERROR: Could not read independent variable properly." << std::endl;
		exit(-1);
		}
		retval = sscanf(cur_token,"%lf",&temp_double); cur_token = strtok (NULL," \t");
		(*indValue)[curLine][ind_vars]=temp_double;
		}  
	
		curLine++;
	
		} 
		else 
		{
		std::cerr << "WARNING: Could not read group token. Skipping line" << std::endl;
		}
		}
	datafile.getline(line,MAXLINE);
		
	}
	datafile.close();
	}
	
	if ( computeScaleFactorFromVolumes )
	{
	std::cout << std::endl;
	std::cout << "--------------------------------------------------------------------" << std::endl;
	std::cout << "WARNING: Reinterpreting scaling column as volumes." << std::endl;
	std::cout << "... computing the scaling factor from the volumes." << std::endl;
	std::cout << "This is different from the old file format." << std::endl;
	std::cout << "MAKE SURE YOU UNDERSTAND WHAT YOU ENABLED WITH --computeScaleFactorFromVolumes" << std::endl;
	std::cout << "--------------------------------------------------------------------" << std::endl;
	std::cout << std::endl;
	
	// reinterpret the scale column in the file as
	// volumes and compute the scale factor from them
	
	// first compute the average volume
	
	double dAverageVolume = 0;
	
	for ( unsigned int iI=0; iI<numSubjects; iI++ )
	{
	dAverageVolume+=(*scaleFactor)[iI][0];
	}
	dAverageVolume/=numSubjects;
	
	std::cout << "Average volume is = " << dAverageVolume << std::endl;
	std::cout << "Computing scale factor by (volume/averageVolume)^(1/3)" << std::endl;
	
	double one_third=1.0/3;
	
	for ( unsigned int iI=0; iI<numSubjects; iI++ )
	{
	(*scaleFactor)[iI][0] = pow(((*scaleFactor)[iI][0])/dAverageVolume,one_third);
	}
	
	}
	
	std::cout << std::endl << std::endl;
	
	// possibly relabel the group associations
	// TODO: Check if this really should be done for the MANCOVA case
	// TODO: It may only be a historic remnant from the old
	// TODO: StatNonParam program
	
	
	int preLabelA, preLabelB;
	unsigned int  i;
	numA=0, numB=0;
	
	for (unsigned int group_type=0;group_type<numGroupTypes; ++group_type)
	{
	// change labels to the predefined ones, this will fail if there are more than 2 labels in the file
	preLabelA = (*groupLabel)[0][group_type]; 
	preLabelB = (*groupLabel)[0][group_type];
	
	for (i = 0; i < numSubjects; i++) {
		if (preLabelA != (*groupLabel)[i][group_type]) {
		if (preLabelB != preLabelA && preLabelB  != (*groupLabel)[i][group_type]) { 
		std::cout << "Error: more than 2 labels in file" << std::endl;
		} else {
		preLabelB = (*groupLabel)[i][group_type];
		}
		}
	}
	for (i = 0; i < numSubjects; i++) 
		{
		if (preLabelA == (*groupLabel)[i][group_type]) 
		{
		(*groupLabel)[i][group_type] = GROUP_A_LABEL ;
		++numA;
		} 
		else if (preLabelB == (*groupLabel)[i][group_type]) 
		{
		(*groupLabel)[i][group_type] = GROUP_B_LABEL ;
		++numB;
		}
	}
	if (debug) {
		std::cout << "data in group_type "<<group_type <<" has been relabeled: " <<  preLabelA << " --> group A = " << GROUP_A_LABEL
			<< " ; " << preLabelB << " --> group B = " << GROUP_B_LABEL << std::endl;
		std::cout << "#(A)= " << numA << "; #(B)= " << numB << std::endl; 
	}
	}// END of For-loop
	
	// TODO Remove leading and ending " from meshFileNames if present

	if(meshFileNames[1].compare(1,1,"/") == 0)
	{
	std::vector< string > meshFileNamesCopy;
	for (unsigned int i = 0; i < numSubjects; i++) 
	{meshFileNamesCopy.push_back(".");}

	if (debug) 
	{
	for (unsigned int i = 0; i < numSubjects; i++) 
	{
	meshFileNamesCopy.at(i).assign(meshFileNames[i],1,meshFileNames[i].size()-2);
	std::cout << "reading Mesh.Copy " << meshFileNamesCopy.at(i) <<std::endl;
	}
	}
	for (unsigned int i = 0; i < numSubjects; i++) 
	{
		meshFileNames[i]=meshFileNamesCopy[i];
	}
        }		
	// debug info
	if (debug) 
	{
	for (unsigned int i = 0; i < numSubjects; i++) 
		{
		std::cout << meshFileNames[i] << " " << (*groupLabel)[i][0] << " " << (*scaleFactor)[i][0] << std::endl;	
		}
	} 
	
}
	
	
void load_Meshes( bool scaleOn,  unsigned int numSubjects, 
			unsigned int numIndependent, 
			vnl_matrix<double>* &scaleFactor,  vnl_matrix<double> * &indValue,
			std::string* &meshFileNames, unsigned int & numFeatures, vnl_matrix<double>  * &featureValue, MeshType::Pointer & surfaceMesh, MeshSpatialObjectType::Pointer & SOMesh )
{
	// Read the meshes
	featureValue = new vnl_matrix<double>;
	
	MeshReaderType::Pointer reader = MeshReaderType::New();
	PointsContainerPointer points;
	for (unsigned int index = 0; index < numSubjects; index++) {
	
	//if the surface file is a vtk file	
	char *extension =const_cast<char*>(strrchr(meshFileNames[index].c_str(),'.'));
	std::cout << "reading Mesh " << meshFileNames[index].c_str() << " - " << extension << std::endl;

	if(!strcmp(extension,".vtk")) //TODO 
	{
		//read vtk file
		vtkPolyDataReader *mesh = vtkPolyDataReader::New();
		mesh->SetFileName(meshFileNames[index].c_str());

		//std::cout << "converting vtk to itk" << std::endl;

		try{
			mesh->Update();		
		}
		catch(itk::ExceptionObject ex)
		{
			std::cout<< "Error reading VTK meshfile:  "<< meshFileNames[index] << std::endl << "ITK error: " << ex.GetDescription()<< std::endl;
		exit(-3);
		}
		
		vtkPolyData *polyData = mesh->GetOutput();
		
		//convert polydata to itk mesh
		int length=strlen(meshFileNames[index].c_str())-4;

		char *meshFile = new char[length+6];
		strcpy(meshFile,meshFileNames[index].c_str());
		meshFile[length]='\0';
		strcat(meshFile,".meta");
		
		vtkPolyDataToitkMesh *vtkItkConverter = new vtkPolyDataToitkMesh () ;
		vtkItkConverter->SetInput ( polyData ) ;
		
		  // write out the itk meta mesh file
		MeshConverterType * itkConverter = new MeshConverterType () ;
		itkMeshSOType::Pointer meshSO = itkMeshSOType::New () ;
		meshSO->SetMesh ( vtkItkConverter->GetOutput () ) ;
		itkConverter->WriteMeta ( meshSO, meshFile ) ;
	
		mesh->Delete () ;
		delete ( itkConverter ) ;
		
		//read mesh file
		try{
			reader->SetFileName(meshFile);
			reader->Update();
		}
		catch(itk::ExceptionObject ex)
		{
			std::cout<< "Error reading META meshfile:  "<< meshFileNames[index] << std::endl << "ITK error: " << ex.GetDescription()<< std::endl;	
		}
		
		
	}
	
	else
	{
		try{
			reader->SetFileName(meshFileNames[index].c_str());
			reader->Update();
		}
		catch(itk::ExceptionObject ex)
		{
			std::cout<< "Error reading META meshfile:  "<< meshFileNames[index] << std::endl << "ITK error: " << ex.GetDescription()<< std::endl;	
		}
	}
	
		MeshReaderType::SceneType::Pointer scene = reader->GetScene();  
		MeshReaderType::SceneType::ObjectListType * objList =  scene->GetObjects(1,NULL);
		
		MeshReaderType::SceneType::ObjectListType::iterator it = objList->begin();
	
		itk::SpatialObject<3> * curObj = *it;
		SOMesh = dynamic_cast<MeshSpatialObjectType*> (curObj);
		surfaceMesh = SOMesh->GetMesh();
		points = surfaceMesh->GetPoints();
		
		if (index == 0) {
		numFeatures = points->Size();
		featureValue->set_size(numSubjects, numFeatures*3 + numIndependent); // Each subjects data is flattened, this data will go into X
		}

		float xCOG=0;float yCOG=0;float zCOG=0;
		for (unsigned int pointID = 0; pointID < numFeatures; pointID++) {
			PointType curPoint =  points->GetElement(pointID);	
			xCOG=xCOG+curPoint[0];
			yCOG=yCOG+curPoint[1];
			zCOG=zCOG+curPoint[2];
		}
		xCOG=xCOG/numFeatures;
		yCOG=yCOG/numFeatures;
		zCOG=zCOG/numFeatures;
		float COG[3];COG[0]=xCOG;COG[1]=yCOG;COG[2]=zCOG;


		for (unsigned int pointID = 0; pointID < numFeatures; pointID++) {
			PointType curPoint =  points->GetElement(pointID);

			for (unsigned int dim = 0; dim < 3; dim++) {
				if (scaleOn) 
				{
				(*featureValue)[index][pointID*3 +dim] = (curPoint[dim]-COG[dim]) / (*scaleFactor)[index][0] + COG[dim];
				} 
				else 
				{
				(*featureValue)[index][pointID*3 +dim] = curPoint[dim];
				}
			}
		}
		
	
	
	
	
	// Store Independent variables after tupel data:
	for (unsigned int ind_vars=0; ind_vars<numIndependent; ++ind_vars)
	featureValue->set_column(numFeatures*3+ind_vars, indValue->get_column(ind_vars));
	}
	
}
	
	//
	//Function added by bp2009
	//Added to include the longitudinal analysis in the study...
	//
void load_KWMreadableInputFile( bool scaleOn,  unsigned int numSubjects, 
			unsigned int numIndependent, 
			vnl_matrix<double>* &scaleFactor,  vnl_matrix<double> * &indValue,
			std::string* &meshFileNames, unsigned int & numFeatures, vnl_matrix<double>  * &featureValue, MeshType::Pointer & surfaceMesh, MeshSpatialObjectType::Pointer & SOMesh )
{
	// Read the KWMfile
	// By now, and only for testing purposes, I am going to ignore a lot of the input parameters to the function... I copied directly from load_Meshes, and there are a lot of things I dont need
	
	featureValue = new vnl_matrix<double>;
	
	for (unsigned int index = 0; index < numSubjects; index++) {
	try
	{
		//std::cout<< "Working on file...  "<< meshFileNames[index] << std::endl << std::endl;
	
		char line[70]; 	
		ifstream inputFile;
		char *aux;
		vnl_vector<double> curPoint(3,0.0);
		int NPoints, pointID, pointCont;
		
	
		inputFile.open(meshFileNames[index].c_str(), std::ios::in);
		inputFile.getline(line,70,'\n');
		aux=strtok(line, " = ");
		aux=strtok(NULL, " = ");
		NPoints=atoi(aux);
		inputFile.getline(line,70,'\n');
		inputFile.getline(line,70,'\n');
	
		pointID=0; pointCont=0;
		//Start reading the Input File
		while(!inputFile.getline(line,70,'\n').eof())
		{ 
			aux=strtok(line, " ");
			while (aux != NULL)
			{
				curPoint[pointCont]=atof(aux);
				//std::cout<< aux << std::endl;
				aux = strtok (NULL, " ");
				pointCont++;
			}
			pointCont=0;
			//std::cout<< "" "Point " << pointID << " values " << curPoint[0] << " " << curPoint[1] << " " << curPoint[2] << std::endl;
	
			if (index == 0) 
			{
				numFeatures = NPoints;
				featureValue->set_size(numSubjects, numFeatures*3 + numIndependent); // Each subjects data is flattened, this data will go into X
			}
	
			for (unsigned int dim = 0; dim < 3; dim++) 
			{
				if (scaleOn) 
				{
					//std::cout<< "Aqui..." << std::endl;
					(*featureValue)[index][pointID*3 +dim] = curPoint[dim] / (*scaleFactor)[index][0];
				} 
				else 
				{
					//std::cout<< "Aqui no scale..." << std::endl;
					(*featureValue)[index][pointID*3 +dim] = curPoint[dim];
				}
			}
	
			//std::cout<< "File id " << index << " Point " << pointID << " values " << curPoint[0] << " " << curPoint[1] << " " << curPoint[2] << std::endl;
			pointID++;
		}
		inputFile.close();
		//End reading the Input File
	
	} //end FOR index
	catch(itk::ExceptionObject ex)
	{
		std::cout<< "Error reading KWMfile:  "<< meshFileNames[index] << std::endl << "ITK error: " << ex.GetDescription()<< std::endl;
		exit(-3);
	}    
	
	//std::cout<< "Finished reading file id " << index << std::endl;
	
	// Store Independent variables after tupel data:
	for (unsigned int ind_vars=0; ind_vars<numIndependent; ++ind_vars)
	featureValue->set_column(numFeatures*3+ind_vars, indValue->get_column(ind_vars));
	}
	
	std::cout<< "KWM Files Readed..." << std::endl;
	
}
	
	//
	//Function modified by bp2009
	//Included new input argument KWMreadableInputFile...
	//Controlling unallowed operations wrt to Mesh Operations (no Meshes appear in a featural analysis)
	//
void write_SurfaceProperties(  std::string outbase,  PointsContainerPointer & meanPoints,  PointsContainerPointer & meanPointsA,  PointsContainerPointer & meanPointsB,   vnl_matrix<double>& diffVectors,  vnl_vector<double>& normProjections, vnl_vector<double>& normDistProjections, vnl_matrix<double>& zScoresProjected, vnl_matrix<double>& zScores, bool writeOutZScores, MeshType::Pointer & surfaceMesh, MeshSpatialObjectType::Pointer & SOMesh, std::string* & meshFileNames, int KWMreadableInputFile)
{
	
	//std::cout<< "Starting to write surface properties..." << std::endl;
	
	if (KWMreadableInputFile == 0)
	{
		
	// if the input file is a csv file
// 	char *extension=strchr(outbase.c_str(),'.');
// 	if(!strcmp(extension,".csv"))
// 	{
// 		int index=outbase.find(".csv",0);
// 		outbase=outbase.erase(index,outbase.size());
// 	}
	// Write out total mean:
	surfaceMesh->SetPoints(meanPoints); 
	SOMesh->SetMesh(surfaceMesh);
	MeshWriterType::Pointer writer = MeshWriterType::New();
	writer->SetInput(SOMesh);
	std::string FilenameAverage(outbase);
	
	FilenameAverage = FilenameAverage + std::string("_meanAll_uncorrected.meta");
	writer->SetFileName(FilenameAverage.c_str());
	writer->Update();
	
	// Write out Group A mean:
	surfaceMesh->SetPoints(meanPointsA); 
	SOMesh->SetMesh(surfaceMesh);
	writer->SetInput(SOMesh);
	std::string FilenameAverageA(outbase);
	FilenameAverageA = FilenameAverageA + std::string("_meanA.meta");
	writer->SetFileName(FilenameAverageA.c_str());
	writer->Update();
	
	// Write out group B mean:
	surfaceMesh->SetPoints(meanPointsB); 
	SOMesh->SetMesh(surfaceMesh);
	writer->SetInput(SOMesh);
	std::string FilenameAverageB(outbase);
	FilenameAverageB = FilenameAverageB + std::string("_meanB.meta");
	writer->SetFileName(FilenameAverageB.c_str());
	writer->Update();
	} //end if (KWMreadableInputFile == 0)
	
	// Create and write out difference vector file:
	
	output_matrix(diffVectors, outbase, std::string("_diffMesh.txt"));
	
	// normal projections
	
	output_vector(normProjections, outbase, std::string("_normProjections.txt"));
	output_vector(normDistProjections, outbase, std::string("_normDistProjections.txt"));
	
	// write out z-scores
	
	if ( writeOutZScores )
	{
	
	unsigned int numSubjects = zScores.columns();
	
	for ( unsigned sub=0; sub<numSubjects; sub++ )
	{
	boost::filesystem::path my_path( meshFileNames[sub] );
	
	std::ostringstream fileNameSuffixP;
	fileNameSuffixP << "-" << my_path.stem() << "-zScore-projected" << ".txt";
	output_vector(zScoresProjected.get_column(sub),outbase,fileNameSuffixP.str());
	std::ostringstream fileNameSuffix;
	fileNameSuffix << "-" << my_path.stem() << "-zScore-mahalanobis" << ".txt";
	output_vector(zScores.get_column(sub),outbase,fileNameSuffix.str());
	
	}
	
	}
	
}
	
	//
	//Function modified by bp2009
	//Included new input argument KWMreadableInputFile...
	//Controlling unallowed operations wrt to Mesh Operations (no Meshes appear in a featural analysis)
	//
void compute_SurfaceProperties( PointsContainerPointer & meanPoints, PointsContainerPointer & meanPointsA, PointsContainerPointer & meanPointsB, vnl_matrix<double> & diffVectors, vnl_vector<double>& normProjections, vnl_vector<double>& normDistProjections, vnl_matrix<double>& zScores, vnl_matrix<double>& zScoresProjected, unsigned int numSubjects, unsigned int numA, unsigned int numB, unsigned int numFeatures, vnl_matrix<int> * &groupLabel, vnl_matrix<double>* &featureValue, vtkPolyDataNormals *&MeshNormals, MeshType::Pointer & surfaceMesh, std::string outbase, std::string* &meshFileNames, int KWMreadableInputFile)
{
	double tmp_pnt[dimension], tmpA[dimension], tmpB[dimension];
	PointType A_pnt, B_pnt;
	unsigned int dim;
	
	// need to make sure that the overall mean
	// is computed correctly
	// let's do this by averaging the mean of A and B
	// rather than by computing the overall mean
	// which would bias towards the group with more subjects
	
	for (unsigned int pointID = 0; pointID < numFeatures; pointID++) 
	{
	for (dim = 0; dim < dimension; dim++)
		{
		tmpA[dim]=0;tmpB[dim]=0;
		}
	
	for (unsigned int sub = 0; sub < numSubjects; ++sub)
		{
		if ((*groupLabel)[sub][0]==GROUP_A_LABEL)
		for (dim=0;dim<dimension;++dim)
		tmpA[dim]+=(*featureValue)[sub][pointID*dimension+dim];
		else
		for (dim=0;dim<dimension;++dim)
		{
			tmpB[dim]+=(*featureValue)[sub][pointID*dimension+dim]; 
		}
		}
	
	// normalize to compute the mean shapes
	// for the individual groups as well as
	// an overall mean shape
	
	if ( numA>0 ) // only if there are elements in group A
		{
		for (dim = 0; dim < dimension; dim++) // Normalize
		{
		tmpA[dim]/=numA; 
		}
		}
	
	if ( numB>0 ) // only if there are elements in group B
		{
		for (dim = 0; dim < dimension; dim++) // Normalize
		{
		tmpB[dim]/=numB;
		}
		}
	
	if ( numA>0 && numB>0 )
		{
		for (dim = 0; dim < dimension; dim++) // Normalize
		{
		tmp_pnt[dim] = 0.5*(tmpA[dim]+tmpB[dim]);
		}
		}
	else if ( numA>0 )
		{
		for (dim = 0; dim < dimension; dim++) // Normalize
		{
		tmp_pnt[dim] = tmpA[dim];
		}
		}
	else if ( numB>0 )
		{
		for (dim = 0; dim < dimension; dim++) // Normalize
		{
		tmp_pnt[dim] = tmpB[dim];
		}
		}
	else
		{
		for (dim = 0; dim < dimension; dim++) // Normalize
		{
		tmp_pnt[dim] = 0;
		}
		}
	
	meanPoints->InsertElement(pointID, PointType(tmp_pnt));
	meanPointsA->InsertElement(pointID, PointType(tmpA));
	meanPointsB->InsertElement(pointID, PointType(tmpB));
	
	if ( numA>0 && numB>0 )
		{
		for (dim=0;dim<dimension;++dim)// Initialize
		diffVectors[pointID][dim]=tmpA[dim]-tmpB[dim];
		}
	else  // one of the groups is empty, so we cannot compute a
		// difference and therefore set everything to zero
		{
		for (dim=0;dim<dimension;++dim)// Initialize
		diffVectors[pointID][dim]= 0 ;
		}
	
	}
	
	
	// compute the mean surface normals
	// as well as the signed distances by projecting the surface
	// difference onto the surface normals
	
		surfaceMesh->SetPoints(meanPoints); 
		
		itkMeshTovtkPolyData *convertMeshToVTK = new itkMeshTovtkPolyData();
		convertMeshToVTK->SetInput(surfaceMesh);
		vtkPolyData * vtkMesh = convertMeshToVTK->GetOutput();
		
		MeshNormals->SetComputePointNormals(1);
		MeshNormals->SetComputeCellNormals(0);
		MeshNormals->SetSplitting(0);
		MeshNormals->AutoOrientNormalsOn();
		MeshNormals->ConsistencyOn();
		MeshNormals->FlipNormalsOff();  // all normals are outward pointing
		MeshNormals->SetInput(vtkMesh);
		MeshNormals->Update();
		
		vtkPolyData * vtkMeshNormals = MeshNormals->GetOutput();
		vtkMeshNormals->Update();
		
		vtkPointData * NormalPoints = vtkMeshNormals->GetPointData();
		vtkDataArray * ArrayNormal = NormalPoints->GetNormals(); 
	
	// now compute the signed distances
	
	double *tempNorm;
	
	for (unsigned int feat=0; feat<numFeatures;++feat)
	{
		if (KWMreadableInputFile == 0)
		{
			tempNorm=ArrayNormal->GetTuple(feat);  //here is crashing... why
			
			// all of this should automatically become zero if there is only
			// one group, because in this case diffVectors will be all zero
			
			normProjections(feat)= diffVectors(feat,0)*tempNorm[0] +diffVectors(feat,1)*tempNorm[1] +diffVectors(feat,2)*tempNorm[2] ;
		} // end if (KWMreadableInputFile == 0) 
	normDistProjections(feat)= sqrt(diffVectors(feat,0)*diffVectors(feat,0) +diffVectors(feat,1)*diffVectors(feat,1) +diffVectors(feat,2)*diffVectors(feat,2) );
	if (normProjections(feat)<0)// Negative, means pointing in
		normDistProjections(feat)=-normDistProjections(feat); 
	}
	
	// Computing the z-scores ...
	
	// first the projected case
	// let's compute the mean projected distance of A wrt. B and the
	// other way around with respect to the mean surface
	
	if (KWMreadableInputFile == 0)
	{
	
	for ( unsigned int feat=0; feat<numFeatures;feat++ )
	{
	double *tempNorm=ArrayNormal->GetTuple(feat);
	PointType currentMeanPoint =  meanPoints->GetElement(feat);
	vnl_vector<double> dProj(numSubjects,0);
	double mu_aP=0;
	double mu_bP=0;
	
	for (unsigned int sub=0; sub<numSubjects; sub++)
	{
	for (unsigned int tup=0; tup<dimension; ++tup)
		{
		// compute the projection along the mean normal
		dProj(sub) += ((*featureValue)[sub][feat*dimension+tup]-currentMeanPoint[tup])*tempNorm[tup];
		}
	// now switch depending on type
	
	if ((*groupLabel)[sub][0]==GROUP_A_LABEL)
		{
		mu_aP+=dProj(sub);
		}
	else
		{
		mu_bP+=dProj(sub);
		}
	}
	// now compute the mean
	if ( numA>0 )
	{
	mu_aP/=numA;
	}
	if ( numB>0 )
	{
	mu_bP/=numB;
	}
	
	// now compute the variance
	
	double sigma_aP=0;
	double sigma_bP=0;
	
	for (unsigned int sub=0; sub<numSubjects; sub++)
	{
	// now switch depending on type
	
	if ((*groupLabel)[sub][0]==GROUP_A_LABEL)
		{
		sigma_aP+=(dProj(sub)-mu_aP)*(dProj(sub)-mu_aP);
		}
	else
		{
		sigma_bP+=(dProj(sub)-mu_bP)*(dProj(sub)-mu_bP);
		}
	}
	// now compute the standard deviation
	if ( numA>0 )
	{
	sigma_aP/=numA;
	sigma_aP=sqrt(sigma_aP);
	}
	if ( numB>0 )
	{
	sigma_bP/=numB;
	sigma_bP=sqrt(sigma_bP);
	}
	
	// and finally the z-scores
	
	for (unsigned int sub=0; sub<numSubjects; sub++)
	{
	
	// now switch depending on type
	
	if ((*groupLabel)[sub][0]==GROUP_A_LABEL)
		{
		zScoresProjected(feat,sub) = (dProj(sub)-mu_bP)/sigma_bP;
		}
	else
		{
		zScoresProjected(feat,sub) = (dProj(sub)-mu_aP)/sigma_aP;
		}
	}
	}
	
	}  // end of if KWMreadableInputFile 
	
	// now do all of the z-score business using the Mahalanobis distance
	
	for ( unsigned int feat=0; feat<numFeatures;feat++ )
	{
	vnl_vector<double> mu_b(dimension,0);
	vnl_vector<double> mu_a(dimension,0);
	
	vnl_vector<double> currentVec(dimension,0);
	
	for (unsigned int sub=0; sub<numSubjects; sub++)
	{
	for (unsigned int tup=0; tup<dimension; ++tup)
		{
		currentVec(tup) = (*featureValue)[sub][feat*dimension+tup];
		}
	// now switch depending on type
	
	if ((*groupLabel)[sub][0]==GROUP_A_LABEL)
		{
		mu_a+=currentVec;
		}
	else
		{
		mu_b+=currentVec;
		}
	}
	// now compute the mean
	if ( numA>0 )
	{
	mu_a/=numA;
	}
	if ( numB>0 )
	{
	mu_b/=numB;
	}
	
	// now compute the variances
	
	vnl_matrix<double> var_b(dimension,dimension,0);
	vnl_matrix<double> var_a(dimension,dimension,0);
	
	for (unsigned int sub=0; sub<numSubjects; sub++)
	{
	
	for (unsigned int tup=0; tup<dimension; ++tup)
		{
		currentVec(tup) = (*featureValue)[sub][feat*dimension+tup];
		}
	
	// now switch depending on type
	
	if ((*groupLabel)[sub][0]==GROUP_A_LABEL)
		{
		var_a+=outer_product<double>(currentVec-mu_a,currentVec-mu_a);
		}
	else
		{
		var_b+=outer_product<double>(currentVec-mu_b,currentVec-mu_b);
		}
	}
	// now compute the standard deviation
	if ( numA>0 )
	{
	var_a/=numA;
	}
	if ( numB>0 )
	{
	var_b/=numB;
	}
	
	// and finally compute the z-scores by computing
	// the Mahalanobis distance
	
	// compute the inverse of the variance
	
	vnl_matrix<double> var_bInv(dimension,dimension);
	vnl_matrix<double> var_aInv(dimension,dimension);
	
	var_aInv = vnl_matrix_inverse<double>(var_a);
	var_bInv = vnl_matrix_inverse<double>(var_b);
	
	for (unsigned int sub=0; sub<numSubjects; sub++)
	{
	
	for (unsigned int tup=0; tup<dimension; ++tup)
		{
		currentVec(tup) = (*featureValue)[sub][feat*dimension+tup];
		}
	
	// now switch depending on type
	
	if ((*groupLabel)[sub][0]==GROUP_A_LABEL)
		{
		zScores(feat,sub) = sqrt( dot_product(currentVec-mu_b,var_bInv*(currentVec-mu_b)) );
		}
	else
		{
		zScores(feat,sub) = sqrt( dot_product(currentVec-mu_a,var_aInv*(currentVec-mu_a)) );
		}
	}
	
	}
}
	
	
vnl_vector<double> fdrCorrection( vnl_vector<double>& rawP, double fdrLevel, double &fdrThresh )
{
	// to determine the adjusted level, we first need to sort the raw p
	// values
	
	vnl_vector<double> pSorted( rawP );
	std::sort( pSorted.begin(), pSorted.end() );
	
	fdrThresh = 0;
	
	unsigned int numFeatures = rawP.size();
	
	for (int pointID=(int)(numFeatures-1);pointID>=0;--pointID)
	{
		if (pSorted(pointID)<=(pointID+1)*fdrLevel/numFeatures) 
		{
			fdrThresh=pSorted(pointID);
		break;
		}
	}
	
	vnl_vector<double> fdrP( rawP );
	double fdrFactor = 0;
	if ( fdrThresh>0 ) 
	{
		fdrFactor = fdrLevel/fdrThresh;
	}
	if ( fdrFactor < 1.0 ) 
	{
		fdrFactor = 1.0;
	}
	
	for (unsigned int pointID=0;pointID<numFeatures;++pointID) 
	{
		if ( fdrThresh>0 ) 
		{
			fdrP(pointID)=rawP[pointID] * fdrFactor;
		} 
		else 
		{
			fdrP(pointID) = 1.0;  // there is nothing signicant (sorry), so let's set the p-value to 1
		}
		if (fdrP(pointID) > 1.0) 
		{
			fdrP(pointID) = 1.0; // bound the p-value at 1.0, p-values over 1.0 do not make sense
		}
	}
		
	return fdrP;  
	
}
	
vnl_vector<double> bonferroniCorrection( vnl_vector<double>& rawP )
{
	vnl_vector<double> bP( rawP );
	unsigned int sz = rawP.size();
	
	for (unsigned int iI=0; iI<sz; iI++ )
	{
	bP(iI)=rawP(iI)*sz;
	if ( bP(iI)>1 )
	bP(iI) = 1;
	}
	return bP;
}


void Meta2VTK(char * infile, char* outfile)
{
  typedef itk::DefaultDynamicMeshTraits < double, 3, 3, double, double > MeshTraitsType ; 
  typedef itk::Mesh < double, 3, MeshTraitsType > itkMeshType ;
  typedef itk::MeshSpatialObject < itkMeshType > itkMeshSOType ;
  typedef itk::MetaMeshConverter < 3, double, MeshTraitsType > MeshConverterType ;
  
  // read the data in meta format
  MeshConverterType * itkConverter = new MeshConverterType() ;
  itkMeshSOType::Pointer meshSO = itkConverter->ReadMeta (infile) ;
  itkMeshType::Pointer mesh = meshSO->GetMesh() ;
  //delete (itkConverter);
  
  // convert to vtk format
  itkMeshTovtkPolyData * ITKVTKConverter = new itkMeshTovtkPolyData;
  ITKVTKConverter->SetInput ( mesh ) ;
  
  // write out the vtk mesh
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New () ;
  writer->SetInput ( ITKVTKConverter->GetOutput () ) ;
  writer->SetFileName ( outfile ) ;
  writer->Update () ;
  
  writer->Delete () ;
  //delete (ITKVTKConverter);

}


void write_ColorMap(std::string outbase,bool interactionTest,double significanceLevel)
{
int end;
if(interactionTest){end=8;}
else{end=3;}   
char InputMetaFile[512];
 char InputVTKFile[512];
  char TextFile[512];
char OutputFile[512];
strcpy(OutputFile,outbase.c_str());
strcat(OutputFile,"_meanAll_uncorrected.vtk");
  strcpy(InputMetaFile,outbase.c_str());
  strcpy(InputVTKFile,outbase.c_str());
  strcat(InputVTKFile,"_meanAll_uncorrected.vtk");
  strcat(InputMetaFile,"_meanAll_uncorrected.meta");
  Meta2VTK(InputMetaFile,InputVTKFile);

  for(int i=0;i<end;i++)
  {

  std::vector<const char*> args;  
  char* data = NULL;
  int length;
  double timeout = 0.05;
  int result;
  //double PvalueColorMapNb =0.05;

std::ostringstream tmp ;
tmp << significanceLevel;
std::string PvalueColorMapString(tmp.str());
//tmp=tmp.str();

    switch(i)
    {
      case 0:
	if(interactionTest){
	strcpy(TextFile,outbase.c_str());
      strcat(TextFile,"_normProjectionsSpearman.txt");  //TODO fichier txt avec les info pour les overlay (activescalar)
        args.push_back("MeshMath"); //PROG appele
      args.push_back(OutputFile); //origine 
      args.push_back(OutputFile); // sortie (ou serotn stoquees les info de l overlay
      args.push_back("-KWMtoPolyData");
      args.push_back(TextFile);
      args.push_back("normProjectionsSpearman");} //non de ton active scalar
	else{
      strcpy(TextFile,outbase.c_str());
      strcat(TextFile,"_mancovaRawP.txt");
      args.push_back("MeshMath");
      args.push_back(OutputFile);
      args.push_back(OutputFile);
      args.push_back("-KWMtoPolyData");
      args.push_back(TextFile);
      args.push_back("RawP");
	}

      break;

    case 1:
	if(interactionTest){
	strcpy(TextFile,outbase.c_str());
      strcat(TextFile,"_normProjectionsSpearmanPval.txt");
        args.push_back("MeshMath");
      args.push_back(OutputFile);
      args.push_back(OutputFile);
      args.push_back("-KWMtoPolyData");
      args.push_back(TextFile);
      args.push_back("normProjectionsSpearmanPval");
      args.push_back("-significanceLevel");
args.push_back(PvalueColorMapString.c_str() );}
	else{
	strcpy(TextFile,outbase.c_str());
      strcat(TextFile,"_mancovaFDRP.txt");
  
      args.push_back("MeshMath");

      args.push_back(OutputFile);
      args.push_back(OutputFile);
      args.push_back("-KWMtoPolyData");
      args.push_back(TextFile);
      args.push_back("FDRP");
      args.push_back("-significanceLevel");
args.push_back(PvalueColorMapString.c_str());}

      break;

      case 2:
	if(interactionTest){
	strcpy(TextFile,outbase.c_str());
      strcat(TextFile,"_normProjectionsPearsonPval.txt");
        args.push_back("MeshMath");
      args.push_back(OutputFile);
      args.push_back(OutputFile);
      args.push_back("-KWMtoPolyData");
      args.push_back(TextFile);
      args.push_back("normProjectionsPearsonPval");
      args.push_back("-significanceLevel");
args.push_back(PvalueColorMapString.c_str());}
	else{
      strcpy(TextFile,outbase.c_str());
      strcat(TextFile,"_normDistProjections.txt");
        args.push_back("MeshMath");
      args.push_back(OutputFile);
      args.push_back(OutputFile);
      args.push_back("-KWMtoPolyData");
      args.push_back(TextFile);
      args.push_back("normDistProjections");}

      break;

	case 3:
      strcpy(TextFile,outbase.c_str());
      strcat(TextFile,"_normDistProjectionsSpearman.txt");
      args.push_back("MeshMath");
      args.push_back(OutputFile);  
      args.push_back(OutputFile);
      args.push_back("-KWMtoPolyData");
      args.push_back(TextFile);
      args.push_back("normDistProjectionsSpearman");
      break;

	case 4:
	strcpy(TextFile,outbase.c_str());
      strcat(TextFile,"_normDistProjectionsPearson.txt");
      args.push_back("MeshMath");
      args.push_back(OutputFile);  
      args.push_back(OutputFile);
      args.push_back("-KWMtoPolyData");
      args.push_back(TextFile);
      args.push_back("normDistProjectionsPearson");
      break;

	case 5:
     strcpy(TextFile,outbase.c_str());
      strcat(TextFile,"_normDistProjectionsSpearmanPval.txt");
      args.push_back("MeshMath");
      args.push_back(OutputFile);  
      args.push_back(OutputFile);
      args.push_back("-KWMtoPolyData");
      args.push_back(TextFile);
      args.push_back("normDistProjectionsSpearmanPval");
      args.push_back("-significanceLevel");
	args.push_back(PvalueColorMapString.c_str());

      break;

	case 6:
	strcpy(TextFile,outbase.c_str());
      strcat(TextFile,"_normDistProjectionsPearsonPval.txt");
      args.push_back("MeshMath");
      args.push_back(OutputFile);  
      args.push_back(OutputFile);
      args.push_back("-KWMtoPolyData");
      args.push_back(TextFile);
      args.push_back("normDistProjectionsPearsonPval");
  args.push_back("-significanceLevel");
args.push_back(PvalueColorMapString.c_str());


      break;

	case 7:
	strcpy(TextFile,outbase.c_str());
      strcat(TextFile,"_normProjectionsPearson.txt");
      args.push_back("MeshMath");
      args.push_back(OutputFile);
      args.push_back(OutputFile);
      args.push_back("-KWMtoPolyData");
      args.push_back(TextFile);
      args.push_back("normProjectionsPearson");

      break;


	}

  
  args.push_back(0);
	

   // Run the application

  itksysProcess* gp = itksysProcess_New();
  itksysProcess_SetCommand(gp, &*args.begin());
  itksysProcess_SetOption(gp,itksysProcess_Option_HideWindow,1);
  itksysProcess_Execute(gp);



  while(int Value = itksysProcess_WaitForData(gp,&data,&length,&timeout)) // wait for 1s
  {
    if ( ((Value == itksysProcess_Pipe_STDOUT) || (Value == itksysProcess_Pipe_STDERR)) && data[0]=='D' )
    {
	strstream st;
      	for(int i=0;i<length;i++) 	
	{
		st<<data[i];
	}
	string dim=st.str();
    }
    	timeout = 0.05;   	
  }

  itksysProcess_WaitForExit(gp, 0);
  
  result = 1;
  switch(itksysProcess_GetState(gp))
  {
  case itksysProcess_State_Exited:
  {
    result = itksysProcess_GetExitValue(gp);
  } break;
  case itksysProcess_State_Error:
  {
    std::cerr<<"Error: Could not run " << args[0]<<":\n";
    std::cerr<<itksysProcess_GetErrorString(gp)<<"\n";
    std::cout<<"Error: Could not run " << args[0]<<":\n";
    std::cout<<itksysProcess_GetErrorString(gp)<<"\n";
  } break;
  case itksysProcess_State_Exception:
  {
    std::cerr<<"Error: "<<args[0]<<" terminated with an exception: "<<itksysProcess_GetExceptionString(gp)<<"\n";
    std::cout<<"Error: "<<args[0]<<" terminated with an exception: "<<itksysProcess_GetExceptionString(gp)<<"\n";
  } break;
  case itksysProcess_State_Starting:
  case itksysProcess_State_Executing:
  case itksysProcess_State_Expired:
  case itksysProcess_State_Killed:
  {
    // Should not get here.
    std::cerr<<"Unexpected ending state after running "<<args[0]<<std::endl;
    std::cout<<"Unexpected ending state after running "<<args[0]<<std::endl;
  } break;
  }
  itksysProcess_Delete(gp);  


}


}

int minmax (const std::string & filename, double *min, double* max) {
    std::ifstream fichier(filename.c_str(),std::ios::in);
    std::string s;

if(fichier){
       	*min = 0;
	*max=0;
	double i; 
        while(std::getline(fichier,s)) 
	{	
		const char *x = s.c_str();
	i = strtod (x,NULL);

		if(i<*min) {*min=i;}
		else if (i>*max){*max=i;}
	}

    }
else{
        std::cout << "can not open the file " << filename << std::endl;
    }

    fichier.close();
    return 0;
}



void write_MRMLScene(std::string outbase,bool interactionTest)
{

if(interactionTest)
	{	
	cout<<"outbase: "<<outbase<<endl;

	char file[512];
	std::strcpy (file,outbase.c_str());
	std::strcat(file,"_MRMLscene.mrml");

	char nameVTK[512];

	std::strcpy (nameVTK,outbase.c_str());
	std::strcat(nameVTK,"_meanAll_uncorrected.vtk");

	std::ofstream MRMLFile(file);


	MRMLFile<<"<MRML userTags=\"\">"<<std::endl;
	MRMLFile<<"<Selection"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLSelectionNode1\"  name=\"vtkMRMLSelectionNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  activeVolumeID=\"NULL\"  secondaryVolumeID=\"NULL\"  activeLabelVolumeID=\"NULL\"  activeFiducialListID=\"vtkMRMLFiducialListNode1\"  activeROIListID=\"NULL\"  activeCameraID=\"NULL\"  activeViewID=\"NULL\"  activeLayoutID=\"vtkMRMLLayoutNode1\"></Selection>"<<std::endl;
	MRMLFile<<"<Interaction"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLInteractionNode1\"  name=\"vtkMRMLInteractionNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\" currentInteractionMode=\"ViewTransform\"  lastInteractionMode=\"ViewTransform\" ></Interaction>"<<std::endl;
	MRMLFile<<"<Layout"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLLayoutNode1\"  name=\"vtkMRMLLayoutNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  currentViewArrangement=\"4\"  guiPanelVisibility=\"1\"  bottomPanelVisibility =\"0\"  guiPanelLR=\"0\"  collapseSliceControllers=\"0\" numberOfCompareViewRows=\"1\"  numberOfCompareViewColumns=\"1\"  numberOfLightboxRows=\"1\"  numberOfLightboxColumns=\"1\"  mainPanelSize=\"400\"  secondaryPanelSize=\"400\"  selectedModule=\"Fiducials\" ></Layout>"<<std::endl;
	MRMLFile<<"<TGParameters"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLChangeTrackerNode1\"  name=\"vtkMRMLChangeTrackerNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  ROIMin=\"-1 -1 -1\"  ROIMax=\"-1 -1 -1\"  SegmentThresholdMin=\"-1\"  SegmentThresholdMax=\"-1\"  Analysis_Intensity_Flag=\"0\"  Analysis_Deformable_Flag=\"0\"  UseITK=\"1\"  RegistrationChoice=\"3\"  ROIRegistration=\"1\"  ResampleChoice=\"3\"  ResampleConst=\"0.5\" ></TGParameters>"<<std::endl;
	MRMLFile<<"<Crosshair"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLCrosshairNode1\"  name=\"vtkMRMLCrosshairNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  crosshairMode=\"NoCrosshair\"  navigation=\"true\"  crosshairBehavior=\"Normal\"  crosshairThickness=\"Fine\"  crosshairRAS=\"0 0 0\" ></Crosshair>"<<std::endl;
	MRMLFile<<"<ClipModels"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLClipModelsNode1\"  name=\"vtkMRMLClipModelsNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  clipType=\"0\"  redSliceClipState=\"0\"  yellowSliceClipState=\"0\"  greenSliceClipState=\"0\" ></ClipModels>"<<std::endl;
	MRMLFile<<"<ScriptedModule"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLScriptedModuleNode1\"  name=\"vtkMRMLScriptedModuleNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\" ModuleName =\"Editor\" parameter0=\"label 1\"></ScriptedModule>"<<std::endl;
	MRMLFile<<"<ScriptedModule"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLScriptedModuleNode2\"  name=\"vtkMRMLScriptedModuleNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\" ModuleName =\"Editor\" parameter0=\"label 1\"></ScriptedModule>"<<std::endl;
	MRMLFile<<"<ScriptedModule"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLScriptedModuleNode3\"  name=\"vtkMRMLScriptedModuleNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\" ModuleName =\"Editor\" parameter0=\"label 1\"></ScriptedModule>"<<std::endl;
	MRMLFile<<"<ScriptedModule"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLScriptedModuleNode4\"  name=\"vtkMRMLScriptedModuleNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\" ModuleName =\"Editor\" parameter0=\"label 1\"></ScriptedModule>"<<std::endl;
	MRMLFile<<"<ScriptedModule"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLScriptedModuleNode5\"  name=\"vtkMRMLScriptedModuleNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\" ModuleName =\"Editor\" parameter0=\"label 1\"></ScriptedModule>"<<std::endl;
	MRMLFile<<"<ScriptedModule"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLScriptedModuleNode6\"  name=\"vtkMRMLScriptedModuleNode2\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\" ModuleName =\"Editor\" parameter0=\"label 1\"></ScriptedModule>"<<std::endl;
	MRMLFile<<"<View"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLViewNode1\"  name=\"vtkMRMLViewNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  active=\"true\"  visibility=\"true\"  fieldOfView=\"200\"  letterSize=\"0.05\"  boxVisible=\"false\"  fiducialsVisible=\"true\"  fiducialLabelsVisible=\"true\"  axisLabelsVisible=\"false\"  backgroundColor=\"0.701961 0.701961 0.905882\"  animationMode=\"Off\"  viewAxisMode=\"LookFrom\"  spinDegrees=\"2\"  spinMs=\"5\"  spinDirection=\"YawLeft\"  rotateDegrees=\"5\"  rockLength=\"200\"  rockCount=\"0\"  stereoType=\"NoStereo\"  renderMode=\"Perspective\" ></View>"<<std::endl;
	MRMLFile<<"<Camera"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLCameraNode1\"  name=\"vtkMRMLCameraNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  position=\"-270.792 183.224 665.108\"  focalPoint=\"-129.448 93.8171 81.265\"  viewUp=\"0.971364 -0.0134305 0.237216\"  parallelProjection=\"false\"  parallelScale=\"1\"  activetag=\"vtkMRMLViewNode1\" ></Camera>"<<std::endl;
	MRMLFile<<"<ModelStorage"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelStorageNode1\"  name=\"vtkMRMLModelStorageNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\" fileName=\""<<nameVTK<<"\"   useCompression=\"1\"  readState=\"0\"  writeState=\"4\" ></ModelStorage>"<<std::endl;
	MRMLFile<<"<ModelDisplay"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelDisplayNode2\"  name=\"vtkMRMLModelDisplayNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  color=\"0.5 0.5 0.5\"  selectedColor=\"1 0 0\"  selectedAmbient=\"0.4\"  ambient=\"0\"  diffuse=\"1\"  selectedSpecular=\"0.5\"  specular=\"0\"  power=\"1\"  opacity=\"1\"  visibility=\"true\"  clipping=\"false\"  sliceIntersectionVisibility=\"false\"  backfaceCulling=\"true\"  scalarVisibility=\"true\"  vectorVisibility=\"false\"  tensorVisibility=\"false\"  autoScalarRange=\"true\"  scalarRange=\"0 100\"  colorNodeRef=\"vtkMRMLColorTableNodeFilecustomLUT_Pval.txt\"  activeScalarName=\"normProjectionsPearsonPval\"  ></ModelDisplay>"<<std::endl;
	MRMLFile<<" <Model"<<std::endl;
	MRMLFile<<" id=\"vtkMRMLModelNode2\"  name=\"normProjectionsPearsonPval\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  transformNodeRef=\"vtkMRMLLinearTransformNode3\"  storageNodeRef=\"vtkMRMLModelStorageNode1\"  userTags=\"\"  displayNodeRef=\"vtkMRMLModelDisplayNode2\" ></Model>"<<std::endl;
	MRMLFile<<" <ModelStorage"<<std::endl;
	MRMLFile<<" id=\"vtkMRMLModelStorageNode2\"  name=\"vtkMRMLModelStorageNode2\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fileName=\""<<nameVTK<<"\"  useCompression=\"1\"  readState=\"0\"  writeState=\"4\" ></ModelStorage>"<<std::endl;
	MRMLFile<<" <ModelDisplay"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLModelDisplayNode3\"  name=\"vtkMRMLModelDisplayNode2\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  color=\"1 0 0\"  selectedColor=\"1 0 0\"  selectedAmbient=\"0.4\"  ambient=\"0\"  diffuse=\"1\"  selectedSpecular=\"0.5\"  specular=\"0\"  power=\"1\"  opacity=\"1\"  visibility=\"true\"  clipping=\"false\"  sliceIntersectionVisibility=\"false\"  backfaceCulling=\"true\"  scalarVisibility=\"true\"  vectorVisibility=\"false\"  tensorVisibility=\"false\"  autoScalarRange=\"true\"  scalarRange=\"0 100\"  colorNodeRef=\"vtkMRMLColorTableNodeFilecustomLUT_Pval.txt\"  activeScalarName=\"normDistProjectionsSpearmanPval\"  ></ModelDisplay>"<<std::endl;
	MRMLFile<<" <Model"<<std::endl;
	MRMLFile<<" id=\"vtkMRMLModelNode3\"  name=\"normDistProjectionSpearmanPval\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  transformNodeRef=\"vtkMRMLLinearTransformNode2\"  storageNodeRef=\"vtkMRMLModelStorageNode2\"  userTags=\"\"  displayNodeRef=\"vtkMRMLModelDisplayNode3\" ></Model>"<<std::endl;
	MRMLFile<<"<FiducialListStorage"<<std::endl;
	MRMLFile<<" id=\"vtkMRMLFiducialListStorageNode1\"  name=\"vtkMRMLFiducialListStorageNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fileName=\"L1.fcsv\"  useCompression=\"1\"  readState=\"4\"  writeState=\"4\" ></FiducialListStorage>"<<std::endl;
	MRMLFile<<"<FiducialList"<<std::endl;
	MRMLFile<<" id=\"vtkMRMLFiducialListNode1\"  name=\"L1\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLFiducialListStorageNode1\"  userTags=\"\" symbolScale=\"0\" symbolType=\"12\" textScale=\"2.5\" visibility=\"1\" color=\"0.4 1 1\" selectedcolor=\"1 0.5 0.5\" ambient=\"0\" diffuse=\"1\" specular=\"0\" power=\"1\" locked=\"0\"  opacity=\"1\" fiducials=\" "<<std::endl;
	MRMLFile<<"id normDistProjectionSpearmanPval labeltext normDistProjectionSpearmanPval xyz -258.896 193.705 118.184 orientationwxyz 0 0 0 1 selected 0 visibility 1"<<std::endl;
	MRMLFile<<"id normProjectionPearsonPval labeltext normProjectionPearsonPval xyz -95.1053 -22.4888 151.176 orientationwxyz 0 0 0 1 selected 0 visibility 1"<<std::endl;
	MRMLFile<<"id normDistProjectionPearsonPval labeltext normDistProjectionPearsonPval xyz -246.652 -22.4514 158.584 orientationwxyz 0 0 0 1 selected 0 visibility 1"<<std::endl;
	MRMLFile<<"id normProjectionDistPearson labeltext normProjectionDistPearson xyz -170.672 -22.2909 162.53 orientationwxyz 0 0 0 1 selected 0 visibility 1"<<std::endl;
	MRMLFile<<"id normProjectionPearson labeltext normProjectionPearson xyz -23.9077 -19.4522 157.244 orientationwxyz 0 0 0 1 selected 0 visibility 1"<<std::endl;
	MRMLFile<<"id normDistProjectionSpearman labeltext normDistProjectionSpearman xyz -184.426 195.646 114.445 orientationwxyz 0 0 0 1 selected 0 visibility 1"<<std::endl;
	MRMLFile<<"id normProjectionSpearmanPval labeltext normProjectionSpearmanPval xyz -106.315 197.193 111.159 orientationwxyz 0 0 0 1 selected 0 visibility 1"<<std::endl;
	MRMLFile<<"id normProjectionSpearman labeltext normProjectionSpearman xyz -32.4964 189.799 117.304 orientationwxyz 0 0 0 1 selected 0 visibility 1\" ></FiducialList>"<<std::endl;
	MRMLFile<<"<ModelStorage"<<std::endl;
	MRMLFile<<" id=\"vtkMRMLModelStorageNode7\"  name=\"vtkMRMLModelStorageNode7\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fileName=\""<<nameVTK<<"\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></ModelStorage>"<<std::endl;
	MRMLFile<<"<ModelDisplay"<<std::endl;
	MRMLFile<<" id=\"vtkMRMLModelDisplayNode10\"  name=\"vtkMRMLModelDisplayNode10\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  color=\"0.5 0.5 0.5\"  selectedColor=\"1 0 0\"  selectedAmbient=\"0.4\"  ambient=\"0\"  diffuse=\"1\"  selectedSpecular=\"0.5\"  specular=\"0\"  power=\"1\"  opacity=\"1\"  visibility=\"true\"  clipping=\"false\"  sliceIntersectionVisibility=\"false\"  backfaceCulling=\"true\"  scalarVisibility=\"true\"  vectorVisibility=\"false\"  tensorVisibility=\"false\"  autoScalarRange=\"true\"  scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFilecustomLUT_normDistProjectionsSpearman.txt\" activeScalarName=\"normDistProjectionsSpearman\"  ></ModelDisplay>"<<std::endl;
	MRMLFile<<"<Model"<<std::endl;
	MRMLFile<<" id=\"vtkMRMLModelNode10\"  name=\"normDistProjectionsSpearman\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\" transformNodeRef=\"vtkMRMLLinearTransformNode4\"  storageNodeRef=\"vtkMRMLModelStorageNode7\"  userTags=\"\"  displayNodeRef=\"vtkMRMLModelDisplayNode10\" ></Model>"<<std::endl;
	MRMLFile<<" <ModelStorage"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLModelStorageNode8\"  name=\"vtkMRMLModelStorageNode8\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fileName=\""<<nameVTK<<"\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></ModelStorage>"<<std::endl;
	MRMLFile<<" <ModelDisplay"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLModelDisplayNode11\"  name=\"vtkMRMLModelDisplayNode11\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  color=\"0.5 0.5 0.5\"  selectedColor=\"1 0 0\"  selectedAmbient=\"0.4\"  ambient=\"0\"  diffuse=\"1\"  selectedSpecular=\"0.5\"  specular=\"0\"  power=\"1\"  opacity=\"1\"  visibility=\"true\"  clipping=\"false\"  sliceIntersectionVisibility=\"false\"  backfaceCulling=\"true\"  scalarVisibility=\"true\"  vectorVisibility=\"false\"  tensorVisibility=\"false\"  autoScalarRange=\"true\"  scalarRange=\"0 100\"  colorNodeRef=\"vtkMRMLColorTableNodeFilecustomLUT_normDistProjectionsPearson.txt\" activeScalarName=\"normDistProjectionsPearson\"  ></ModelDisplay>"<<std::endl;
	MRMLFile<<" <Model"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLModelNode11\"  name=\"normDistProjectionsPearson\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  transformNodeRef=\"vtkMRMLLinearTransformNode5\"  storageNodeRef=\"vtkMRMLModelStorageNode8\"  userTags=\"\"  displayNodeRef=\"vtkMRMLModelDisplayNode11\" ></Model>"<<std::endl;
	MRMLFile<<"<ModelStorage"<<std::endl;
	MRMLFile<<" id=\"vtkMRMLModelStorageNode10\"  name=\"vtkMRMLModelStorageNode10\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fileName=\""<<nameVTK<<"\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></ModelStorage>"<<std::endl;
	MRMLFile<<" <ModelDisplay"<<std::endl;
	MRMLFile<<" id=\"vtkMRMLModelDisplayNode13\"  name=\"vtkMRMLModelDisplayNode13\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  color=\"0.5 0.5 0.5\"  selectedColor=\"1 0 0\"  selectedAmbient=\"0.4\"  ambient=\"0\"  diffuse=\"1\"  selectedSpecular=\"0.5\"  specular=\"0\"  power=\"1\"  opacity=\"1\"  visibility=\"true\"  clipping=\"false\"  sliceIntersectionVisibility=\"false\"  backfaceCulling=\"true\"  scalarVisibility=\"true\"  vectorVisibility=\"false\"  tensorVisibility=\"false\"  autoScalarRange=\"true\"  scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFilecustomLUT_Pval.txt\" activeScalarName=\"normDistProjectionsPearsonPval\"  ></ModelDisplay>"<<std::endl;
	MRMLFile<<" <Model"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelNode13\"  name=\"normDistProjectionsPearsonPval\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  transformNodeRef=\"vtkMRMLLinearTransformNode6\"  storageNodeRef=\"vtkMRMLModelStorageNode10\"  userTags=\"\"  displayNodeRef=\"vtkMRMLModelDisplayNode13\" ></Model>"<<std::endl;
	MRMLFile<<" <ModelStorage"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLModelStorageNode11\"  name=\"vtkMRMLModelStorageNode11\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fileName=\""<<nameVTK<<"\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></ModelStorage>"<<std::endl;
	MRMLFile<<" <ModelDisplay"<<std::endl;    
	MRMLFile<<" id=\"vtkMRMLModelDisplayNode14\"  name=\"vtkMRMLModelDisplayNode14\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  color=\"0.5 0.5 0.5\"  selectedColor=\"1 0 0\"  selectedAmbient=\"0.4\"  ambient=\"0\"  diffuse=\"1\"  selectedSpecular=\"0.5\"  specular=\"0\"  power=\"1\"  opacity=\"1\"  visibility=\"true\"  clipping=\"false\"  sliceIntersectionVisibility=\"false\"  backfaceCulling=\"true\"  scalarVisibility=\"true\"  vectorVisibility=\"false\"  tensorVisibility=\"false\"  autoScalarRange=\"true\"  scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFilecustomLUT_normProjectionsPearson.txt\" activeScalarName=\"normProjectionsPearson\"  ></ModelDisplay>"<<std::endl;
	MRMLFile<<"<Model"<<std::endl;
	MRMLFile<<" id=\"vtkMRMLModelNode14\"  name=\"normProjectionsPearson\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  transformNodeRef=\"vtkMRMLLinearTransformNode8\"  storageNodeRef=\"vtkMRMLModelStorageNode11\"  userTags=\"\"  displayNodeRef=\"vtkMRMLModelDisplayNode14\" ></Model>"<<std::endl;
	MRMLFile<<"<ModelStorage"<<std::endl;
	MRMLFile<<" id=\"vtkMRMLModelStorageNode12\"  name=\"vtkMRMLModelStorageNode12\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fileName=\""<<nameVTK<<"\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></ModelStorage>"<<std::endl;
	MRMLFile<<"<ModelDisplay"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLModelDisplayNode15\"  name=\"vtkMRMLModelDisplayNode15\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  color=\"0.5 0.5 0.5\"  selectedColor=\"1 0 0\"  selectedAmbient=\"0.4\"  ambient=\"0\"  diffuse=\"1\"  selectedSpecular=\"0.5\"  specular=\"0\"  power=\"1\"  opacity=\"1\"  visibility=\"true\"  clipping=\"false\"  sliceIntersectionVisibility=\"false\"  backfaceCulling=\"true\"  scalarVisibility=\"true\"  vectorVisibility=\"false\"  tensorVisibility=\"false\"  autoScalarRange=\"true\"  scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFilecustomLUT_normProjectionsSpearman.txt\" activeScalarName=\"normProjectionsSpearman\"  ></ModelDisplay>"<<std::endl;
	MRMLFile<<" <Model"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLModelNode15\"  name=\"normProjectionsSpearman\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  transformNodeRef=\"vtkMRMLLinearTransformNode7\"  storageNodeRef=\"vtkMRMLModelStorageNode12\"  userTags=\"\"  displayNodeRef=\"vtkMRMLModelDisplayNode15\" ></Model>"<<std::endl;
	MRMLFile<<"<ModelStorage"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLModelStorageNode14\"  name=\"vtkMRMLModelStorageNode14\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fileName=\""<<nameVTK<<"\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></ModelStorage>"<<std::endl;
	MRMLFile<<" <ModelDisplay"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLModelDisplayNode17\"  name=\"vtkMRMLModelDisplayNode17\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  color=\"0.5 0.5 0.5\"  selectedColor=\"1 0 0\"  selectedAmbient=\"0.4\"  ambient=\"0\"  diffuse=\"1\"  selectedSpecular=\"0.5\"  specular=\"0\"  power=\"1\"  opacity=\"1\"  visibility=\"true\"  clipping=\"false\"  sliceIntersectionVisibility=\"false\"  backfaceCulling=\"true\"  scalarVisibility=\"true\"  vectorVisibility=\"false\"  tensorVisibility=\"false\"  autoScalarRange=\"true\"  scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFilecustomLUT_Pval.txt\" activeScalarName=\"normProjectionsSpearmanPval\"  ></ModelDisplay>"<<std::endl;
	MRMLFile<<"<Model"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLModelNode17\"  name=\"normProjectionsSpearmanPval\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  transformNodeRef=\"vtkMRMLLinearTransformNode1\"  storageNodeRef=\"vtkMRMLModelStorageNode14\"  userTags=\"\"  displayNodeRef=\"vtkMRMLModelDisplayNode17\" ></Model>"<<std::endl;
	MRMLFile<<" <TransformStorage"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLTransformStorageNode2\"  name=\"vtkMRMLTransformStorageNode2\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fileName=\"LinearTransform2.tfm\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></TransformStorage>"<<std::endl;
	MRMLFile<<" <LinearTransform"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLLinearTransformNode6\"  name=\"LinearTransform8\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLTransformStorageNode3\"  userTags=\"\"  matrixTransformToParent=\"1 0 0 -99  0 1 0 -117  0 0 1 0  0 0 0 1\" ></LinearTransform>"<<std::endl;
	MRMLFile<<" <LinearTransform"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLLinearTransformNode7\"  name=\"LinearTransform3\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLTransformStorageNode4\"  userTags=\"\"  matrixTransformToParent=\"1 0 0 126  0 1 0 71  0 0 1 0  0 0 0 1\" ></LinearTransform>"<<std::endl;
	MRMLFile<<" <LinearTransform"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLLinearTransformNode8\"  name=\"LinearTransform4\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLTransformStorageNode5\"  userTags=\"\"  matrixTransformToParent=\"1 0 0 126  0 1 0 -117  0 0 1 0  0 0 0 1\" ></LinearTransform>"<<std::endl;
	MRMLFile<<" <LinearTransform"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLLinearTransformNode4\"  name=\"LinearTransform5\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLTransformStorageNode6\"  userTags=\"\"  matrixTransformToParent=\"1 0 0 -24  0 1 0 71  0 0 1 0  0 0 0 1\" ></LinearTransform>"<<std::endl;
	MRMLFile<<" <LinearTransform"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLLinearTransformNode5\"  name=\"LinearTransform6\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLTransformStorageNode7\"  userTags=\"\"  matrixTransformToParent=\"1 0 0 -24  0 1 0 -117  0 0 1 0  0 0 0 1\" ></LinearTransform>"<<std::endl;
	MRMLFile<<" <LinearTransform"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLLinearTransformNode1\"  name=\"LinearTransform7\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLTransformStorageNode8\"  userTags=\"\"  matrixTransformToParent=\"1 0 0 51  0 1 0 71  0 0 1 0  0 0 0 1\" ></LinearTransform>"<<std::endl;
	MRMLFile<<" <LinearTransform"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLLinearTransformNode2\"  name=\"LinearTransform1\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLTransformStorageNode1\"  userTags=\"\"  matrixTransformToParent=\"1 0 0 -99  0 1 0 71  0 0 1 0  0 0 0 1\" ></LinearTransform>"<<std::endl;
	MRMLFile<<" <LinearTransform"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLLinearTransformNode3\"  name=\"LinearTransform2\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLTransformStorageNode2\"  userTags=\"\"  matrixTransformToParent=\"1 0 0 51  0 1 0 -117  0 0 1 0  0 0 0 1\" ></LinearTransform>"<<std::endl;
	MRMLFile<<" <TransformStorage"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLTransformStorageNode3\"  name=\"vtkMRMLTransformStorageNode3\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fileName=\"LinearTransform8.tfm\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></TransformStorage>"<<std::endl;
	MRMLFile<<"<TransformStorage"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLTransformStorageNode4\"  name=\"vtkMRMLTransformStorageNode4\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fileName=\"LinearTransform3.tfm\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></TransformStorage>"<<std::endl;
	MRMLFile<<" <TransformStorage"<<std::endl;
	MRMLFile<<" id=\"vtkMRMLTransformStorageNode5\"  name=\"vtkMRMLTransformStorageNode5\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fileName=\"LinearTransform4.tfm\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></TransformStorage>"<<std::endl;
	MRMLFile<<"<TransformStorage"<<std::endl;
	MRMLFile<<" id=\"vtkMRMLTransformStorageNode6\"  name=\"vtkMRMLTransformStorageNode6\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fileName=\"LinearTransform5.tfm\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></TransformStorage>"<<std::endl;
	MRMLFile<<" <TransformStorage"<<std::endl;
	MRMLFile<<"  id=\"vtkMRMLTransformStorageNode7\"  name=\"vtkMRMLTransformStorageNode7\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fileName=\"LinearTransform6.tfm\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></TransformStorage>"<<std::endl;
	MRMLFile<<" <TransformStorage"<<std::endl;
	MRMLFile<<" id=\"vtkMRMLTransformStorageNode8\"  name=\"vtkMRMLTransformStorageNode8\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fileName=\"LinearTransform7.tfm\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></TransformStorage>"<<std::endl;
	MRMLFile<<" <Slice"<<std::endl;
	MRMLFile<<" id=\"vtkMRMLSliceNode1\"  name=\"Green\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fieldOfView=\"389.844 250 1\"  dimensions=\"499 320 1\"  activeSlice=\"0\"  layoutGridRows=\"1\"  layoutGridColumns=\"1\"  sliceToRAS=\"-1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1\"  layoutName=\"Green\"  orientation=\"Coronal\"  jumpMode=\"1\"  sliceVisibility=\"false\"  widgetVisibility=\"false\"  useLabelOutline=\"false\"  sliceSpacingMode=\"0\"  prescribedSliceSpacing=\"1 1 1\" ></Slice>"<<std::endl;
	MRMLFile<<"<SliceComposite"<<std::endl;
	MRMLFile<<" id=\"vtkMRMLSliceCompositeNode1\"  name=\"vtkMRMLSliceCompositeNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  backgroundVolumeID=\"\"  foregroundVolumeID=\"\"  labelVolumeID=\"\"  compositing=\"0\"  labelOpacity=\"1\"  linkedControl=\"0\"  foregroundGrid=\"0\"  backgroundGrid=\"0\"  labelGrid=\"1\"  fiducialVisibility=\"1\"  fiducialLabelVisibility=\"1\"  sliceIntersectionVisibility=\"0\"  layoutName=\"Green\"  annotationMode=\"All\"  doPropagateVolumeSelection=\"1\" ></SliceComposite>"<<std::endl;
	MRMLFile<<"<Slice"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLSliceNode2\"  name=\"Red\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fieldOfView=\"389.062 250 1\"  dimensions=\"498 320 1\"  activeSlice=\"0\"  layoutGridRows=\"1\"  layoutGridColumns=\"1\"  sliceToRAS=\"-1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1\"  layoutName=\"Red\"  orientation=\"Axial\"  jumpMode=\"1\"  sliceVisibility=\"false\"  widgetVisibility=\"false\"  useLabelOutline=\"false\"  sliceSpacingMode=\"0\"  prescribedSliceSpacing=\"1 1 1\" ></Slice>"<<std::endl;
	MRMLFile<<"<SliceComposite"<<std::endl;
	MRMLFile<<" id=\"vtkMRMLSliceCompositeNode2\"  name=\"vtkMRMLSliceCompositeNode2\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  backgroundVolumeID=\"\"  foregroundVolumeID=\"\"  labelVolumeID=\"\"  compositing=\"0\"  labelOpacity=\"1\"  linkedControl=\"0\"  foregroundGrid=\"0\"  backgroundGrid=\"0\"  labelGrid=\"1\"  fiducialVisibility=\"1\"  fiducialLabelVisibility=\"1\"  sliceIntersectionVisibility=\"0\"  layoutName=\"Red\"  annotationMode=\"All\"  doPropagateVolumeSelection=\"1\" ></SliceComposite>"<<std::endl;
	MRMLFile<<"<Slice"<<std::endl;
	MRMLFile<<" id=\"vtkMRMLSliceNode3\"  name=\"Yellow\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fieldOfView=\"389.844 250 1\"  dimensions=\"499 320 1\"  activeSlice=\"0\"  layoutGridRows=\"1\"  layoutGridColumns=\"1\"  sliceToRAS=\"0 0 1 0 -1 0 0 0 0 1 0 0 0 0 0 1\"  layoutName=\"Yellow\"  orientation=\"Sagittal\"  jumpMode=\"1\"  sliceVisibility=\"false\"  widgetVisibility=\"false\"  useLabelOutline=\"false\"  sliceSpacingMode=\"0\"  prescribedSliceSpacing=\"1 1 1\" ></Slice>"<<std::endl;
	MRMLFile<<" <SliceComposite"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLSliceCompositeNode3\"  name=\"vtkMRMLSliceCompositeNode3\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  backgroundVolumeID=\"\"  foregroundVolumeID=\"\"  labelVolumeID=\"\"  compositing=\"0\"  labelOpacity=\"1\"  linkedControl=\"0\"  foregroundGrid=\"0\"  backgroundGrid=\"0\"  labelGrid=\"1\"  fiducialVisibility=\"1\"  fiducialLabelVisibility=\"1\"  sliceIntersectionVisibility=\"0\"  layoutName=\"Yellow\"  annotationMode=\"All\"  doPropagateVolumeSelection=\"1\" ></SliceComposite>"<<std::endl;
	MRMLFile<<"<ColorTableStorage"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableStorageNode1\"  name=\"vtkMRMLColorTableStorageNode1\"  hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\"customLUT_normDistProjectionsSpearman.txt\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></ColorTableStorage>"<<std::endl;
	MRMLFile<<"<ColorTable"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableNodeFilecustomLUT_normDistProjectionsSpearman.txt\"  name=\"customLUT_normDistProjectionsSpearman.txt\"  selected=\"false\"  storageNodeRef=\"vtkMRMLColorTableStorageNode1\"  userTags=\"\" type=\"14\" numcolors=\"256\"></ColorTable>"<<std::endl;
	/*MRMLFile<<"<ColorTableStorage"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableStorageNode2\"  name=\"vtkMRMLColorTableStorageNode2\"  hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\"customLUT_normDistProjectionsSpearmanPval.txt\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></ColorTableStorage>"<<std::endl;
 	MRMLFile<<"<ColorTable"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableNodeFilecustomLUT_normDistProjectionsSpearmanPval.txt\" name=\"customLUT_normDistProjectionsSpearmanPval.txt\"   selected=\"false\"  storageNodeRef=\"vtkMRMLColorTableStorageNode2\"  userTags=\"\" type=\"14\" numcolors=\"256\"></ColorTable>"<<std::endl;*/
	MRMLFile<<"<ColorTableStorage"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableStorageNode3\"  name=\"vtkMRMLColorTableStorageNode3\"  hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\"customLUT_normDistProjectionsPearson.txt\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></ColorTableStorage>"<<std::endl;
 	MRMLFile<<"<ColorTable"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableNodeFilecustomLUT_normDistProjectionsPearson.txt\"  name=\"customLUT_normDistProjectionsPearson.txt\"   hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLColorTableStorageNode3\"  userTags=\"\" type=\"14\" numcolors=\"256\"></ColorTable>"<<std::endl;
	/*MRMLFile<<"<ColorTableStorage"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableStorageNode4\"  name=\"vtkMRMLColorTableStorageNode4\"  hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\"customLUT_normDistProjectionsPearsonPval.txt\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></ColorTableStorage>"<<std::endl;
 	MRMLFile<<"<ColorTable"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableNodeFilecustomLUT_normDistProjectionsPearsonPval.txt\"  name=\"customLUT_normDistProjectionsPearsonPval.txt\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLColorTableStorageNode4\"  userTags=\"\" type=\"14\" numcolors=\"256\"></ColorTable>"<<std::endl;*/
	MRMLFile<<"<ColorTableStorage"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableStorageNode5\"  name=\"vtkMRMLColorTableStorageNode5\"  hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\"customLUT_normProjectionsPearson.txt\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></ColorTableStorage>"<<std::endl;
 	MRMLFile<<"<ColorTable"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableNodeFilecustomLUT_normProjectionsPearson.txt\"  name=\"customLUT_normProjectionsPearson.txt\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLColorTableStorageNode5\"  userTags=\"\" type=\"14\" numcolors=\"256\"></ColorTable>"<<std::endl;
	/*MRMLFile<<"<ColorTableStorage"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableStorageNode6\"  name=\"vtkMRMLColorTableStorageNode6\"  hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\"customLUT_normProjectionsPearsonPval.txt\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></ColorTableStorage>"<<std::endl;
 	MRMLFile<<"<ColorTable"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableNodeFilecustomLUT_normProjectionsPearsonPval.txt\"  name=\"customLUT_normProjectionsPearsonPval.txt\"   hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLColorTableStorageNode6\"  userTags=\"\" type=\"14\" numcolors=\"256\"></ColorTable>"<<std::endl;*/
	MRMLFile<<"<ColorTableStorage"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableStorageNode7\"  name=\"vtkMRMLColorTableStorageNode7\"  hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\"customLUT_normProjectionsSpearman.txt\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></ColorTableStorage>"<<std::endl;
 	MRMLFile<<"<ColorTable"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableNodeFilecustomLUT_normProjectionsSpearman.txt\"  name=\"customLUT_normProjectionsSpearman.txt\"    hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLColorTableStorageNode7\"  userTags=\"\" type=\"14\" numcolors=\"256\"></ColorTable>"<<std::endl;
	MRMLFile<<"<ColorTableStorage"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableStorageNode8\"  name=\"vtkMRMLColorTableStorageNode8\"  hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\"customLUT_normProjectionsSpearmanPval.txt\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></ColorTableStorage>"<<std::endl;
 	MRMLFile<<"<ColorTable"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableNodeFilecustomLUT_Pval.txt\"  name=\"customLUT_Pval.txt\"   hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLColorTableStorageNode8\"  userTags=\"\" type=\"14\" numcolors=\"256\"></ColorTable>"<<std::endl;
	MRMLFile<<"</MRML>"<<std::endl;
	}







	
	else
	{	
	cout<<"outbase: "<<outbase<<endl;

	char minmaxf[512];
	std::strcpy (minmaxf,outbase.c_str());
	std::strcat(minmaxf,"_normDistProjections.txt");
	double min, max;
	minmax(minmaxf,&min,&max);

	char file[512];
	std::strcpy (file,outbase.c_str());
	std::strcat(file,"_MRMLscene.mrml");


	char meanA[512];
	char meanB[512];
	//char RawP[512];
	//char FDRP[512];
	//char diffMeshVect[512];
	//char normDistProjections[512];

	std::strcpy (meanA,outbase.c_str());
	std::strcat(meanA,"_meanA.meta");

	std::strcpy (meanB,outbase.c_str());
	std::strcat(meanB,"_meanB.meta");

/*	std::strcpy (RawP,outbase.c_str());
	std::strcat(RawP,"_meanAll_uncorrected_RawP.vtk");

	std::strcpy (FDRP,outbase.c_str());
	std::strcat(FDRP,"_meanAll_uncorrected_FDRP.vtk");

	std::strcpy (diffMeshVect,outbase.c_str());
	std::strcat(diffMeshVect,"_meanAll_uncorrected_diffMesh.vtk");

	std::strcpy (normDistProjections,outbase.c_str());
	std::strcat(normDistProjections,"_meanAll_uncorrected_normDistProjections.vtk");*/


	char nameVTK[512];

	std::strcpy (nameVTK,outbase.c_str());
	std::strcat(nameVTK,"_meanAll_uncorrected.vtk");


	
	std::ofstream MRMLFile(file);
	MRMLFile<<"<MRML userTags=\"\">"<<std::endl;

	MRMLFile<<"<Selection"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLSelectionNode1\" name=\"vtkMRMLSelectionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" activeVolumeID=\"NULL\" secondaryVolumeID=\"NULL\" activeLabelVolumeID=\"NULL\" activeFiducialListID=\"vtkMRMLFiducialListNode5\" activeROIListID=\"NULL\" activeCameraID=\"NULL\" activeViewID=\"NULL\" activeLayoutID=\"vtkMRMLLayoutNode1\"></Selection>"<<std::endl;

	MRMLFile<<"<Interaction"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLInteractionNode1\" name=\"vtkMRMLInteractionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" currentInteractionMode=\"ViewTransform\" lastInteractionMode=\"ViewTransform\"></Interaction>"<<std::endl;
	MRMLFile<<"<Layout"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLLayoutNode1\" name=\"vtkMRMLLayoutNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" currentViewArrangement=\"4\" guiPanelVisibility=\"1\" bottomPanelVisibility =\"0\" guiPanelLR=\"0\" numberOfCompareViewRows=\"0\" numberOfCompareViewColumns=\"0\" numberOfLightboxRows=\"1\" numberOfLightboxColumns=\"1\" mainPanelSize=\"400\" secondaryPanelSize=\"400\"></Layout>"<<std::endl;
	MRMLFile<<"<View"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLViewNode1\" name=\"vtkMRMLViewNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" active=\"false\" fieldOfView=\"200\" letterSize=\"0.05\" boxVisible=\"false\" fiducialsVisible=\"true\" fiducialLabelsVisible=\"true\" axisLabelsVisible=\"false\" backgroundColor=\"0.701961 0.701961 0.905882\" animationMode=\"Off\" viewAxisMode=\"LookFrom\" spinDegrees=\"2\" spinMs=\"5\" spinDirection=\"YawLeft\" rotateDegrees=\"5\" rockLength=\"200\" rockCount=\"0\" stereoType=\"NoStereo\" renderMode=\"Perspective\"></View>"<<std::endl;
	MRMLFile<<"<Camera"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLCameraNode1\" name=\"vtkMRMLCameraNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" position=\"-126.044 36.3633 -530.443\" focalPoint=\"-92.6637 22.8341 -33.6216\" viewUp=\"-0.99766 -0.0154332 0.0666092\" parallelProjection=\"false\" parallelScale=\"1\" active=\"false\"></Camera>"<<std::endl;
	MRMLFile<<"<TGParameters"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLChangeTrackerNode1\" name=\"vtkMRMLChangeTrackerNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" ROIMin=\"-1 -1 -1\" ROIMax=\"-1 -1 -1\" SegmentThresholdMin=\"-1\" SegmentThresholdMax=\"-1\" Analysis_Intensity_Flag=\"0\" Analysis_Deformable_Flag=\"0\" UseITK=\"1\"></TGParameters>"<<std::endl;
	MRMLFile<<"<VolumeRenderingSelection"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLVolumeRenderingSelectionNode1\" name=\"vtkMRMLVolumeRenderingSelectionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" activeVolumeID=\"NULL\" activeVolumeRenderingID=\"NULL\"></VolumeRenderingSelection>\\"<<std::endl;
	MRMLFile<<"<Crosshair"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLCrosshairNode1\" name=\"vtkMRMLCrosshairNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" crosshairMode=\"NoCrosshair\" crosshairBehavior=\"Normal\" crosshairThickness=\"Fine\" crosshairRAS=\"-17.719 127.24 -1.90859\"></Crosshair>"<<std::endl;
	MRMLFile<<"<VolumeRenderingSelection"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLVolumeRenderingParametersNode1\" name=\"vtkMRMLVolumeRenderingSelectionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" volumeNodeID=\"NULL\" croppingEnabled=\"0\" ROINodeID=\"NULL\" volumePropertyNodeID=\"NULL\"></VolumeRenderingSelection>"<<std::endl;
	MRMLFile<<"<VolumeRenderingSelection"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLVolumeRenderingParametersNode2\" name=\"vtkMRMLVolumeRenderingSelectionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" volumeNodeID=\"NULL\" croppingEnabled=\"0\" ROINodeID=\"NULL\" volumePropertyNodeID=\"NULL\"></VolumeRenderingSelection>"<<std::endl;
	MRMLFile<<"<ClipModels"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLClipModelsNode1\" name=\"vtkMRMLClipModelsNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" clipType=\"0\" redSliceClipState=\"0\" yellowSliceClipState=\"0\" greenSliceClipState=\"0\"></ClipModels>\\"<<std::endl;
	MRMLFile<<"<VolumeRenderingSelection"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLVolumeRenderingParametersNode3\" name=\"vtkMRMLVolumeRenderingSelectionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" volumeNodeID=\"NULL\" croppingEnabled=\"0\" ROINodeID=\"NULL\" volumePropertyNodeID=\"NULL\"></VolumeRenderingSelection>"<<std::endl;
	MRMLFile<<"<VolumeRenderingSelection"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLVolumeRenderingParametersNode4\" name=\"vtkMRMLVolumeRenderingSelectionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" volumeNodeID=\"NULL\" croppingEnabled=\"0\" ROINodeID=\"NULL\" volumePropertyNodeID=\"NULL\"></VolumeRenderingSelection>"<<std::endl;
	MRMLFile<<"<ScriptedModule"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLScriptedModuleNode1\" name=\"vtkMRMLScriptedModuleNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" ModuleName =\"Editor\" parameter0= \"label 1\"></ScriptedModule>\\"<<std::endl;
	MRMLFile<<"<ModelStorage"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelStorageNode1\" name=\"vtkMRMLModelStorageNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\""<<meanA<<"\" useCompression=\"1\" readState=\"0\" writeState=\"4\"></ModelStorage>"<<std::endl;
	MRMLFile<<"<ModelDisplay"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelDisplayNode1\" name=\"vtkMRMLModelDisplayNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"1 0 0\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0\" diffuse=\"1\" selectedSpecular=\"0.5\" specular=\"0\" power=\"1\" opacity=\"0.6\" visibility=\"true\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"false\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFullRainbow\" activeScalarName=\"\" ></ModelDisplay>"<<std::endl;
	MRMLFile<<"<Model"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelNode1\" name=\"ShapeAnalysisModule_OutputFile_meanA\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" transformNodeRef=\"vtkMRMLLinearTransformNode9\" storageNodeRef=\"vtkMRMLModelStorageNode1\" userTags=\"\" displayNodeRef=\"vtkMRMLModelDisplayNode1\"></Model>"<<std::endl;
	MRMLFile<<"<ModelStorage"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelStorageNode2\" name=\"vtkMRMLModelStorageNode2\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\""<<meanB<<"\" useCompression=\"1\" readState=\"0\" writeState=\"4\"></ModelStorage>"<<std::endl;
	MRMLFile<<"<ModelDisplay"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelDisplayNode2\" name=\"vtkMRMLModelDisplayNode2\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"0 0.501961 1\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0\" diffuse=\"1\" selectedSpecular=\"0.5\" specular=\"0\" power=\"1\" opacity=\"0.4\" visibility=\"true\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"false\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" activeScalarName=\"\" ></ModelDisplay>"<<std::endl;
	MRMLFile<<"<Model"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelNode2\" name=\"ShapeAnalysisModule_OutputFile_meanB\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" transformNodeRef=\"vtkMRMLLinearTransformNode9\" storageNodeRef=\"vtkMRMLModelStorageNode2\" userTags=\"\" displayNodeRef=\"vtkMRMLModelDisplayNode2\"></Model>"<<std::endl;
	MRMLFile<<"<ModelStorage"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelStorageNode7\" name=\"vtkMRMLModelStorageNode7\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\""<<nameVTK<<"\" useCompression=\"1\" readState=\"0\" writeState=\"4\"></ModelStorage>"<<std::endl;
	MRMLFile<<"<ModelDisplay"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelDisplayNode10\" name=\"vtkMRMLModelDisplayNode10\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"0.5 0.5 0.5\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0\" diffuse=\"1\" selectedSpecular=\"0.5\" specular=\"0\" power=\"1\" opacity=\"1\" visibility=\"true\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"true\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFilecustomLUT_normDistProjections.txt\" activeScalarName=\"normDistProjections\" ></ModelDisplay>"<<std::endl;
	MRMLFile<<"<Model"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelNode10\" name=\"ShapeAnalysisModule_OutputFile_normDistProjections\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLModelStorageNode7\" userTags=\"\" displayNodeRef=\"vtkMRMLModelDisplayNode10\"></Model>"<<std::endl;
	MRMLFile<<"<LinearTransform"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLLinearTransformNode4\" name=\"LinearTransform3\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLTransformStorageNode3\" userTags=\"\" matrixTransformToParent=\"1 -0 0 -139  -0 1 -0 -12  0 -0 1 -0  -0 0 -0 1\"></LinearTransform>"<<std::endl;
	MRMLFile<<"<LinearTransform"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLLinearTransformNode5\" name=\"LinearTransform4\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLTransformStorageNode4\" userTags=\"\" matrixTransformToParent=\"1 -0 0 -200  -0 1 -0 -12  0 -0 1 -0  -0 0 -0 1\"></LinearTransform>"<<std::endl;

	MRMLFile<<"<ModelStorage"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelStorageNode8\" name=\"vtkMRMLModelStorageNode8\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\""<<nameVTK<<"\" useCompression=\"1\" readState=\"0\" writeState=\"4\"></ModelStorage>\\"<<std::endl;
	MRMLFile<<"<ModelDisplay"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelDisplayNode11\" name=\"vtkMRMLModelDisplayNode11\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"0.5 0.5 0.5\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0\" diffuse=\"1\" selectedSpecular=\"0.5\" specular=\"0\" power=\"1\" opacity=\"1\" visibility=\"true\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"true\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFilecustomLUT_RawP.txt\" activeScalarName=\"RawP\" ></ModelDisplay>"<<std::endl;
	MRMLFile<<"<Model"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelNode11\" name=\"ShapeAnalysisModule_OutputFile_RawP\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" transformNodeRef=\"vtkMRMLLinearTransformNode5\" storageNodeRef=\"vtkMRMLModelStorageNode8\" userTags=\"\" displayNodeRef=\"vtkMRMLModelDisplayNode11\"></Model>"<<std::endl;
	MRMLFile<<"<ModelStorage"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelStorageNode9\" name=\"vtkMRMLModelStorageNode9\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\""<<nameVTK<<"\" useCompression=\"1\" readState=\"0\" writeState=\"0\"></ModelStorage>"<<std::endl;
	MRMLFile<<"<LinearTransform"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLLinearTransformNode9\" name=\"LinearTransform5\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLTransformStorageNode5\" userTags=\"\" matrixTransformToParent=\"1 -0 0 -109  -0 1 -0 -12  0 -0 1 -0  -0 0 -0 1\"></LinearTransform>"<<std::endl;
	MRMLFile<<"<ModelStorage"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelStorageNode15\" name=\"vtkMRMLModelStorageNode15\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\""<<nameVTK<<"\" useCompression=\"1\" readState=\"0\" writeState=\"4\"></ModelStorage>\\n"<<std::endl;
	MRMLFile<<"<ModelDisplay"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelDisplayNode18\" name=\"vtkMRMLModelDisplayNode18\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"0.5 0.5 0.5\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0\" diffuse=\"1\" selectedSpecular=\"0.5\" specular=\"0\" power=\"1\" opacity=\"1\" visibility=\"true\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"true\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFilecustomLUT_FDRP.txt\"  activeScalarName=\"FDRP\" ></ModelDisplay>"<<std::endl;
	MRMLFile<<"<Model"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelNode18\" name=\"ShapeAnalysisModule_OutputFile_FDRP\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" transformNodeRef=\"vtkMRMLLinearTransformNode14\" storageNodeRef=\"vtkMRMLModelStorageNode15\" userTags=\"\" displayNodeRef=\"vtkMRMLModelDisplayNode18\"></Model>"<<std::endl;
  	MRMLFile<<"<LinearTransform"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLLinearTransformNode14\" name=\"LinearTransform10\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLTransformStorageNode6\" userTags=\"\" matrixTransformToParent=\"1 -0 0 -170  -0 1 -0 -12  0 -0 1 -0  -0 0 -0 1\"></LinearTransform>"<<std::endl;
	MRMLFile<<"<ModelDisplay"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelDisplayNode21\" name=\"vtkMRMLModelDisplayNode21\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"0.5 0.5 0.5\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0\" diffuse=\"1\" selectedSpecular=\"0.5\" specular=\"0\" power=\"1\" opacity=\"1\" visibility=\"false\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"true\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFilecustomLUT.txt\" activeScalarName=\"norm_distance_projections\" ></ModelDisplay>"<<std::endl;
	
	MRMLFile<<"<ModelStorage"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelStorageNode21\" name=\"vtkMRMLModelStorageNode16\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\""<<meanA<<"\" useCompression=\"1\" readState=\"0\" writeState=\"4\"></ModelStorage>"<<std::endl;
	MRMLFile<<"<ModelDisplay"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelDisplayNode24\" name=\"vtkMRMLModelDisplayNode24\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"1 0 0\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0\" diffuse=\"1\" selectedSpecular=\"0.6\" specular=\"0\" power=\"1\" opacity=\"0.4\" visibility=\"true\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"false\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFullRainbow\" activeScalarName=\"\" ></ModelDisplay>"<<std::endl;
	MRMLFile<<"<Model"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelNode24\" name=\"ShapeAnalysisModule_OutputFile_meanA\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" transformNodeRef=\"vtkMRMLLinearTransformNode21\" storageNodeRef=\"vtkMRMLModelStorageNode21\" userTags=\"\" displayNodeRef=\"vtkMRMLModelDisplayNode24\"></Model>"<<std::endl;
	
	MRMLFile<<"<LinearTransform"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLLinearTransformNode21\" name=\"LinearTransform13\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLTransformStorageNode9\" userTags=\"\" matrixTransformToParent=\"1 -0 0 -110  -0 1 -0 60  0 -0 1 -0  -0 0 -0 1\"></LinearTransform>"<<std::endl;

	MRMLFile<<"<ModelStorage"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelStorageNode16\" name=\"vtkMRMLModelStorageNode15\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\""<<meanA<<"\" useCompression=\"1\" readState=\"0\" writeState=\"4\"></ModelStorage>"<<std::endl;
	MRMLFile<<"<ModelDisplay"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelDisplayNode19\" name=\"vtkMRMLModelDisplayNode19\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"1 0 0\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0\" diffuse=\"1\" selectedSpecular=\"0.5\" specular=\"0\" power=\"1\" opacity=\"1\" visibility=\"true\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"false\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFullRainbow\" activeScalarName=\"\" ></ModelDisplay>"<<std::endl;
	MRMLFile<<"<Model"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelNode19\" name=\"ShapeAnalysisModule_OutputFile_meanA\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" transformNodeRef=\"vtkMRMLLinearTransformNode15\" storageNodeRef=\"vtkMRMLModelStorageNode15\" userTags=\"\" displayNodeRef=\"vtkMRMLModelDisplayNode19\"></Model>"<<std::endl;
	
	MRMLFile<<"<ModelStorage"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelStorageNode17\" name=\"vtkMRMLModelStorageNode17\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\""<<meanB<<"\" useCompression=\"1\" readState=\"0\" writeState=\"4\"></ModelStorage>"<<std::endl;
	MRMLFile<<"<ModelDisplay"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelDisplayNode20\" name=\"vtkMRMLModelDisplayNode20\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"0 0.501961 1\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0\" diffuse=\"1\" selectedSpecular=\"0.5\" specular=\"0\" power=\"1\" opacity=\"1\" visibility=\"true\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"false\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFullRainbow\" activeScalarName=\"\" ></ModelDisplay>"<<std::endl;
	MRMLFile<<"<Model"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelNode20\" name=\"ShapeAnalysisModule_OutputFile_meanB\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" transformNodeRef=\"vtkMRMLLinearTransformNode16\" storageNodeRef=\"vtkMRMLModelStorageNode17\" userTags=\"\" displayNodeRef=\"vtkMRMLModelDisplayNode20\"></Model>"<<std::endl;
	
	MRMLFile<<"<ModelDisplay"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelDisplayNode22\" name=\"vtkMRMLModelDisplayNode22\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"0.5 0.5 0.5\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0\" diffuse=\"1\" selectedSpecular=\"0.5\" specular=\"0\" power=\"1\" opacity=\"1\" visibility=\"false\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"true\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFilecustomLUT.txt\" activeScalarName=\"norm_distance_projections\" ></ModelDisplay>"<<std::endl;
	MRMLFile<<"<ModelStorage"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelStorageNode22\" name=\"vtkMRMLModelStorageNode17\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\""<<meanB<<"\" useCompression=\"1\" readState=\"0\" writeState=\"4\"></ModelStorage>"<<std::endl;
	MRMLFile<<"<ModelDisplay"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelDisplayNode25\" name=\"vtkMRMLModelDisplayNode25\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"0 0.501961 1\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0\" diffuse=\"1\" selectedSpecular=\"0.5\" specular=\"0\" power=\"1\" opacity=\"0.6\" visibility=\"true\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"false\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFullRainbow\" activeScalarName=\"\" ></ModelDisplay>"<<std::endl;
	MRMLFile<<"<Model"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLModelNode25\" name=\"ShapeAnalysisModule_OutputFile_meanB\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" transformNodeRef=\"vtkMRMLLinearTransformNode21\" storageNodeRef=\"vtkMRMLModelStorageNode22\" userTags=\"\" displayNodeRef=\"vtkMRMLModelDisplayNode25\"></Model>"<<std::endl;
	MRMLFile<<"<LinearTransform"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLLinearTransformNode15\" name=\"LinearTransform11\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLTransformStorageNode7\" userTags=\"\" matrixTransformToParent=\"1 -0 0 -61  -0 1 -0 65  0 -0 1 -0  -0 0 -0 1\"></LinearTransform>"<<std::endl;
	MRMLFile<<"<LinearTransform"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLLinearTransformNode16\" name=\"LinearTransform12\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLTransformStorageNode8\" userTags=\"\" matrixTransformToParent=\"1 -0 0 -61  -0 1 -0 -10  0 -0 1 -0  -0 0 -0 1\"></LinearTransform>"<<std::endl;
	MRMLFile<<"<FiducialList"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLFiducialListNode6\" name=\"FiducialListMeanB\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLFiducialListStorageNode4\" userTags=\"\" symbolScale=\"0\" symbolType=\"11\" textScale=\"6\" visibility=\"1\" color=\"0.4 1 1\" selectedcolor=\"0 0 0\" ambient=\"0\" diffuse=\"1\" specular=\"0\" power=\"1\" locked=\"0\" opacity=\"1\" fiducials=\""<<std::endl;
	MRMLFile<<"id Mean labeltext Mean B xyz -20.881771 19.564474 36.170582 orientationwxyz 0 0 0 1 selected 1 visibility 1\"></FiducialList>"<<std::endl;
	MRMLFile<<"<FiducialList"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLFiducialListNode4\" name=\"FiducialListMeanA\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLFiducialListStorageNode4\" userTags=\"\" symbolScale=\"0\" symbolType=\"11\" textScale=\"6\" visibility=\"1\" color=\"0.4 1 1\" selectedcolor=\"0 0 0\" ambient=\"0\" diffuse=\"1\" specular=\"0\" power=\"1\" locked=\"0\" opacity=\"1\" fiducials=\""<<std::endl;
	MRMLFile<<"id Means labeltext Mean A xyz -22.31890730 83.159927 18.968515 orientationwxyz 0 0 0 1 selected 1 visibility 1\"></FiducialList>"<<std::endl;
	MRMLFile<<"<FiducialList"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLFiducialListNode7\" name=\"FiducialListOverlays\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLFiducialListStorageNode4\" userTags=\"\" symbolScale=\"0\" symbolType=\"11\" textScale=\"6\" visibility=\"1\" color=\"0.4 1 1\" selectedcolor=\"0 0 0\" ambient=\"0\" diffuse=\"1\" specular=\"0\" power=\"1\" locked=\"0\" opacity=\"1\" fiducials=\""<<std::endl;
	MRMLFile<<"id Overlays labeltext Overlays(Mean A and B) xyz -69.488144 72.463066 17.649809 orientationwxyz 0 0 0 1 selected 1 visibility 1\"></FiducialList>"<<std::endl;
 	MRMLFile<<"<FiducialList"<<std::endl;
 	MRMLFile<<"id=\"vtkMRMLFiducialListNode2\" name=\"FiducialListOthers\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLFiducialListStorageNode1\" userTags=\"\" symbolScale=\"0\" symbolType=\"11\" textScale=\"6\" visibility=\"1\" color=\"0.4 1 1\" selectedcolor=\"0 0 0.25098\" ambient=\"0\" diffuse=\"1\" specular=\"0\" power=\"1\" locked=\"0\" opacity=\"1\" fiducials=\""<<std::endl;
 	MRMLFile<<"id FDR labeltext FDR-P xyz -129.448 5.60482 8.68362 orientationwxyz 0 0 0 1 selected 1 visibility 1"<<std::endl;
	MRMLFile<<"id Mean labeltext Mean Difference ( Min: "<<min<<" Max: "<<max<<") xyz 45.56 23.31 6.65557 orientationwxyz 0 0 0 1 selected 1 visibility 1"<<std::endl;
	MRMLFile<<"id Raw  labeltext Raw P xyz -162.785 -9.31485 -4.46765 orientationwxyz 0 0 0 1 selected 1 visibility 1\"></FiducialList>"<<std::endl;
	MRMLFile<<"<Slice"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLSliceNode1\" name=\"Green\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fieldOfView=\"304.687 250 1\" dimensions=\"390 320 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"-1 0 0 0 0 -0 1 127.24 0 1 0 0 0 0 0 1\" layoutName=\"Green\" orientation=\"Coronal\" jumpMode=\"1\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>"<<std::endl;
	MRMLFile<<"<Slice"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLSliceNode2\" name=\"Red\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fieldOfView=\"304.687 250 1\" dimensions=\"390 320 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"-1 0 0 0 0 1 0 0 0 0 1 -1.9086 0 0 0 1\" layoutName=\"Red\" orientation=\"Axial\" jumpMode=\"1\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>"<<std::endl;
	MRMLFile<<"<Slice"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLSliceNode3\" name=\"Yellow\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fieldOfView=\"304.687 250 1\" dimensions=\"390 320 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"-0 0 1 -17.719 -1 -0 0 0 0 1 -0 0 0 0 0 1\" layoutName=\"Yellow\" orientation=\"Sagittal\" jumpMode=\"1\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>"<<std::endl;
	MRMLFile<<"<SliceComposite"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLSliceCompositeNode3\" name=\"vtkMRMLSliceCompositeNode3\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" backgroundVolumeID=\"\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"0\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Yellow\" annotationMode=\"All\" doPropagateVolumeSelection=\"1\"></SliceComposite>"<<std::endl;
	MRMLFile<<"<SliceComposite"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLSliceCompositeNode2\" name=\"vtkMRMLSliceCompositeNode2\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" backgroundVolumeID=\"\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"0\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Red\" annotationMode=\"All\" doPropagateVolumeSelection=\"1\"></SliceComposite>"<<std::endl;
	MRMLFile<<"<SliceComposite"<<std::endl;
	MRMLFile<<"id=\"vtkMRMLSliceCompositeNode1\" name=\"vtkMRMLSliceCompositeNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" backgroundVolumeID=\"\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"0\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Green\" annotationMode=\"All\" doPropagateVolumeSelection=\"1\"></SliceComposite>"<<std::endl;
	

	MRMLFile<<"<ColorTableStorage"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableStorageNode11\"  name=\"vtkMRMLColorTableStorageNode11\"  hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\"customLUT_FDRP.txt\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></ColorTableStorage>"<<std::endl;
 	MRMLFile<<"<ColorTable"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableNodeFilecustomLUT_FDRP.txt\"  name=\"customLUT_FDRP.txt\"  description=\"A color table read in from a text file, each line of the format: IntegerLabel  Name  R  G  B  Alpha\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLColorTableStorageNode11\"  userTags=\"\" type=\"14\" numcolors=\"256\"></ColorTable>"<<std::endl;

	MRMLFile<<"<ColorTableStorage"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableStorageNode13\"  name=\"vtkMRMLColorTableStorageNode13\"  hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\"customLUT_RawP.txt\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></ColorTableStorage>"<<std::endl;
 	MRMLFile<<"<ColorTable"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableNodeFilecustomLUT_RawP.txt\"  name=\"customLUT_RawP.txt\"  description=\"A color table read in from a text file, each line of the format: IntegerLabel  Name  R  G  B  Alpha\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLColorTableStorageNode13\"  userTags=\"\" type=\"14\" numcolors=\"256\"></ColorTable>"<<std::endl;

	MRMLFile<<"<ColorTableStorage"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableStorageNode14\"  name=\"vtkMRMLColorTableStorageNode14\"  hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\"customLUT_normDistProjections.txt\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></ColorTableStorage>"<<std::endl;
 	MRMLFile<<"<ColorTable"<<std::endl;
  	MRMLFile<<"id=\"vtkMRMLColorTableNodeFilecustomLUT_normDistProjections.txt\"  name=\"customLUT_normDistProjections.txt\"  description=\"A color table read in from a text file, each line of the format: IntegerLabel  Name  R  G  B  Alpha\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLColorTableStorageNode14\"  userTags=\"\" type=\"14\" numcolors=\"256\"></ColorTable>"<<std::endl;


	MRMLFile<<"</MRML>"<<std::endl;
	}


}

void write_commandline_txt(std::string outbase)
{

	

}

	


	
