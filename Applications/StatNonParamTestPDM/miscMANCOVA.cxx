	#include "miscMANCOVA.h"
	
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
			unsigned int numGroupTypes,vector<int>  groupTypeColumns,vector<int> independentColumns, bool scaleOn, int scaleColumn, unsigned int &numSubjects, unsigned int &numA, unsigned int &numB,
			vnl_matrix<int>* &groupLabel, vnl_matrix<double>* &scaleFactor, 
			std::string* &meshFileNames,  vnl_matrix<double> * &indValue, bool computeScaleFactorFromVolumes )
	{
	
	
	cout<<"filename: "<<filename<<" numIndependent: "<<numIndependent<<endl;
	//bool interactionTest =false; // Replace
	//unsigned int testColumn=1;
	const int MAXLINE  = 5000; 
	static char line [ MAXLINE ];
	
	char *extension=strrchr(filename.c_str(),'.');
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
	char *extension =strrchr(meshFileNames[index].c_str(),'.');
	if(!strcmp(extension,".vtk"))
	{
		//read vtk file
		vtkPolyDataReader *mesh = vtkPolyDataReader::New();
		mesh->SetFileName(meshFileNames[index].c_str());
		
		try{
			mesh->Update();		
		}
		catch(itk::ExceptionObject ex)
		{
			std::cout<< "Error reading meshfile:  "<< meshFileNames[index] << std::endl << "ITK error: " << ex.GetDescription()<< std::endl;
		exit(-3);
		}
		
		vtkPolyData *polyData= vtkPolyData::New();
		polyData=mesh->GetOutput();
		
		//convert polydata to itk mesh
		char *meshFile;
		int length=strlen(meshFileNames[index].c_str())-4;
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
			std::cout<< "Error reading meshfile:  "<< meshFileNames[index] << std::endl << "ITK error: " << ex.GetDescription()<< std::endl;	
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
			std::cout<< "Error reading meshfile:  "<< meshFileNames[index] << std::endl << "ITK error: " << ex.GetDescription()<< std::endl;	
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
		
		for (unsigned int pointID = 0; pointID < numFeatures; pointID++) {
		PointType curPoint =  points->GetElement(pointID);
		for (unsigned int dim = 0; dim < 3; dim++) {
		if (scaleOn) 
		{
		(*featureValue)[index][pointID*3 +dim] = curPoint[dim] / (*scaleFactor)[index][0];
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
	
	for (unsigned int pointID=0;pointID<numFeatures;++pointID) 
	{
	if ( fdrThresh>0 ) 
	{
	fdrP(pointID)=(rawP[pointID]/fdrThresh)*fdrLevel;
	} 
	else 
	{
	fdrP(pointID) = 1.0;  // there is nothing signicant (sorry), so let's set the p-value to 1
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
