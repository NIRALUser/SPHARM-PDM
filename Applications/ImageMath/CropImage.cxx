/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: CropImage.cxx,v $
Language:  C++
Date:      $Date: 2009/10/23 14:11:14 $
Version:   $Revision: 1.9 $

Copyright (c) 2002 Insight Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include <itkImage.h>
#include <itkImageFileReader.h> 
#include <itkCastImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageFileWriter.h>
#include <itkMetaImageIO.h>

#include <itkExtractImageFilter.h>
#include <itkConstantPadImageFilter.h>

#include "argio.h"

using namespace itk;
using namespace std;

static int debug = 1;

int main(int argc, const char* argv[])
{


  if (argc <= 1 || ipExistsArgument(argv,"-h ") || ipExistsArgument(argv,"-help") || 
      ipExistsArgument(argv,"-usage"))
    {
      cout << "usage: " << argv[0] << endl
	   << "       infile [-o outfile] [-v] " << endl
	   << "       [-border n | -region px,py,pz,w,h,d | -autocrop <thresh>" << endl
	   << " -v           verbose mode" << endl
	   << " -border n    crop a border of n voxels uniformly form x,y,z" << endl
	   << " -region px,py,pz,w,h,d   crop image including starting point px,py,pz" << endl
	   << "              (startindex is 0) with dimensions width(w),height(h) " << endl
	   << "              and depth(d)" << endl
	   << " -embed px,py,pz,w,h,d   uncrop/embed/pad into an image at position px,py,pz" << endl
	   << "              (startindex is 0) and enlarged image dimensions are width(w),height(h) " << endl
	   << "              and depth(d)" << endl
	   << " -fEmbed <textfilename> same use as the embed option, but px,py,pz,w,h,d are written in the textfile, in that order" << endl
	   << " -autocrop  <thresh>  crop away background areas (below thresh)" << endl
	   << " -autoborder n  additional option to specify safety border around autocrop" << endl
	   << " -zpacks <#regions>  Uniformly subdivide z-axis into #regions images" << endl
	   << " -interl <#im> <axis> <filename> Create #im images taking 1 slice every #im along " << endl
	   << "                the \"axis\" axis in the original image (X=0, Y=1, Z=2(def))" << endl
	   << endl;
      exit(0) ;
    }

  // argument processing
  string fileName(argv[1]);

  debug =  ipExistsArgument(argv,"-v"); 
  string outfileName(ipGetStringArgument(argv, "-o", ""));

  if (outfileName.empty()) {
    outfileName = "output.mha";
    cout << "no outputname specified using " << outfileName << endl;
  } 
  if (debug) cout << "Output filename: " << outfileName << endl;
  bool borderOn = ipExistsArgument(argv,"-border"); 
  bool regionOn = ipExistsArgument(argv,"-region"); 
  bool embedOn = ipExistsArgument(argv,"-embed"); 
  bool fEmbedOn = ipExistsArgument(argv,"-fEmbed"); 
  bool autoCropOn = ipExistsArgument(argv,"-autocrop");
  bool zpacksOn = ipExistsArgument(argv,"-zpacks");
  bool interlOn = ipExistsArgument(argv,"-interl");

  int border = 0;
  const int numCropParam = 6;
  int cropParam[numCropParam];
  if (borderOn) 
    border = ipGetIntArgument(argv,"-border",0); 
  if (regionOn) {  
    char *tmp_str    = ipGetStringArgument(argv, "-region", NULL);
    int numDim       = ipExtractIntTokens(cropParam, tmp_str, numCropParam);
    if (numDim != numCropParam) {              
      cerr << argv[0] << ": region needs "<< numCropParam << " parameters.\n";
      exit(-1);
    }
    free(tmp_str);
  } 
  if (embedOn) {  
    char *tmp_str    = ipGetStringArgument(argv, "-embed", NULL);
    int numDim       = ipExtractIntTokens(cropParam, tmp_str, numCropParam);
    if (numDim != numCropParam) {              
      cerr << argv[0] << ": region needs "<< numCropParam << " parameters.\n";
      exit(-1);
    }
    free(tmp_str);
  } 
  if ((regionOn && embedOn)||(regionOn && fEmbedOn)) {
    cerr << "embed and region cannot be used at the same time" << endl;
    exit(-1);
  }  
  if (fEmbedOn) {
    string EmbedTextFileName = ipGetStringArgument(argv, "-fEmbed", "");
    ifstream file(EmbedTextFileName.c_str(),ios::in);
    if (file)
    {
	    int px,py,pz,w,h,d;
	    file >> px >> py >> pz >> w >> h >> d;
	    cropParam[0] = px;
	    cropParam[1] = py;
	    cropParam[2] = pz;
	    cropParam[3] = w;
	    cropParam[4] = h;
	    cropParam[5] = d;
    }
    else
    {
	    cerr << "opening file error" << endl;
    }
  }

  int nbinter = 0;
  int axisinter = 0;
  std::string coilid_filename;
  if(interlOn){
    int maxnbarg = 3;
    char** values = new char*[3];
    for(int i = 0 ; i < 3 ; i++)
      values[i] = new char[50];
    //Get the arguments associated with interl option
    ipGetStringMultipArgument(argv,"-interl",values,maxnbarg);
    //Give the number of pair ID/name
    nbinter = atoi(values[0]);
    axisinter = atoi(values[1]);
    coilid_filename = values[2];    

    if(axisinter < 0 || axisinter > 2) {
      std::cerr << "the axis parameter is either 0,1 or 2" << std::endl;
      return -1;
    }
    if(nbinter < 1) {
      std::cerr << "The # of images has to be a postive number greater than 0" << std::endl;
      return -1;
    }    
    if(coilid_filename.empty()) {
      std::cerr << "File name missing" << std::endl;
      return -1;
    }
  }

  const int autoCropThresh = ipGetIntArgument(argv,"-autocrop",127);
  const int autoBorder = ipGetIntArgument(argv,"-autoborder",1);
  const int zpackNumber = ipGetIntArgument(argv,"-zpacks",1);

  // argument processing done, do real work

  // Defitions
  const int Dimension = 3;
  typedef short ImagePixelType;

  typedef itk::Point< double, Dimension > PointType;

  typedef Image<ImagePixelType,Dimension>  ImageType;
  typedef ImageType::RegionType ImageRegionType;

  typedef ImageFileReader< ImageType > VolumeReaderType;
  typedef ImageFileWriter<ImageType> VolumeWriterType;

  typedef ExtractImageFilter<ImageType, ImageType> CropFilterType;
  typedef ConstantPadImageFilter<ImageType, ImageType> EmbedFilterType;
  
  typedef ImageRegionConstIterator< ImageType > ConstIterator;
  typedef ImageRegionIterator< ImageType > Iterator;

  // read file
  VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
  if (debug) cout << "Loading file " << endl;
  try
    {
      imageReader->SetFileName(fileName.c_str()) ;
      imageReader->Update() ;
    }
  catch (ExceptionObject e)
    {
      e.Print(cout) ;
      exit(0) ;
    }

  ImageRegionType extractRegion, fullRegion;
  if (zpacksOn) {
    fullRegion = (imageReader->GetOutput())->GetLargestPossibleRegion();
    if (debug) cout << "size of the original image " << fullRegion.GetSize(0) 
		    << "," << fullRegion.GetSize(1) << "," << fullRegion.GetSize(2) << endl; 

    extractRegion.SetIndex(0,0);
    extractRegion.SetIndex(1,0);
    extractRegion.SetSize(0,fullRegion.GetSize(0));
    extractRegion.SetSize(1,fullRegion.GetSize(1));

    int NumSlices = fullRegion.GetSize(2) / zpackNumber;

    if (debug) cout << "cropping zpacks: " << zpackNumber <<"*" << NumSlices <<   endl;
    char * origName = strdup(outfileName.c_str());
    char * dotPtr = strchr(origName,'.');
    dotPtr[0] = '\0'; dotPtr++;

    for (int pack = 0; pack < zpackNumber; pack++) {
      extractRegion.SetIndex(2, pack * NumSlices);
      extractRegion.SetSize(2,NumSlices);
      CropFilterType::Pointer cropFilter = CropFilterType::New();
      cropFilter->SetInput(imageReader->GetOutput());
      cropFilter->SetExtractionRegion(extractRegion);
      cropFilter->Update();
      
      if (debug) cout << "writing output data" << endl;
      VolumeWriterType::Pointer writer = VolumeWriterType::New();
      
      char outfileP[1024];
      sprintf(outfileP, "%s_P%03d.%s",origName,pack, dotPtr);
      
      writer->SetFileName(outfileP); 
      writer->SetInput(cropFilter->GetOutput());
      writer->Write();
    }
  }else if(interlOn){

    //Check if the input image  number of slices along the specified axis is a multiple of nbinterl
    int nbslices; 
    int dim[Dimension];
    ImageType::RegionType imageRegion = imageReader->GetOutput()->GetLargestPossibleRegion(); 
    dim[0] = imageRegion.GetSize(0);
    dim[1] = imageRegion.GetSize(1);
    dim[2] = imageRegion.GetSize(2);
    
    nbslices = dim[axisinter];
    if(debug) std::cout << "nbslice: " << nbslices << " nbinter: " << nbinter << std::endl;
    float remainer = std::fmod(static_cast<float>(nbslices), static_cast<float>(nbinter));

    //If the input number is not a mutliple of the slice number
    if(remainer != 0) {
      std::cerr << "The number of slice (" << nbslices << ") along this axis is not a multiple of " << nbinter << std::endl;
    }else{

      //Read the coilid file
      std::vector<int> coilid_list;
      std::ifstream coilid_file;
      coilid_file.open(coilid_filename.c_str(),std::ifstream::in);
      while(coilid_file.good())
	{
	  char lineval[30];
	  coilid_file.getline(lineval,30);
	  std::string strline = lineval; 
	  //make sure the line is not empty
	  if(strline.size())
	    {
	      int id = -1;
	      int num = 0;
	      sscanf(strline.c_str(),"%d %*s %d",&num,&id);
	      coilid_list.push_back(id);
	    }
	}
      coilid_file.close();

      if(nbslices != coilid_list.size())
	{
	  std::cerr << "Coil ID list has a different number of line than the number of slices (" << coilid_list.size() << " Vs. " << nbslices << ")" << std::endl;
	  return -1;
	}

      //Compute the number of slice for each output image
      int sizealongaxis = nbslices / nbinter;
      std::cout  << "Div: " << sizealongaxis << " = " << nbslices << " / " << nbinter << std::endl;
      //Set the output image
      ImageType::RegionType outputRegion;
      ImageType::RegionType::IndexType outputStart;
      outputStart[0] = 0;
      outputStart[1] = 0;
      outputStart[2] = 0;	
      ImageType::SizeType outputSize;
      for(int d = 0 ; d < Dimension ; d++) {
	if(d == axisinter)
	  outputSize[d] = sizealongaxis;
	else
	  outputSize[d] = dim[d];	   
      }

      outputRegion.SetIndex(outputStart);
      outputRegion.SetSize(outputSize);
      
      PointType outputOrigin;
      outputOrigin[0] = 0;
      outputOrigin[1] = 0;
      outputOrigin[2] = 0;
      
      ImageType::SpacingType outputSpacing;
      outputSpacing[0] = imageReader->GetOutput()->GetSpacing()[0];
      outputSpacing[1] = imageReader->GetOutput()->GetSpacing()[1];
      outputSpacing[2] = imageReader->GetOutput()->GetSpacing()[2];

      //Create the output images
      std::vector<ImageType::Pointer> outimages_list;           
      for(int im = 0 ; im < nbinter+1 ; im++) {

	ImageType::Pointer curr_image = ImageType::New();
	ImageType::Pointer outputImage = ImageType::New();
	curr_image->SetRegions(outputRegion);
	curr_image->SetSpacing(outputSpacing);
	curr_image->SetOrigin(outputOrigin);
	curr_image->Allocate();
	
	outimages_list.push_back(curr_image);
      }

      if(debug) std::cout << "Output images, initialized" << std::endl;
      //Create the global image from all the slices. Initialization
      //Iterator init 
      ImageType::Pointer global_outimage = outimages_list.at(nbinter);
      Iterator globoutputIt(global_outimage, global_outimage->GetLargestPossibleRegion());
      globoutputIt.GoToBegin();
      //Individual slice init
      int ind = 0;
      int vecdim[2];
      for(int d = 0 ; d < Dimension ; d++) {
	if(d != axisinter) {
	  vecdim[ind] = dim[d];
	  ind++;
	  }	
      }
      std::vector<double> globalslice;

      ImageRegionType extractRegion;
      //Now that the output images are created, we fill them up.
      //For each slice of the original image
      for(int slice = 0 ; slice < dim[axisinter] ; slice++) {
	
	//ID of the slice in the output volume
	int curr_slice = static_cast<int>(floor(slice/nbinter));

	int outimageID = static_cast<int>(fmod(static_cast<float>(slice), static_cast<float>(nbinter)));
	
	//reset the global image pixel value
	if(outimageID == 0){
	  globalslice.clear();
	  for(int i = 0 ; i < vecdim[0] * vecdim[1] ; i++)
	    globalslice.push_back(0);
	}

	//	ImageType::Pointer curr_outimage = outimages_list.at(outimageID);
	//We use the sequence of slices from the coilid file (listed in the coilid_list vector)
	ImageType::Pointer curr_outimage = outimages_list.at(coilid_list.at(slice)-1);

	if(debug) std::cout << "Slice: " << slice << std::endl;

	//Set the ID of the slice to extract
	for(int d = 0 ; d < Dimension ;  d++) {
	  if(d == axisinter) {
	    extractRegion.SetIndex(d,slice);
	    extractRegion.SetSize(d,1);
	  } else {
	    extractRegion.SetIndex(d,0);
	    extractRegion.SetSize(d,dim[d]);
	  }
	}

	//Slice extraction
	CropFilterType::Pointer cropFilter = CropFilterType::New();
	cropFilter->SetInput(imageReader->GetOutput());
	cropFilter->SetExtractionRegion(extractRegion);
	cropFilter->Update();	 	

	//Set the iterator over the extracted slice
	ConstIterator inputIt(cropFilter->GetOutput(),cropFilter->GetOutput()->GetLargestPossibleRegion());
	Iterator outputIt(curr_outimage, curr_outimage->GetLargestPossibleRegion());

	//Iterator of the global image
	std::vector<double>::iterator globalpix;
	globalpix = globalslice.begin();

	//Copy the voxels
	for(inputIt.GoToBegin() ; !inputIt.IsAtEnd() ; ++inputIt) {       
	  ImageRegionType::IndexType outIndex;
	  
	  for(int d = 0 ; d < Dimension ; d++) {	    
	    if(d == axisinter)
	      outIndex[d] = curr_slice;
	    else
	      outIndex[d] = inputIt.GetIndex()[d];
	  }
	  outputIt.SetIndex(outIndex);	    
	  outputIt.Set(inputIt.Get());
	  (*globalpix) += static_cast<double>(inputIt.Get() * inputIt.Get());
	  
	  //If it's the last occurence of the sequence of nbinter slices, we copy the voxel values
	  //here we save the image which is the average of the nbinter images (Pi = sqrt(p1^2 + p2^2 +...) )
	  if(outimageID == nbinter-1)
	    {
	      globoutputIt.SetIndex(outIndex);
	      globoutputIt.Set(static_cast<ImagePixelType>(sqrt((*globalpix))));		
	    }
	  ++globalpix;
	}	  	
      }   

      //Save the output images 
      //Get the extension:
      bool zipfile = false;
      char *zip = new char[10];
      strcpy(zip,"");
      std::string::size_type dotpos = outfileName.find_last_of(".",outfileName.size());
      std::string extension = outfileName.substr(dotpos+1, outfileName.size() - dotpos);
      //If the extension contains a gz, then we get rid of it to get the real extension, and add it when creating the file name
      if(extension == "gz") {
	outfileName = outfileName.substr(0,dotpos);
	zipfile = true;
	strcpy(zip,".gz");
      }

      char * origName = strdup(outfileName.c_str());
      char * dotPtr = strrchr(origName,'.');
      dotPtr[0] = '\0'; dotPtr++;

      for(int im = 0 ; im < nbinter+1 ; im++) {
	//Create writer pointer
	VolumeWriterType::Pointer writer = VolumeWriterType::New();
	
	char outfileP[1024];
	if(im != nbinter)
	  sprintf(outfileP, "%s_P%03d.%s%s",origName,im, dotPtr,zip);
	else
	  sprintf(outfileP, "%s_Paverage.%s%s",origName, dotPtr,zip);

	writer->SetFileName(outfileP);
	writer->SetInput(outimages_list.at(im));
	writer->Write();
      }
      delete []zip;
   
    } //endif remainer

  }else if ((embedOn)||(fEmbedOn)) {
    int dim[Dimension];
    ImageType::RegionType imageRegion = imageReader->GetOutput()->GetLargestPossibleRegion();
    dim[0] = imageRegion.GetSize(0);
    dim[1] = imageRegion.GetSize(1);
    dim[2] = imageRegion.GetSize(2);

    unsigned long upperfactors[3];
    unsigned long lowerfactors[3];
    lowerfactors[0] = cropParam[0]; lowerfactors[1] = cropParam[1]; lowerfactors[2] = cropParam[2];
    upperfactors[0] = cropParam[3]-dim[0]-lowerfactors[0]; 
    upperfactors[1] = cropParam[4]-dim[1]-lowerfactors[1]; 
    upperfactors[2] = cropParam[5]-dim[2]-lowerfactors[2]; 
    if (debug) cout << "embedding " << endl;

    EmbedFilterType::Pointer embedFilter = EmbedFilterType::New();
    embedFilter->SetInput(imageReader->GetOutput());
    embedFilter->SetConstant(0);  
    embedFilter->SetPadLowerBound(lowerfactors);
    embedFilter->SetPadUpperBound(upperfactors);
    embedFilter->UpdateLargestPossibleRegion();
    
    if (debug) cout << "writing output data " << outfileName << endl;
    VolumeWriterType::Pointer writer = VolumeWriterType::New();
    writer->SetFileName(outfileName.c_str()); 
    writer->SetInput(embedFilter->GetOutput());
    writer->Write();

  } else {
    if (regionOn) {
      extractRegion.SetIndex(0,cropParam[0]);
      extractRegion.SetIndex(1,cropParam[1]);
      extractRegion.SetIndex(2,cropParam[2]);
      extractRegion.SetSize(0,cropParam[3]);
      extractRegion.SetSize(1,cropParam[4]);
      extractRegion.SetSize(2,cropParam[5]);
    } else if (borderOn) {
      extractRegion = (imageReader->GetOutput())->GetLargestPossibleRegion();
      extractRegion.SetIndex(0,extractRegion.GetIndex(0) + border);
      extractRegion.SetIndex(1,extractRegion.GetIndex(1) + border);
      extractRegion.SetIndex(2,extractRegion.GetIndex(2) + border);
      extractRegion.SetSize(0, extractRegion.GetSize(0) - 2 * border);
      extractRegion.SetSize(1, extractRegion.GetSize(1) - 2 * border);
      extractRegion.SetSize(2, extractRegion.GetSize(2) - 2 * border);
    } else if (autoCropOn) {
      //autoCropThresh;
      unsigned int minx, miny, minz, maxx, maxy, maxz;
      minx = miny = minz = 100000;
      maxx = maxy = maxz = 0;
      
      // first grow for all seed (constrained) and then do connected component labeling 
      Iterator iter (imageReader->GetOutput(), imageReader->GetOutput()->GetLargestPossibleRegion());
      while ( !iter.IsAtEnd() )
	{
	  Iterator::IndexType ind = iter.GetIndex();
	  ImagePixelType value =  iter.Get();
	  if (value >= autoCropThresh) {
	    if (ind[0] < minx) minx = ind[0];
	    if (ind[1] < miny) miny = ind[1];
	    if (ind[2] < minz) minz = ind[2];
	    if (ind[0] > maxx) maxx = ind[0];
	    if (ind[1] > maxy) maxy = ind[1];
	    if (ind[2] > maxz) maxz = ind[2];
	  }
	  ++iter;
	}
      if (!maxx && !maxy && !maxz) {
	cerr << "image has no pixels larger than " << autoCropThresh << endl;
	exit(-1);
      }
      // enlarge by one voxel border
      ImageType::RegionType imageRegion = imageReader->GetOutput()->GetLargestPossibleRegion();
      if (minx >= autoBorder) minx -= autoBorder;
      else minx = 0;
      if (miny >= autoBorder) miny -= autoBorder;
      else miny = 0;
      if (minz >= autoBorder) minz -= autoBorder;
      else minz = 0;
      if (maxx < imageRegion.GetSize(0)-autoBorder) maxx += autoBorder;
      else maxx = imageRegion.GetSize(0) - 1;
      if (maxy < imageRegion.GetSize(1)-autoBorder) maxy += autoBorder;
      else maxy = imageRegion.GetSize(1) - 1;
      if (maxz < imageRegion.GetSize(2)-autoBorder) maxz += autoBorder;
      else maxz = imageRegion.GetSize(2) - 1;
      
      extractRegion.SetIndex(0,minx);
      extractRegion.SetIndex(1,miny);
      extractRegion.SetIndex(2,minz);
      extractRegion.SetSize(0,maxx-minx+1);
      extractRegion.SetSize(1,maxy-miny+1);
      extractRegion.SetSize(2,maxz-minz+1);
    }
    
    int dim[Dimension];
    ImageType::RegionType imageRegion = imageReader->GetOutput()->GetLargestPossibleRegion();
    dim[0] = imageRegion.GetSize(0);
    dim[1] = imageRegion.GetSize(1);
    dim[2] = imageRegion.GetSize(2);
    if (debug) cout << "size of the original image " << dim[0] 
		    << "," << dim[1] << "," << dim[2]<< endl; 
    if (debug) cout << "cropping (x,y,z,w,h,d) " << extractRegion.GetIndex(0) << ","
		    << extractRegion.GetIndex(1) << "," 
		    << extractRegion.GetIndex(2) << "," 
		    << extractRegion.GetSize(0) << "," 
		    << extractRegion.GetSize(1) << "," 
		    << extractRegion.GetSize(2) <<  endl;
    
    CropFilterType::Pointer cropFilter = CropFilterType::New();
    cropFilter->SetInput(imageReader->GetOutput());
    cropFilter->SetExtractionRegion(extractRegion);
    //cropFilter->ReleaseDataFlagOn();
    cropFilter->Update();
    
    if (debug) cout << "writing output data" << outfileName << endl;
    VolumeWriterType::Pointer writer = VolumeWriterType::New();
    writer->SetFileName(outfileName.c_str()); 
    writer->SetInput(cropFilter->GetOutput());
    writer->Write();
  }

  return 0 ;
}

