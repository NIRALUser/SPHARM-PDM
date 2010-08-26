#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>

#include <itkImage.h>
#include <itkImageFileReader.h> 
#include <itkCastImageFilter.h>
#include <itkHistogram.h>
#include <itkImageRegionIterator.h>
#include <itkThresholdImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h> 

#include <itkCurvatureFlowImageFilter.h> 
#include <itkDiscreteGaussianImageFilter.h> 
#include <itkBinaryBallStructuringElement.h> 
#include <itkBinaryCrossStructuringElement.h>
#include <itkBinaryDilateImageFilter.h> 
#include <itkBinaryErodeImageFilter.h> 
#include <itkBinaryThresholdImageFilter.h>

#include <itkResampleImageFilter.h>

#include <itkMinimumMaximumImageCalculator.h> 
#include <itkBSplineInterpolateImageFunction.h> 
#include <itkAffineTransform.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

#include <itkCastImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkAntiAliasBinaryImageFilter.h>

#include <itkImageFileWriter.h>
#include <itkMetaImageIO.h>

#include <itkAddImageFilter.h>

#include <itkImageRegionIteratorWithIndex.h>
#include <itkIndex.h>

#include <math.h>
#include <algorithm>
#include "itkListSample.h"
#include "itkVector.h"
#include "itkScalarImageToListAdaptor.h"
#include "itkMembershipSample.h"
#include "itkCovarianceCalculator.h"

#include "SegPostProcessCLPCLP.h"

using namespace itk;
using namespace std;

static int debug = 1;
static void draw_fill_inside_image(unsigned short *image, int *dim, int new_label);
static int NoDiagConnect (unsigned short *image, int *dim);
static void clear_edge(unsigned short *image, int *dims, int clear_label);
bool searchList(int label, vector<int> list);



int main( int argc, char * argv[] )
{

  PARSE_ARGS;
  //if (debug) cerr<<"begins"<<endl;
  const int Dimension = 3;
  const int NO_LABEL = 0;
  
  double variance[Dimension];
  for(int i=0; i<Dimension; i++) variance[i]=variance_vect[i];
  
  double  spacing[Dimension];
  if ( scaleOn == true ) {
    for(int i=0; i<Dimension;i++) spacing[i]=spacing_vect[i];
    if (isotropicOn) {
      isotropicOn = false;
      cerr << "--space option overrides -isotropic option" << endl;
    }
  }
  
  
  //if (debug) cerr<<"typedef"<<endl;
// types used in this routine
  typedef double  CoordRepType;
  typedef unsigned short ImagePixelType;
  typedef Image<ImagePixelType,Dimension>  ImageType;
  typedef float SmoothImagePixelType;
  typedef ImageType::SizeType ImageSizeType;
  typedef Image<SmoothImagePixelType,Dimension> SmoothImageType;
  typedef ImageFileReader<ImageType> VolumeReaderType;
  typedef ImageFileWriter<ImageType> VolumeWriterType;

  typedef ResampleImageFilter<ImageType , ImageType > ResamplerType;
  typedef AffineTransform<CoordRepType, Dimension> TransformType;

  typedef BSplineInterpolateImageFunction<ImageType, CoordRepType> SplineInterpolFunctionType;
  typedef LinearInterpolateImageFunction<ImageType, CoordRepType> LinearInterpolFunctionType;
  typedef NearestNeighborInterpolateImageFunction<ImageType, CoordRepType> NNInterpolFunctionType;

  typedef ThresholdImageFilter< ImageType > threshFilterType;
  typedef BinaryThresholdImageFilter< ImageType, ImageType > binThreshFilterType;
  typedef BinaryBallStructuringElement<ImagePixelType,Dimension> StructuringElementType;
  typedef BinaryCrossStructuringElement<ImagePixelType,Dimension> AsymStructuringElementType;

  typedef BinaryDilateImageFilter<ImageType, ImageType, StructuringElementType> dilateFilterType;
  typedef BinaryErodeImageFilter<ImageType, ImageType, StructuringElementType> erodeFilterType;
  typedef BinaryErodeImageFilter<ImageType, ImageType, AsymStructuringElementType> erodeAsymFilterType;
  typedef ConnectedComponentImageFilter<ImageType, ImageType> CCLFilterType;
  typedef RelabelComponentImageFilter<ImageType, ImageType> CCLRelabelFilterType;

  typedef CastImageFilter<ImageType,SmoothImageType> castInputFilterType;
  typedef CastImageFilter<SmoothImageType,ImageType> castOutputFilterType;
  typedef DiscreteGaussianImageFilter<SmoothImageType, SmoothImageType> gaussFilterType;

  typedef CastImageFilter<ImageType ,SmoothImageType > CastToRealFilterType;
  typedef AntiAliasBinaryImageFilter<SmoothImageType, SmoothImageType> AntiAliasFilterType;
  typedef RescaleIntensityImageFilter<SmoothImageType, ImageType > RescaleFilter;

  typedef ImageRegionIteratorWithIndex< ImageType > IteratorType;
  typedef Index<Dimension> IndexType;

  typedef short PixelType;
  typedef ImageFileReader< ImageType >          VolumeReaderType;
  typedef ImageRegionIterator< ImageType >      Iterator;
  typedef ImageType::Pointer                    ImagePointer;
  typedef itk::Statistics::ScalarImageToListAdaptor<ImageType>	ImageSampleType;
  typedef ImageSampleType::Iterator				ImageSampleIterator;
  typedef itk::Statistics::MembershipSample< ImageSampleType >	MembershipSampleType;
  typedef itk::Statistics::CovarianceCalculator< MembershipSampleType::ClassSampleType >	CovarianceAlgorithmType;
  typedef itk::Vector< float, 3 > VectorType;
  
  
  typedef AddImageFilter< ImageType, ImageType,  ImageType > addFilterType;
  typedef MaskImageFilter< ImageType, ImageType, ImageType >  maskFilterType;

  ImageType::Pointer image ;
  VolumeReaderType::Pointer labelReader = VolumeReaderType::New();
  ImageType::Pointer maskImage ;
  ImageType::Pointer labelImage ;
  ImageType::Pointer EMSImage ;
  ImageType::Pointer strippedImage ;
  
  //bool skullstrippingOn=0 , savemaskOn=0, deleteVesselsOn=0;
  bool skullstrippingOn=0 , savemaskOn=0;
 /* if(!delVesselCorrectedImName.empty())
  deleteVesselsOn=1;*/
  
     //threshold for skullstripping
  if(!notstrippedfileName.empty())
    skullstrippingOn=1;
  if(skullstrippingOn)
  {
    if(!maskfileName.empty())
      savemaskOn=1;
    if(!fileName.empty()) 
    {
    if (debug) cout << "Loading labelfile " << endl;
    try
      {
	labelReader->SetFileName(fileName.c_str()) ;
	labelReader->Update() ;
	labelImage = labelReader->GetOutput() ;
      }
    catch (ExceptionObject e)
      {
	e.Print(cout) ;
	exit(0) ;
      }
    if (debug) cout << "Labelfile loaded " << endl;
    if (debug) cout << "Thresholding " << endl;
    //extract gray matter
    binThreshFilterType::Pointer binthreshFilterGM = binThreshFilterType::New();
    binThreshFilterType::Pointer binthreshFilterWM = binThreshFilterType::New();
    binThreshFilterType::Pointer binthreshFilterMWM = binThreshFilterType::New();
    binThreshFilterType::Pointer binthreshFilterCSF = binThreshFilterType::New();
    binthreshFilterGM->SetInput(labelImage);
    binthreshFilterGM->SetLowerThreshold(gm);
    binthreshFilterGM->SetUpperThreshold(gm);
    binthreshFilterGM->SetOutsideValue (0);
    binthreshFilterGM->SetInsideValue (1);
    try {
      binthreshFilterGM->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }        
    //extract white matter
    binthreshFilterWM->SetInput(labelImage);
    binthreshFilterWM->SetLowerThreshold(wm);
    binthreshFilterWM->SetUpperThreshold(wm);
    binthreshFilterWM->SetOutsideValue (0);
    binthreshFilterWM->SetInsideValue (1);
    try {
      binthreshFilterWM->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }        
    //extract csf
    binthreshFilterCSF->SetInput(labelImage);
    binthreshFilterCSF->SetLowerThreshold(csf);
    binthreshFilterCSF->SetUpperThreshold(csf);
    binthreshFilterCSF->SetOutsideValue (0);
    binthreshFilterCSF->SetInsideValue (1);
    try {
      binthreshFilterCSF->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
      }

    if (mwm > 0)
      {
	// extract myelinated WM
	binthreshFilterMWM->SetInput(labelImage);
	binthreshFilterMWM->SetLowerThreshold(mwm);
	binthreshFilterMWM->SetUpperThreshold(mwm);
	binthreshFilterMWM->SetOutsideValue (0);
	binthreshFilterMWM->SetInsideValue (1);
	try {
	  binthreshFilterMWM->Update();
	}
	catch (ExceptionObject & err) {
	  cerr << "ExceptionObject caught!" << endl;
	  cerr << err << endl;
	  return EXIT_FAILURE;	
	}
    
      }
    //add gm, wm and csf to have the binary mask
    addFilterType::Pointer addFilter1 = addFilterType::New();
    addFilterType::Pointer addFilter2 = addFilterType::New();
    addFilter1->SetInput1(binthreshFilterGM->GetOutput());
    addFilter1->SetInput2(binthreshFilterWM->GetOutput());
      try {
	addFilter1->Update();
      }
      catch (ExceptionObject & err) {
	cerr << "ExceptionObject caught!" << endl;
	cerr << err << endl;
	return EXIT_FAILURE;	
      }   
      
      addFilter2->SetInput1(binthreshFilterCSF->GetOutput());
      addFilter2->SetInput2(addFilter1->GetOutput());
      try {
	addFilter2->Update();
      }
      catch (ExceptionObject & err) {
	cerr << "ExceptionObject caught!" << endl;
	cerr << err << endl;
	return EXIT_FAILURE;	
      }   

      if (mwm > 0)
	{
	  // if myelinated white matter
	  addFilterType::Pointer addFilter3 = addFilterType::New();
	  
	  addFilter3->SetInput1(binthreshFilterMWM->GetOutput());
	  addFilter3->SetInput2(addFilter2->GetOutput());
	  try {
	    addFilter3->Update();
	  }
	  catch (ExceptionObject & err) {
	    cerr << "ExceptionObject caught!" << endl;
	    cerr << err << endl;
	    return EXIT_FAILURE;	
	  }
	  maskImage=addFilter3->GetOutput();
	}
      else
	maskImage=addFilter2->GetOutput();	 
    }
    if (debug) cout << "Threshold finished " << endl;
    
    //dilate option
    if(dilateOn){      
      StructuringElementType structuringElement;
      structuringElement.SetRadius( 1 ); 
      structuringElement.CreateStructuringElement( );
    
      dilateFilterType::Pointer dilateFilter = dilateFilterType::New(); 
      dilateFilter->SetInput(maskImage);
      dilateFilter->SetDilateValue (1);
      dilateFilter->SetKernel( structuringElement );
      try {
	dilateFilter->Update();
      }
      catch (ExceptionObject & err) {
	cerr << "ExceptionObject caught!" << endl;
	cerr << err << endl;
    	return EXIT_FAILURE;	
      }   
      
      maskImage = dilateFilter->GetOutput();
      
    }
    VolumeReaderType::Pointer EMSReader = VolumeReaderType::New();
    if (debug) cout << "Loading EMSfile " << endl;
    try{
	     EMSReader->SetFileName(notstrippedfileName.c_str()) ;
	     EMSReader->Update() ;
	     EMSImage = EMSReader->GetOutput() ;
     }
     catch (ExceptionObject e){
	     e.Print(cout) ;
	     exit(0) ;
     }
     if (debug) cout << "EMSfile loaded " << endl;
    }
  
  if (debug) cout << "Loading file " << fileName << endl;
  try
    {
      if(skullstrippingOn)
	image=maskImage;
      else{
	VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
	imageReader->SetFileName(fileName.c_str()) ;
	imageReader->Update() ;
	image = imageReader->GetOutput() ;
      }
    }
  catch (ExceptionObject e)
    {
      e.Print(cout) ;
      exit(0) ;
    }
  
  int dim[Dimension];

  // HERE the processing starts
  
  // prepare for my own routines
  ImageType::RegionType imageRegion = image->GetBufferedRegion();
  dim[0] = imageRegion.GetSize(0);
  dim[1] = imageRegion.GetSize(1);
  dim[2] = imageRegion.GetSize(2);
  ImageType::IndexType nullIndex;
  nullIndex[0] = 0;
  nullIndex[1] = 0;
  nullIndex[2] = 0;
  if (debug) cout << "size of the image " << dim[0] << "," << dim[1] << "," << dim[2]<< endl; 

  const int elemSize = 1;
  const int LABEL_VAL = 255;
  ImageType::Pointer procImage ;

  if (label != NO_LABEL ) {
    // Thresholding at label -> set to 0 otherwise
    if (debug) cout << "extracting object " << label << endl; 
    
    threshFilterType::Pointer threshFilter = threshFilterType::New();
    threshFilter->SetInput(image);
    threshFilter->ThresholdOutside(label, label);
    threshFilter->SetOutsideValue (0);
    threshFilter->Update();
    image = threshFilter->GetOutput();
    
  }
  {
    // Set object to 255 -> Tresholding all but label
    if (debug) cout << "set object to " << LABEL_VAL  << endl; 
    threshFilterType::Pointer threshFilter = threshFilterType::New();
    threshFilter->SetInput(image);
    threshFilter->ThresholdAbove(0);
    threshFilter->SetOutsideValue (LABEL_VAL);
    threshFilter->Update();
    procImage = threshFilter->GetOutput();
  }

  // clear edge and fill inside 
  if (debug) cout << "fill inside " << endl; 

  imageRegion = procImage->GetBufferedRegion();
  dim[0] = imageRegion.GetSize(0);
  dim[1] = imageRegion.GetSize(1);
  dim[2] = imageRegion.GetSize(2);
  ImagePixelType *data = &((*procImage)[nullIndex]);
  clear_edge(data, dim, 0);
  draw_fill_inside_image(data, dim, LABEL_VAL);


  // Rescaling the data ?
  if (isotropicOn) {
    const itk::Vector<double,3> inSpacing = procImage->GetSpacing();
    double minSpacing = inSpacing[0];
    for (int i = 1; i < Dimension; i++) if (minSpacing > inSpacing[i]) minSpacing = inSpacing[i];
    for (int i = 0; i < Dimension; i++) spacing[i] = minSpacing;
    scaleOn = true;
  }

  if (scaleOn) {
    if (linearOn) {
      if (debug) {
	cout << " resampling (Linear interpol) the data to voxelsize ";
	for (int i = 0; i < Dimension; i++) cout << spacing[i] << ", ";
	cout << endl;
      }
      ResamplerType::Pointer resampleFilter = ResamplerType::New();
      TransformType::Pointer transform = TransformType::New();
      //SplineInterpolFunctionType::Pointer interpolFunction = SplineInterpolFunctionType::New();
      LinearInterpolFunctionType::Pointer interpolFunction = LinearInterpolFunctionType::New();

      transform->SetIdentity();
      interpolFunction->SetInputImage(procImage);
      //interpolFunction->SetSplineOrder(3);

      ResamplerType::SizeType size;
      imageRegion = procImage->GetBufferedRegion();
      dim[0] = imageRegion.GetSize(0);
      dim[1] = imageRegion.GetSize(1);
      dim[2] = imageRegion.GetSize(2);
      const itk::Vector<double,3> inSpacing = image->GetSpacing();
      for (int i = 0; i < Dimension; i++) {
	size[i] = (unsigned int) ((double) inSpacing[i] * dim[i] / spacing[i] );
      }
      resampleFilter->SetInput(procImage); 
      resampleFilter->SetSize(size);
      resampleFilter->SetOutputSpacing(spacing);
      resampleFilter->SetOutputOrigin(image->GetOrigin());
      resampleFilter->SetTransform(transform.GetPointer());
      resampleFilter->SetInterpolator(interpolFunction.GetPointer());
      resampleFilter->Update();

      if (debug) cout << "Tresholding at 50% " << endl;
      binThreshFilterType::Pointer binTreshFilter = binThreshFilterType::New();
      binTreshFilter->SetInput(resampleFilter->GetOutput());
      binTreshFilter->SetUpperThreshold(LABEL_VAL);
      binTreshFilter->SetLowerThreshold(LABEL_VAL/2); 
      binTreshFilter->SetOutsideValue (0);
      binTreshFilter->SetInsideValue (LABEL_VAL);
      binTreshFilter->Update();

      procImage = binTreshFilter->GetOutput();
    } else {

      if (debug) {
	cout << " resampling (NN interpol) the data to voxelsize ";
	for (int i = 0; i < Dimension; i++) cout << spacing[i] << ", ";
	cout << endl;
      }
      ResamplerType::Pointer resampleFilter = ResamplerType::New();
      TransformType::Pointer transform = TransformType::New();
      NNInterpolFunctionType::Pointer interpolFunction = NNInterpolFunctionType::New();

      transform->SetIdentity();
      interpolFunction->SetInputImage(procImage);
    
      ResamplerType::SizeType size;
      imageRegion = procImage->GetBufferedRegion();
      dim[0] = imageRegion.GetSize(0);
      dim[1] = imageRegion.GetSize(1);
      dim[2] = imageRegion.GetSize(2);
      const itk::Vector<double,3> inSpacing = image->GetSpacing();
      for (int i = 0; i < Dimension; i++) {
	size[i] = (unsigned int) ((double) inSpacing[i] * dim[i] / spacing[i] );
      }
      resampleFilter->SetInput(procImage); 
      resampleFilter->SetSize(size);
      resampleFilter->SetOutputSpacing(spacing);
      resampleFilter->SetOutputOrigin(image->GetOrigin());
      resampleFilter->SetTransform(transform.GetPointer());
      resampleFilter->SetInterpolator(interpolFunction.GetPointer());
      resampleFilter->Update();

      procImage = resampleFilter->GetOutput();
    }
    imageRegion = procImage->GetBufferedRegion();
    dim[0] = imageRegion.GetSize(0);
    dim[1] = imageRegion.GetSize(1);
    dim[2] = imageRegion.GetSize(2);
    if (debug) cout << "size of the resized image " << dim[0] << "," << dim[1] << "," << dim[2]<< endl; 
  }

  // Closing 1
  if (!noClosing)
  {
    if (asymCloseOn) {
      if (debug) cout << "Closing 1" << endl;
      StructuringElementType structuringElement;
      AsymStructuringElementType asymStructuringElement;
      
      dilateFilterType::Pointer dilateFilter = dilateFilterType::New();  
      erodeAsymFilterType::Pointer erodeFilter = erodeAsymFilterType::New();  
      
      dilateFilter->SetInput(procImage);
      erodeFilter->SetInput(dilateFilter->GetOutput());
      dilateFilter->SetDilateValue (LABEL_VAL);
      erodeFilter->SetErodeValue (LABEL_VAL);
      
      structuringElement.SetRadius( elemSize );  // 3x3x3 structuring element
      asymStructuringElement.SetRadius( elemSize );
      structuringElement.CreateStructuringElement( );
      asymStructuringElement.CreateStructuringElement( );
      
      dilateFilter->SetKernel( structuringElement );
      erodeFilter->SetKernel( asymStructuringElement );
    
      erodeFilter->Update();
      
      procImage = erodeFilter->GetOutput();
    } else {
      if (debug) cout << "Closing 1" << endl;
      StructuringElementType structuringElement;
      
      dilateFilterType::Pointer dilateFilter = dilateFilterType::New();  
      erodeFilterType::Pointer erodeFilter = erodeFilterType::New();  
      
      dilateFilter->SetInput(procImage);
      erodeFilter->SetInput(dilateFilter->GetOutput());
      dilateFilter->SetDilateValue (LABEL_VAL);
      erodeFilter->SetErodeValue (LABEL_VAL);
      
      structuringElement.SetRadius( elemSize );  // 3x3x3 structuring element
      structuringElement.CreateStructuringElement( );
      
      dilateFilter->SetKernel( structuringElement );
      erodeFilter->SetKernel( structuringElement );
    
      erodeFilter->Update();
      
      procImage = erodeFilter->GetOutput();
    }
  }

  // CC 1 & threshold
  if (!noCCL) {
    if (debug) cout << "Connected Component Labeling 1"  << endl;
    CCLFilterType::Pointer CCLFilter = CCLFilterType::New();
    CCLRelabelFilterType::Pointer RelabelFilter = CCLRelabelFilterType::New();
    binThreshFilterType::Pointer binTreshFilter = binThreshFilterType::New();
    CCLFilter->SetInput(procImage);
    RelabelFilter->SetInput(CCLFilter->GetOutput());
    binTreshFilter->SetInput(RelabelFilter->GetOutput());
    binTreshFilter->SetUpperThreshold(1);
    binTreshFilter->SetLowerThreshold(1);// only largest object (labels are sorted largest to smallest)
    binTreshFilter->SetOutsideValue (0);
    binTreshFilter->SetInsideValue (LABEL_VAL);
    binTreshFilter->Update();
    
    procImage = binTreshFilter->GetOutput();
  }

  // Conversion to float & gauss filter & Conversion to unsigned short & Tresholding at LABEL_VAL/2

  if (gaussianOn) {

    if (debug) cout << "conversions, smoothing: variance " << variance[0] << "," << variance[1] << "," << variance[2] << endl;
    {
      castInputFilterType::Pointer convInFilter = castInputFilterType::New();
      convInFilter->SetInput(procImage);
      
      gaussFilterType::Pointer smoothFilter = gaussFilterType::New();  
      
      smoothFilter->SetInput(convInFilter->GetOutput());
      smoothFilter->SetVariance(variance);
      
      castOutputFilterType::Pointer convOutFilter = castOutputFilterType::New();
      
      convOutFilter->SetInput(smoothFilter->GetOutput());
      convOutFilter->Update();
      
      
      if (debug) cout << "Tresholding at 50% " << endl;
      binThreshFilterType::Pointer binTreshFilter = binThreshFilterType::New();
      binTreshFilter->SetInput(convOutFilter->GetOutput());
      binTreshFilter->SetUpperThreshold(LABEL_VAL);
      binTreshFilter->SetLowerThreshold(LABEL_VAL/2); 
      binTreshFilter->SetOutsideValue (0);
      binTreshFilter->SetInsideValue (LABEL_VAL);
      binTreshFilter->Update();
      
      if (debug) cout << "Closing 2 " << endl;
      if (!noClosing)
	{
	  if (asymCloseOn) {
	    StructuringElementType structuringElement;
	    AsymStructuringElementType asymStructuringElement;
	    dilateFilterType::Pointer dilateFilter = dilateFilterType::New();  
	    erodeAsymFilterType::Pointer erodeFilter = erodeAsymFilterType::New();  
	    
	    dilateFilter->SetInput(binTreshFilter->GetOutput());
	    erodeFilter->SetInput(dilateFilter->GetOutput());
	    dilateFilter->SetDilateValue (LABEL_VAL);
	    erodeFilter->SetErodeValue (LABEL_VAL);
	    structuringElement.SetRadius( elemSize );  // 3x3x3 structuring element
	    asymStructuringElement.SetRadius( elemSize );
	    asymStructuringElement.CreateStructuringElement( );
	    structuringElement.CreateStructuringElement( );
	    dilateFilter->SetKernel( structuringElement );
	    erodeFilter->SetKernel( asymStructuringElement );
	    
	    erodeFilter->Update();
	    
	    procImage = erodeFilter->GetOutput();
	  } else {
	    StructuringElementType structuringElement;
	    dilateFilterType::Pointer dilateFilter = dilateFilterType::New();  
	    erodeFilterType::Pointer erodeFilter = erodeFilterType::New();  
	    
	    dilateFilter->SetInput(binTreshFilter->GetOutput());
	    erodeFilter->SetInput(dilateFilter->GetOutput());
	    dilateFilter->SetDilateValue (LABEL_VAL);
	    erodeFilter->SetErodeValue (LABEL_VAL);
	    structuringElement.SetRadius( elemSize );  // 3x3x3 structuring element
	    structuringElement.CreateStructuringElement( );
	    dilateFilter->SetKernel( structuringElement );
	    erodeFilter->SetKernel( structuringElement );
	    
	    erodeFilter->Update();
	    
	    procImage = erodeFilter->GetOutput();
	  }
	}
    }

    // CC 2 & treshold
    if (!noCCL) {
      if (debug) cout << "Connected Component Labeling 1"  << endl;
      CCLFilterType::Pointer CCLFilter = CCLFilterType::New();
      CCLRelabelFilterType::Pointer RelabelFilter = CCLRelabelFilterType::New();
      binThreshFilterType::Pointer binTreshFilter = binThreshFilterType::New();
      CCLFilter->SetInput(procImage);
      RelabelFilter->SetInput(CCLFilter->GetOutput());
      binTreshFilter->SetInput(RelabelFilter->GetOutput());
      binTreshFilter->SetUpperThreshold(1);
      binTreshFilter->SetLowerThreshold(1);// only largest object (labels are sorted largest to smallest)
      binTreshFilter->SetOutsideValue (0);
      binTreshFilter->SetInsideValue (LABEL_VAL);
      binTreshFilter->Update();
    
      procImage = binTreshFilter->GetOutput();
    } 
  }

  if (!noLS) { // do level set smoothing!
    if (debug) cout << "Level Set Smoothing"  << endl;
    CastToRealFilterType::Pointer toReal = CastToRealFilterType::New();
    toReal->SetInput(procImage );
    
    AntiAliasFilterType::Pointer antiAliasFilter = AntiAliasFilterType::New();
    antiAliasFilter->SetInput( toReal->GetOutput() );
    antiAliasFilter->SetMaximumRMSError( maximumRMSError );
    antiAliasFilter->SetNumberOfIterations( numberOfIterations );
    antiAliasFilter->SetNumberOfLayers( 2 );
    antiAliasFilter->Update();
  
    RescaleFilter::Pointer rescale = RescaleFilter::New();
    rescale->SetInput( antiAliasFilter->GetOutput() );
    rescale->SetOutputMinimum(   0 );
    rescale->SetOutputMaximum( 255 );

    binThreshFilterType::Pointer binTreshFilter = binThreshFilterType::New();
    binTreshFilter->SetInput(rescale->GetOutput());
    binTreshFilter->SetUpperThreshold(255);
    binTreshFilter->SetLowerThreshold(127); 
    binTreshFilter->SetOutsideValue (0);
    binTreshFilter->SetInsideValue (1);
    binTreshFilter->Update();
    procImage = binTreshFilter->GetOutput();
  }

  // connect everything by 6 connectedness, so that no diagonals in data
  if (debug) cout << "enforce strict 6 connectedness" << endl; 
  data = &((*procImage)[nullIndex]);
  imageRegion = procImage->GetBufferedRegion();
  dim[0] = imageRegion.GetSize(0);
  dim[1] = imageRegion.GetSize(1);
  dim[2] = imageRegion.GetSize(2);
  clear_edge(data, dim, 0);
  NoDiagConnect(data,dim);

  if(deleteVessels){
	  cout << "1" << endl;
  /*  if(!MaskfileName.empty()){
      savemaskOn=1;
    }
    if (debug) cout << "Loading labelfile " << endl;
    try
    {
      labelReader->SetFileName(fileName.c_str()) ;
      labelReader->Update() ;
      labelImage = labelReader->GetOutput() ;
    }
    catch (ExceptionObject e)
    {
      e.Print(cout) ;
      exit(0) ;
    }
    
    if (debug) cout << "Loading corrected T1 image" << endl;
    try
    {
      VolumeReaderType::Pointer correctedImageReader = VolumeReaderType::New();
      correctedImageReader->SetFileName(delVesselCorrectedImName.c_str());
      correctedImageReader->Update() ;
      correctedImage = correctedImageReader->GetOutput() ;
    }
    catch (ExceptionObject e)
    {
      e.Print(cout) ;
      exit(0) ;
  }*/
    
    ImageSampleType::Pointer imageSample;
    ImageSampleType::Pointer labelSample;
    MembershipSampleType::Pointer membershipSample;
    CovarianceAlgorithmType::Pointer covarianceAlgorithm;
    MembershipSampleType::ClassSampleType::ConstPointer classSample;
    vector<int> labelList;
    
    
    cout << "2" << endl;
    
    
    imageSample = ImageSampleType::New();
    labelSample = ImageSampleType::New();
    imageSample->SetImage(EMSImage);
    labelSample->SetImage(labelImage);

    membershipSample = MembershipSampleType::New();
    membershipSample->SetSample(imageSample);

    covarianceAlgorithm = CovarianceAlgorithmType::New();
    
    //Fill labelList vector with the labels IDs
    PixelType label ;
    ImageSampleIterator iterLabel = labelSample->Begin() ;
    while( iterLabel != labelSample->End() )
    {
      label = iterLabel.GetMeasurementVector()[0];
      if (label != 0 && !searchList(label,labelList)) labelList.push_back( label );
      ++iterLabel;
    }

    sort(labelList.begin(), labelList.end());
    
    
    cout << "3" << endl;
    
    
    int nbLabels = labelList.size();
    membershipSample->SetNumberOfClasses(nbLabels);

    ImageSampleIterator iter = imageSample->Begin() ;
    iterLabel = labelSample->Begin() ;  //set iterLabel Iterator to the beginning of labelSample
    while( iter != imageSample->End() )
    {
      label = iterLabel.GetMeasurementVector()[0];
      if (label != 0 ) {
        membershipSample->AddInstance(label,iter.GetInstanceIdentifier());
      }
      ++iter ;
      ++iterLabel;
    }
	
    int l;
    l=labelList[0]; //the WM is the first in the list
    classSample = membershipSample->GetClassSample(l);
    covarianceAlgorithm->SetInputSample( classSample );
    covarianceAlgorithm->Update();
	
    double meanInt, standartDev, thr;
    meanInt = (*covarianceAlgorithm->GetMean())[0];
    standartDev = sqrt((*covarianceAlgorithm->GetOutput())(0,0));
    thr = meanInt + standartDev;
    
    
    cout << "4" << endl;
    
    
    ImageType::Pointer labelImageSmooth ;
    
    castInputFilterType::Pointer convInFilter = castInputFilterType::New();
    convInFilter->SetInput(EMSImage);
    castOutputFilterType::Pointer convOutFilter = castOutputFilterType::New();
    double smoothSize = 2.0;
    int numIter = 1;
    cout << "smoothing (Gaussian): size " << smoothSize << ", iterations " << numIter << endl;
    gaussFilterType::Pointer smoothFilter = gaussFilterType::New();  
    for ( int i=0; i < numIter; i++) {
      smoothFilter->SetInput(convInFilter->GetOutput());
      smoothFilter->SetVariance(smoothSize);
      convOutFilter->SetInput(smoothFilter->GetOutput());
      try {
        convOutFilter->Update();
      }
      catch (ExceptionObject & err) {
        cerr << "ExceptionObject caught!" << endl;
        cerr << err << endl;
        return EXIT_FAILURE;	
      }
    }   
     
    
    cout << "5" << endl;
    
    
    labelImageSmooth = convOutFilter->GetOutput();
    
    IteratorType it_labelIm(labelImage, labelImage->GetRequestedRegion());
    IteratorType it_procIm(procImage, procImage->GetRequestedRegion());
    IteratorType it_EMSIm(EMSImage, EMSImage->GetRequestedRegion());
    IteratorType it_labelSmoothIm(labelImageSmooth, labelImageSmooth->GetRequestedRegion());
    
    ImageType::IndexType index_labelIm, index_procIm, index_EMSIm, index_labelSmoothIm;
    
    it_labelIm.GoToBegin();
    it_procIm.GoToBegin();
    it_EMSIm.GoToBegin();
    it_labelSmoothIm.GoToBegin();
    
    while((!it_labelIm.IsAtEnd()) && (!it_procIm.IsAtEnd()) && (!it_EMSIm.IsAtEnd()))
    {
      index_labelIm = it_labelIm.GetIndex();
      index_procIm = it_procIm.GetIndex();
      index_EMSIm = it_EMSIm.GetIndex();
      index_labelSmoothIm = it_labelSmoothIm.GetIndex();
      
      int subtractionIm = ((procImage->GetPixel(index_procIm)) - (labelImage->GetPixel(index_labelIm)));

      if((subtractionIm == 1) && (EMSImage->GetPixel(index_EMSIm) > thr))
      {
	      it_EMSIm.Set(it_labelSmoothIm.Get());
	      it_procIm.Set(2);
      }
      ++it_labelIm;
      ++it_procIm;
      ++it_EMSIm;
      ++it_labelSmoothIm;
    }
  }
  
  
  cout << "6" << endl;
  
  
  //write the volume
  try
    {
      if (debug && (!skullstrippingOn || savemaskOn)) cout << "writing processed data" << endl;
      VolumeWriterType::Pointer writer = VolumeWriterType::New();
	writer->UseCompressionOn();
      if (savemaskOn){
	writer->SetFileName(maskfileName.c_str());
	writer->SetInput(procImage);
	writer->Write();
      }
      else{
	if(!skullstrippingOn){
	  writer->SetFileName(outfileName.c_str());
	  writer->SetInput(procImage);
	  writer->Write();
	}
      }
      /*if(deleteVesselsOn){
        writer->SetFileName(outfileName.c_str());
        writer->SetInput(correctedImage);
        writer->Write();
      }*/
    }
    catch (ExceptionObject e)
    {
      e.Print(cout) ;
      exit(0) ;
    }
  
    
  //second part of skullstripping : using the binary mask
  if(skullstrippingOn){
   // ImageType::Pointer EMSImage ;
   // ImageType::Pointer strippedImage ;
   /* VolumeReaderType::Pointer EMSReader = VolumeReaderType::New();
    if (debug) cout << "Loading EMSfile " << endl;
    try{
      EMSReader->SetFileName(notstrippedfileName.c_str()) ;
      EMSReader->Update() ;
      EMSImage = EMSReader->GetOutput() ;
    }
    catch (ExceptionObject e){
      e.Print(cout) ;
      exit(0) ;
    }
	  if (debug) cout << "EMSfile loaded " << endl;*/
    if (debug) cout << "masking images  " << endl;

    maskFilterType::Pointer maskFilter = maskFilterType::New();
    maskFilter->SetInput1(EMSImage);
    maskFilter->SetInput2(procImage);
    try {
      maskFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    } 
    strippedImage = maskFilter->GetOutput();

    VolumeWriterType::Pointer writer2 = VolumeWriterType::New();
    writer2->SetFileName(outfileName.c_str());
    writer2->UseCompressionOn();
    writer2->SetInput(strippedImage);
    writer2->Write();

  }
  
  cout<<"end of prog"<<endl;  
  return 0 ;
}


static int NoDiagConnect (unsigned short *image, int *dim) 
  // does not allow connection via diagonals only, enforces strict 6 connectedness
  // image has to be of type unsigned short
{

  //z axis
  int dimx = dim[0];
  int dimy = dim[1];
  int dimz = dim[2];
  bool correctionNeeded = true;
  int cnt = 0;

  while (correctionNeeded) {
    cnt++;
    if (debug) cout << "NoDiag scan " << cnt << endl; 
    correctionNeeded = false;
    int dy = dimx*dimy;
    int dx = dimx;

   for (int i = 1; i < dimx - 1; i++) {
    for (int j = 1; j < dimy - 1; j++) {
     for (int k = 1; k < dimz - 1; k++) {
       unsigned short val = image[i + j * dimx + k * dy];
       if (val != 0) {
         // x,y 
         if ((image[i-1+j*dx+k*dy] == 0) && (image[i+(j-1)*dx+k*dy] == 0) && (image[i-1+(j-1)*dx+k*dy] != 0)) {
           correctionNeeded = true;
           image[i-1+j*dx+k*dy] = val;
         }
         if ((image[i+1+j*dx+k*dy] == 0) && (image[i+(j+1)*dx+k*dy] == 0) && (image[i+1+(j+1)*dx+k*dy] != 0)) {
           correctionNeeded = true;
           image[i+1+j*dx+k*dy] = val;
         }
         if ((image[i+1+j*dx+k*dy] == 0) && (image[i+(j-1)*dx+k*dy] == 0) && (image[i+1+(j-1)*dx+k*dy] != 0)) {
           correctionNeeded = true;
           image[i+1+j*dx+k*dy] = val;
         }
         if ((image[i-1+j*dx+k*dy] == 0) && (image[i+(j+1)*dx+k*dy] == 0) && (image[i-1+(j+1)*dx+k*dy] != 0)) {
           correctionNeeded = true;
           image[i-1+j*dx+k*dy] = val;
         }
         
         // xz
         if ((image[i-1+j*dx+k*dy] == 0) && (image[i+j*dx+(k-1)*dy] == 0) && (image[i-1+j*dx+(k-1)*dy] != 0)) {
           correctionNeeded = true;
           image[i-1+j*dx+k*dy] = val;
         }
         if ((image[i+1+j*dx+k*dy] == 0) && (image[i+j*dx+(k+1)*dy] == 0) && (image[i+1+j*dx+(k+1)*dy] != 0)) {
           correctionNeeded = true;
           image[i+1+j*dx+k*dy] = val;
         }
         if ((image[i+1+j*dx+k*dy] == 0) && (image[i+j*dx+(k-1)*dy] == 0) && (image[i+1+j*dx+(k-1)*dy] != 0)) {
           correctionNeeded = true;
           image[i+1+j*dx+k*dy] = val;
         }
         if ((image[i-1+j*dx+k*dy] == 0) && (image[i+j*dx+(k+1)*dy] == 0) && (image[i-1+j*dx+(k+1)*dy] != 0)) {
           correctionNeeded = true;
           image[i-1+j*dx+k*dy] = val;
         }
         
         // yz
         if ((image[i+(j-1)*dx+k*dy] == 0) && (image[i+j*dx+(k-1)*dy] == 0) && (image[i+(j-1)*dx+(k-1)*dy] != 0)) {
           correctionNeeded = true;
           image[i+(j-1)*dx+k*dy] = val;
         }
         if ((image[i+(j+1)*dx+k*dy] == 0) && (image[i+j*dx+(k+1)*dy] == 0) && (image[i+(j+1)*dx+(k+1)*dy] != 0)) {
           correctionNeeded = true;
           image[i+(j+1)*dx+k*dy] = val;
         }
         if ((image[i+(j+1)*dx+k*dy] == 0) && (image[i+j*dx+(k-1)*dy] == 0) && (image[i+(j+1)*dx+(k-1)*dy] != 0)) {
           correctionNeeded = true;
           image[i+(j+1)*dx+k*dy] = val;
         }
         if ((image[i+(j-1)*dx+k*dy] == 0) && (image[i+j*dx+(k+1)*dy] == 0) && (image[i+(j-1)*dx+(k+1)*dy] != 0)) {
           correctionNeeded = true;
           image[i+(j-1)*dx+k*dy] = val;
         }
       }
     }
      }
    }
  }

  return 1;
}

static void clear_edge(unsigned short *image, int *dims, int clear_label)
  // clears the edge of the image
{
  int size_plane = dims[0]*dims[1];
  int size_line = dims[0];

  for (int z = 0; z < dims[2]; z++) {
    for (int y = 0; y < dims[1]; y++) {
      if ( (y == 0) || (y == dims[1]-1) ||
        (z == 0) || (z == dims[2]-1) ) { // draw whole plane
        for (int x = 0; x < dims[0] ; x++) 
          image[x +  size_line * y + size_plane * z] = clear_label;
        } else { // draw edges of x
        image[0 +  size_line * y + size_plane * z] = clear_label;
        image[size_line - 1 +  size_line * y + size_plane * z] = clear_label;
      }
    }
  }

}


static void 
draw_fill_inside_image(unsigned short *image, int *dims, int new_label)
  // Fill the 'inside' part of an image (closed objects that do not touch the
  // image edges)
{
  int size_plane = dims[0]*dims[1];
  int size_line = dims[0];

  unsigned short label;
  if (new_label > 1) label = new_label-1;
  else label = new_label+1;
  
  unsigned short background = 0;

  // fill image edges -> won't work if object touches the edge !!
  for (int z = 0; z < dims[2]; z++) {
    for (int y = 0; y < dims[1]; y++) {
      if ( (y == 0) || (y == dims[1]-1) ||
        (z == 0) || (z == dims[2]-1) ) { // draw whole plane
     for (int x = 0; x < dims[0] ; x++) 
       image[x +  size_line * y + size_plane * z] = label;
      } else { // draw edges of x
     image[0 +  size_line * y + size_plane * z] = label;
     image[size_line - 1 +  size_line * y + size_plane * z] = label;
      }
    }
  }
  
  // forward propagation
  for (int z = 1; z < dims[2]-1; z++) {
    for (int y = 1; y < dims[1]-1; y++) {
      for (int x = 1; x < dims[0]-1; x++) {
     int index = x +  size_line * y + size_plane * z;
     if (image[index] == background &&    // check past neighborhood
         (image[index - 1] == label || 
          image[index + size_line] == label || 
          image[index - size_plane] == label 
          )) {
       image[index] = label;
     }
      }
    }
  }

  // backward propagation
  for (int z = dims[2]-2; z > 0; z--) {
    for (int y = dims[1]-2; y > 0; y--) {
      for (int x = dims[0]-2; x > 0; x--) {
     int index = x +  size_line * y + size_plane * z;
     if (image[index] == background &&    // check past neighborhood
         (image[index + 1] == label || 
          image[index + size_line] == label || 
          image[index + size_plane] == label 
          )) {
       image[index] = label;
     }
      }
    }
  }

  // reassign labels
  for (int i = 0; i < dims[2]*dims[1]*dims[0]; i++) {
    if (image[i] == label) image[i] = background;
    else if (image[i] == background) image[i] = new_label;
  }
}

bool searchList(int label, vector<int> list)
{
	int size = list.size();
	bool found = 0;
	for (int i=0; (i<size) && (found == 0) ; i++) if (list[i]==label) found = 1;
	return found;
}
