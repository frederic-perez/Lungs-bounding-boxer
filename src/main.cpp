// --

#include <iostream>
#include <limits>
#include <string>

#include "aux-raw-compiler-warnings-off++begin.h"

	#include <itkBinaryImageToShapeLabelMapFilter.h>
	#include <itkBinaryThresholdImageFilter.h>
	#include <itkConnectedComponentImageFilter.h>
	#include <itkConnectedThresholdImageFilter.h>
	#include <itkImage.h>
	#include <itkImageFileReader.h>
	#include <itkImageFileWriter.h>
	#include <itkLabelMapToLabelImageFilter.h>
	#include <itkShapeKeepNObjectsLabelMapFilter.h>
	#include <itkShapeOpeningLabelMapFilter.h>

#include "aux-raw-compiler-warnings-off++end.h"

using ImageCT = itk::Image<double, 3>;
using ImageBinary = itk::Image<unsigned char,3>;
using ImageLabels = itk::Image<size_t, 3>;

bool
FindLungs(const ImageCT& a_imgCT, ImageBinary::Pointer& a_imgLabels)
{
	// Thresholding the input image
	//
	std::clog << __func__ << ": Thresholding the input image..." << std::endl;

	using BinaryThresholdType =
		itk::BinaryThresholdImageFilter<ImageCT, ImageBinary>;
	BinaryThresholdType::Pointer thresholdFilter = BinaryThresholdType::New();

	thresholdFilter->SetInput(&a_imgCT);
	const ImageCT::ValueType lowerThreshold = -1000.;
	const ImageCT::ValueType upperThreshold = -400.;
	thresholdFilter->SetUpperThreshold(upperThreshold);
	thresholdFilter->SetLowerThreshold(lowerThreshold);
	const ImageLabels::ValueType backgroundValue = 0;
	const ImageLabels::ValueType foregroundValue = 1;
	thresholdFilter->SetOutsideValue(backgroundValue);
	thresholdFilter->SetInsideValue(foregroundValue);

	try {
		thresholdFilter->Update();
	}	catch (itk::ExceptionObject& e) {
		std::cerr << __func__
			<< ": thresholdFilter->Update: ITK Exception caught: " << e << std::endl;
		return false;
	}


	// Remove small blobs
	//
	std::clog << __func__ << ": Removing small blobs..." << std::endl;

	using LabelObj3DType = itk::ShapeLabelObject<size_t, 3>;
	using LabelMap3DType = itk::LabelMap<LabelObj3DType>;

	using ShapeLabelFilter =
		itk::BinaryImageToShapeLabelMapFilter<ImageBinary, LabelMap3DType>;
	ShapeLabelFilter::Pointer featureExtractor = ShapeLabelFilter::New();
	featureExtractor->SetInput(thresholdFilter->GetOutput());
	featureExtractor->SetInputForegroundValue(foregroundValue);
	featureExtractor->SetOutputBackgroundValue(backgroundValue);
	try {
		featureExtractor->Update();
	}	catch (itk::ExceptionObject& e) {
		std::cerr << __func__
			<< ": featureExtractor->Update: ITK Exception caught: " << e << std::endl;
		return false;
	}

	LabelMap3DType::Pointer labelMap = featureExtractor->GetOutput();
	size_t numLabels = labelMap->GetNumberOfLabelObjects();

	std::clog << __func__ << ": #labels after featureExtractor: " << numLabels
		<< std::endl;

	// Remove labels with volume smaller than lambda ml
	//
	const double lambda = 200000.;
	std::clog << __func__
		<< ": Removing labels with volume smaller than lambda (" << lambda
		<< ") mm^3..." << std::endl;

	using LabelSelector = itk::ShapeOpeningLabelMapFilter<LabelMap3DType>;
	LabelSelector::Pointer filterLowVolume = LabelSelector::New();
	filterLowVolume->SetInput(labelMap);
	filterLowVolume->SetAttribute("PhysicalSize");
	filterLowVolume->SetLambda(lambda);
	filterLowVolume->SetReverseOrdering(false);
	try {
		filterLowVolume->Update();
	} catch (itk::ExceptionObject& e) {
		std::cerr << __func__
			<< ": filterLowVolume->Update: ITK Exception caught: " << e << std::endl;
		return false;
	}

	numLabels = filterLowVolume->GetOutput()->GetNumberOfLabelObjects();
	std::clog << __func__ << ": #labels after small blobs removal: "
		<< numLabels << std::endl;

#undef _LABEL_DEBUGGING_
#ifdef _LABEL_DEBUGGING_
	using ImageWriterLabels = itk::ImageFileWriter<ImageLabels>;
	if (numLabels < std::numeric_limits<ImageWriterLabels::InputImagePixelType>::max()) {
		std::clog << __func__
			<< ": Writting Intermediate Label Image" << std::endl;
		using LabelMap2LabelImage = itk::LabelMapToLabelImageFilter<LabelMap3DType, ImageLabels>;
		LabelMap2LabelImage::Pointer labelMapToLabelImageFilter = LabelMap2LabelImage::New();
		labelMapToLabelImageFilter->SetInput(filterLowVolume->GetOutput());
		//labelMapToLabelImageFilter->Update();
		ImageWriterLabels::Pointer writer = ImageWriterLabels::New();
		writer->SetInput(labelMapToLabelImageFilter->GetOutput());
		writer->SetFileName("LabelsAfterLargeBlogsRemoved.mha");
		try {
			writer->Update();
		} catch(itk::ExceptionObject& e) {
			std::cerr << __func__
				<< ": writer->Update: ITK Exception caught: " << e << std::endl;
			return false;
		}
	}
#endif

	// Keep just one object
	//
	std::clog << __func__ << ": Keeping just one object..." << std::endl;

	using KeepNObjectsFilter =
		itk::ShapeKeepNObjectsLabelMapFilter<LabelMap3DType>;
	KeepNObjectsFilter::Pointer keepNObjectsFilter = KeepNObjectsFilter::New();
	keepNObjectsFilter->SetInput(filterLowVolume->GetOutput());
	keepNObjectsFilter->SetReverseOrdering(true);
	const itk::SizeValueType maxNumObjects = 1;
	keepNObjectsFilter->SetNumberOfObjects(maxNumObjects);
	keepNObjectsFilter->SetAttribute(
		KeepNObjectsFilter::LabelObjectType::NUMBER_OF_PIXELS_ON_BORDER);
	try {
		keepNObjectsFilter->Update();
	}	catch (itk::ExceptionObject& e) {
		std::cerr << __func__
			<< ": keepNObjectsFilter->Update: ITK Exception caught: " << e
			<< std::endl;
		return false;
	}

	std::clog << __func__ << ": Number of labels found: "
		<< keepNObjectsFilter->GetNumberOfObjects() << std::endl;

	using LabelMapToLabelImageFilterType =
		itk::LabelMapToLabelImageFilter<LabelMap3DType, ImageBinary>;
	LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter
		= LabelMapToLabelImageFilterType::New();
	labelMapToLabelImageFilter->SetInput(keepNObjectsFilter->GetOutput());
	try {
		labelMapToLabelImageFilter->Update();
	}	catch (itk::ExceptionObject& e) {
		std::cerr << __func__
			<< ": labelMapToLabelImageFilter->Update: ITK Exception caught: " << e
			<< std::endl;
		return false;
	}

	// Connected component filtering
	//
	std::clog << __func__ << ": Connected component filtering..." << std::endl;

	using ConnectedComponentImageFilterType =
		itk::ConnectedComponentImageFilter <ImageBinary, ImageBinary>;

	ConnectedComponentImageFilterType::Pointer connected =
		ConnectedComponentImageFilterType::New();
	connected->SetInput(labelMapToLabelImageFilter->GetOutput());
	try {
		connected->Update();
	}	catch (itk::ExceptionObject& e) {
		std::cerr << __func__
			<< ": connected->Update: ITK Exception caught: " << e << std::endl;
		return false;
	}

	a_imgLabels = connected->GetOutput();
	a_imgLabels->Register();

	return true;
}

bool
RemoveExternalStuffByFloodFilling(
	ImageCT::Pointer& a_imgCT,
	const std::string& a_filenameOutputLabelsSpyHalfWay)
{
	// Getting the CT vertices (indices)
	//
	std::clog << __func__ << ": Computing the 8 vertices (indices) of the CT, "
		<< "seeds for region growing filter..." << std::endl;

	using IteratorType = itk::ImageRegionConstIteratorWithIndex<ImageCT>;
	IteratorType it(a_imgCT, a_imgCT->GetRequestedRegion());

	ImageCT::IndexValueType iMin = 90000000;
	ImageCT::IndexValueType jMin = iMin, kMin = iMin;
	ImageCT::IndexValueType iMax = -2000000;
	ImageCT::IndexValueType jMax = iMax, kMax = iMax;
	it.GoToBegin();
	for (; !it.IsAtEnd(); ++it) {
		const ImageCT::IndexType idx = it.GetIndex();
		iMin = std::min<ImageCT::IndexValueType>(iMin, idx[0]);
		jMin = std::min<ImageCT::IndexValueType>(jMin, idx[1]);
		kMin = std::min<ImageCT::IndexValueType>(kMin, idx[2]);
		iMax = std::max<ImageCT::IndexValueType>(iMax, idx[0]);
		jMax = std::max<ImageCT::IndexValueType>(jMax, idx[1]);
		kMax = std::max<ImageCT::IndexValueType>(kMax, idx[2]);
	}

	std::clog << __func__ << ": Bounding box: ["
		<< iMin << ", " << iMax << "] x ["
		<< jMin << ", " << jMax << "] x ["
		<< kMin << ", " << kMax << "]" << std::endl;

	std::list<ImageCT::IndexType> regionGrowerSeeds;
	const ImageCT::IndexType corner1 = { { iMin, jMin, kMax } };
	const ImageCT::IndexType corner2 = { { iMin, jMax, kMax } };
	const ImageCT::IndexType corner3 = { { iMax, jMax, kMax } };
	const ImageCT::IndexType corner4 = { { iMax, jMin, kMax } };
	const ImageCT::IndexType corner5 = { { iMin, jMin, kMin } };
	const ImageCT::IndexType corner6 = { { iMin, jMax, kMin } };
	const ImageCT::IndexType corner7 = { { iMax, jMax, kMin } };
	const ImageCT::IndexType corner8 = { { iMax, jMin, kMin } };

	regionGrowerSeeds.push_back(corner1);
	regionGrowerSeeds.push_back(corner2);
	regionGrowerSeeds.push_back(corner3);
	regionGrowerSeeds.push_back(corner4);
	regionGrowerSeeds.push_back(corner5);
	regionGrowerSeeds.push_back(corner6);
	regionGrowerSeeds.push_back(corner7);
	regionGrowerSeeds.push_back(corner8);

	// Region growing
	//
	std::clog << __func__
		<< ": Setting and executing the region growing filter..." << std::endl;

	using RegionGrower =
		itk::ConnectedThresholdImageFilter<ImageCT, ImageCT>;
	RegionGrower::Pointer regionGrower = RegionGrower::New();
	regionGrower->ReleaseDataFlagOn();
	regionGrower->SetInput(a_imgCT);
	const ImageCT::PixelType lowerThreshold = -1024;
	const ImageCT::PixelType upperThreshold = -200;
	regionGrower->SetLower(lowerThreshold);
	regionGrower->SetUpper(upperThreshold);
	const ImageCT::PixelType replaceValue = 666;
	regionGrower->SetReplaceValue(replaceValue);
	for (const auto& seed : regionGrowerSeeds)
		regionGrower->AddSeed(seed);
	try {
		regionGrower->Update();
	} catch (itk::ExceptionObject& e) {
		std::cerr << __func__
			<< ": regionGrower->Update: ITK Exception caught: " << e << std::endl;
		return false;
	}

	// Mask filter
	//
	std::clog << __func__
		<< ": Setting and executing a mask filter..." << std::endl;

	using MaskFilterType = itk::MaskImageFilter<ImageCT, ImageCT>;
	MaskFilterType::Pointer maskFilter = MaskFilterType::New();
	maskFilter->SetInput(a_imgCT);
	maskFilter->SetMaskImage(regionGrower->GetOutput());
	maskFilter->SetMaskingValue(666);
	try {
		maskFilter->Update();
	} catch (itk::ExceptionObject& e) {
		std::cerr << __func__
			<< ": maskFilter->Update: ITK Exception caught: " << e << std::endl;
		return false;
	}

	// Update a_imgCT
	//
	std::clog << __func__
		<< ": Setting the input/output CT image as the result of the mask filter..."
		<< std::endl;

	a_imgCT = maskFilter->GetOutput();
	a_imgCT->Register();

	// Write output (spy) file if requested
	//
	if (!a_filenameOutputLabelsSpyHalfWay.empty()) {

		std::clog << __func__ << ": Writing the output (spy) file as requested..."
			<< std::endl;

		using WriterType = itk::ImageFileWriter<ImageCT>;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName(a_filenameOutputLabelsSpyHalfWay);
		writer->SetInput(a_imgCT);
		try {
			writer->Update();
		} catch (itk::ExceptionObject& e) {
			std::cerr << __func__
				<< ": writer->Update: ITK Exception caught: " << e << std::endl;
			return false;
		}
	}

	return true;
}

void
ComputeAndOutputBoundingBox(const ImageBinary& a_imgLabels)
{
	using IteratorType = itk::ImageRegionConstIteratorWithIndex<ImageBinary>;
	IteratorType it(&a_imgLabels, a_imgLabels.GetRequestedRegion());

	ImageBinary::IndexValueType iMin = 90000000;
	ImageBinary::IndexValueType jMin = iMin, kMin = iMin;
	ImageBinary::IndexValueType iMax = -2000000;
	ImageBinary::IndexValueType jMax = iMax, kMax = iMax;
	it.GoToBegin();

	for (; !it.IsAtEnd(); ++it) {
		const ImageLabels::PixelType value = it.Get();
		if (value != 0) {
			const ImageLabels::IndexType idx = it.GetIndex();
			iMin = std::min<ImageBinary::IndexValueType>(iMin, idx[0]);
			jMin = std::min<ImageBinary::IndexValueType>(jMin, idx[1]);
			kMin = std::min<ImageBinary::IndexValueType>(kMin, idx[2]);
			iMax = std::max<ImageBinary::IndexValueType>(iMax, idx[0]);
			jMax = std::max<ImageBinary::IndexValueType>(jMax, idx[1]);
			kMax = std::max<ImageBinary::IndexValueType>(kMax, idx[2]);
		}
	}

	std::clog << __func__ << ": Bounding box of non-zero labels: ["
		<< iMin << ", " << iMax << "] x ["
		<< jMin << ", " << jMax << "] x ["
		<< kMin << ", " << kMax << "]" << std::endl;
}

int
main(int argc, char* argv[])
{
	// Parsing command line
	//
	if (argc != 3 && argc != 4) {
		const std::string folderWindows = " E:\\7ickle\\Lungs-bounding-boxer\\";
		const std::string folderLinux = " /mnt/hgfs/E/7ickle/Lungs-bounding-boxer/";
		std::cerr << "Usage: " << argv[0]
			<< " <filenameInputCT> <filenameOutputLabels>"
			<< " [<filenameOutputLabelsSpyHalfWay>]\n"
			<< "Example call (Windows): " << argv[0]
				<< folderWindows << "inputCT.mha"
				<< folderWindows << "outputLabels-Windows.mha"
				<< folderWindows << "outputLabelsSpyHalfWay-Windows.mha\n"
			<< "Example call (Linux): " << argv[0]
				<< folderLinux << "inputCT.mha"
				<< folderLinux << "outputLabels-Linux.mha"
				<< folderLinux << "outputLabelsSpyHalfWay-Linux.mha"
			<< std::endl;
		return EXIT_FAILURE;
	}

	const std::string filenameInputCT = argv[1];
	const std::string filenameOutputLabels = argv[2];
	const std::string filenameOutputLabelsSpyHalfWay = argc == 3 ? "" : argv[3];

	std::clog << __func__ << ": Parsed parameters:\n"
		<< " <filenameInputCT> = " << filenameInputCT << '\n'
		<< " <filenameOutputLabels> = " << filenameOutputLabels << '\n'
		<< " <filenameOutputLabelsSpyHalfWay> = "
		<< (filenameOutputLabelsSpyHalfWay.empty() ?
			"[empty]" : filenameOutputLabelsSpyHalfWay)
		<< std::endl;

	// Read the input image
	//
	std::clog << __func__ << ": Reading the input image..." << std::endl;

	using ReaderType = itk::ImageFileReader<ImageCT>;

	ReaderType::Pointer inputImage = ReaderType::New();
	inputImage->SetFileName(filenameInputCT);
	try {
		inputImage->Update();
	} catch (itk::ExceptionObject& e) {
		std::cerr << __func__
			<< ": inputImage->Update: ITK Exception caught: " << e << std::endl;
		return EXIT_FAILURE;
	}

	ImageCT::Pointer imageCT = inputImage->GetOutput();

	// Create the image of labels
	//
	std::clog << __func__ << ": Creating the image of labels..." << std::endl;

	ImageBinary::Pointer imgLabels = ImageBinary::New();
	imgLabels->SetRegions(imageCT->GetLargestPossibleRegion());
	imgLabels->CopyInformation(imageCT);
	imgLabels->Allocate();
	imgLabels->FillBuffer(0);

	// Get rid of external stuff (modifies imageCT)
	//
	bool succeeded =
		RemoveExternalStuffByFloodFilling(imageCT, filenameOutputLabelsSpyHalfWay);
	if (!succeeded) {
		std::cerr << __func__
			<< ": Error: RemoveExternalByFloodFilling did not succeed. Aborting...\n";
		return EXIT_FAILURE;
	}

	// Find lungs (modifies imgLabels)
	//
	succeeded =	FindLungs(*imageCT, imgLabels);
	if (!succeeded) {
		std::cerr << __func__
			<< ": Error: FindLungs did not succeed. Aborting...\n";
		return EXIT_FAILURE;
	}

	// Final serious part of the application
	//
	ComputeAndOutputBoundingBox(*imgLabels);

	// Write the output image
	//
	std::clog << __func__ << ": Writing the output image..." << std::endl;

	using WriterType = itk::ImageFileWriter<ImageBinary>;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(filenameOutputLabels);
	writer->SetInput(imgLabels);
	try {
		writer->Update();
	} catch (itk::ExceptionObject& e) {
		std::cerr << __func__
			<< ": writer->Update(): ITK Exception caught: " << e << std::endl;
		return EXIT_FAILURE;
	}

	std::clog << __func__ << " finished. Exiting..." << std::endl;

	return EXIT_SUCCESS;
}

// -- eof
