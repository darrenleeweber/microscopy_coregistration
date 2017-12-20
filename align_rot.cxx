#include "itkImage.h"

#include "itkImageFileReader.h"

#include "itkRegionOfInterestImageFilter.h"

#include "itkRigid3DTransform.h"
#include "itkEuler3DTransform.h"

#include "itkImageRegistrationMethod.h"

#include "itkNormalizedCorrelationImageToImageMetric.h"

#include "itkLinearInterpolateImageFunction.h"

#include "itkAmoebaOptimizer.h"

#include "itkResampleImageFilter.h"

#include "itkImageFileWriter.h"

#include "itkCommand.h"


class CommandIterationUpdateAmoeba : public itk::Command 
    {
public:
    typedef  CommandIterationUpdateAmoeba   Self;
    typedef  itk::Command             Superclass;
    typedef itk::SmartPointer<Self>  Pointer;
    itkNewMacro( Self );
protected:
    CommandIterationUpdateAmoeba() 
        {
        m_IterationNumber=0;
        }
public:
    typedef itk::AmoebaOptimizer         OptimizerType;
    typedef   const OptimizerType   *    OptimizerPointer;

    void Execute(itk::Object *caller, const itk::EventObject & event)
        {
        Execute( (const itk::Object *)caller, event);
        }

    void Execute(const itk::Object * object, const itk::EventObject & event)
        {
        OptimizerPointer optimizer = 
            dynamic_cast< OptimizerPointer >( object );
        if( itk::FunctionEvaluationIterationEvent().CheckEvent( &event ) )
            {
            std::cout << m_IterationNumber++ << "   ";
            std::cout << optimizer->GetCachedValue() << "   ";
            std::cout << optimizer->GetCachedCurrentPosition() << std::endl;
            }
        }
private:
    unsigned long m_IterationNumber;
    };


int main( int argc, char *argv[] )
    {

    if( argc < 10 )
        {
        std::cerr << "Missing Parameters " << std::endl;
        std::cerr << "Usage: " << argv[0] << std::endl;
        std::cerr << " fixed_file  moving_file " << std::endl;
        std::cerr << " out_file init_trans_step " << std::endl;
        std::cerr << " init_rot_step min_step" << std::endl;
        std::cerr << " trans_x trans_y trans_z " << std::endl;
        return EXIT_FAILURE;
        }

    char * fixed_file      =       argv[1];
    char * moving_file     =       argv[2];
    char * out_file        =       argv[3];
    double init_trans_step = atof( argv[4] );
    double init_rot_step   = atof( argv[5] );
    double min_step        = atof( argv[6] );
    double trans_x         = atof( argv[7] );
    double trans_y         = atof( argv[8] );
    double trans_z         = atof( argv[9] );

    const    unsigned int    Dimension = 3;
    typedef  unsigned char   PixelType;

    typedef itk::Image< PixelType, Dimension >  ImageType;

    typedef itk::Euler3DTransform
        < double > TransformType;

    typedef  itk::AmoebaOptimizer  OptimizerType;

    typedef itk::NormalizedCorrelationImageToImageMetric
        < ImageType, ImageType >      MetricType;

    typedef itk::LinearInterpolateImageFunction
        < ImageType, double >              InterpolatorType;	

    typedef itk::ImageRegistrationMethod
        < ImageType, ImageType >      RegistrationType;

    MetricType::Pointer         metric        = MetricType::New();
    OptimizerType::Pointer      optimizer     = OptimizerType::New();
    InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
    RegistrationType::Pointer   registration  = RegistrationType::New();

    metric->ComputeGradientOff();

    registration->SetMetric(        metric        );
    registration->SetOptimizer(     optimizer     );
    registration->SetInterpolator(  interpolator  );

    TransformType::Pointer  transform = TransformType::New();
    registration->SetTransform( transform );

    typedef itk::ImageFileReader< ImageType  > ImageReaderType;
    ImageReaderType::Pointer  fixedImageReader  = ImageReaderType::New();
    ImageReaderType::Pointer movingImageReader  = ImageReaderType::New();

    fixedImageReader->SetFileName(  fixed_file );
    movingImageReader->SetFileName( moving_file );

    fixedImageReader->Update();

    ImageType::RegionType fixed_region = fixedImageReader->GetOutput()->GetLargestPossibleRegion();
    ImageType::SizeType fixed_size = fixed_region.GetSize();
    ImageType::IndexType fixed_index = fixed_region.GetIndex();
    ImageType::SpacingType fixed_spacing = fixedImageReader->GetOutput()->GetSpacing();
    ImageType::PointType fixed_origin = fixedImageReader->GetOutput()->GetOrigin();

    std::cout << "Fixed x size = " << fixed_size[0] << std::endl;
    std::cout << "Fixed y size = " << fixed_size[1] << std::endl;
    std::cout << "Fixed z size = " << fixed_size[2] << std::endl;
    std::cout << "Fixed x origin = " << fixed_origin[0] << std::endl;
    std::cout << "Fixed y origin = " << fixed_origin[1] << std::endl;
    std::cout << "Fixed z origin = " << fixed_origin[2] << std::endl;

    fixed_size[0] -= (long unsigned int)
        ( fabs( trans_x - fixed_origin[0] ) / fixed_spacing[0] );
    fixed_size[1] -= (long unsigned int)
        ( fabs( trans_y - fixed_origin[1] ) / fixed_spacing[1] );
    fixed_size[2] -= (long unsigned int)
        ( fabs( trans_z - fixed_origin[2] ) / fixed_spacing[2] );

    std::cout << "Fixed x overlap = " << fixed_size[0] << std::endl;
    std::cout << "Fixed y overlap = " << fixed_size[1] << std::endl;
    std::cout << "Fixed z overlap = " << fixed_size[2] << std::endl;

    if ( trans_x > 0.0 )
        fixed_index[0] += (long int)( ( trans_x - fixed_origin[0] ) / fixed_spacing[0] );
    if ( trans_y > 0.0 )
        fixed_index[1] += (long int)( ( trans_y - fixed_origin[1] ) / fixed_spacing[1] );
    if ( trans_z > 0.0 )
        fixed_index[2] += (long int)( ( trans_z - fixed_origin[2] ) / fixed_spacing[2] );

    std::cout << "Fixed x start = " << fixed_index[0] << std::endl;
    std::cout << "Fixed y start = " << fixed_index[1] << std::endl;
    std::cout << "Fixed z start = " << fixed_index[2] << std::endl;

    fixed_region.SetSize( fixed_size );
    fixed_region.SetIndex( fixed_index );

    typedef itk::RegionOfInterestImageFilter
        < ImageType, ImageType > ROIFilterType;
    ROIFilterType::Pointer  fixed_roi_filter = ROIFilterType::New();
    fixed_roi_filter->SetRegionOfInterest(   fixed_region );
    fixed_roi_filter->SetInput(   fixedImageReader->GetOutput() );
    fixed_roi_filter->Update();
    ImageType::Pointer fixed_roi = fixed_roi_filter->GetOutput();

    std::cout << "Deleting used fixed filters." << std::endl;

    fixed_roi_filter = 0; // Delete
    fixedImageReader = 0; // Delete

    movingImageReader->Update();

    ImageType::RegionType moving_region 
        = movingImageReader->GetOutput()->GetLargestPossibleRegion();
    ImageType::SizeType moving_size = moving_region.GetSize();
    ImageType::IndexType moving_index = moving_region.GetIndex();
    ImageType::SpacingType moving_spacing 
        = movingImageReader->GetOutput()->GetSpacing();
    ImageType::PointType moving_origin 
        = movingImageReader->GetOutput()->GetOrigin();

    std::cout << "Moving x size = " << moving_size[0] << std::endl;
    std::cout << "Moving y size = " << moving_size[1] << std::endl;
    std::cout << "Moving z size = " << moving_size[2] << std::endl;

    moving_size[0] -= (long unsigned int)
        ( fabs( trans_x - fixed_origin[0] ) / moving_spacing[0] );
    moving_size[1] -= (long unsigned int)
        ( fabs( trans_y - fixed_origin[1] ) / moving_spacing[1] );
    moving_size[2] -= (long unsigned int)
        ( fabs( trans_z - fixed_origin[2] ) / moving_spacing[2] );

    std::cout << "Moving x overlap = " << moving_size[0] << std::endl;
    std::cout << "Moving y overlap = " << moving_size[1] << std::endl;
    std::cout << "Moving z overlap = " << moving_size[2] << std::endl;

    if ( trans_x < 0.0 )
        moving_index[0] -= (long int)
            ( ( trans_x - fixed_origin[0] ) / moving_spacing[0] );
    if ( trans_y < 0.0 )
        moving_index[1] -= (long int)
            ( ( trans_y - fixed_origin[1] ) / moving_spacing[1] );
    if ( trans_z < 0.0 )
        moving_index[2] -= (long int)
            ( ( trans_z - fixed_origin[2] ) / moving_spacing[2] );

    std::cout << "Moving x start = " << moving_index[0] << std::endl;
    std::cout << "Moving y start = " << moving_index[1] << std::endl;
    std::cout << "Moving z start = " << moving_index[2] << std::endl;

    moving_region.SetSize( moving_size );
    moving_region.SetIndex( moving_index );


    ROIFilterType::Pointer moving_roi_filter = ROIFilterType::New();
    moving_roi_filter->SetRegionOfInterest( moving_region );
    moving_roi_filter->SetInput( movingImageReader->GetOutput() );
    moving_roi_filter->Update();
    ImageType::Pointer moving_roi = moving_roi_filter->GetOutput();

    std::cout << "Deleting used moving filters." << std::endl;

    moving_roi_filter = 0; // Delete
    movingImageReader = 0; // Delete


    registration->SetFixedImageRegion
        ( fixed_roi->GetBufferedRegion() );

    registration->SetFixedImage(   fixed_roi );
    registration->SetMovingImage( moving_roi );

    // Setup Transform
    transform->SetIdentity();
    TransformType::InputPointType transform_center;
    transform_center[0] = 0.0;
    transform_center[1] = 0.0;
    transform_center[2] = 0.0;
    transform->SetCenter( transform_center );
    TransformType::OutputVectorType transform_offset;
    transform_offset[0] = -trans_x;
    transform_offset[1] = -trans_y;
    transform_offset[2] = -trans_z;
    transform->SetOffset( transform_offset );

    TransformType::ParametersType initialParameters = transform->GetParameters();
    registration->SetInitialTransformParameters( initialParameters );

    std::cout << " Transform has " << transform->GetNumberOfParameters() 
        << " parameters." << std::endl;

    typedef OptimizerType::ParametersType AmoebaParametersType;
    AmoebaParametersType initialStep
        ( transform->GetNumberOfParameters() );
    initialParameters[0] = init_rot_step;
    initialParameters[1] = init_rot_step;
    initialParameters[2] = init_rot_step;
    initialParameters[3] = init_trans_step;
    initialParameters[4] = init_trans_step;
    initialParameters[5] = init_trans_step;
    optimizer->SetInitialSimplexDelta( initialParameters );

    optimizer->SetParametersConvergenceTolerance( min_step );
    optimizer->AutomaticInitialSimplexOff();
    optimizer->SetMaximumNumberOfIterations( 2000 );


    CommandIterationUpdateAmoeba::Pointer observer = 
        CommandIterationUpdateAmoeba::New();
    optimizer->AddObserver( itk::FunctionEvaluationIterationEvent(), observer );

    try 
        { 
        registration->StartRegistration(); 
        std::cout << "Restarting registration" <<std::endl;

        TransformType::ParametersType intermediateParameters 
            = registration->GetLastTransformParameters();

        registration->SetInitialTransformParameters( intermediateParameters );

        registration->StartRegistration();
        } 
    catch( itk::ExceptionObject & err ) 
        { 
        std::cerr << "ExceptionObject caught !" << std::endl; 
        std::cerr << err << std::endl; 
        return EXIT_FAILURE;
        } 

    OptimizerType::ParametersType finalParameters = 
        registration->GetLastTransformParameters();

    std::cout << "Final parameters : " << std::endl
        << finalParameters << std::endl;

    TransformType::Pointer final_transform = TransformType::New();
    final_transform->SetParameters( finalParameters );

    TransformType::Pointer inverse_transform = TransformType::New();
    inverse_transform->SetCenter( transform_center );
    final_transform->GetInverse( inverse_transform );

    ImageReaderType::Pointer movingImageReReader  
        = ImageReaderType::New();
    movingImageReReader->SetFileName( moving_file );
    movingImageReReader->Update();

    ImageType::Pointer reread_image 
        = movingImageReReader->GetOutput();
    ImageType::RegionType reread_region 
        = reread_image->GetLargestPossibleRegion();
    ImageType::RegionType::SizeType reread_size 
        = reread_region.GetSize();
    ImageType::SpacingType reread_spacing 
        = reread_image->GetSpacing();
    ImageType::PointType reread_origin 
        = reread_image->GetOrigin();

    ImageType::PointType transformed_origin 
        = inverse_transform->TransformPoint( reread_origin );

    typedef itk::ResampleImageFilter
        < ImageType, ImageType > ResampleFilterType;
    ResampleFilterType::Pointer resample = ResampleFilterType::New();
    resample->SetDefaultPixelValue( 0 );
    resample->SetTransform( final_transform );
    resample->SetOutputOrigin( transformed_origin );
    resample->SetOutputSpacing( reread_spacing );
    resample->SetSize( reread_size );
    resample->SetInput( movingImageReReader->GetOutput() );


    typedef itk::ImageFileWriter< ImageType > MovingWriterType;
    MovingWriterType::Pointer moving_writer = MovingWriterType::New();
    moving_writer->SetFileName( out_file );
    moving_writer->SetInput( resample->GetOutput() );

    std::cout << "Writing output image...";

    try 
        { 
        moving_writer->Update();
        } 
    catch( itk::ExceptionObject & err ) 
        { 
        std::cerr << "ExceptionObject caught !" << std::endl; 
        std::cerr << err << std::endl; 
        return EXIT_FAILURE;
        } 

    std::cout << "Done" << std::endl;



    moving_writer = 0; // Delete
    reread_image = 0; // Delete
    movingImageReReader = 0; // Delete
    observer = 0; // Delete
    moving_roi = 0; // Delete
    fixed_roi = 0; // Delete
    transform = 0; // Delete
    registration = 0; // Delete
    interpolator = 0; // Delete
    optimizer = 0; // Delete
    metric = 0; // Delete


    return 0;
    }

