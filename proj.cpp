#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkLookupTable.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <vtkShortArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkCutter.h>
#include <vtkPlane.h>
#include <vtkPointData.h>
#include <vtkGraphicsFactory.h>
#include <vtkPNGWriter.h>
#include <vtkImageData.h>
#include <vtkCamera.h>
#include <vtkContourFilter.h>
#include <vtkProperty.h>
#include <vtkVectorNorm.h>
#include <ctime>
#include <cmath>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>

#include "config.h"

#define BIG

//*/
#define FICHIER MY_MESHES_PATH "/Frog.raw_256_256_44_CHAR.raw"

int gridSize = 256;
int YgridSize = 256;
int ZgridSize = 44;

#define CHAR

int startexploreval = 13;
int endexploreval = 20;

//*/

/*/
#define FICHIER MY_MESHES_PATH "/Mystere5_2048_2048_756_SHORT.raw"
 
 int gridSize = 2048;
 int YgridSize = 2048;
 int ZgridSize = 756;

 #define SHORT

int startexploreval=1;
int endexploreval=65000;

 //*/

/*/
#define FICHIER MY_MESHES_PATH "/Mystere6_1118_2046_694_CHAR.raw"

int gridSize = 1118;
int YgridSize = 2046;
int ZgridSize = 694;

#define CHAR

int startexploreval=1;
int endexploreval=255;

//*/


/*
#define FICHIER MY_MESHES_PATH "/Mystere1_512_512_134_SHORT.raw"
int gridSize = 512;
int YgridSize = 512;
int ZgridSize = 134;
 
 #define SHORT

 int startexploreval=1;
 int endexploreval=65000;
//*/


/*
#define FICHIER MY_MESHES_PATH "/Mystere4_512_512_322_SHORT.raw"
int gridSize = 512;
int YgridSize = 512;
int ZgridSize = 322;
 
 #define SHORT

 int startexploreval=1;
 int endexploreval=65000;
//*/


const char *location = FICHIER;

int winSize = 500;
int numPasses = 10;
int nbimages = 10;

const char *prefix = "";
int passNum = 0;

using std::cerr;
using std::endl;

// Function prototypes
vtkRectilinearGrid *ReadGrid(int zStart, int zEnd);
void WriteImage(const char *name, const float *rgba, int width, int height);
bool ComposeImageZbuffer(float *rgba_out, float *zbuffer, int image_width, int image_height);


int main(int argc, char *argv[]) {
	srand(time(nullptr));
//    int zStart = 0;
//    int zEnd = ZgridSize;

	int npixels = winSize * winSize;
	vtkRectilinearGrid *reader = nullptr;


	vtkLookupTable *lut = vtkLookupTable::New();

	lut->SetHueRange(0.1, 0.0);
	lut->SetSaturationRange(0.0, 1.0);
	lut->SetValueRange(1.0, 255.0);
	lut->SetNumberOfColors(100);
	lut->Build();


	vtkRenderer *ren = vtkRenderer::New();
	double bounds[6] = {0.00001, 1 - 0.00001, 0.00001, 1 - 0.00001, 0.00001, 1 - 0.00001};

	ren->ResetCamera(bounds);

	bool once = true;
	vtkRenderWindow *renwin = vtkRenderWindow::New();
	renwin->SetSize(winSize, winSize);

	renwin->AddRenderer(ren);


	int zStart = 0;


	int zEnd = ZgridSize;


	reader = ReadGrid(zStart, zEnd);

	vtkContourFilter *cf = vtkContourFilter::New();
	cf->SetNumberOfContours(1);
	// cf->SetValue(0, 20.0);
	int valcont = startexploreval;
	cf->SetValue(1, valcont);

	cf->SetInputData(reader);

	int maxsize = std::max(gridSize, std::max(YgridSize, ZgridSize));
	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	transform->Scale(gridSize / (float) maxsize, YgridSize / (float) maxsize, ZgridSize / (float) maxsize);
	vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
	transformFilter->SetInputConnection(cf->GetOutputPort());
	transformFilter->SetTransform(transform);

	vtkDataSetMapper *mapper = vtkDataSetMapper::New();
	//     mapper->SetInputData( reader );
	mapper->SetInputConnection(transformFilter->GetOutputPort());

	vtkActor *actor = vtkActor::New();
	actor->SetMapper(mapper);

	float startval = valcont;
	float endval = endexploreval;

	mapper->SetScalarRange(startexploreval, endexploreval);

	mapper->SetLookupTable(lut);

	ren->AddActor(actor);
	ren->SetViewport(0, 0, 1, 1);


	// vtkCamera *cam = ren->GetActiveCamera();
	vtkCamera *cam;
//            if(once){
//                once=false;
//                cam= ren->GetActiveCamera();
//             cam->SetFocalPoint(0.5, 0.5, 0.5);
//                cam->SetPosition(0., .0, 3.);
//                cam->SetViewUp(0., -1.0, 0.0);
//                bounds=ren->ComputeVisiblePropBounds();
//                std::cout<<"bounds : "<<bounds[0]<<" ; "<<bounds[1]<<" ; "<<bounds[2]<<" ; "<<bounds[3]<<" ; "<<bounds[4]<<" ; "<<bounds[5]<<" ; "<<endl<<fflush;
//            }
//            else{ren->SetActiveCamera(cam);
//               
//            }

	// cam= ren->GetActiveCamera();
	//cam->SetFocalPoint(0.5, 0.5, 0.5);
	//cam->SetPosition(-0., .0, 3.);
	//cam->SetViewUp(0., -1.0, 0.0);
	//bounds=ren->ComputeVisiblePropBounds();
	// std::cout<<"bounds : "<<bounds[0]<<" ; "<<bounds[1]<<" ; "<<bounds[2]<<" ; "<<bounds[3]<<" ; "<<bounds[4]<<" ; "<<bounds[5]<<" ; "<<endl<<fflush;

	//       ren->SetActiveCamera (cam);
	//           cam= ren->GetActiveCamera();
	//   cam->SetPosition(0.5, 3.0, 0.5);
	//  cam->SetViewUp(0., 0.0, 1.0);
	//      cam->SetFocalPoint(0., 0.0, 0.);

//        

	//*
	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renwin);
	renwin->Render();
	iren->Start();
	//*/

	reader->Delete();
	mapper->Delete();
	cf->Delete();
	//          lut->Delete();
	ren->RemoveActor(actor);
	actor->Delete();

	ren->Delete();
	renwin->Delete();
}

// You should not need to modify these routines.

vtkRectilinearGrid *
ReadGrid(int zStart, int zEnd) {
	int i;

	/*   if (zStart < 0 || zEnd < 0 || zStart >= gridSize || zEnd >= gridSize || zStart > zEnd)
	 {
	 cerr << prefix << "Invalid range: " << zStart << "-" << zEnd << endl;
	 return NULL;
	 }
	 */
	std::ifstream ifile(location);
	if (ifile.fail()) {
		cerr << prefix << "Unable to open file: " << location << "!!" << endl;
		throw std::runtime_error("can't find the file!! Check the name and the path of this file? ");
	}

	cerr << prefix << "Reading from " << zStart << " to " << zEnd << endl;

	vtkRectilinearGrid *rg = vtkRectilinearGrid::New();
	vtkFloatArray *X = vtkFloatArray::New();
	X->SetNumberOfTuples(gridSize);
	for (i = 0; i < gridSize; i++)
		X->SetTuple1(i, i / (gridSize - 1.0));
	rg->SetXCoordinates(X);
	X->Delete();
	vtkFloatArray *Y = vtkFloatArray::New();
	Y->SetNumberOfTuples(YgridSize);
	for (i = 0; i < YgridSize; i++)
		Y->SetTuple1(i, i / (YgridSize - 1.0));
	rg->SetYCoordinates(Y);
	Y->Delete();
	vtkFloatArray *Z = vtkFloatArray::New();
	int numSlicesToRead = zEnd - zStart + 1;
	Z->SetNumberOfTuples(numSlicesToRead);
	for (i = zStart; i <= zEnd; i++)
		Z->SetTuple1(i - zStart, i / (ZgridSize - 1.0));
	rg->SetZCoordinates(Z);
	Z->Delete();

	rg->SetDimensions(gridSize, YgridSize, numSlicesToRead);

	unsigned int valuesPerSlice = gridSize * YgridSize;

#if defined(SHORT)
	unsigned int bytesPerSlice   = sizeof(short)*valuesPerSlice;

#elif defined(CHAR)
	unsigned int bytesPerSlice = sizeof(char) * valuesPerSlice;

#elif  defined(FLOAT)
	unsigned int bytesPerSlice   = sizeof(float)*valuesPerSlice;

#else
#error Unsupported choice setting
#endif


#if defined(SMALL)
	unsigned int offset          = (unsigned int)zStart * (unsigned int)bytesPerSlice;
   unsigned int bytesToRead     = bytesPerSlice*numSlicesToRead;
   unsigned int valuesToRead    = valuesPerSlice*numSlicesToRead;
#elif defined(BIG)
	unsigned long long offset = (unsigned long long) zStart * bytesPerSlice;
	unsigned long long bytesToRead = (unsigned long long) bytesPerSlice * numSlicesToRead;
	unsigned int valuesToRead = (unsigned int) valuesPerSlice * numSlicesToRead;
#else
#error Unsupported choice setting
#endif


#if defined(SHORT)
	vtkUnsignedShortArray *scalars = vtkUnsignedShortArray::New();
	scalars->SetNumberOfTuples(valuesToRead);
	unsigned short *arr = scalars->GetPointer(0);

#elif defined(CHAR)
	vtkUnsignedCharArray *scalars = vtkUnsignedCharArray::New();
	scalars->SetNumberOfTuples(valuesToRead);
	unsigned char *arr = scalars->GetPointer(0);

#elif  defined(FLOAT)
	vtkFloatArray *scalars = vtkFloatArray::New();
	scalars->SetNumberOfTuples(valuesToRead);
	float *arr = scalars->GetPointer(0);
#else
#error Unsupported choice setting
#endif


	ifile.seekg(offset, ios::beg);
	ifile.read((char *) arr, bytesToRead);
	ifile.close();

	int min = +2147483647;
	int max = 0;


	for (int i = 0; i < valuesToRead; i++) {
		if (min > (scalars->GetPointer(0))[i]) min = (scalars->GetPointer(0))[i];
		if (max < (scalars->GetPointer(0))[i]) max = (scalars->GetPointer(0))[i];

		if (rand() % (valuesToRead / 20) == 0) {
#if defined(SHORT)
			std::cout<<(unsigned short)(scalars->GetPointer(0))[i]<<" ";

#elif defined(CHAR)
			std::cout << +(unsigned char) (scalars->GetPointer(0))[i] << " ";

#elif  defined(FLOAT)
			std::cout<<(float)(scalars->GetPointer(0))[i]<<" ";

#else
#error Unsupported choice setting
#endif


		}
	}


	std::cout << endl << fflush;
	std::cout << "min value read: " << min << endl << fflush;
	std::cout << "max value read: " << max << endl << fflush;


	scalars->SetName("entropy");
	rg->GetPointData()->AddArray(scalars);
	scalars->Delete();

	vtkFloatArray *pr = vtkFloatArray::New();
	pr->SetNumberOfTuples(valuesToRead);
	for (int i = 0; i < valuesToRead; i++)
		pr->SetTuple1(i, passNum);

	pr->SetName("pass_num");
	rg->GetPointData()->AddArray(pr);
	pr->Delete();

	rg->GetPointData()->SetActiveScalars("entropy");

	cerr << prefix << " Done reading" << endl;
	return rg;
}


void WriteImage(const char *name, const float *rgba, int width, int height) {
	vtkImageData *img = vtkImageData::New();
	img->SetDimensions(width, height, 1);
#if VTK_MAJOR_VERSION <= 5
	img->SetNumberOfScalarComponents(3);
	img->SetScalarTypeToUnsignedChar();
#else
	img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
#endif

	for (int i = 0; i < width; i++)
		for (int j = 0; j < height; j++) {
			unsigned char *ptr = (unsigned char *) img->GetScalarPointer(i, j, 0);
			int idx = j * width + i;
			ptr[0] = (unsigned char) (255 * rgba[4 * idx + 0]);
			ptr[1] = (unsigned char) (255 * rgba[4 * idx + 1]);
			ptr[2] = (unsigned char) (255 * rgba[4 * idx + 2]);
		}

	vtkPNGWriter *writer = vtkPNGWriter::New();
	writer->SetInputData(img);
	writer->SetFileName(name);
	writer->Write();

	img->Delete();
	writer->Delete();
}

bool ComposeImageZbuffer(float *rgba_out, float *zbuffer, int image_width, int image_height) {
	int npixels = image_width * image_height;

	float min = 1;
	float max = 0;
	for (int i = 0; i < npixels; i++) {
		if (zbuffer[i] < min) min = zbuffer[i];
		if (zbuffer[i] > max) max = zbuffer[i];

	}
	std::cout << "min:" << min;
	std::cout << "max:" << max << "  ";

	float coef = 1.f / ((max - min));

	std::cout << "coef:" << coef << "  ";

	for (int i = 0; i < npixels; i++) {
		rgba_out[i * 4] = (zbuffer[i] == 1.0 ? 0 : 1 - coef * (zbuffer[i] - min));
		rgba_out[i * 4 + 1] = (zbuffer[i] == 1.0 ? 0 : 1 - coef * (zbuffer[i] - min));
		rgba_out[i * 4 + 2] = (zbuffer[i] == 1.0 ? 0 : 1 - coef * (zbuffer[i] - min));
		rgba_out[i * 4 + 3] = 0.0;
	}
	return false;
}

