#include <string>
#include <iostream>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkLookupTable.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <vtkCutter.h>
#include <vtkPlane.h>
#include <vtkPointData.h>
#include <vtkImageData.h>
#include <vtkCamera.h>
#include <vtkContourFilter.h>
#include <vtkProperty.h>
#include <ctime>
#include <cmath>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkPNGWriter.h>
#include <vtkUnsignedShortArray.h>
#include <mpi.h>

#include "config.h"

using Time = std::chrono::time_point<std::chrono::system_clock>;
#define NOW std::chrono::system_clock::now()

#define BIG

#define FICHIER MY_MESHES_PATH "/Mystere4_512_512_322_SHORT.raw"

int gridSize = 512;
int YgridSize = 512;
int ZgridSize = 322;

//#define CHAR
#define SHORT

int startexploreval = 54000;
int endexploreval = 58000;

const char *location = FICHIER;

int winSize = 500;
int numPasses = 20;
int passNum = 0;

const char *prefix = "";

int pid = -1, nprocs = -1;

// Function prototypes
vtkRectilinearGrid *ReadGrid(int zStart, int zEnd);
void WriteImage(const char *name, const float *rgba, int width, int height);
vtkRectilinearGrid* ParallelReadGrid();
void CompositeImage(const float *rgba_in, float *zbuffer, float *rgba_out, int image_width, int image_height);

int main(int argc, char *argv[]) {

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	Time startTime = NOW;

	srand(time(nullptr));

	vtkLookupTable *lut = vtkLookupTable::New();
	lut->SetHueRange(0.1, 0.0);
	lut->SetSaturationRange(0.0, 1.0);
	lut->SetValueRange(1.0, 255.0);
	lut->SetNumberOfColors(100);
	lut->Build();

	vtkRenderer *ren = vtkRenderer::New();
	double bounds[6] = {0.00001, 1 - 0.00001, 0.00001, 1 - 0.00001, 0.00001, 1 - 0.00001};
	ren->ResetCamera(bounds);

	vtkRenderWindow *renwin = vtkRenderWindow::New();
	renwin->SetOffScreenRendering(true);
	renwin->SetSize(winSize, winSize);
	renwin->AddRenderer(ren);

	vtkRectilinearGrid *reader = ParallelReadGrid();

	vtkContourFilter *cf = vtkContourFilter::New();
	cf->SetNumberOfContours(1);
	int valcont = startexploreval;
	cf->SetValue(1, valcont);
	cf->SetInputData(reader);

	int maxsize = std::max(gridSize, std::max(YgridSize, ZgridSize));
	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	transform->Scale(
			(float) gridSize / (float) maxsize,
			(float) YgridSize / (float) maxsize,
			(float) ZgridSize / (float) maxsize
	);
	//transform->Scale(1.f, 1.f, 1.f);
	vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
	transformFilter->SetInputConnection(cf->GetOutputPort());
	transformFilter->SetTransform(transform);

	vtkDataSetMapper *mapper = vtkDataSetMapper::New();
	mapper->SetInputConnection(transformFilter->GetOutputPort());

	vtkActor *actor = vtkActor::New();
	actor->SetMapper(mapper);

	mapper->SetScalarRange(startexploreval, endexploreval);
	mapper->SetLookupTable(lut);

	ren->AddActor(actor);
	ren->SetViewport(0, 0, 1, 1);

	vtkCamera *cam;
	cam = ren->GetActiveCamera();
	cam->SetFocalPoint(0.5, 0.5, 0.5);
	cam->SetPosition(-0., .0, 3.);
	cam->SetViewUp(0., -1.0, 0.0);
	//bounds = ren->ComputeVisiblePropBounds();

	float *rgba, *zbuffer;

	renwin->Render();
	rgba = renwin->GetRGBAPixelData(0, 0, (winSize - 1), (winSize - 1), 1);
	zbuffer = renwin->GetZbufferData(0, 0, (winSize - 1), (winSize - 1));

	auto *new_rgba = new float[4 * winSize * winSize];

	CompositeImage(rgba, zbuffer, new_rgba, winSize, winSize);
	if (pid == 0) WriteImage("final.png", new_rgba, winSize, winSize);

	free(rgba);
	free(zbuffer);
	free(new_rgba);

	//reader = ReadGrid(0, ZgridSize);
	//cf->SetInputData(reader);
	//cf->Update();
	//
	//vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	//iren->SetRenderWindow(renwin);
	//renwin->SetOffScreenRendering(false);
	//renwin->Render();
	//iren->Start();

	std::chrono::duration<double> elapsed_seconds = NOW - startTime;
	std::clog << "Time: " << elapsed_seconds.count() << "s" << std::endl;

	//iren->Delete();
	reader->Delete();
	mapper->Delete();
	cf->Delete();
	lut->Delete();
	ren->RemoveActor(actor);
	actor->Delete();
	ren->Delete();
	renwin->Delete();

	MPI_Finalize();
}

vtkRectilinearGrid *ReadGrid(int zStart, int zEnd) {
	int i;

	std::ifstream ifile(location);
	if (ifile.fail()) {
		//std::clog << prefix << "Unable to open file: " << location << "!!" << std::endl;
		throw std::runtime_error("can't find the file!! Check the name and the path of this file? ");
	}

	std::clog << prefix << "Reading from " << zStart << " to " << zEnd << std::endl;

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
	unsigned int bytesPerSlice = sizeof(short) * valuesPerSlice;

#elif defined(CHAR)
	unsigned int bytesPerSlice = sizeof(char) * valuesPerSlice;

#elif  defined(FLOAT)
	unsigned int bytesPerSlice   = sizeof(float)*valuesPerSlice;

#else
#error Unsupported choice setting
#endif


#if defined(SMALL)
	unsigned int offset      = (unsigned int) zStart * (unsigned int) bytesPerSlice;
   unsigned int bytesToRead  = bytesPerSlice*numSlicesToRead;
   unsigned int valuesToRead = valuesPerSlice*numSlicesToRead;
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
			//std::cout << (unsigned short) (scalars->GetPointer(0))[i] << " ";
#elif defined(CHAR)
			//std::cout << +(unsigned char) (scalars->GetPointer(0))[i] << " ";
#elif  defined(FLOAT)
			std::cout<<(float)(scalars->GetPointer(0))[i]<<" ";
#else
#error Unsupported choice setting
#endif


		}
	}

	//std::cout << endl << fflush;
	//std::cout << "min value read: " << min << endl << fflush;
	//std::cout << "max value read: " << max << endl << fflush;

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

	//std::clog << prefix << " Done reading" << std::endl;
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
			auto *ptr = (unsigned char *) img->GetScalarPointer(i, j, 0);
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

vtkRectilinearGrid* ParallelReadGrid() {
	int zCount = pid < ZgridSize % nprocs
				 ? (ZgridSize / nprocs) + 1
				 : ZgridSize / nprocs;
	int zStart = pid < ZgridSize % nprocs
				 ? ((ZgridSize / nprocs) + 1) * pid
				 : ((ZgridSize / nprocs) + 1) * (ZgridSize % nprocs) + (ZgridSize / nprocs) * (pid - (ZgridSize % nprocs));
	int zEnd = zStart + zCount;

	if (zEnd == ZgridSize) --zEnd;

	return ReadGrid(zStart, zEnd);
}

void CompositeImage(const float *rgba_in, float *zbuffer, float *rgba_out, int image_width, int image_height) {
	int npixels = image_width * image_height;

	auto *zbuffer_min = new float[npixels];
	MPI_Allreduce(
			zbuffer, zbuffer_min,
			npixels, MPI_FLOAT, MPI_MIN,
			MPI_COMM_WORLD
	);

	auto *rgba_tmp = new float[4 * npixels];
	for (int i = 0; i < npixels; ++i) {
		if (zbuffer[i] <= zbuffer_min[i]) {
			rgba_tmp[i * 4] = rgba_in[i * 4];
			rgba_tmp[i * 4 + 1] = rgba_in[i * 4 + 1];
			rgba_tmp[i * 4 + 2] = rgba_in[i * 4 + 2];
			rgba_tmp[i * 4 + 3] = rgba_in[i * 4 + 3];
		}
	}
	MPI_Reduce(rgba_tmp, rgba_out, 4 * npixels, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

	delete[] zbuffer_min;
	delete[] rgba_tmp;
}