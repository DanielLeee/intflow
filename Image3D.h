#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include "Array3D.h"
#include "dbh.h"
#include "stdlib.h"
using namespace std;

extern int NumScales;
extern T_input GrayAvgBaseLine; // 4.0 ADNI: 55.0 
extern T_input SIFTAvgBaseLine; // 8.0 ADNI: 1.0
extern const int OriSiftDim; // 48
extern const int SiftDim;

/*
static int NumScales = 1;
const int OriSiftDim = 48; // 48
const int SiftDim = OriSiftDim * NumScales + 1;
static T_input GrayAvgBaseLine = 4.0; // 4.0 ADNI: 55.0 
static T_input SIFTAvgBaseLine = 8.0; // 8.0 ADNI: 1.0
*/

class Image3D
{
public:
	Array3D_float *data;
	Array3D_vect *sift;
	Array3D_int *mask;
	string dataFileName;
	int Layer, Height, Width;
	struct dsr dataHDR;

	void ParseHDR(string fileName, struct dsr& hdr);
	void ReadFloatImg(Array3D_float* target, string fileName, int dataType);
	void ReadIntImg(Array3D_int* target, string fileName, int dataType);
	void WriteHDR(string fileName, struct dsr& hdr);
	void WriteIMG(string fileName, Array3D_float& data, int dataType);
	void WriteIMG(string fileName, Array3D_int& data, int dataType);
	void OutputDataImage(string fileName);
	void ChangeSIFTByCof(T_input cof);
	void ChangeGrayByCof(T_input cof);
	bool InBound(int x, int y, int z);
	void ComputeSIFT(Array3D_float *imsrc, Array3D_vect *imsift, int arrOffset = 0, int cellSize = 1, int stepSize = 1, bool isBoundaryIncluded = true, int nBins = 6, int alpha = 1, int winSize = 1);
	void ComputeDerivatives(Array3D_float *imsrc, Array3D_float *dx, Array3D_float *dy, Array3D_float *dz, bool isAdvancedFilter);
	void GetDirections(T_input dirs[][3], int nBins);
	void ThreeDimensionalFiltering(Array3D_vect *imsrc, Array3D_vect *imdst, int fSize);

	Image3D(string _dataFileName);
	Image3D();
	~Image3D();
};

class Image3DWithSeg : public Image3D
{
public:
	Array3D_int *seg;
	string segFileName;
	struct dsr segHDR;

	void OutputSegImage(string fileName);
	void OutputDataSeg(string fileName0, string fileName1);
    
    void ApplyFlow(string fileName, Image3DWithSeg *refImg);
    void PrintFlowDetail(string flowFilePath, string outTxtPath);
    void OutputDeformationAmount(string flowFilePath, string outTxtPath);

	Image3DWithSeg(string dataFileName, string _segFileName);
	Image3DWithSeg();
	~Image3DWithSeg();
};

