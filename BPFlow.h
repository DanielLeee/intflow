#pragma once
#include "Array3D.h"
#include "memory.h"
#include <string>
#include <cstring>
#include "stdio.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "stdlib.h"
#include "float.h"
using namespace std;

const int nNeighbors = 6;

//-----------------------------------------------------------------------------------
// the class for BP flow
//-----------------------------------------------------------------------------------
class BPFlow
{
public:
	int Height, Width, Slice, Volume, nChannels, Height2, Width2, Slice2;
	int winSize;
	T_input s, d, gamma;

	Array3D_vect *pIm1, *pIm2;
	Array3D_vect *Im_s, *Im_d;
	Array3D_vect *foo;
	
	Array3D_int *pOffset[3];
	Array3D_int *pX[3];
	Array3D_vect *pDataTerm;
	Array3D_vect *pRangeTerm[3];
	Array3D_vect *pSpatialMessage[3][nNeighbors];
	Array3D_vect *pDualMessage[3];
	Array3D_vect *pBelief[3];

	BPFlow();
	~BPFlow();

	void setPara(T_input _s, T_input _d);
	void LoadImages(Array3D_vect* pImage1, Array3D_vect* pImage2);
	void LoadOffset(Array3D_int *vx, Array3D_int *vy, Array3D_int *vz);
	void LoadWinSize(int wsize);

	void ComputeDataTerm();
	void ComputeRangeTerm(T_input _gamma);
	void ComputeBelief();
	void ComputeVelocity(Array3D_vect *mFlow);
	void FindOptimalSolution();
	void AllocateMessage();
	void AllocateBuffer();

	void MessagePassing(int nIterations, int nHierarchy, T_input &energy, bool onTop);
	void BP_S(int count);
	
	template<class T>
	bool InsideImage(T x,T y,T z);

	int EnforceRange(const int x, const int MaxValue);

	void UpdateSpatialMessage(int x,int y,int z,int plane,int direction);
	void UpdateDualMessage(int x,int y,int z,int plane);
	void Add2Message(T_input *m1, T_input *m2, int len, T_input cof = 1);

	T_input GetEnergy();

	//------------------------------------------------------------------------
	// multi-grid belief propagation
	void generateCoarserLevel(BPFlow& bp);
	void propagateFinerLevel(BPFlow& bp);	
	void Smoothing(Array3D_vect *im1, Array3D_vect *im2);
	void ImageFilter(Array3D_vect *src, Array3D_vect *dst, T_input *filter, int fsize);
	void ResizeImageHalf(Array3D_vect *src, Array3D_vect *dst);
	void ReduceImageHalf(Array3D_vect *src, Array3D_vect *dst, bool computeAvg = true);
	void ReduceImageHalf(Array3D_int *src, Array3D_int *dst, bool computeAve = true);
	T_input FindVectorMin(T_input *v, int len);
};
