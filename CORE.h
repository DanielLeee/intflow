#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include "Image3D.h"
#include "BPFlow.h"
#include "CORE.h"

struct argumentArray
{
	T_input alpha;
	T_input d;
	T_input gamma;
	int nIterations;
	int nHierarchy;
	int wsize;
	argumentArray(){}
	argumentArray(T_input alpha, T_input d, T_input gamma, int nIterations, int nHierarchy, int wsize)
	{
		this->alpha = alpha;
		this->d = d;
		this->gamma = gamma;
		this->nIterations = nIterations;
		this->nHierarchy = nHierarchy;
		this->wsize = wsize;
	}
};

void mexDiscreteFlow(Array3D_vect *im1, Array3D_vect *im2, argumentArray args, Array3D_int *vx, Array3D_int *vy, Array3D_int *vz, T_input &energy, Array3D_vect *mFlow)
{
	BPFlow bpflow;

	bpflow.LoadImages(im1, im2);
	bpflow.setPara(args.alpha,args.d);
	bpflow.LoadOffset(vx,vy,vz);
	bpflow.LoadWinSize(args.wsize);
	bpflow.ComputeDataTerm();
	bpflow.ComputeRangeTerm(args.gamma);
	//--------------Debug Only--------------Begin
	cout << "BPFLOW " << "MessagePassing" << endl;
	//--------------Debug Only--------------End
	bpflow.MessagePassing(args.nIterations, args.nHierarchy, energy, true);
	bpflow.ComputeVelocity(mFlow);

	return;
}
