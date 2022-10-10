#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include "Image3D.h"
#include "BPFlow.h"
#include "CORE.h"
#include "stdlib.h"
#include <algorithm>

static T_input alpha = 3; // 3   ADNI: 2
static T_input d = alpha * 20; // alpha * 20
static T_input gammaPara = 0.005; // 0.005
static int nlevels = 6; // 6   ADNI: 3
static int wsize = 3; // 3   ADNI: 2
static int nIterations = 60; // 60

struct Pyramid
{
    Array3D_vect* im1;
    Array3D_vect* im2;
    Array3D_int* xx, *yy, *zz;
} *pyrd = new Pyramid[nlevels + 1];

void SIFTflowc2f(Array3D_vect* im1, Array3D_vect* im2, Array3D_int*& vx, Array3D_int*& vy, Array3D_int*& vz, T_input &energy);

int iround(T_input x)
{
    return int(round(x));
}

void OutputFlow(string outputFlowName, Array3D_int* vx, Array3D_int* vy, Array3D_int* vz)
{
    int N = vx->Volume;

    ofstream fout(outputFlowName.c_str(), ios::binary);

    fout.write((char *)vx->data, sizeof(int) * N);
    fout.write((char *)vy->data, sizeof(int) * N);
    fout.write((char *)vz->data, sizeof(int) * N);

    fout.close();

    return;
}

void SIFT_FLOW_ALGO(Image3DWithSeg* fixImg, Image3DWithSeg* movImg, string outputFlowName, T_input &energy)
{
    Array3D_vect* im1 = fixImg->sift;
    Array3D_vect* im2 = movImg->sift;
    Array3D_int *vx, *vy, *vz;

    vx = new Array3D_int(fixImg->Layer, fixImg->Height, fixImg->Width);
    vy = new Array3D_int(fixImg->Layer, fixImg->Height, fixImg->Width);
    vz = new Array3D_int(fixImg->Layer, fixImg->Height, fixImg->Width);

    SIFTflowc2f(im1, im2, vx, vy, vz, energy);

    OutputFlow(outputFlowName, vx, vy, vz);

    delete vx;
    delete vy;
    delete vz;

    return;
}

Array3D_vect* imresize3D(Array3D_vect* im)
{
    int Layer = im->Layer;
    int Height = im->Height;
    int Width = im->Width;
    int Layer2 = (Layer + 1) / 2;
    int Height2 = (Height + 1) / 2;
    int Width2 = (Width + 1) / 2;

    Array3D_vect *ret = new Array3D_vect(Layer2, Height2, Width2, SiftDim);
    for(int i = 0; i < Layer2; i++)
    for(int j = 0; j < Height2; j++)
    for(int k = 0; k < Width2; k++)
    for (int u = 0; u < SiftDim; u++)
    {
        T_input v = 0;
        int temsum = 0;
        for(int di = 0; di < 2; di++)
        for(int dj = 0; dj < 2; dj++)
        for (int dk = 0; dk < 2; dk++)
        {
            int ii = 2 * i + di;
            int jj = 2 * j + dj;
            int kk = 2 * k + dk;
            if (0 <= ii && ii < Layer && 0 <= jj && jj < Height && 0 <= kk && kk < Width)
            {
                v += im->at(ii, jj, kk)[u];
                temsum++;
            }
        }
        if (temsum == 0)
            ret->at(i, j, k)[u] = v;
        else
            ret->at(i, j, k)[u] = v / temsum;
    }
    return ret;
}

void SIFTflowc2f(Array3D_vect* im1, Array3D_vect* im2, Array3D_int*& vx, Array3D_int*& vy, Array3D_int*& vz, T_input &energy)
{

    pyrd[1].im1 = im1;
    pyrd[1].im2 = im2;

    for(int i = 2; i <= nlevels; i++)
    {
        pyrd[i].im1 = imresize3D(pyrd[i - 1].im1);
        pyrd[i].im2 = imresize3D(pyrd[i - 1].im2);
    }

    for(int i = 1; i <= nlevels; i++)
    {
        int Layer = pyrd[i].im1->Layer;
        int Height = pyrd[i].im1->Height;
        int Width = pyrd[i].im1->Width;
        int Layer2 = pyrd[i].im2->Layer;
        int Height2 = pyrd[i].im2->Height;
        int Width2 = pyrd[i].im2->Width;
        
        pyrd[i].xx = new Array3D_int(Layer, Height, Width);
        pyrd[i].yy = new Array3D_int(Layer, Height, Width);
        pyrd[i].zz = new Array3D_int(Layer, Height, Width);

        // Be careful with the Layer-Height-Width Order
        for(int x = 0; x < Layer; x++)
        for(int y = 0; y < Height; y++)
        for(int z = 0; z < Width; z++)
        {
            pyrd[i].xx->at(x, y, z) = iround(z * (Width2 - 1) / (Width - 1) - z);
            pyrd[i].yy->at(x, y, z) = iround(y * (Height2 - 1) / (Height - 1) - y);
            pyrd[i].zz->at(x, y, z) = iround(x * (Layer2 - 1) / (Layer - 1)- x);
        }
    }

    for (int i = nlevels; i >= 1; i--)
    {
        //--------------Debug Only--------------Begin
        cout << "Pyramid " << i << endl;
        //--------------Debug Only--------------End

        int Layer = pyrd[i].im1->Layer;
        int Height = pyrd[i].im1->Height;
        int Width = pyrd[i].im1->Width;

        Array3D_vect *flow;

        Array3D_vect* Im1 = pyrd[i].im1;
        Array3D_vect* Im2 = pyrd[i].im2;

        argumentArray arg;

        if(i == nlevels) // top level
        {
            vx = pyrd[i].xx;
            vy = pyrd[i].yy;
            vz = pyrd[i].zz;
        }
        else
        {
            Array3D_int* old_vx = vx;
            Array3D_int* old_vy = vy;
            Array3D_int* old_vz = vz;

            vx = new Array3D_int(Layer, Height, Width);
            vy = new Array3D_int(Layer, Height, Width);
            vz = new Array3D_int(Layer, Height, Width);

            for(int x = 0; x < Layer; x++)
            for(int y = 0; y < Height; y++)
            for(int z = 0; z < Width; z++)
            {
                vx->at(x, y, z) = pyrd[i].xx->at(x, y, z) + (old_vx->at(x/2, y/2, z/2) - pyrd[i+1].xx->at(x/2, y/2, z/2)) * 2;
                vy->at(x, y, z) = pyrd[i].yy->at(x, y, z) + (old_vy->at(x/2, y/2, z/2) - pyrd[i+1].yy->at(x/2, y/2, z/2)) * 2;
                vz->at(x, y, z) = pyrd[i].zz->at(x, y, z) + (old_vz->at(x/2, y/2, z/2) - pyrd[i+1].zz->at(x/2, y/2, z/2)) * 2;
            }

            delete old_vx;
            delete old_vy;
            delete old_vz;
        }

        arg = argumentArray(alpha, d, gammaPara*pow(2, i - 1), nIterations, max(2, nlevels - i), wsize + i - 1);

        flow = new Array3D_vect(Layer, Height, Width, 3);
        mexDiscreteFlow(Im1, Im2, arg, vx, vy, vz, energy, flow);

        for (int x = 0; x < Layer; x++)
        for (int y = 0; y < Height; y++)
        for (int z = 0; z < Width; z++)
        {
            vx->at(x, y, z) = (int)(flow->at(x, y, z)[0]);
            vy->at(x, y, z) = (int)(flow->at(x, y, z)[1]);
            vz->at(x, y, z) = (int)(flow->at(x, y, z)[2]);
        }

        delete flow;

        // For only output flow in the future, you should delete this block
        if (i > 1)
        {
            delete pyrd[i].im1;
            delete pyrd[i].im2;
        }
    }

    delete[] pyrd;

    return;
}
