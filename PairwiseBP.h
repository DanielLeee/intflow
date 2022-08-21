#pragma once
#include <iostream>
#include <fstream>
#include "Array3D.h"
#include "float.h"
#include <algorithm>
using namespace std;

int dx[] = { 1, 0, 0, -1, 0, 0 };
int dy[] = { 0, 1, 0, 0, -1, 0 };
int dz[] = { 0, 0, 1, 0, 0, -1 };

class PairwiseBP
{
public:
	// Input Var
	int Layer, Height, Width, nCategory;
	Array3D_vect* pDataTerm;
	Array3D_float* pSmoothTerm[3];

	// State Var
	Array3D_vect* pMessage[6];
	Array3D_vect* pBelief;

	// Output Var
	Array3D_int* result;

	// Functions

	void Add2Message(T_input *m1, T_input *m2, int len)
	{
		for (int i = 0; i<len; i++)
			m1[i] += m2[i];

		return;
	}

	void CopyMessage(T_input *m1, T_input *m2, int len)
	{
		memcpy(m1, m2, sizeof(T_input) * len);

		return;
	}

	T_input GetMinValue(T_input *m, int len)
	{
		T_input minv = FLT_MAX;

		for (int i = 0; i < len; i++)
			minv = min(minv, m[i]);

		return minv;
	}

	int GetMinIndex(T_input *m, int len)
	{
		int mink = 0;

		for (int i = 1; i < len; i++)
			if (m[mink] > m[i])
				mink = i;

		return mink;
	}

	void ModifyByMinThreshold(T_input *m, T_input threshold, int len)
	{
		for (int i = 0; i < len; i++)
			m[i] = min(m[i], threshold);

		return;
	}

	void MinusValue(T_input *m, T_input v, int len)
	{
		for (int i = 0; i < len; i++)
			m[i] -= v;

		return;
	}

	void AllocateBuffer(int _Layer, int _Height, int _Width, int _nCategory)
	{
		Layer = _Layer, Height = _Height, Width = _Width, nCategory = _nCategory;
		pDataTerm = new Array3D_vect(Layer, Height, Width, nCategory);
		for (int i = 0; i < 3; i++)
			pSmoothTerm[i] = new Array3D_float(Layer, Height, Width);
		for (int i = 0; i < 6; i++)
			pMessage[i] = new Array3D_vect(Layer, Height, Width, nCategory);
		pBelief = new Array3D_vect(Layer, Height, Width, nCategory);
		result = new Array3D_int(Layer, Height, Width);

		return;
	}

	void ReleaseBuffer()
	{
		delete pDataTerm;
		for (int i = 0; i < 3; i++)
			delete pSmoothTerm[i];
		for (int i = 0; i < 6; i++)
			delete pMessage[i];
		delete pBelief;
	}

	void LoadDataTerm(Array3D_vect* input_pDataTerm)
	{
		AllocateBuffer(input_pDataTerm->Layer, input_pDataTerm->Height, input_pDataTerm->Width, input_pDataTerm->nDim);
		memcpy(pDataTerm->data, input_pDataTerm->data, sizeof(T_input) * input_pDataTerm->nElements);

		return;
	}

	void LoadSmoothness(Array3D_float* input_pSmoothTerm[3])
	{
		for (int d = 0; d < 3; d++)
			memcpy(pSmoothTerm[d]->data, input_pSmoothTerm[d]->data, sizeof(T_input) * input_pSmoothTerm[d]->Volume);

		return;
	}

	PairwiseBP()
	{

	}

	PairwiseBP(int _Layer, int _Height, int _Width, int _nCategory)
	{
		AllocateBuffer(_Layer, _Height, _Width, _nCategory);
	}

	PairwiseBP(PairwiseBP *oriBP) // generateCoarserLevel
	{
		AllocateBuffer((oriBP->Layer + 1) / 2, (oriBP->Height + 1) / 2, (oriBP->Width + 1) / 2, oriBP->nCategory);
		for (int c = 0; c < nCategory; c++)
		for (int i = 0; i < Layer; i++)
		for (int j = 0; j < Height; j++)
		for (int k = 0; k < Width; k++)
		{
			pDataTerm->at(i, j, k)[c] = 0;
			for (int di = 0; di < 2; di++)
			for (int dj = 0; dj < 2; dj++)
			for (int dk = 0; dk < 2; dk++)
			{
				int ni = 2 * i + di;
				int nj = 2 * j + dj;
				int nk = 2 * k + dk;
				if (oriBP->validPosition(ni, nj, nk))
					pDataTerm->at(i, j, k)[c] += oriBP->pDataTerm->at(ni, nj, nk)[c];
			}
		}

		// I think we can merge Smooth in a better way.
		for (int dir = 0; dir < 3; dir++)
		for (int i = 0; i < Layer; i++)
		for (int j = 0; j < Height; j++)
		for (int k = 0; k < Width; k++)
		{
			pSmoothTerm[dir]->at(i, j, k) = 0;
			for (int di = 0; di < 2; di++)
			for (int dj = 0; dj < 2; dj++)
			for (int dk = 0; dk < 2; dk++)
			{
				int ni = 2 * i + di;
				int nj = 2 * j + dj;
				int nk = 2 * k + dk;
				if (oriBP->validPosition(ni, nj, nk))
					pSmoothTerm[dir]->at(i, j, k) += oriBP->pSmoothTerm[dir]->at(ni, nj, nk);
			}
		}
	}

	void propagateFinerLevel(PairwiseBP *toBP) // from me (coarse) to a finer
	{
		for (int i = 0; i < toBP->Layer; i++)
		for (int j = 0; j < toBP->Height; j++)
		for (int k = 0; k < toBP->Width; k++)
		{
			int ni = i / 2;
			int nj = j / 2;
			int nk = k / 2;
			for (int dir = 0; dir < 6; dir++)
			for (int c = 0; c < nCategory; c++)
				toBP->pMessage[dir]->at(i, j, k)[c] = pMessage[dir]->at(ni, nj, nk)[c];
		}
	}

	~PairwiseBP()
	{
		ReleaseBuffer();
	}

	bool validPosition(int x, int y, int z)
	{
		if (0 <= x && x < Layer)
		if (0 <= y && y < Height)
		if (0 <= z && z < Width)
			return true;
		return false;
	}

	void ComputeBelief()
	{
		T_input *message = new T_input[nCategory];

		for (int i = 0; i < Layer; i++)
		for (int j = 0; j < Height; j++)
		for (int k = 0; k < Width; k++)
		{
			CopyMessage(message, pDataTerm->at(i, j, k), nCategory);

			for (int dir = 0; dir < 6; dir++)
			{
				int ni = i - dx[dir];
				int nj = j - dy[dir];
				int nk = k - dz[dir];
				if (validPosition(ni, nj, nk))
					Add2Message(message, pMessage[dir]->at(i, j, k), nCategory);
			}
			CopyMessage(pBelief->at(i, j, k), message, nCategory);
			result->at(i, j, k) = GetMinIndex(message, nCategory);
		}

		delete[] message;
	}

	int opposite(int dir){ return (dir + 3) % 6; }

	void UpdateMessage(int x, int y, int z, int dir)
	{
		int nx = x + dx[dir];
		int ny = y + dy[dir];
		int nz = z + dz[dir];

		if (!validPosition(nx, ny, nz))
			return;

		T_input d;
		if (dir < 3)
			d = pSmoothTerm[dir]->at(x, y, z);
		else
			d = pSmoothTerm[dir - 3]->at(nx, ny, nz);

		T_input *message = new T_input[nCategory];
		CopyMessage(message, pDataTerm->at(x, y, z), nCategory);

		for (int fromDir = 0; fromDir < 6; fromDir++)
		{
			if (fromDir == opposite(dir))
				continue;
			int orix = x - dx[fromDir];
			int oriy = y - dy[fromDir];
			int oriz = z - dz[fromDir];
			if (!validPosition(orix, oriy, oriz))
				continue;
			Add2Message(message, pMessage[fromDir]->at(x, y, z), nCategory);
		}

		ModifyByMinThreshold(message, GetMinValue(message, nCategory) + d, nCategory);
		MinusValue(message, GetMinValue(message, nCategory), nCategory);
		CopyMessage(pMessage[dir]->at(nx, ny, nz), message, nCategory);

		delete[] message;
	}

	void BP_S(int count)
	{
		if (count % 2 == 0)
		{
			for (int i = 0; i < Layer; i++)
			for (int j = 0; j < Height; j++)
			for (int k = 0; k < Width; k++)
			for (int dir = 0; dir < 3; dir++)
				UpdateMessage(i, j, k, dir);
		}
		else
		{
			for (int i = Layer - 1; i >= 0; i--)
			for (int j = Height - 1; j >= 0; j--)
			for (int k = Width - 1; k >= 0; k--)
			for (int dir = 0; dir < 3; dir++)
				UpdateMessage(i, j, k, dir + 3);
		}
	}

	void MessagePassing(int nIterations, int nHierarchy, bool isTopLevel = false)
	{
		if (nHierarchy > 0)
		{
			PairwiseBP son(this);
			son.MessagePassing(6, nHierarchy - 1);
			son.propagateFinerLevel(this);

			delete son.result;
		}

		for (int count = 0; count < nIterations; count++)
		{
			BP_S(count);
		}

		if (isTopLevel)
			ComputeBelief();

		return;
	}

};
