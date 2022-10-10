#pragma once
#include <cstring>
#define T_input float

class Array3D_int
{
public:
    int* data;
    int Layer, Height, Width, Volume;
    Array3D_int(int _Layer, int _Height, int _Width)
    {
        Layer = _Layer;
        Height = _Height;
        Width = _Width;
        Volume = Layer * Height * Width;
        data = new int[Volume];
        memset(data, 0, sizeof(int)* Volume);
    }
    ~Array3D_int()
    {
        delete []data;
    }
    int& at(int l, int h, int w)
    {
        int offset = 0;
        offset += l * Height * Width;
        offset += h * Width;
        offset += w;
        return data[offset];
    }
    int& at(int index)
    {
        return data[index];
    }
    void setVal(int val)
    {
        for (int i = 0; i < Volume; ++i)
            data[i] = val;
    }
};


class Array3D_float
{
public:
    T_input* data;
    int Layer, Height, Width, Volume;
    Array3D_float(int _Layer, int _Height, int _Width)
    {
        Layer = _Layer;
        Height = _Height;
        Width = _Width;
        Volume = Layer * Height * Width;
        data = new T_input[Volume];
        memset(data, 0, sizeof(T_input)* Volume);
    }
    ~Array3D_float()
    {
        delete[] data;
    }
    T_input& at(int l, int h, int w)
    {
        int offset = 0;
        offset += l * Height * Width;
        offset += h * Width;
        offset += w;
        return data[offset];
    }
    T_input& at(int index)
    {
        return data[index];
    }
};

class Array3D_vect
{
public:
    T_input* data;
    int Layer, Height, Width, nDim, Volume;
    size_t nElements;
    Array3D_vect(Array3D_vect *old)
    {
        Layer = old->Layer;
        Height = old->Height;
        Width = old->Width;
        nDim = old->nDim;
        Volume = old->Volume;
        nElements = old->nElements;
        data = new T_input[nElements];
        memcpy(data, old->data, sizeof(T_input) * old->nElements);
    }
    Array3D_vect(int _Layer, int _Height, int _Width, int _nDim)
    {
        Layer = _Layer;
        Height = _Height;
        Width = _Width;
        nDim = _nDim;
        Volume = Layer * Height * Width;
        nElements = (size_t)Volume * nDim;
        data = new T_input[nElements];
        memset(data, 0, sizeof(T_input)* nElements);
    }
    ~Array3D_vect()
    {
        delete[] data;
    }
    T_input& at(size_t index)
    {
        return data[index];
    }
    T_input& at(int index, int dimindex)
    {
        size_t offset = (size_t) index * nDim + dimindex;

        return data[offset];
    }
    T_input* at(int l, int h, int w)
    {
        size_t offset = (size_t)l * Height * Width * nDim + (size_t)h * Width * nDim + (size_t)w * nDim;

        return &(data[offset]);
    }
    T_input& at(int l, int h, int w, int dimindex)
    {
        return at(l, h, w)[dimindex];
    }
    void setValue(T_input d)
    {
        for (size_t i = 0; i < nElements; i++)
            data[i] = d;
        return;
    }
};
