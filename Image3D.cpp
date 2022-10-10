#include "Image3D.h"

int NumScales = 1;
const int OriSiftDim = 48; // 48
const int SiftDim = OriSiftDim * NumScales + 1;
T_input GrayAvgBaseLine = 4.0; // 4.0 ADNI: 55.0
T_input SIFTAvgBaseLine = 8.0; // 8.0 ADNI: 1.0

T_input sqr(T_input x){ return x*x; }
int sqr(int x){ return x*x; }

int enforceRange(int x, int minval, int maxval)
{
    return max(minval, min(x, maxval));
}

int FloatComp(const void *a, const void *b)
{
    T_input fa = *(T_input *)a;
    T_input fb = *(T_input *)b;

    if (fa > fb) return 1;
    else if (fa < fb) return -1;
    else return 0;
}


Image3DWithSeg::Image3DWithSeg(string dataFileName, string _segFileName) : Image3D(dataFileName)
{
    seg = NULL;
    segFileName = _segFileName;

    if (segFileName != "")
    {
        ParseHDR(segFileName, segHDR);
        seg = new Array3D_int(Layer, Height, Width);
        cout << "reading " << segFileName << '\n';
        ReadIntImg(seg, segFileName, segHDR.dime.datatype);
    }

}

Image3DWithSeg::Image3DWithSeg()
{
}

Image3DWithSeg::~Image3DWithSeg()
{
    delete seg;
}

Image3D::Image3D(string _dataFileName)
{
    data = NULL;
    sift = NULL;
    mask = NULL;

    dataFileName = _dataFileName;

    ParseHDR(dataFileName, dataHDR);
    data = new Array3D_float(Layer, Height, Width);
    ReadFloatImg(data, dataFileName, dataHDR.dime.datatype);

    sift = new Array3D_vect(Layer, Height, Width, SiftDim);
    
    int *cellSizeArr = new int[NumScales];
    for (int i = 0; i < NumScales; ++i)
        cellSizeArr[i] = i + 1;

    for (int i = 0; i < NumScales; ++i)
        ComputeSIFT(data, sift, i * OriSiftDim, cellSizeArr[i]);

    for (int i = 0; i < sift->Volume; ++i)
        sift->at(i, SiftDim - 1) = data->at(i);

    ChangeSIFTByCof(SIFTAvgBaseLine);
    ChangeGrayByCof(GrayAvgBaseLine);

    mask = new Array3D_int(Layer, Height, Width);
    mask->setVal(1);

    delete[] cellSizeArr;
}

Image3D::Image3D()
{
}

Image3D::~Image3D()
{
    delete data;
    delete sift;
    delete mask;
}


void Image3D::ThreeDimensionalFiltering(Array3D_vect *imsrc, Array3D_vect *imdst, int fSize)
{
    int Layer = imsrc->Layer;
    int Height = imsrc->Height;
    int Width = imsrc->Width;
    int nChannels = imsrc->nDim;

    T_input *filter = new T_input[2 * fSize + 1];

    filter[0] = filter[fSize + 1] = 0.25;
    for (int i = 1; i < fSize + 1; ++i)
        filter[i] = 1;
    for (int i = fSize + 2; i < fSize * 2 + 1; ++i)
        filter[i] = 0;

    for (int di = -fSize; di <= fSize; ++di)
    for (int dj = -fSize; dj <= fSize; ++dj)
    for (int dk = -fSize; dk <= fSize; ++dk)
    {
        T_input cof = filter[di + fSize] * filter[dj + fSize] * filter[dk + fSize];

        if (cof == 0)
            continue;

        for (int i = 0; i < Layer; ++i)
        for (int j = 0; j < Height; ++j)
        for (int k = 0; k < Width; ++k)
        for (int c = 0; c < nChannels; ++c)
            imdst->at(i, j, k, c) += imsrc->at(enforceRange(i + di, 0, Layer - 1), enforceRange(j + dj, 0, Height - 1), enforceRange(k + dk, 0, Width - 1), c) * cof;
    }

    delete[] filter;

    return;
}

void Image3D::GetDirections(T_input dirs[][3], int nBins)
{
    T_input dir6[6][3] = { { 0, 0, 1 }, { 0, 1, 0 }, { 1, 0, 0 }, { 0, 0, -1 }, { 0, -1, 0 }, { -1, 0, 0 } };

    if (nBins == 6)
    {
        for (int i = 0; i < nBins; ++i)
        for (int j = 0; j < 3; ++j)
            dirs[i][j] = dir6[i][j];
    }
    else if (nBins == 8)
    {
        for (int i = -1, cnt = 0; i <= 1; i += 2)
        for (int j = -1; j <= 1; j += 2)
        for (int k = -1; k <= 1; k += 2)
            dirs[cnt][0] = i, dirs[cnt][1] = j, dirs[cnt][2] = k, ++cnt;
    }
    else if (nBins == 26)
    {
        for (int i = -1, cnt = 0; i <= 1; ++i)
        for (int j = -1; j <= 1; ++j)
        for (int k = -1; k <= 1; ++k)
            if (i != 0 || j != 0 || k != 0)
                dirs[cnt][0] = i, dirs[cnt][1] = j, dirs[cnt][2] = k, ++cnt;
    }


    return;
}

void Image3D::ComputeDerivatives(Array3D_float *imsrc, Array3D_float *dx, Array3D_float *dy, Array3D_float *dz, bool isAdvancedFilter)
{
    if (!isAdvancedFilter)
    {
        for (int i = 0; i < Layer; ++i)
        for (int j = 0; j < Height; ++j)
        for (int k = 0; k < Width; ++k)
            dx->at(i, j, k) = imsrc->at(enforceRange(i + 1, 0, Layer - 1), j, k) - imsrc->at(i, j, k),
            dy->at(i, j, k) = imsrc->at(i, enforceRange(j + 1, 0, Height - 1), k) - imsrc->at(i, j, k),
            dz->at(i, j, k) = imsrc->at(i, j, enforceRange(k + 1, 0, Width - 1)) - imsrc->at(i, j, k);
    }
    else
    {
        const int fSize = 2;
        T_input filter[fSize * 2 + 1] = { 1, -8, 0, 8, -1 };

        for (int i = 0; i < fSize * 2 + 1; ++i)
            filter[i] /= 12;

        for (int i = 0; i < Layer; ++i)
        for (int j = 0; j < Height; ++j)
        for (int k = 0; k < Width; ++k)
        for (int d = -fSize; d <= fSize; ++d)
            dx->at(i, j, k) += imsrc->at(enforceRange(i + d, 0, Layer - 1), j, k) * filter[d + fSize],
            dy->at(i, j, k) += imsrc->at(i, enforceRange(j + d, 0, Height - 1), k) * filter[d + fSize],
            dz->at(i, j, k) += imsrc->at(i, j, enforceRange(k + d, 0, Width - 1)) * filter[d + fSize];

        return;
    }

    return;
}

void Image3D::ComputeSIFT(Array3D_float *imsrc, Array3D_vect *imsift, int arrOffset, int cellSize, int stepSize, bool isBoundaryIncluded, int nBins, int alpha, int winSize)
{
    int nPixels = Layer * Height * Width;
    int siftdim = nBins * (int)pow(2 * winSize, 3);
    int sift_layer = Layer / stepSize;
    int sift_height = Height / stepSize;
    int sift_width = Width / stepSize;
    int x_shift = 0, y_shift = 0, z_shift = 0;
    T_input dirs[26][3];
    T_input smooth = 0.000001;

    Array3D_float *imdx = new Array3D_float(Layer, Height, Width);
    Array3D_float *imdy = new Array3D_float(Layer, Height, Width);
    Array3D_float *imdz = new Array3D_float(Layer, Height, Width);
    Array3D_float *mag = new Array3D_float(Layer, Height, Width);
    Array3D_vect *gradient = new Array3D_vect(Layer, Height, Width, 3);
    Array3D_vect *imband = new Array3D_vect(Layer, Height, Width, nBins);
    Array3D_vect *imband_cell = new Array3D_vect(Layer, Height, Width, nBins);

    if (!isBoundaryIncluded)
    {
        sift_layer = (Layer - 4 * cellSize) / stepSize;
        sift_height = (Height - 4 * cellSize) / stepSize;
        sift_width = (Width - 4 * cellSize) / stepSize;
        x_shift = y_shift = z_shift = 2 * cellSize;
    }

    ComputeDerivatives(data, imdx, imdy, imdz, false);
    GetDirections(dirs, nBins);

    for (int i = 0; i < nPixels; ++i)
    {
        mag->at(i) = sqrt(sqr(imdx->at(i)) + sqr(imdy->at(i)) + sqr(imdz->at(i)));
        gradient->at(i * 3 + 0) = imdx->at(i) / (mag->at(i) + smooth);
        gradient->at(i * 3 + 1) = imdy->at(i) / (mag->at(i) + smooth);
        gradient->at(i * 3 + 2) = imdz->at(i) / (mag->at(i) + smooth);
    }

    for (int i = 0; i < nPixels; ++i)
    for (int k = 0; k < nBins; ++k)
    {
        T_input temp;

        temp = max(gradient->at(i * 3 + 0) * dirs[k][0] + gradient->at(i * 3 + 1) * dirs[k][1] + gradient->at(i * 3 + 2) * dirs[k][2], (T_input)0);
        temp = pow(temp, alpha);
        imband->at(i, k) = temp * mag->at(i);
    }

    ThreeDimensionalFiltering(imband, imband_cell, cellSize);

    for (int i = 0; i < sift_layer; ++i)
    for (int j = 0; j < sift_height; ++j)
    for (int k = 0; k < sift_width; ++k)
    {
        int count = 0;
        T_input mag = 0;

        for (int ii = -winSize; ii < winSize; ++ii)
        for (int jj = -winSize; jj < winSize; ++jj)
        for (int kk = -winSize; kk < winSize; ++kk)
        {
            int x = enforceRange(x_shift + i * stepSize + ii * cellSize, 0, Layer - 1);
            int y = enforceRange(y_shift + j * stepSize + jj * cellSize, 0, Height - 1);
            int z = enforceRange(z_shift + k * stepSize + kk * cellSize, 0, Width - 1);

            memcpy(imsift->at(i, j, k) + arrOffset + count * nBins, imband_cell->at(x, y, z), sizeof(T_input) * nBins);
            ++count;
        }
        for (int x = 0; x < siftdim; ++x)
            mag += sqr(imsift->at(i, j, k, arrOffset + x));
        mag = sqrt(mag);
        for (int x = 0; x < siftdim; ++x)
            imsift->at(i, j, k, arrOffset + x) /= (mag + smooth);
    }

    delete imdx;
    delete imdy;
    delete imdz;
    delete mag;
    delete gradient;
    delete imband;
    delete imband_cell;

    return;
}

bool Image3D::InBound(int x, int y, int z)
{
    if (0 <= x && x < Layer && 0 <= y && y < Height && 0 <= z && z < Width)
        return true;
    else
        return false;
}

void Image3D::ChangeGrayByCof(T_input grayAvgBaseLine)
{
    T_input graySum = 0;
    T_input cof;
    
    for (int i = 0; i < Layer; i++)
    for (int j = 0; j < Height; j++)
    for (int k = 0; k < Width; k++)
        graySum += sift->at(i, j, k)[SiftDim - 1];

    cout << dataFileName << " Gray Average: " << (graySum / (Layer * Height * Width)) << endl;

    cof = Layer * Height * Width * grayAvgBaseLine / graySum;

    for (int i = 0; i < Layer; i++)
    for (int j = 0; j < Height; j++)
    for (int k = 0; k < Width; k++)
    {
        sift->at(i, j, k)[SiftDim - 1] *= cof;
        data->at(i, j, k) *= cof;
    }

    return;
}

void Image3D::ChangeSIFTByCof(T_input siftAvgBaseLine)
{
    T_input siftSum = 0;
    T_input cof;

    for (int i = 0; i < Layer; i++)
    for (int j = 0; j < Height; j++)
    for (int k = 0; k < Width; k++)
    for (int id = 0; id < SiftDim - 1; id++)
        siftSum += sift->at(i, j, k)[id];

    cout << dataFileName << " SIFT Average: " << (siftSum / (Layer * Height * Width)) << endl;

    cof = Layer * Height * Width * siftAvgBaseLine / siftSum;

    for (int i = 0; i < Layer; i++)
    for (int j = 0; j < Height; j++)
    for (int k = 0; k < Width; k++)
    for (int id = 0; id < SiftDim - 1; id++)
        sift->at(i, j, k)[id] *= cof;

    return;
}

void Image3D::ParseHDR(string fileName, struct dsr& hdr)
{
    ifstream fp(fileName + ".hdr", ios::binary);
    if (!fp)
        exit(1);
    fp.read((char *)&hdr, sizeof(struct dsr));
    fp.close();
    
    Width = hdr.dime.dim[1];
    Height = hdr.dime.dim[2];
    Layer = hdr.dime.dim[3];

    return;
}

void Image3D::ReadIntImg(Array3D_int* target, string fileName, int dataType)
{
    int N = target->Volume;

    ifstream fin((fileName + ".img").c_str(), ios::binary);

    if (dataType == DT_DOUBLE)
    {
        double *tem = new double[N];
        fin.read((char *)tem, sizeof(double)* N);
        for (int i = 0; i < N; ++i)
            target->data[i] = (int)tem[i];
        delete[] tem;
    }
    else if (dataType == DT_FLOAT)
    {
        float *tem = new float[N];
        fin.read((char *)tem, sizeof(float)* N);
        for (int i = 0; i < N; ++i)
            target->data[i] = (int)tem[i];
        delete[] tem;
    }
    else if (dataType == DT_SIGNED_INT)
    {
        int* tem = new int[N];
        fin.read((char*)tem, sizeof(int) * N);
        for (int i = 0; i < N; ++i)
            target->data[i] = (int)tem[i];
        delete[] tem;
    }
    else if (dataType == DT_SIGNED_SHORT)
    {
        short *tem = new short[N];
        fin.read((char *)tem, sizeof(short)* N);
        for (int i = 0; i < N; ++i)
            target->data[i] = (int)tem[i];
        delete[] tem;
    }
    else if (dataType == DT_UNSIGNED_CHAR)
    {
        unsigned char *tem = new unsigned char[N];
        fin.read((char *)tem, sizeof(unsigned char)* N);
        for (int i = 0; i < N; ++i)
            target->data[i] = (int)tem[i];
        delete[] tem;
    }
    else
    {
        fprintf(stderr, "0 length not support!\n");
        exit(1);
    }

    fin.close();

    return;
}

void Image3D::ReadFloatImg(Array3D_float* target, string fileName, int dataType)
{
    int N = target->Volume;

    ifstream fin((fileName + ".img").c_str(), ios::binary);

    if (dataType == DT_DOUBLE)
    {
        double* tem = new double[N];
        fin.read((char*)tem, sizeof(double) * N);
        for (int i = 0; i < N; ++i)
            target->data[i] = (T_input)tem[i];
        delete[] tem;
    }
    else if (dataType == DT_FLOAT)
    {
        float* tem = new float[N];
        fin.read((char*)tem, sizeof(float) * N);
        for (int i = 0; i < N; ++i)
            target->data[i] = (T_input)tem[i];
        delete[] tem;
    }
    else if (dataType == DT_SIGNED_INT)
    {
        int* tem = new int[N];
        fin.read((char*)tem, sizeof(int) * N);
        for (int i = 0; i < N; ++i)
            target->data[i] = (T_input)tem[i];
        delete[] tem;
    }
    else if (dataType == DT_SIGNED_SHORT)
    {
        short* tem = new short[N];
        fin.read((char*)tem, sizeof(short) * N);
        for (int i = 0; i < N; ++i)
            target->data[i] = (T_input)tem[i];
        delete[] tem;
    }
    else if (dataType == DT_UNSIGNED_CHAR)
    {
        unsigned char* tem = new unsigned char[N];
        fin.read((char*)tem, sizeof(unsigned char) * N);
        for (int i = 0; i < N; ++i)
            target->data[i] = (T_input)tem[i];
        delete[] tem;
    }
    else
    {
        fprintf(stderr, "1 length not support!\n");
        exit(1);
    }

    fin.close();

    return;
}

void Image3D::WriteHDR(string fileName, struct dsr& hdr)
{
    ofstream fout((fileName + ".hdr").c_str(), ios::binary);
    fout.write((char *)&hdr, sizeof(struct dsr));
    fout.close();

    return;
}

void Image3D::WriteIMG(string fileName, Array3D_float& data, int dataType)
{
    int N = data.Volume;
    ofstream fout((fileName + ".img").c_str(), ios::binary);

    if (dataType == DT_DOUBLE)
    {
        double *tem = new double[N];
        for (int i = 0; i < N; ++i)
            tem[i] = (double)data.data[i];
        fout.write((char *)tem, sizeof(double)* N);
        delete[] tem;
    }
    else if (dataType == DT_FLOAT)
    {
        float *tem = new float[N];
        for (int i = 0; i < N; ++i)
            tem[i] = (float)data.data[i];
        fout.write((char *)tem, sizeof(float)* N);
        delete[] tem;
    }
    else if (dataType == DT_SIGNED_INT)
    {
        int* tem = new int[N];
        for (int i = 0; i < N; ++i)
            tem[i] = (int)data.data[i];
        fout.write((char*)tem, sizeof(int) * N);
        delete[] tem;
    }
    else if (dataType == DT_SIGNED_SHORT)
    {
        short *tem = new short[N];
        for (int i = 0; i < N; ++i)
            tem[i] = (short)data.data[i];
        fout.write((char *)tem, sizeof(short)* N);
        delete[] tem;
    }
    else if (dataType == DT_UNSIGNED_CHAR)
    {
        unsigned char *tem = new unsigned char[N];
        for (int i = 0; i < N; ++i)
            tem[i] = (unsigned char)data.data[i];
        fout.write((char *)tem, sizeof(unsigned char)* N);
        delete[] tem;
    }
    else
    {
        fprintf(stderr, "2 length not support!\n");
        exit(1);
    }

    fout.close();

    return;
}

void Image3D::WriteIMG(string fileName, Array3D_int& data, int dataType)
{
    int N = data.Volume;
    ofstream fout((fileName + ".img").c_str(), ios::binary);

    if (dataType == DT_DOUBLE)
    {
        double* tem = new double[N];
        for (int i = 0; i < N; ++i)
            tem[i] = (double)data.data[i];
        fout.write((char*)tem, sizeof(double) * N);
        delete[] tem;
    }
    else if (dataType == DT_FLOAT)
    {
        float* tem = new float[N];
        for (int i = 0; i < N; ++i)
            tem[i] = (float)data.data[i];
        fout.write((char*)tem, sizeof(float) * N);
        delete[] tem;
    }
    else if (dataType == DT_SIGNED_INT)
    {
        int* tem = new int[N];
        for (int i = 0; i < N; ++i)
            tem[i] = (int)data.data[i];
        fout.write((char*)tem, sizeof(int) * N);
        delete[] tem;
    }
    else if (dataType == DT_SIGNED_SHORT)
    {
        short* tem = new short[N];
        for (int i = 0; i < N; ++i)
            tem[i] = (short)data.data[i];
        fout.write((char*)tem, sizeof(short) * N);
        delete[] tem;
    }
    else if (dataType == DT_UNSIGNED_CHAR)
    {
        unsigned char* tem = new unsigned char[N];
        for (int i = 0; i < N; ++i)
            tem[i] = (unsigned char)data.data[i];
        fout.write((char*)tem, sizeof(unsigned char) * N);
        delete[] tem;
    }
    else
    {
        fprintf(stderr, "3 length not support!\n");
        exit(1);
    }

    fout.close();

    return;
}

void Image3D::OutputDataImage(string fileName)
{
    if (data == NULL)
        return;

    WriteHDR(fileName, dataHDR);
    WriteIMG(fileName, *data, dataHDR.dime.datatype);

    return;
}


void Image3DWithSeg::OutputSegImage(string fileName)
{
    if (seg == NULL)
        return;

    WriteHDR(fileName, segHDR);
    WriteIMG(fileName, *seg, segHDR.dime.datatype);

    return;
}


void Image3DWithSeg::OutputDataSeg(string fileName0, string fileName1)
{
    OutputDataImage(fileName0);
    OutputSegImage(fileName1);

    return;
}

void Image3DWithSeg::ApplyFlow(string fileName, Image3DWithSeg *refImg)
{
    ifstream fin(fileName.c_str(), ios::binary);

    if (!fin)
        return;

    int newLayer = refImg->Layer;
    int newHeight = refImg->Height;
    int newWidth = refImg->Width;
    int N = newLayer * newHeight * newWidth;

    Array3D_float* newdata = new Array3D_float(newLayer, newHeight, newWidth);
    Array3D_int* newseg = new Array3D_int(newLayer, newHeight, newWidth);
    Array3D_vect* newsift = new Array3D_vect(newLayer, newHeight, newWidth, SiftDim);
    Array3D_int* newmask = new Array3D_int(newLayer, newHeight, newWidth);
    Array3D_int* vvx = new Array3D_int(newLayer, newHeight, newWidth);
    Array3D_int* vvy = new Array3D_int(newLayer, newHeight, newWidth);
    Array3D_int* vvz = new Array3D_int(newLayer, newHeight, newWidth);

    fin.read((char *)vvx->data, sizeof(int)* N);
    fin.read((char *)vvy->data, sizeof(int)* N);
    fin.read((char *)vvz->data, sizeof(int)* N);

    for (int i = 0; i < newLayer; i++)
    for (int j = 0; j < newHeight; j++)
    for (int k = 0; k < newWidth; k++)
    {
        int vx = vvx->at(i, j, k);
        int vy = vvy->at(i, j, k);
        int vz = vvz->at(i, j, k);
        int I = i + vz;
        int J = j + vy;
        int K = k + vx;

        if (I < 0 || J < 0 || K < 0 || I >= Layer || J >= Height || K >= Width)
            newmask->at(i, j, k) = 0;
        else
            newmask->at(i, j, k) = 1;

        if (I < 0) I = 0;
        if (J < 0) J = 0;
        if (K < 0) K = 0;
        if (I >= Layer) I = Layer - 1;
        if (J >= Height) J = Height - 1;
        if (K >= Width) K = Width - 1;

        newdata->at(i, j, k) = data->at(I, J, K);
        newseg->at(i, j, k) = seg->at(I, J, K);
        for (int u = 0; u < SiftDim; u++)
            newsift->at(i, j, k)[u] = sift->at(I, J, K)[u];
    }

    fin.close();
    
    delete vvx;
    delete vvy;
    delete vvz;

    delete data;
    delete seg;
    delete sift;

    data = newdata;
    seg = newseg;
    sift = newsift;
    mask = newmask;

    Layer = newLayer;
    Height = newHeight;
    Width = newWidth;

    dataHDR = refImg->dataHDR;
    segHDR = refImg->dataHDR;
    dataHDR.dime.bitpix = 32;
    dataHDR.dime.datatype = 16;
    segHDR.dime.bitpix = 8;
    segHDR.dime.datatype = 2;

    return;
}

void Image3DWithSeg::PrintFlowDetail(string flowFilePath, string outTxtPath)
{
    ifstream fin(flowFilePath.c_str(), ios::binary);
    ofstream fout(outTxtPath.c_str());

    if (!fin)
        return;

    int Layer = data->Layer;
    int Height = data->Height;
    int Width = data->Width;
    int N = Layer * Height * Width;

    Array3D_float* newdata = new Array3D_float(Layer, Height, Width);
    Array3D_int* newseg = new Array3D_int(Layer, Height, Width);
    Array3D_int* vvx = new Array3D_int(Layer, Height, Width);
    Array3D_int* vvy = new Array3D_int(Layer, Height, Width);
    Array3D_int* vvz = new Array3D_int(Layer, Height, Width);

    fin.read((char *)vvx->data, sizeof(int)* N);
    fin.read((char *)vvy->data, sizeof(int)* N);
    fin.read((char *)vvz->data, sizeof(int)* N);

    for (int i = 0; i < Layer; i++)
    for (int j = 0; j < Height; j++)
    for (int k = 0; k < Width; k++)
    {
        int vx, vy, vz;

        vx = vvx->at(i, j, k);
        vy = vvy->at(i, j, k);
        vz = vvz->at(i, j, k);


        int I = i + vz;
        int J = j + vy;
        int K = k + vx;


        if (I < 0) I = 0;
        if (J < 0) J = 0;
        if (K < 0) K = 0;
        if (I >= Layer) I = Layer - 1;
        if (J >= Height) J = Height - 1;
        if (K >= Width) K = Width - 1;

        fout << k + 1 << "," << j + 1 << "," << i + 1 << " -> " << K + 1 << "," << J + 1 << "," << I + 1 << "   " << data->at(I, J, K) << " " << seg->at(I, J, K) << endl;

        newdata->at(i, j, k) = data->at(I, J, K);
        newseg->at(i, j, k) = seg->at(I, J, K);
    }


    fin.close();
    fout.close();

    delete vvx;
    delete vvy;
    delete vvz;

    delete data;
    delete seg;
    delete sift;

    data = newdata;
    seg = newseg;

    return;
}


void Image3DWithSeg::OutputDeformationAmount(string flowFilePath, string outTxtPath)
{
    ifstream fin(flowFilePath.c_str(), ios::binary);
    ofstream fout(outTxtPath.c_str());

    if (!fin)
        return;

    int N = Layer * Height * Width;
    double sum = 0;

    Array3D_int* vvx = new Array3D_int(Layer, Height, Width);
    Array3D_int* vvy = new Array3D_int(Layer, Height, Width);
    Array3D_int* vvz = new Array3D_int(Layer, Height, Width);

    fin.read((char *)vvx->data, sizeof(int)* N);
    fin.read((char *)vvy->data, sizeof(int)* N);
    fin.read((char *)vvz->data, sizeof(int)* N);

    for (int i = 0; i < Layer; i++)
    for (int j = 0; j < Height; j++)
    for (int k = 0; k < Width; k++)
    {
        int vx, vy, vz;

        vx = vvx->at(i, j, k);
        vy = vvy->at(i, j, k);
        vz = vvz->at(i, j, k);

        sum += sqrt(sqr(vx) + sqr(vy) + sqr(vz));
    }

    fout << sum << endl;

    fin.close();
    fout.close();

    delete vvx;
    delete vvy;
    delete vvz;

    return;
}
