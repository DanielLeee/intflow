#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <map>
#include <cstring>
#include "Image3D.h"
#include "PairwiseBP.h"
#include <algorithm>
using namespace std;

int Layer, Height, Width;

// algorithm parameters
static T_input epsilon = 0.333; // 0.333
static T_input beta = 0.01; // 0.01
static T_input prioreps = 0.2; // 0.2

// image parameters
static int nCategory = 208;
static T_input worst = 100;

int nTraining;
bool noFlow = false;
vector <string> TrainingImgNameList;
vector <string> TrainingSegNameList;
vector <string> FlowNameList;
string TestImgName;
Image3DWithSeg* TrainingImg;
Image3DWithSeg* TestImg;

Array3D_vect* pDataTerm; // From AfterWarp
Array3D_float* pSmoothTerm[3]; // From Test[i]
Array3D_vect* objPrior; // From Train[1..n]
Array3D_vect* graydiff;
Array3D_vect* siftdiff;

T_input *objPriorMax;
T_input *satuPara;

extern T_input sqr(T_input x);
extern int sqr(int x);
extern int enforceRange(int x, int minval, int maxval);

bool atBoard(int x, int y, int z)
{
    if(x == 0 || y == 0 || z == 0)
        return true;
    if(x == Layer-1 || y == Height-1 || z == Width-1)
        return true;
    return false;
}

void readTrainingNameList(char** argv, int argc, int start)
{
    int i;
    for (i = start; i < argc && argv[i][0] != '-'; i += 2)
    {
        string st;
        st = argv[i];
        TrainingImgNameList.push_back(st);
        st = argv[i + 1];
        TrainingSegNameList.push_back(st);
    }
        
    nTraining = (i - start) / 2;
}

void readFlowName(char** argv, int argc, int start)
{
    int i;
    for (i = start; i < argc && argv[i][0] != '-'; i ++)
    {
        string st;
        st = argv[i];
        FlowNameList.push_back(st);
    }

    return;
}

void readTestName(char** argv, int argc, int start)
{
    string st;
    st = argv[start];
    TestImgName = st;

    return;
}

void calcSmoothTerm()
{
    int pixelNum = Layer * Height * Width;
    
    Array3D_float t(Layer, Height, Width);

    memcpy(t.data, TestImg->data, sizeof(T_input) * pixelNum);

    for(int i = 0; i < 3; i++)
        pSmoothTerm[i] = new Array3D_float(Layer, Height, Width);
    

    T_input sums[3];
    sums[0] = sums[1] = sums[2] = 0;


    for(int i = 0; i < Layer; i++)
    for(int j = 0; j < Height; j++)
    for(int k = 0; k < Width; k++)
    for(int dir = 0; dir < 3; dir ++)
    {
        if(dir == 0 && i+1 == Layer)
            continue;
        if(dir == 1 && j+1 == Height)
            continue;
        if(dir == 2 && k+1 == Width)
            continue;
        T_input diff = 0;
        if(dir == 0)
            diff = t.at(i, j, k) - t.at(i+1, j, k);
        if(dir == 1)
            diff = t.at(i, j, k) - t.at(i, j+1, k);
        if(dir == 2)
            diff = t.at(i, j, k) - t.at(i, j, k+1);
            sums[dir] += sqr(diff);
    }

    for(int i = 0; i < 3; i++)
        sums[i] = 2 * sums[i] / pixelNum;

    for(int i = 0; i < Layer; i++)
    for(int j = 0; j < Height; j++)
    for(int k = 0; k < Width; k++)
    for(int dir = 0; dir < 3; dir ++)
    {
        T_input diff;
        if ((dir == 0 && i + 1 == Layer) ||
            (dir == 1 && j + 1 == Height) || 
            (dir == 2 && k + 1 == Width))
            diff = 0;
        else
        {
            if (dir == 0)
                diff = t.at(i, j, k) - t.at(i + 1, j, k);
            if (dir == 1)
                diff = t.at(i, j, k) - t.at(i, j + 1, k);
            if (dir == 2)
                diff = t.at(i, j, k) - t.at(i, j, k + 1);
        }
        T_input t = (epsilon + exp(-sqr(diff) / sums[dir])) / (epsilon + 1);
        pSmoothTerm[dir]->at(i, j, k) = pow(t, 1);
    }


    for (int i = 0; i < Layer; i++)
    for (int j = 0; j < Height; j++)
    for (int k = 0; k < Width; k++)
    for (int dir = 0; dir < 3; dir++)
    {
        pSmoothTerm[dir]->at(i, j, k) *= beta;
    }

    return;
}


void calcLikelihood()
{

    for(int i = 0; i < Layer; i++)
    for(int j = 0; j < Height; j++)
    for(int k = 0; k < Width; k++)
        if (TrainingImg->mask->at(i, j, k) == 1)
        {
            T_input sum_error;
            int c = TrainingImg->seg->at(i, j, k);
                
            sum_error = 0;
            for (int x = 0; x < OriSiftDim; x++)
                sum_error += (T_input)abs(TestImg->sift->at(i, j, k)[x] - TrainingImg->sift->at(i, j, k)[x]);
            sum_error /= OriSiftDim;
            // siftdiff->at(i, j, k)[c] = __min(siftdiff->at(i, j, k)[c], sum_error);
            siftdiff->at(i, j, k)[c] += sum_error;

            int ws = 3;
            int stepSize = 1;
            sum_error = 0;
            for (int di = -ws; di <= ws; di += stepSize)
            for (int dj = -ws; dj <= ws; dj += stepSize)
            for (int dk = -ws; dk <= ws; dk += stepSize)
            {
                int ii = enforceRange(i + di, 0, Layer - 1);
                int jj = enforceRange(j + dj, 0, Height - 1);
                int kk = enforceRange(k + dk, 0, Width - 1);

                sum_error += (T_input)abs(TestImg->data->at(ii, jj, kk) - TrainingImg->data->at(ii, jj, kk));
            }
            sum_error /= pow(2 * ws / stepSize + 1, 3);
            // graydiff->at(i, j, k)[c] = __min(graydiff->at(i, j, k)[c], sum_error);
            graydiff->at(i, j, k)[c] += sum_error;
        }

    return;
}

void calcDataTerm()
{
    pDataTerm = new Array3D_vect(Layer, Height, Width, nCategory);

    for (int i = 0; i < Layer; i++)
    for (int j = 0; j < Height; j++)
    for (int k = 0; k < Width; k++)
    for (int c = 0; c < nCategory; c++)
    {
        pDataTerm->at(i, j, k)[c] = 10 * siftdiff->at(i, j, k)[c] + graydiff->at(i, j, k)[c] + 1 * objPrior->at(i, j, k)[c];
        // pDataTerm->at(i, j, k)[c] = (siftdiff->at(i, j, k)[c] + 0.001) * (graydiff->at(i, j, k)[c] + 0.001) * (objPrior->at(i, j, k)[c] + 0.001);
    }

    return;
}

void calcObjPrior_Part1()
{

    for (int i = 0; i < Layer; i++)
    for (int j = 0; j < Height; j++)
    for (int k = 0; k < Width; k++)
    {
        int c = TrainingImg->seg->at(i, j, k);
        
        ++objPrior->at(i, j, k)[c];
        objPriorMax[c] = max(objPriorMax[c], (T_input)(objPrior->at(i, j, k)[c]));
    }

    return;
}

void calcObjPrior_Part2()
{
    for (int i = 0; i < Layer; i++)
    for (int j = 0; j < Height; j++)
    for (int k = 0; k < Width; k++)
    for (int c = 0; c < nCategory; c++)
    {
        // mean is better
        
        if (objPrior->at(i, j, k)[c] == 0)
            siftdiff->at(i, j, k)[c] = graydiff->at(i, j, k)[c] = worst;
        else
            siftdiff->at(i, j, k)[c] /= objPrior->at(i, j, k)[c], graydiff->at(i, j, k)[c] /= objPrior->at(i, j, k)[c];
        
        objPrior->at(i, j, k)[c] = (objPrior->at(i, j, k)[c] + prioreps) / (objPriorMax[c] + prioreps);
        objPrior->at(i, j, k)[c] *= satuPara[c];
        objPrior->at(i, j, k)[c] = (-log(objPrior->at(i, j, k)[c])) / (-log(prioreps / (nTraining + prioreps)));
        // objPrior->at(i, j, k)[c] = 1 - objPrior->at(i, j, k)[c];
    }

    return;
}

void printPredInfo(Array3D_int* pred, Array3D_int* target, PairwiseBP* bp)
{
    ofstream fout("..\\data\\voxelinfo");

    for (int i = 0; i < Layer; i++)
    for (int j = 0; j < Height; j++)
    for (int k = 0; k < Width; k++)
        if (pred->at(i, j, k) != target->at(i, j, k))
        {
            fout << i << "\t" << j << "\t" << k << "\t" << pred->at(i, j, k) << "\t" << target->at(i, j, k);
            fout << "\t" << siftdiff->at(i, j, k)[pred->at(i, j, k)] << "\t" << graydiff->at(i, j, k)[pred->at(i, j, k)] << "\t" << objPrior->at(i, j, k)[pred->at(i, j, k)] << "\t" << pDataTerm->at(i, j, k)[pred->at(i, j, k)];
            fout << "\t" << siftdiff->at(i, j, k)[target->at(i, j, k)] << "\t" << graydiff->at(i, j, k)[target->at(i, j, k)] << "\t" << objPrior->at(i, j, k)[target->at(i, j, k)] << "\t" << pDataTerm->at(i, j, k)[target->at(i, j, k)];
            fout << endl;
        }

    fout.close();

    return;
}

void runTest(string outputsegname)
{
    TestImg = new Image3DWithSeg(TestImgName, "");

    Layer = TestImg->data->Layer;
    Height = TestImg->data->Height;
    Width = TestImg->data->Width;

    objPrior = new Array3D_vect(Layer, Height, Width, nCategory);
    graydiff = new Array3D_vect(Layer, Height, Width, nCategory);
    siftdiff = new Array3D_vect(Layer, Height, Width, nCategory);
    objPriorMax = new T_input[nCategory];
    satuPara = new T_input[nCategory];
    for (int i = 0; i < nCategory; i++)
        objPriorMax[i] = 0, satuPara[i] = 1;

    for (int i = 0; i < nTraining; i++)
    {
        cout << "Reading Training " << i << endl;
        TrainingImg = new Image3DWithSeg(TrainingImgNameList[i], TrainingSegNameList[i]);
        if (!noFlow)
            TrainingImg->ApplyFlow(FlowNameList[i], TestImg);

        calcObjPrior_Part1();
        calcLikelihood();

        delete TrainingImg;
    }

    calcObjPrior_Part2();

    calcDataTerm();
    calcSmoothTerm();


    // Message Passing
    PairwiseBP* bp = new PairwiseBP();
    bp->LoadDataTerm(pDataTerm);
    bp->LoadSmoothness(pSmoothTerm);
    bp->MessagePassing(12, 3, true);

    Array3D_int* mySeg = bp->result;

    // Output Seg
    TestImg->seg = mySeg;
    TestImg->segHDR = TestImg->dataHDR;
    TestImg->segHDR.dime.bitpix = 8;
    TestImg->segHDR.dime.datatype = DT_UNSIGNED_CHAR;
    TestImg->OutputSegImage(outputsegname);

    // printPredInfo(mySeg, TestImg->seg, bp);

    delete objPrior;
    delete graydiff;
    delete siftdiff;
    delete bp;
    delete TestImg;
    delete[] objPriorMax;
    delete[] satuPara;

    return;
}

int findPara(const char *pattern, int argc, char** argv)
{
    for(int i = 0; i < argc; ++i)
        if (strcmp(argv[i], pattern) == 0)
            return i;

    return -1;
}

int main(int argc, char** argv)
{

    // -train -flow -test -output -noflow -eps -beta -prioreps -nlabels -worst
    int traini, flowi, testi, outputi, noflowi, epsi, betai, priorepsi, nlabelsi, worsti, grayavgi, siftavgi;

    traini = findPara("-train", argc, argv);
    flowi = findPara("-flow", argc, argv);
    testi = findPara("-test", argc, argv);
    outputi = findPara("-output", argc, argv);
    noflowi = findPara("-noflow", argc, argv);
    epsi = findPara("-eps", argc, argv);
    betai = findPara("-beta", argc, argv);
    priorepsi = findPara("-prioreps", argc, argv);
    nlabelsi = findPara("-nlabels", argc, argv);
    worsti = findPara("-worst", argc, argv);
    grayavgi = findPara("-grayavg", argc, argv);
    siftavgi = findPara("-siftavg", argc, argv);

    if (noflowi != -1)
        noFlow = true;
    if (epsi != -1)
        epsilon = atof(argv[epsi + 1]);
    if (betai != -1)
        beta = atof(argv[betai + 1]);
    if (priorepsi != -1)
        prioreps = atof(argv[priorepsi + 1]);
    if (nlabelsi != -1)
        nCategory = atoi(argv[nlabelsi + 1]);
    if (worsti != -1)
        worst = atof(argv[worsti + 1]);
    if (grayavgi != -1) 
        GrayAvgBaseLine = atof(argv[grayavgi + 1]);
    if (siftavgi != -1) 
        SIFTAvgBaseLine = atof(argv[siftavgi + 1]);

    // Read names
    readTrainingNameList(argv, argc, traini + 1);
    if (!noFlow)
        readFlowName(argv, argc, flowi + 1);
    readTestName(argv, argc, testi + 1);
    string outputsegname = argv[outputi + 1];

    // Run
    runTest(outputsegname);

    return 0;
}

