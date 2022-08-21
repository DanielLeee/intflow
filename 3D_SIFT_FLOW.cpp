#include <iostream>
#include <fstream>
#include "Image3D.h"
#include "SIFT_FLOW_ALGO.h"
using namespace std;

void OutputEnergyFile(string fileName, T_input energy)
{
	ofstream fout(fileName);

	fout << energy << endl;

	fout.close();

	return;
}

void CalcSIFTFlowAndWarp(Image3DWithSeg* fix, Image3DWithSeg* mov, string outputFlowName, string outputEnergyName)
{
	T_input energy;

	SIFT_FLOW_ALGO(fix, mov, outputFlowName, energy);

	OutputEnergyFile(outputEnergyName, energy);

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
	// -fix -mov -o -alpha -gamma -nlevels -wsize -niter -nscale -grayavg -siftavg
	int fixi, movi, outputi, alphai, gammai, nlevelsi, wsizei, niteri, nscalei, grayavgi, siftavgi;
	
    fixi = findPara("-fix", argc, argv);
    movi = findPara("-mov", argc, argv);
    outputi = findPara("-o", argc, argv);
    alphai = findPara("-alpha", argc, argv);
    gammai = findPara("-gamma", argc, argv);
    nlevelsi = findPara("-nlevels", argc, argv);
    wsizei = findPara("-wsize", argc, argv);
    niteri = findPara("-niter", argc, argv);
    nscalei = findPara("-nscale", argc, argv);
    grayavgi = findPara("-grayavg", argc, argv);
    siftavgi = findPara("-siftavg", argc, argv);

    if (alphai != -1)
        alpha = atof(argv[alphai + 1]), d = alpha * 20;
    if (gammai != -1)
        gammaPara = atof(argv[gammai + 1]);
    if (nlevelsi != -1)
        nlevels = atoi(argv[nlevelsi + 1]);
    if (wsizei != -1)
        wsize = atoi(argv[wsizei + 1]);
    if (niteri != -1)
        nIterations = atoi(argv[niteri + 1]);
    if (nscalei != -1)
        NumScales = atoi(argv[nscalei + 1]);
    if (grayavgi != -1)
        GrayAvgBaseLine = atof(argv[grayavgi + 1]);
    if (siftavgi != -1)
        SIFTAvgBaseLine = atof(argv[siftavgi + 1]);

	Image3DWithSeg* testIMG = NULL;
	Image3DWithSeg* trainIMG = NULL;

	// Read Images
	testIMG = new Image3DWithSeg(argv[fixi + 1], "");
	trainIMG = new Image3DWithSeg(argv[movi + 1], "");

	// SIFT Flow
	cout << "SIFT Flow: " << testIMG->dataFileName << " -> " << trainIMG->dataFileName << endl;
	CalcSIFTFlowAndWarp(testIMG, trainIMG, argv[outputi + 1], argv[outputi + 2]);

	delete trainIMG;
	delete testIMG;

	return 0;
}
