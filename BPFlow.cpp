#include "BPFlow.h"

BPFlow::BPFlow(void)
{

}

BPFlow::~BPFlow(void)
{
    for (int i = 0; i < 3; i++)
    {
        delete pX[i];
        delete pRangeTerm[i];
        delete pDualMessage[i];
        delete pBelief[i];

        for (int j = 0; j < nNeighbors; j++)
            delete pSpatialMessage[i][j];
    }

    delete foo;
    delete Im_s;
    delete Im_d;
    delete pDataTerm;
}

void BPFlow::LoadOffset(Array3D_int *vx, Array3D_int *vy, Array3D_int *vz)
{
    pOffset[0] = vx;
    pOffset[1] = vy;
    pOffset[2] = vz;

    return;
}

void BPFlow::LoadWinSize(int wsize)
{
    winSize = wsize;

    return;
}

int BPFlow::EnforceRange(const int x, const int MaxValue)
{
    return min(max(x, 0), MaxValue - 1);
}

//---------------------------------------------------------------------------------
// set the parameter of the model
// the regularization is a truncated L1 norm: min(s|v(x,y)-v(x+1,y)|,d)
// sigma is used to penalize large displacement
//---------------------------------------------------------------------------------

void BPFlow::setPara(T_input _s, T_input _d)
{
    s=_s;
    d=_d;

    if(Width>0)
    {
        foo = new Array3D_vect(Slice, Height, Width, 3);
        Im_s = new Array3D_vect(Slice, Height, Width, 3);
        Im_d = new Array3D_vect(Slice, Height, Width, 3);
        Im_s->setValue(s);
        Im_d->setValue(d);
    }
    else
        printf("The image dimension has not been specified! Call LoadImages() first\n");
}

//----------------------------------------------------------------------
// function to load images
//----------------------------------------------------------------------

void BPFlow::LoadImages(Array3D_vect* pImage1, Array3D_vect* pImage2)
{
    Width=pImage1->Width;
    Height=pImage1->Height;
    Slice=pImage1->Layer;
    Volume=Width*Height*Slice;
    nChannels=pImage1->nDim;
    Width2=pImage2->Width;
    Height2=pImage2->Height;
    Slice2=pImage2->Layer;

    pIm1 = pImage1;
    pIm2 = pImage2;

    return;
}


//------------------------------------------------------------------------------------------------
// function to verify whether a poit is inside image boundary or not
//------------------------------------------------------------------------------------------------

template <class T>
bool BPFlow::InsideImage(T x,T y,T z)
{
    if(x >= 0 && x < Width2 && y >= 0 && y < Height2 && z >= 0 && z < Slice2)
        return true;
    else
        return false;
}

//------------------------------------------------------------------------------------------------
// function to compute range term
//------------------------------------------------------------------------------------------------

void BPFlow::ComputeRangeTerm(T_input _gamma)
{
    gamma=_gamma;

    for (int i = 0; i < 3; i++)
        pRangeTerm[i] = new Array3D_vect(Slice, Height, Width, 2 * winSize + 1);

    for (int t = 0; t < Slice; t++)              //index over z
    for (int i = 0; i < Height; i++)            // index over y
    for (int j = 0; j < Width; j++)        // index over x
    for (int plane = 0; plane < 3; plane++)
    for (int w = -winSize; w <= winSize; w++)
        pRangeTerm[plane]->at(t, i, j)[w + winSize] = gamma*fabs((T_input)w + pOffset[plane]->at(t, i, j));

    return;
}

//------------------------------------------------------------------------------------------------
// function to compute data term
//------------------------------------------------------------------------------------------------

void BPFlow::ComputeDataTerm()
{
    // allocate the buffer for data term
    pDataTerm = new Array3D_vect(Slice, Height, Width, (2 * winSize + 1) * (2 * winSize + 1) * (2 * winSize + 1));

    T_input HistMin, HistMax, HistInterval, HistAverage;
    long long* pHistogramBuffer;
    const int nBins=200000;
    long long total = 0; // total is the total number of plausible matches, used to normalize the histogram
    pHistogramBuffer = new long long[nBins];
    memset(pHistogramBuffer, 0, sizeof(long long)*nBins);
    HistMin=1000000000;
    HistMax=0;
    HistAverage = 0;
    //--------------------------------------------------------------------------------------------------
    // step 1. the first sweep to compute the data term for the visible matches
    //--------------------------------------------------------------------------------------------------
    for(int t = 0; t < Slice; t++)              //index over z
    for(int i = 0; i < Height; i++)            // index over y
    for(int j = 0; j < Width; j++)        // index over x
    {
        int WinLength=winSize*2+1;

        // loop over a local window
        for(int m=-winSize;m<=winSize;m++)  // index over z
        for(int k=-winSize;k<=winSize;k++)  // index over y
        for(int l=-winSize;l<=winSize;l++)  // index over x
        {
            int x = EnforceRange(j+pOffset[0]->at(t, i, j)+l, Width2);
            int y = EnforceRange(i+pOffset[1]->at(t, i, j)+k, Height2); 
            int z = EnforceRange(t+pOffset[2]->at(t, i, j)+m, Slice2);

            // if the point is outside the image boundary then continue
            // if(!InsideImage(x,y,z))
            //    continue;

            T_input foo = 0;
            for(int n=0;n<nChannels;n++)
                foo+=fabs(pIm1->at(t, i, j)[n]-pIm2->at(z, y, x)[n]); // L1 norm

            HistAverage += foo;

            pDataTerm->at(t, i, j)[(m + winSize)*WinLength*WinLength + (k + winSize)*WinLength + l + winSize] = foo;
            HistMin=min(HistMin,foo);
            HistMax=max(HistMax,foo);
            ++total;
        }
    }

    delete[] pHistogramBuffer;
    return;

    // compute the histogram info
    HistInterval=(HistMax-HistMin)/nBins;
    HistAverage /= total;
    printf("HistAverage: %f\n", HistAverage);
    //--------------------------------------------------------------------------------------------------
    // step 2. get the histogram of the matching
    //--------------------------------------------------------------------------------------------------
    for(int t=0;t<Slice;t++)           // index over z
    for(int i=0;i<Height;i++)            // index over y
    for(int j=0;j<Width;j++)        // index over x
    {
        int WinLength = winSize * 2 + 1;

        // loop over a local window
        for(int m=-winSize;m<=winSize;m++)  // index over z
        for(int k=-winSize;k<=winSize;k++)  // index over y
        for(int l=-winSize;l<=winSize;l++)  // index over x
        {
            int x = j + pOffset[0]->at(t, i, j) + l;
            int y = i + pOffset[1]->at(t, i, j) + k;
            int z = t + pOffset[2]->at(t, i, j) + m;
            
            // if the point is outside the image boundary then continue
            if(!InsideImage(x,y,z))
                continue;

            int foo = (int)min((pDataTerm->at(t, i, j)[(m + winSize)*WinLength*WinLength + (k + winSize)*WinLength + l + winSize] - HistMin) / HistInterval, (T_input)(nBins - 1));
            pHistogramBuffer[foo]++;
        }
    }

    T_input DefaultMatchingScore;
    long long preSum = 0;
    for(int i=0;i<nBins;i++)
    {
        preSum += pHistogramBuffer[i];
        if (preSum * 2 >= total) // find the matching score
        {
            DefaultMatchingScore=max(i,1)*HistInterval+HistMin; 
            break;
        }
    }


    //--------------------------------------------------------------------------------------------------
    // step 3. assigning the default matching score to the outside matches
    //--------------------------------------------------------------------------------------------------
    for(int t=0;t<Slice;t++)           // index over z
    for(int i=0;i<Height;i++)            // index over y
    for(int j=0;j<Width;j++)        // index over x
    {
        int WinLength = winSize * 2 + 1;

        // loop over a local window
        for(int m=-winSize;m<=winSize;m++)  // index over z
        for(int k=-winSize;k<=winSize;k++)  // index over y
        for(int l=-winSize;l<=winSize;l++)  // index over x
        {
            int x = j + pOffset[0]->at(t, i, j) + l;
            int y = i + pOffset[1]->at(t, i, j) + k;
            int z = t + pOffset[2]->at(t, i, j) + m;
            
            int _ptr = (m + winSize)*WinLength*WinLength + (k + winSize)*WinLength + l + winSize;
            // if the point is outside the image boundary then continue
            if(!InsideImage(x,y,z))
                pDataTerm->at(t, i, j)[_ptr] = DefaultMatchingScore;
        }
    }

    delete[] pHistogramBuffer;
}

//------------------------------------------------------------------------------------------------
//    function to allocate buffer
//------------------------------------------------------------------------------------------------

void BPFlow::AllocateBuffer()
{
    for (int i = 0; i<3; i++)
        pOffset[i] = new Array3D_int(Slice, Height, Width);
    pDataTerm = new Array3D_vect(Slice, Height, Width, (2 * winSize + 1) * (2 * winSize + 1) * (2 * winSize + 1));

    return;
}

//------------------------------------------------------------------------------------------------
//    function to allocate buffer for the messages
//------------------------------------------------------------------------------------------------

void BPFlow::AllocateMessage()
{
    // allocate the buffers for the messages
    int WinLength = 2 * winSize + 1;
    for (int i=0;i<3;i++)
    {
        for (int j = 0; j < nNeighbors; j ++)
            pSpatialMessage[i][j] = new Array3D_vect(Slice, Height, Width, WinLength);
        pDualMessage[i] = new Array3D_vect(Slice, Height, Width, WinLength);
        pBelief[i] = new Array3D_vect(Slice, Height, Width, WinLength);
    }

    return;
}

//------------------------------------------------------------------------------------------------
// function for belief propagation
//------------------------------------------------------------------------------------------------

void BPFlow::MessagePassing(int nIterations, int nHierarchy, T_input &energy, bool onTop)
{
    AllocateMessage();

    if(nHierarchy>0)
    {
        BPFlow bp;
        generateCoarserLevel(bp);
        bp.MessagePassing(24,nHierarchy-1, energy, false);
        bp.propagateFinerLevel(*this);

        for (int i = 0; i < 3; i++)
            delete bp.pOffset[i];
    }

    for (int i = 0; i < 3; i++)
        pX[i] = new Array3D_int(Slice, Height, Width);

    for(int count=0;count<nIterations;count++)
    {
        BP_S(count);
    }

    if (onTop)
    {
        ComputeBelief();
        FindOptimalSolution();
        energy = GetEnergy();
    }

    return;
}

void BPFlow::BP_S(int count)
{
    //  ---- Update Order( Dual, forward-backward Spatial) ----
    int k = count % 3;

    // if (count % 6<3)
    for (int t = 0; t<Slice; t++)
    for (int i = 0; i<Height; i++)
    for (int j = 0; j<Width; j++)
    {
        UpdateDualMessage(j, i, t, k);
    }

    //forward update
    for(int t=0;t<Slice;t++)
    for(int i=0;i<Height;i++)
    for(int j=0;j<Width;j++)
    {
        UpdateSpatialMessage(j,i,t,k,0);
        UpdateSpatialMessage(j,i,t,k,2);
        UpdateSpatialMessage(j,i,t,k,4);
    }
    // backward upate
    for(int t=Slice-1;t>=0;t--)
    for(int i=Height-1;i>=0;i--)
    for(int j=Width-1;j>=0;j--)
    {
        UpdateSpatialMessage(j,i,t,k,1);
        UpdateSpatialMessage(j,i,t,k,3);
        UpdateSpatialMessage(j,i,t,k,5); 
    }
}

//------------------------------------------------------------------------------------------------
//  update the message from (x0,y0,plane) to the neighbors on the same plane
//    the encoding of the direction
//             2
//             |
//             |
//             v
//    0 ------> <------- 1
//             ^
//             |
//             |
//             3
//------------------------------------------------------------------------------------------------

void BPFlow::UpdateSpatialMessage(int x, int y, int z, int plane, int direction)
{
    // eliminate impossible messages
    if (direction==0 && x==Width-1)
        return;
    if (direction==1 && x==0)
        return;
    if (direction==2 && y==Height-1)
        return;
    if (direction==3 && y==0)
        return;
    if (direction==4 && z==Slice-1)
        return;    
    if (direction==5 && z==0)
        return;
        
    int nStates=winSize*2+1;
    int wsize = winSize;

    T_input *message_org = new T_input[nStates];

    memset(message_org, 0, sizeof(T_input)* nStates);

    int x1=x,y1=y,z1=z; // get the destination
    switch(direction)
    {
        case 0:  x1++;  break;
        case 1:  x1--;  break;
        case 2:  y1++;  break;
        case 3:  y1--;  break;
        case 4:  z1++;  break;
        case 5:  z1--;  break;
    }

    s = Im_s->at(z1, y1, x1)[plane];
    d = Im_d->at(z1, y1, x1)[plane];

    T_input* message = pSpatialMessage[plane][direction]->at(z1, y1, x1);

    // initialize the message from the dual plane
    Add2Message(message_org, pDualMessage[plane]->at(z, y, x), nStates);

    // add the range term
    Add2Message(message_org, pRangeTerm[plane]->at(z, y, x), nStates);
    
    // add spatial messages
    if(x>0 && direction!=1) // add left to right
        Add2Message(message_org, pSpatialMessage[plane][0]->at(z, y, x), nStates);
    if(x<Width-1 && direction!=0) // add right to left 
        Add2Message(message_org, pSpatialMessage[plane][1]->at(z, y, x), nStates);
    if(y>0 && direction!=3) // add top down
        Add2Message(message_org, pSpatialMessage[plane][2]->at(z, y, x), nStates);
    if(y<Height-1 && direction!=2) // add bottom up
        Add2Message(message_org, pSpatialMessage[plane][3]->at(z, y, x), nStates);
    if(z>0 && direction!=5) // add bottom up
        Add2Message(message_org, pSpatialMessage[plane][4]->at(z, y, x), nStates);
    if(z<Slice-1 && direction!=4) // add bottom up
        Add2Message(message_org, pSpatialMessage[plane][5]->at(z, y, x), nStates);

    // use distance transform function to impose smoothness compatibility
    T_input Min = FindVectorMin(message_org, nStates) + d;
    for(int l=1;l<nStates;l++)
        message_org[l]=min(message_org[l],message_org[l-1]+s);
    for(int l=nStates-2;l>=0;l--)
        message_org[l]=min(message_org[l],message_org[l+1]+s);


    // transform the compatibility 
    int shift = pOffset[plane]->at(z1, y1, x1) - pOffset[plane]->at(z, y, x);
    if(abs(shift)>wsize+wsize) // the shift is too big that there is no overlap
    {
        if (shift>0)
            for(int l=0;l<nStates;l++)
                message[l]=l*s;
        else
            for(int l=0;l<nStates;l++)
                message[l] = (nStates - 1 - l)*s;
    }
    else
    {
        int start=max(-wsize,shift-wsize);
        int end=min(wsize,shift+wsize);
        for(int i=start;i<=end;i++)
            message[i-shift+wsize]=message_org[i+wsize];
        if(start-shift+wsize>0)
            for(int i=start-shift+wsize-1;i>=0;i--)
                message[i]=message[i+1]+s;
        if(end-shift+wsize<nStates)
            for(int i=end-shift+wsize+1;i<nStates;i++)
                message[i]=message[i-1]+s;
    }

    // put back the threshold
    for(int l=0;l<nStates;l++)
        message[l]=min(message[l],Min);

    // normalize the message by subtracting the minimum value
    Min = FindVectorMin(message, nStates);
    for(int l=0;l<nStates;l++)
        message[l]-=Min;

    delete[] message_org;
}

void BPFlow::Add2Message(T_input *m1, T_input *m2, int len, T_input cof)
{
    for (int i = 0; i<len; i++)
        m1[i] += m2[i] * cof;

    return;
}

//------------------------------------------------------------------------------------------------
// update dual message passing from one plane to the other
//------------------------------------------------------------------------------------------------

void BPFlow::UpdateDualMessage(int x, int y, int z, int plane)
{
    T_input Min;
    int wsize = winSize;
    int nstates = wsize * 2 + 1;

    s = Im_s->at(z, y, x)[plane];
    d = Im_d->at(z, y, x)[plane];

    T_input *message_org_X = new T_input[nstates];
    T_input *message_org_Y = new T_input[nstates];
    T_input *message_org_Z = new T_input[nstates];

    memset(message_org_X, 0, sizeof(T_input)* nstates);
    memset(message_org_Y, 0, sizeof(T_input)* nstates);
    memset(message_org_Z, 0, sizeof(T_input)* nstates);

    // add spatial messages
    if (x>0)  //add left to right
    {
        Add2Message(message_org_X, pSpatialMessage[0][0]->at(z, y, x), nstates);
        Add2Message(message_org_Y, pSpatialMessage[1][0]->at(z, y, x), nstates);
        Add2Message(message_org_Z, pSpatialMessage[2][0]->at(z, y, x), nstates);
    }
    if (x<Width - 1) // add right to left
    {
        Add2Message(message_org_X, pSpatialMessage[0][1]->at(z, y, x), nstates);
        Add2Message(message_org_Y, pSpatialMessage[1][1]->at(z, y, x), nstates);
        Add2Message(message_org_Z, pSpatialMessage[2][1]->at(z, y, x), nstates);
    }
    if (y>0) // add top down
    {
        Add2Message(message_org_X, pSpatialMessage[0][2]->at(z, y, x), nstates);
        Add2Message(message_org_Y, pSpatialMessage[1][2]->at(z, y, x), nstates);
        Add2Message(message_org_Z, pSpatialMessage[2][2]->at(z, y, x), nstates);
    }
    if (y<Height - 1) // add bottom up
    {
        Add2Message(message_org_X, pSpatialMessage[0][3]->at(z, y, x), nstates);
        Add2Message(message_org_Y, pSpatialMessage[1][3]->at(z, y, x), nstates);
        Add2Message(message_org_Z, pSpatialMessage[2][3]->at(z, y, x), nstates);
    }
    if (z>0) // 
    {
        Add2Message(message_org_X, pSpatialMessage[0][4]->at(z, y, x), nstates);
        Add2Message(message_org_Y, pSpatialMessage[1][4]->at(z, y, x), nstates);
        Add2Message(message_org_Z, pSpatialMessage[2][4]->at(z, y, x), nstates);
    }
    if (z<Slice - 1) // 
    {
        Add2Message(message_org_X, pSpatialMessage[0][5]->at(z, y, x), nstates);
        Add2Message(message_org_Y, pSpatialMessage[1][5]->at(z, y, x), nstates);
        Add2Message(message_org_Z, pSpatialMessage[2][5]->at(z, y, x), nstates);
    }

    T_input *message_X = pDualMessage[0]->at(z, y, x);
    T_input *message_Y = pDualMessage[1]->at(z, y, x);
    T_input *message_Z = pDualMessage[2]->at(z, y, x);

    if (plane == 0)
    {
        //  Layer YZ  --->  Layer X
        for (int xx = 0; xx < nstates; xx++)
        {
            Min = FLT_MAX;
            for (int zz = 0; zz < nstates; zz++)
            {
                for (int yy = 0; yy < nstates; yy++)
                {
                    Min = min(Min, pDataTerm->at(z, y, x)[zz * nstates * nstates + yy * nstates + xx] + message_org_Y[yy] + message_org_Z[zz]);
                }
            }

            message_X[xx] = Min;
        }

        Min = FindVectorMin(message_X, nstates);
        for (int l = 0; l<nstates; l++)
            message_X[l] -= Min;
    }
    else if (plane == 1)
    {
        //  Layer XZ  --->  Layer Y
        for (int yy = 0; yy < nstates; yy++)
        {
            Min = FLT_MAX;
            for (int zz = 0; zz < nstates; zz++)
            {
                for (int xx = 0; xx < nstates; xx++)
                {
                    Min = min(Min, pDataTerm->at(z, y, x)[zz * nstates * nstates + yy * nstates + xx] + message_org_X[xx] + message_org_Z[zz]);
                }
            }

            message_Y[yy] = Min;
        }

        Min = FindVectorMin(message_Y, nstates);
        for (int l = 0; l<nstates; l++)
            message_Y[l] -= Min;
    }
    else if (plane == 2)
    {
        //  Layer XY  --->  Layer Z
        for (int zz = 0; zz < nstates; zz++)
        {
            Min = FLT_MAX;
            for (int yy = 0; yy < nstates; yy++)
            {
                for (int xx = 0; xx < nstates; xx++)
                {
                    Min = min(Min, pDataTerm->at(z, y, x)[zz * nstates * nstates + yy * nstates + xx] + message_org_X[xx] + message_org_Y[yy]);
                }
            }

            message_Z[zz] = Min;
        }

        Min = FindVectorMin(message_Z, nstates);
        for (int l = 0; l<nstates; l++)
            message_Z[l] -= Min;
        
    }

    delete[] message_org_X;
    delete[] message_org_Y;
    delete[] message_org_Z;
}

//------------------------------------------------------------------------------------------------
// compute belief
//------------------------------------------------------------------------------------------------

void BPFlow::ComputeBelief()
{
    int nStates = winSize * 2 + 1;

    for(int plane=0;plane<3;plane++)
    for(int t=0;t<Slice;t++)
    for(int i=0;i<Height;i++)
    for(int j=0;j<Width;j++)
    {
        T_input *belief = pBelief[plane]->at(t, i, j);

        // add range term
        Add2Message(belief, pRangeTerm[plane]->at(t, i, j), nStates);

        // add message from the dual layer
        Add2Message(belief, pDualMessage[plane]->at(t, i, j), nStates);

        if(j>0)
            Add2Message(belief, pSpatialMessage[plane][0]->at(t, i, j), nStates);
        if(j<Width-1)
            Add2Message(belief, pSpatialMessage[plane][1]->at(t, i, j), nStates);
        if(i>0)
            Add2Message(belief, pSpatialMessage[plane][2]->at(t, i, j), nStates);
        if(i<Height-1)
            Add2Message(belief, pSpatialMessage[plane][3]->at(t, i, j), nStates);
        if(t>0)
            Add2Message(belief, pSpatialMessage[plane][4]->at(t, i, j), nStates);
        if(t<Slice-1)
            Add2Message(belief, pSpatialMessage[plane][5]->at(t, i, j), nStates);
    }
}

void BPFlow::FindOptimalSolution()
{
    int nStates = winSize * 2 + 1;

    for(int plane=0; plane<3; plane++)
    for (int t = 0; t<Slice; t++)
    for (int i = 0; i<Height; i++)
    for (int j = 0; j<Width; j++)
    {
        T_input Min;
        int index=0;
        Min = FLT_MAX;
        for(int l=0;l<nStates;l++)
        if (Min>pBelief[plane]->at(t, i, j)[l])
        {
            Min = pBelief[plane]->at(t, i, j)[l];
            index = l;
        }
        pX[plane]->at(t, i, j) = index;
    }

    return;
}


//------------------------------------------------------------------------------------------------
// function to get energy
//------------------------------------------------------------------------------------------------

T_input BPFlow::GetEnergy()
{
    T_input energy = 0;
    for(int t=0;t<Slice;t++)
    for(int i=0;i<Height;i++)
    for(int j=0;j<Width;j++)
    {
        for (int k = 0; k < 3; k++)
        {
            int tt = t, ii = i, jj = j;

            s = Im_s->at(t, i, j)[k];
            d = Im_d->at(t, i, j)[k];

            if (j < Width - 1)
            {
                jj++;
                energy += min((T_input)abs(pX[k]->at(t, i, j) - winSize + pOffset[k]->at(t, i, j) - pX[k]->at(tt, ii, jj) + winSize - pOffset[k]->at(tt, ii, jj))*s, d);
                jj--;
            }
            if (i < Height - 1)
            {
                ii++;
                energy += min((T_input)abs(pX[k]->at(t, i, j) - winSize + pOffset[k]->at(t, i, j) - pX[k]->at(tt, ii, jj) + winSize - pOffset[k]->at(tt, ii, jj))*s, d);
                ii--;
            }
            if (t < Slice - 1)
            {
                tt++;
                energy += min((T_input)abs(pX[k]->at(t, i, j) - winSize + pOffset[k]->at(t, i, j) - pX[k]->at(tt, ii, jj) + winSize - pOffset[k]->at(tt, ii, jj))*s, d);
                tt--;
            }
        }

        int vx = pX[0]->at(t, i, j);
        int vy = pX[1]->at(t, i, j);
        int vz = pX[2]->at(t, i, j);
        int nStates=winSize*2+1;
        energy += pDataTerm->at(t, i, j)[vz*nStates*nStates + vy*nStates + vx];
        for(int k=0;k<3;k++)
            energy += pRangeTerm[k]->at(t, i, j)[pX[k]->at(t, i, j)];
    }
    
    return energy;
}

void BPFlow::ComputeVelocity(Array3D_vect *mFlow)
{
    for (int t = 0; t<Slice; t++)   //index over z
    for (int i = 0; i<Height; i++)    // index over y
    for (int j = 0; j<Width; j++)    // index over x
    for (int plane = 0; plane < 3; plane++)
        mFlow->at(t, i, j)[plane] = pX[plane]->at(t, i, j) + pOffset[plane]->at(t, i, j) - winSize;

    return;
}

void BPFlow::ImageFilter(Array3D_vect *src, Array3D_vect *dst, T_input *filter, int fsize)
{
    int wsize = fsize * 2 + 1;
    int nChannels = src->nDim;

    for (int t = 0; t < dst->Layer; t++)
    for (int i = 0; i < dst->Height; i++)
    for (int j = 0; j < dst->Width; j++)
    {
        for (int x = 0; x < nChannels; ++x)
            dst->at(t, i, j)[x] = 0;

        for (int m = -fsize; m <= fsize; m++)
        for (int u = -fsize; u <= fsize; u++)
        for (int v = -fsize; v <= fsize; v++)
        {
            T_input w = filter[(m + fsize)*wsize*wsize + (u + fsize)*wsize + v + fsize];

            int tt = EnforceRange(t + m, src->Layer);
            int ii = EnforceRange(i + u, src->Height);
            int jj = EnforceRange(j + v, src->Width);
            
            for (int k = 0; k < nChannels; k++)
                dst->at(t, i, j)[k] += src->at(tt, ii, jj)[k] * w;
        }
    }

    return;
}

void BPFlow::Smoothing(Array3D_vect *im1, Array3D_vect *im2)
{
    T_input filter3D[27] = { 0.0062, 0.0214, 0.0062, 0.0214, 0.0735, 0.0214, 0.0062, 0.0214, 0.0062,
                            0.0214, 0.0735, 0.0214, 0.0735, 0.2520, 0.0735, 0.0214, 0.0735, 0.0214,
                            0.0062, 0.0124, 0.0062, 0.0214, 0.0735, 0.0214, 0.0062, 0.0214, 0.0062 };

    ImageFilter(im1, im2, filter3D, 1);

    return;
}

void BPFlow::ResizeImageHalf(Array3D_vect *src, Array3D_vect *dst)
{
    int nChannels = dst->nDim;

    memset(dst->data, 0, sizeof(T_input) * dst->nElements);
    for (int t = 0; t<dst->Layer; t++)
    for (int i = 0; i<dst->Height; i++)
    for (int j = 0; j<dst->Width; j++)
    {
        int x = EnforceRange(j * 2 + 1, src->Width);
        int y = EnforceRange(i * 2 + 1, src->Height);
        int z = EnforceRange(t * 2 + 1, src->Layer);

        for (int l = 0; l < nChannels; l++)
            dst->at(t, i, j)[l] = src->at(z, y, x)[l];
    }

    return;
}

void BPFlow::ReduceImageHalf(Array3D_vect *src, Array3D_vect *dst, bool computeAvg)
{
    int nChannels = dst->nDim;

    memset(dst->data, 0, sizeof(T_input) * dst->nElements);
    for (int t = 0; t<dst->Layer; t++)
    for (int i = 0; i<dst->Height; i++)
    for (int j = 0; j<dst->Width; j++)
    {
        int sum = 0;

        for (int tt = 0; tt<2; tt++)
        for (int ii = 0; ii<2; ii++)
        for (int jj = 0; jj<2; jj++)
        {
            int x = j * 2 + jj;
            int y = i * 2 + ii;
            int z = t * 2 + tt;
            if (y<src->Height && x<src->Width && z<src->Layer)
            {
                for (int l = 0; l < nChannels; l++)
                    dst->at(t, i, j)[l] += src->at(z, y, x)[l];
                sum++;
            }
        }

        if (sum != 0 && computeAvg)
            for (int l = 0; l < nChannels; l++)
                dst->at(t, i, j)[l] /= sum;
    }

    return;
}

void BPFlow::ReduceImageHalf(Array3D_int *src, Array3D_int *dst, bool computeAvg)
{
    memset(dst->data, 0, sizeof(int) * dst->Volume);
    for (int t = 0; t < dst->Layer; t++)
    for (int i = 0; i < dst->Height; i++)
    for (int j = 0; j < dst->Width; j++)
    {
        int sum = 0;

        for (int tt = 0; tt < 2; tt++)
        for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
        {
            int x = j * 2 + jj;
            int y = i * 2 + ii;
            int z = t * 2 + tt;
            if (y < src->Height && x < src->Width && z < src->Layer)
            {
                dst->at(t, i, j) += src->at(z, y, x);
                sum++;
            }
        }

        if (sum != 0 && computeAvg)
            dst->at(t, i, j) /= sum;
    }

    return;
}

T_input BPFlow::FindVectorMin(T_input *v, int len)
{
    T_input ret = FLT_MAX;
    for (int i = 0; i < len; i++)
        ret = min(ret, v[i]);

    return ret;
}

//------------------------------------------------------------------------------------------------
// multi-grid belie propagation
//------------------------------------------------------------------------------------------------

void BPFlow::generateCoarserLevel(BPFlow &bp)
{
    //------------------------------------------------------------------------------------------------
    // set the dimensions and parameters
    //------------------------------------------------------------------------------------------------
    bp.Width = (Width + 1) / 2;
    bp.Height = (Height + 1) / 2;    
    bp.Slice = (Slice + 1) / 2;
    bp.Volume = bp.Width * bp.Height * bp.Slice;
    bp.nChannels = nChannels;
    bp.winSize = winSize;
    bp.setPara(s, d);
    bp.AllocateBuffer();

    Smoothing(Im_s, foo);
    ReduceImageHalf(foo, bp.Im_s);

    Smoothing(Im_d, foo);
    ReduceImageHalf(foo, bp.Im_d);

    //------------------------------------------------------------------------------------------------
    // allocate buffers
    //------------------------------------------------------------------------------------------------
    for(int i=0;i<3;i++)
        ReduceImageHalf(pOffset[i], bp.pOffset[i]);
    //------------------------------------------------------------------------------------------------
    // generate data term
    //------------------------------------------------------------------------------------------------
    ReduceImageHalf(pDataTerm, bp.pDataTerm);
    //------------------------------------------------------------------------------------------------
    // generate range term
    //------------------------------------------------------------------------------------------------
    bp.ComputeRangeTerm(gamma / 2);

    return;
}

void BPFlow::propagateFinerLevel(BPFlow &bp)
{
    int WinLength = winSize * 2 + 1;

    for(int t=0;t<bp.Slice;t++)  
    for(int i=0;i<bp.Height;i++)
    for(int j=0;j<bp.Width;j++)
    {
        int z=t/2;
        int y=i/2;
        int x=j/2;

        for (int k = 0; k < 3; k++)
        {
            memcpy(bp.pDualMessage[k]->at(t, i, j), pDualMessage[k]->at(z, y, x), sizeof(T_input) * WinLength);
            for (int l = 0; l < nNeighbors; l++)
                memcpy(bp.pSpatialMessage[k][l]->at(t, i, j), pSpatialMessage[k][l]->at(z, y, x), sizeof(T_input) * WinLength);
        }
    }

    return;
}

