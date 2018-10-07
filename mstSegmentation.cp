/**************************************************
 Author: David Esparza Alba
 Date: Someday in 2011
 Ritsumeikan University
 **************************************************/

#include <highgui.h>
#include <cv.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <algorithm>

using namespace std;

enum CHANNELS
{
    RED = 0,
    GREEN,
    BLUE,
    INTENSITY
};

typedef struct stImage
{
    double val[4];
}Image;

typedef struct stEdge
{
    int v1;
    int v2;
    double w;
}Edge;

IplImage *imgSrc;
IplImage *imgDst;
IplImage *imgGray;

Edge *E;
int *p, *rango, *C, *COL;
int nRows;
int nColumns;
double *Int;
double thao = 10000.0;

Image ** ConvImgToDbl(IplImage *);
void ConvDblToImg(Image **, IplImage *, bool);
void Smooth(Image **, Image **, int);
void MST_Segmentation(Image **, int);
int compare(const void *, const void *);

void make_set(int);
void link(int, int);
int find_set(int);
void union_set(int, int);

int Min(double, double);
int Max(double, double);

//Original main
int main()
{
    int i, j;
    
    Image **X;
    //Image **imgR;
    //Image **imgG;
    //Image **imgB;
    Image **imgI;
    
    srand(time(NULL));
    //freopen("mat4.txt", "w", stdout);
    
    cvNamedWindow("MST");
    cvNamedWindow("MST Gray");
    
    imgSrc = cvLoadImage("florencia.jpg");
    imgDst = cvCreateImage(cvSize(imgSrc->width, imgSrc->height), imgSrc->depth, imgSrc->nChannels);
    imgGray = cvCreateImage(cvSize(imgSrc->width, imgSrc->height), imgSrc->depth, imgSrc->nChannels);
    
    nRows = imgSrc->height;
    nColumns = imgSrc->width;
    
    X = ConvImgToDbl(imgSrc);
    imgI = ConvImgToDbl(imgSrc);
    
    ConvDblToImg(X, imgGray, true);
    
    for(i=0; i<nRows; i++)
    {
        for(j=0; j<nColumns; j++)
        {
            if(j == 0)
                printf("%.0lf", X[i][j].val[INTENSITY]);
            else
                printf(" %.0lf", X[i][j].val[INTENSITY]);
        }
        printf("\n");
    }
    
    
    Smooth(X, imgI, INTENSITY);
    
    for(i=0; i<nRows; i++)
    {
        for(j=0; j<nColumns; j++)
        {
            if(j == 0)
                printf("%.2lf", imgI[i][j].val[INTENSITY]);
            else
                printf(" %.2lf", imgI[i][j].val[INTENSITY]);
        }
        printf("\n");
    }

    
    MST_Segmentation(imgI, INTENSITY);
    
    ConvDblToImg(imgI, imgDst, false);
    
    cvShowImage("MST", imgDst);
    cvSaveImage("mst_tiger_2.jpg", imgDst);
    
    cvShowImage("MST Gray", imgGray);
    //cvSaveImage("mstG3.jpg", imgGray);
    
    
    cvWaitKey(0);
    
    cvReleaseImage(&imgSrc);
    cvReleaseImage(&imgDst);
    cvReleaseImage(&imgGray);
    
    cvDestroyWindow("MST");
    cvDestroyWindow("MST Gray");
    
    return 0;
}

void MST_Segmentation(Image **I, int channel)
{
    int i, j, k, r, c;
    int u, v, set1, set2, newSet;
    int C1, C2;
    int nEdges = 0;
    int r1, c1, r2, c2;
    double Int1, Int2;
    double MInt, cost;
    
    E = new Edge[5*nRows*nColumns];
    p = new int[5*nRows*nColumns];
    rango = new int[5*nRows*nColumns];
    Int = new double[5*nRows*nColumns];
    C = new int[5*nRows*nColumns];
    COL = new int[5*nRows*nColumns];
    
    for(i=0; i<nRows; i++)
    {
        for(j=0; j<nColumns; j++)
        {
            k = i*nColumns + j;
            
            make_set(k);
            
            COL[k] = -1;
            Int[k] = 0.0;
            C[k] = 1;
            
            for(r=i-1; r<=i; r++)
            {
                for(c=j-1; c<=j+1; c++)
                {
                    if(r == i && c == j)
                        continue;
                    
                    if(r < 0 || c < 0 || c >= nColumns || r >= nRows)
                        continue;
                    
                    E[nEdges].v1 = i*nColumns + j;
                    E[nEdges].v2 = r*nColumns + c;
                    E[nEdges].w = fabs(I[i][j].val[channel] - I[r][c].val[channel]);
                    
                    nEdges++;
                 }
            }
        }
    }
    
    qsort(E, nEdges, sizeof(E[0]), compare);
    
    //printf("nEdges: %d\n", nEdges);
    
    for(i=0; i<nEdges; i++)
    {
        u = E[i].v1;
        v = E[i].v2;
        cost = E[i].w;
        
        set1 = find_set(u);
        set2 = find_set(v);
        Int1 = Int[set1];
        Int2 = Int[set2];
        C1 = C[set1];
        C2 = C[set2];
        
        MInt = Min(Int1 + thao/C1, Int2 + thao/C2);
        
        if(set1 != set2 && cost < MInt)
        {
            union_set(u, v);
            
            newSet = find_set(u);
            
            Int[newSet] = Max(Int1, Int2);
            Int[newSet] = Max(Int[newSet], cost);
            C[newSet] = C1 + C2;
        }
    }
    
    for(i=0; i<nRows*nColumns; i++)
    {
        r1 = i/nColumns;
        c1 = i%nColumns;
        
        k = find_set(i);
        
        if(COL[k] == -1)
        {
            COL[k] = i;
            
            I[r1][c1].val[RED] = rand()%255;
            I[r1][c1].val[GREEN] = rand()%255;
            I[r1][c1].val[BLUE] = rand()%255;
            I[r1][c1].val[INTENSITY] = (I[r1][c1].val[RED] + I[r1][c1].val[GREEN] +I[r1][c1].val[BLUE])/3.0;
        }
        else 
        {
            r2 = COL[k]/nColumns;
            c2 = COL[k]%nColumns;
            
            I[r1][c1].val[RED] = I[r2][c2].val[RED];
            I[r1][c1].val[GREEN] = I[r2][c2].val[GREEN];
            I[r1][c1].val[BLUE] = I[r2][c2].val[BLUE];
            I[r1][c1].val[INTENSITY] = I[r2][c2].val[INTENSITY];       
        }
        
    }
}

void Smooth(Image **S, Image **D, int channel)
{
    int i, j;
    int r, c;
    double sum, k;
    
    double M[5][5] = {
        {2, 4, 5, 4, 2}, 
        {4, 9, 12, 9, 4},
        {5, 12, 15, 12, 5},
        {4, 9, 12, 9, 4},
        {2, 4, 5, 4, 2}
    };
    
    for(i=0; i<nRows; i++)
    {
        for(j=0; j<nColumns; j++)
        {
            k = 0.0;
            sum = 0.0;
            
            for(r=i-2; r<=i+2; r++)
            {
                for(c=j-2; c<=j+2; c++)
                {
                    if(r < 0 || r >= nRows || c < 0 || c >= nColumns)
                        continue;
                    
                    k += S[r][c].val[channel]*M[r-i+2][c-j+2];
                    sum += M[r-i+2][c-j+2];
                }
            }
            
            D[i][j].val[channel] = k/sum;
        }
    }
}

Image ** ConvImgToDbl(IplImage *img)
{
    CvScalar color;
    Image **T;
    int i, j;
    
    T = new Image *[img->height];
    for(i=0; i<img->height; i++)
    {
        T[i] = new Image[img->width];
        for(j=0; j<img->width; j++)
        {
            color = cvGet2D(img, i, j);
            
            T[i][j].val[RED] = color.val[RED];
            T[i][j].val[GREEN] = color.val[GREEN];
            T[i][j].val[BLUE] = color.val[2];
            T[i][j].val[INTENSITY] = (color.val[RED] + color.val[GREEN] + color.val[BLUE])/3.0;
        }
    }
    
    return T;
}

void ConvDblToImg(Image **I, IplImage *img, bool grayScale)
{
    int i, j;
    int r, g, b;
    
    for(i=0; i<nRows; i++)
    {
        for(j=0; j<nColumns; j++)
        {
            if(grayScale == true)
                r = g = b = (int)I[i][j].val[INTENSITY];
            else 
            {
                r = (int)I[i][j].val[RED];
                g = (int)I[i][j].val[GREEN];
                b = (int)I[i][j].val[BLUE];
            }
            
            cvSet2D(img, i, j, cvScalar(r, g, b));
        }
    }
}

void make_set(int x)
{
    p[x] = x;
    rango[x] = 0;
}

void link(int x, int y)
{
    if (rango[x] > rango[y])
        p[y] = x;
    else
    {
        p[x] = y;
        if (rango[x] == rango[y])
            rango[y] = rango[y] + 1;
    }
}

int find_set(int x)
{
    if (x != p[x])
        p[x] = find_set(p[x]);
    return p[x];
}

void union_set(int x, int y)
{
    link(find_set(x), find_set(y));
}

int Min(double a, double b)
{
    if(a < b)
        return a;
    else 
        return b;
}

int Max(double a, double b)
{
    if(a > b)
        return a;
    else 
        return b;
}

int compare(const void *a, const void *b)
{
    Edge *sp1 = (Edge *)a;
    Edge *sp2 = (Edge *)b;
    
    if(sp1->w < sp2->w)
        return -1;
    else if(sp1->w > sp2->w)
        return 1;
    
    return 0;
}
