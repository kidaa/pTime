#include <stdio.h>
#include <stdlib.h>

#include "allocRoutine.h"

double *alloc1D(int size)
{
    double *arr;
    arr = (double*)malloc(size*sizeof(double));
    return arr;
}

int *allocInt1D(int size)
{
    int *arr;
    arr = (int*)malloc(size*sizeof(int));
    return arr;
}

double **alloc2D(int dimX, int dimY)
{
    double **arr;
    int i;
    arr = (double**)malloc(dimX*sizeof(double*));
    for(i = 0; i < dimX; i++)
        arr[i] = (double*)malloc(dimY*sizeof(double));
    return arr;
}

char **allocChar2D(int numStr, int lengthStr)
{
    char **arr;
    int i;
    arr = (char**)malloc(numStr*sizeof(char*));
    for(i = 0; i < numStr; i++)
        arr[i] = (char*)malloc(lengthStr*sizeof(char));
    return arr;
}

char *allocChar1D(int size)
{
    char *arr;
    arr = (char*)malloc(size*sizeof(char));
    return arr;
}

void free1D(double **arr)
{
    free(*arr);
}

void freeInt1D(int **arr)
{
    free(*arr);
}

void free2D(double **arr, int dimX)
{
    int i;
    for(i = 0; i < dimX; i++)
        free(arr[i]);
    free(arr);
}

void freeChar2D(char **arr, int dimX)
{
    int i;
    for(i = 0; i < dimX; i++)
        free(arr[i]);
    free(arr);
}

void freeChar1D(char **arr)
{
    free(*arr);
}

void setZero1D(double *arr, int N)
{
    int i;
    for(i = 0; i < N; i++)
        arr[i] = 0.0;
}

void setZero2D(double **arr, int dimX, int dimY)
{
    int i,j;
    for(i = 0; i < dimX; i++)
    {
        for(j = 0; j < dimY; j++)
            arr[i][j] = 0.0;
    }
}
