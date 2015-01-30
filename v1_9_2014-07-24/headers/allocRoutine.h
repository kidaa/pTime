#ifndef ALLOCROUTINE_H
#define ALLOCROUTINE_H

#include <stdio.h>
#include <stdlib.h>

double *alloc1D(int size);
int *allocInt1D(int size);
double **alloc2D(int dimX, int dimY);
char **allocChar2D(int numStr, int lengthStr);
char *allocChar1D(int size);
void free1D(double **arr);
void freeInt1D(int **arr);
void free2D(double **arr, int dimX);
void freeChar2D(char **arr, int dimX);
void freeChar1D(char **arr);
void setZero1D(double *arr, int N);
void setZero2D(double **arr, int dimX, int dimY);

#endif
