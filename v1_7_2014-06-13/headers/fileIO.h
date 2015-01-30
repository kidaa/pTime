#ifndef FILEIO_TEMP_H
#define FILEIO_TEMP_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern char elemName[20][3];
extern double elemMass[20];

void readXYZ(double **coord, double **mass, int **index, int *numPart, const char *path);
void writeNW(double *coord, char **name, int numPart, char *pathNW, char *pathF, char *pathPrefix);				// IF WRITING FROM PTIME.C, CONVERT X FIRST!!
void IDtoName(char ***name, double *mass, int numPart);
void writeTrajt(double *coord, int *ID, double *mass, double *vel, int numPart, const char *path);				// IF WRITING FROM PTIME.C, CONVERT X FIRST!!
void readTrajt(double **coord, int **ID, double **mass, double **vel, int numPart, const char *path);
void readTrajtPTime(double *coord, int *ID, double *mass, double *vel, int numPart, const char *path, int index);
void readForce(double **force, int *ID /*for validation purpose*/, int numPart, const char *path);
void readForcePTime(double *force, int *ID, int numPart, const char *path);
void writeHist(double **rHist, double **aHist, int histNum, int kValue, char *path, int M, int numPart);
void writeTraj(double *XNext, double *mass, int *ID, int IValue, int M, int numPart, char *path);	// EDIT THIS TO "w" option in file writing rather than "a". Updated April 30, 2014
void writeTrajOne(double *x, double *v, double *mass, int *ID, int iVal, int numPart, char *path);
void readHist(double **rHist, double **aHist, int M, int numPart, char *path, int histNum);

#endif
