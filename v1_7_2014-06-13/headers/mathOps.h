#ifndef MATHOPS_TEMP_H
#define MATHOPS_TEMP_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double propagateX(double coord, double vel, double a, double dt);
double propagateV(double vel, double aPrev, double aNext, double dt);
double propagateX_au(double coord, double vel, double a, double dt);
double propagateV_au(double vel, double aPrev, double aNext, double dt);
void calcTilde(double *Fr, double *Fv, double *rTilde, double *vTilde, double **rHist, double **aHist, int M, int numPart, int histNum, int kValue, double dt);
void calcXtR(double *rTilde, double *vTilde, double *Fr, double *Hr, int ind, int numPart, double dt);
void calcXtV(double *vTilde, double *Fv, double *Hr, int ind, int numPart, double dt);
void calcHr(double *Hr, double **rHist, double **aHist, double *rTilde, int numPart, int histNum, int kValue, int ind);
double vecDot(double *v1, double *v2, int size);
double vecNormSq(double *v, int size);
void XNext(double *XNext, double *rTilde, double *vTilde, double *rNormal, double *vNormal, int numPart, int M);

#endif
