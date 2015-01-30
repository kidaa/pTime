#ifndef ARRAYMANIP_TEMP_H
#define ARRAYMANIP_TEMP_H

#include <stdio.h>
#include <stdlib.h>

void sliceX(double *XNext, double *tcoord, double *tvel, int index, int numPart);
void storeHist(double **rHist, double **aHist, double *x, double *x0, double *f, double *f0, double *mass, const double auToKg, const double angToM, int M, int h, int N, int k, char *path);

#endif
