#include <stdio.h>
#include <stdlib.h>

#include "arrayManip.h"

void sliceX(double *XNext, double *tcoord, double *tvel, int index, int numPart)
{
    int i;
    for(i = 0; i < numPart; i++)
    {
        tcoord[3*i] = XNext[(2*index*3*numPart) + (3*i)];
        tcoord[3*i+1] = XNext[(2*index*3*numPart) + (3*i+1)];
        tcoord[3*i+2] = XNext[(2*index*3*numPart) + (3*i+2)];
        tvel[3*i] = XNext[((2*index+1)*3*numPart) + (3*i)];
        tvel[3*i+1] = XNext[((2*index+1)*3*numPart) + (3*i+1)];
        tvel[3*i+2] = XNext[((2*index+1)*3*numPart) + (3*i+2)];
    }
}

void storeHist(double **rHist, double **aHist, double *x, double *x0, double *f, double *f0, double *mass, const double auToKg, const double angToM, int M, int h, int N, int k, char *path)
{
    int j;
    if(k == 1)
    {
        if(h < 2)
        {
            printf("histNum is less than 2: Impossible to implement the algorithm!");
            exit(1);
        }
        else
        {
            // MAY 21, 2014 Experiment...
            // NO f0 STORAGE!! LEAVE IT AS ZERO!!
            for(j = 0; j < 3*M*N; j++)
            {
                // CONVERT EVERYTHING TO SI UNITS RIGHT HERE!!
                aHist[0][j] = f[j]/(mass[(j%(3*N))/3]*auToKg);				// CONVERSION IS NOT NECESSARY!! FORCES ARE ALREADY IN ATOMIC UNIT!!
                //aHist[1][j] = f0[(j>=3*N)?(j%(3*N)):j]/(mass[(j%(3*N))/3]*auToKg);
                rHist[0][j] = x[j]/angToM; 		// [MAY 26, 2014] --> NOW IN BOHR RADIUS FROM ANGSTROMS!!	
                //rHist[1][j] = x0[(j>=3*N)?(j%(3*N)):j]*angToM;
            }
        }
    }
    else
    {
        if(h < 2)
        {
            printf("histNum is less than 2: Impossible to implement the algorithm!");
            exit(1);
        }
        else
        {
            readHist(rHist,aHist,M,N,path,h);		// THE FUNCTION CALL AUTOMATICALLY UPDATES THE HISTORY!!
            for(j = 0; j < 3*M*N; j++)
            {
                aHist[0][j] = f[j]/(mass[(j%(3*N))/3]*auToKg);
                rHist[0][j] = x[j]/angToM;
            }
        }
    }
}
