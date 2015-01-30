#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mathOps.h"
#include "allocRoutine.h"

double propagateX(double coord, double vel, double a, double dt)
{
    double coordNext;
    coordNext = (coord*1.0e-10) + vel*dt + 0.5*dt*dt*a;
    coordNext /= (1.0e-10);
    return coordNext;
}

double propagateV(double vel, double aPrev, double aNext, double dt)
{
    double velNext;
    velNext = vel + 0.5*dt*(aPrev+aNext);
    return velNext;
}

double propagateX_au(double coord, double vel, double a, double dt)
{
    double coordNext;
    coordNext = coord + vel*dt + 0.5*dt*dt*a;
    return coordNext;
}

double propagateV_au(double vel, double aPrev, double aNext, double dt)
{
    double velNext;
    velNext = vel + 0.5*dt*(aPrev+aNext);
    return velNext;
}

void calcTilde(double *Fr, double *Fv, double *rTilde, double *vTilde, double **rHist, double **aHist, int M, int numPart, int histNum, int kValue, double dt)
{
    double *Hr;
    int i, j;

    Hr = alloc1D(3*M*numPart);
    setZero1D(Hr,3*M*numPart);
    
    for(i = 0; i < M; i++)
    {
        calcXtR(rTilde,vTilde,Fr,Hr,i,numPart,dt);
        calcHr(Hr,rHist,aHist,rTilde,numPart,histNum,kValue,i);
        calcXtV(vTilde,Fv,Hr,i,numPart,dt);    
    }

    free1D(&Hr); 
}

void calcXtR(double *rTilde, double *vTilde, double *Fr, double *Hr, int ind, int numPart, double dt)
{
    int i,midIndX,midIndY,midIndZ,midPrevX,midPrevY,midPrevZ;

    if(ind == 0)
    {
        for(i = 0; i < numPart; i++)
        {
            midIndX = 3*ind*numPart+(3*i);
            midIndY = 3*ind*numPart+(3*i+1);
            midIndZ = 3*ind*numPart+(3*i+2);
            rTilde[midIndX] = Fr[midIndX];		// set main unit system for rTilde to SI
            rTilde[midIndY] = Fr[midIndY];
            rTilde[midIndZ] = Fr[midIndZ];
        }
    }
    else
    {
        for(i = 0; i < numPart; i++)
        {
            midIndX = 3*ind*numPart+(3*i);
            midIndY = 3*ind*numPart+(3*i+1);
            midIndZ = 3*ind*numPart+(3*i+2);
            midPrevX = 3*(ind-1)*numPart+(3*i);
            midPrevY = 3*(ind-1)*numPart+(3*i+1);
            midPrevZ = 3*(ind-1)*numPart+(3*i+2);
            rTilde[midIndX] = Fr[midIndX] + rTilde[midPrevX] + vTilde[midPrevX]*dt + 0.5*dt*dt*Hr[midPrevX];
            rTilde[midIndY] = Fr[midIndY] + rTilde[midPrevY] + vTilde[midPrevY]*dt + 0.5*dt*dt*Hr[midPrevY];
            rTilde[midIndZ] = Fr[midIndZ] + rTilde[midPrevZ] + vTilde[midPrevZ]*dt + 0.5*dt*dt*Hr[midPrevZ];
        }
    }
}

void calcXtV(double *vTilde, double *Fv, double *Hr, int ind, int numPart, double dt)
{
    int i,midIndX,midIndY,midIndZ,midPrevX,midPrevY,midPrevZ;

    if(ind == 0)
    {
        for(i = 0; i < numPart; i++)
        {
            midIndX = 3*ind*numPart+(3*i);
            midIndY = 3*ind*numPart+(3*i+1);
            midIndZ = 3*ind*numPart+(3*i+2);
            vTilde[midIndX] = Fv[midIndX] + 0.5*Hr[midIndX]*dt;
            vTilde[midIndY] = Fv[midIndY] + 0.5*Hr[midIndY]*dt;
            vTilde[midIndZ] = Fv[midIndZ] + 0.5*Hr[midIndZ]*dt;
        }
    }
    else
    {
        for(i = 0; i < numPart; i++)
        {
            midIndX = 3*ind*numPart+(3*i);
            midIndY = 3*ind*numPart+(3*i+1);
            midIndZ = 3*ind*numPart+(3*i+2);
            midPrevX = 3*(ind-1)*numPart+(3*i);
            midPrevY = 3*(ind-1)*numPart+(3*i+1);
            midPrevZ = 3*(ind-1)*numPart+(3*i+2);
            vTilde[midIndX] = Fv[midIndX] + vTilde[midPrevX] + 0.5*dt*(Hr[midPrevX]+Hr[midIndX]);
            vTilde[midIndY] = Fv[midIndY] + vTilde[midPrevY] + 0.5*dt*(Hr[midPrevY]+Hr[midIndY]);
            vTilde[midIndZ] = Fv[midIndZ] + vTilde[midPrevZ] + 0.5*dt*(Hr[midPrevZ]+Hr[midIndZ]);
        }
    }
}

void calcHr(double *Hr, double **rHist, double **aHist, double *rTilde, int numPart, int histNum, int kValue, int ind)
{
    // [JUNE 5, 2014] - Add a boolean flag to check whether rho is less than our tolerance called "flag".
    int i,j,midIndX,midIndY,midIndZ;
    double *st,*stPrev,*rDiff,*aDiff;
    double gamma_k,rho_k;
    int lim;
    enum BOOLVAL {FALSE,TRUE} flag;

    // ALLOCATE ARRAYS
    st = alloc1D(3*numPart);
    stPrev = alloc1D(3*numPart);
    rDiff = alloc1D(3*numPart);
    aDiff = alloc1D(3*numPart);

    // MAY 21, 2014
    // CHECK IF THE TERNARY OPERATOR RETURNS THE CORRECT VALUE: If k = 1, should return 1...
    //printf("%d\n",((5 >= histNum-1)?histNum-1:5));
    //exit(1);
    lim = (kValue >= histNum-1)?(histNum-1):kValue;
    printf("%d\n",lim);
    
    // [JUNE 5, 2014] - Boolean "flag" initialized to "FALSE"
    flag = FALSE;
    for(i = 0; i < ((kValue >= histNum-1)?histNum-1:kValue); i++)
    {
        for(j = 0; j < 3*numPart; j++)
        {
            if(!flag)
                st[j] = rTilde[3*ind*numPart+j];
            else
                st[j] = stPrev[j] - (gamma_k/rho_k)*rDiff[j];
            //printf("{%.2e,%.2e} ",aHist[i][3*ind*numPart+j],aHist[i][3*ind*numPart+j]);
            rDiff[j] = rHist[i][3*ind*numPart+j] - rHist[i+1][3*ind*numPart+j]; 
            aDiff[j] = aHist[i][3*ind*numPart+j] - aHist[i+1][3*ind*numPart+j];
            printf("{%.2e,%.2e}, ",rDiff[j],st[j]);
        }
        //exit(1);
        rho_k = vecNormSq(rDiff,3*numPart);
        gamma_k = vecDot(rDiff,st,3*numPart);

        // [JUNE 5, 2014] - Check if rho_k is less than the tolerance and assign the corresponding boolean values.
        if(rho_k > 1.0e-15)
            flag = TRUE;
        else
            flag = FALSE;

        printf("\n{%.2e,%.2e,%d}\n",rho_k,gamma_k,i);
        //
        // [JUNE 5, 2014] - Update Hessian only when rho_k is greater than the tolerance.
        if(flag)
        {
            for(j = 0; j < 3*numPart; j++)
            {
                Hr[ind*3*numPart+j] += (gamma_k/rho_k)*aDiff[j];
                stPrev[j] = st[j];
            }
        }
    }

    // FREE ARRAYS
    free1D(&st);
    free1D(&stPrev);
    free1D(&rDiff);
    free1D(&aDiff);
}

double vecDot(double *v1, double *v2, int size)
{
    int i;
    double sum;
    sum = 0.0;
    for(i = 0; i < size; i++)
        sum += (v1[i]*v2[i]);
    return sum;
}

double vecNormSq(double *v, int size)
{
    int i;
    double sum;
    sum = 0.0;
    for(i = 0; i < size; i++)
        sum += (v[i]*v[i]);
    return sum;
}

void calcXNext(double *XNext, double *rTilde, double *vTilde, double *rNormal, double *vNormal, int numPart, int M)
{
    int i;
    //const double toMe = 0.529;
    //const double v_au = 2.1877e+6;
    // STORE EVERYTHING BACK TO THE SAME OLD FORMAT!!
    for(i = 0; i < 3*M*numPart; i++)
    {
        XNext[2*(i/(3*numPart))*3*numPart+(i%(3*numPart))] = (rNormal[i] - rTilde[i]);
        XNext[(2*(i/(3*numPart))+1)*3*numPart+(i%(3*numPart))] = (vNormal[i] - vTilde[i]);
    }
}
