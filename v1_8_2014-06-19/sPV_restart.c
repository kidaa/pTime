#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "fileIO.h"
#include "allocRoutine.h"
#include "mathOps.h"

#define MAXLEN 120
#define LEN 20

// LAST UPDATED: Apr 18, 2014 1:01PM PDT
// Edition according to my framework in Evernote!!
// Removal of "writeTraj()", and output to step 1 instead of step 0

int main(int argc, char *argv[])
{
    // ASSUME option -N to read for numPart, and argv[2] contains the value
    double dt,massFac,*coord,*vel,*velNext,*fPrev,*fNext,*mass,*massKg,kesum,keprevsum,keave,sum,Tave,TTarget;
    int *ID,numPart,i,ind,optInd,optCount,isLastStep,TOptionIndex,isIni;
    char *pathTrajtIn,*pathTrajtOut,*pathForcePrevIn,*pathForceNextIn,*pathNWOut,*pathForceOut,*tTrajtIn,*tTrajtOut,*tForcePrevIn,*tForceNextIn,*tForceOut,*tNWOut,**name,*pathTempOut;
    char *dirPrefix,*pathTrajOut,*pathPropOut /* SAVED FOR LATER VERSIONS TO OUTPUT PROPERTIES!! */;
    const double toMe = 1.6605e-27/9.1094e-31;	// [MAY 27, 2014] - Atomic unit conversion factor
    const double a0 = 0.592;
    const double v_au = 2.1877e+6;
    const double kb = 3.1668e-6;		// [MAY 27, 2014] - Atomic unit Kb
    FILE *writeTemp;

    // argv[2] = -N $numPart, argv[4] = -P $pathpref, argv[6] = -i $i, argv[8] -K $kesum, argv[10] -T Ttarget
    // CHECKING OPTIONS!!
    optInd = 1;
    optCount = 0;
    while(optInd < argc)
    {
        if(strcmp(argv[optInd],"-N")==0)
        {
            numPart = atoi(argv[optInd+1]);
            optInd += 2;
            optCount++;
        }
        // [JUNE 9, 2014] - Add the temperature option (Y/N)
        else if(strcmp(argv[optInd],"-T")==0)
        {
            TOptionIndex = optInd+1;
            optInd += 2;
            optCount++; 
        }
        // [JUNE 9, 2014] - Readd the target temperature option
        else if(strcmp(argv[optInd],"--TT")==0)
        {
            TTarget = (double)(atof(argv[optInd+1]));
            optInd += 2;
            optCount++;
        }
        else if(strcmp(argv[optInd],"-I")==0)
        {
            ind = atoi(argv[optInd+1]);
            optInd += 2;
            optCount++;
        }
        else if(strcmp(argv[optInd],"-K")==0)
        {
            keprevsum = (double)(atof(argv[optInd+1]));
            optInd += 2;
            optCount++;
        }
        // [MAY 27, 2014] - Add "--dt" option
        else if(strcmp(argv[optInd],"--dt")==0)
        {
            dt = (double)(atof(argv[optInd+1]));
            optInd += 2;
            optCount++;
        }
        // [JUNE 16, 2014] - Add "-b" (boolean) option - takes only 0 or 1
        // 0 means not initial position, 1 means is so
        else if(strcmp(argv[optInd],"-b")==0)
        {
            isIni = atoi(argv[optInd+1]);
            optInd += 2;
            optCount++;
        }
        else if(strcmp(argv[optInd],"-l")==0)
        {
            if((strcmp(argv[optInd+1],"True")==0)||(strcmp(argv[optInd+1],"true")==0))
                isLastStep = 1;
            else if((strcmp(argv[optInd+1],"False")==0)||(strcmp(argv[optInd+1],"false")==0))
                isLastStep = 0;
            else
            {
                printf("Invalid option for -l tag!! Exiting the program...\n");
                exit(1);
            }
            optInd += 2;
            optCount++;
        }
        else if(strcmp(argv[optInd],"-P")==0)
        {
            dirPrefix = allocChar1D(MAXLEN);
            strcpy(dirPrefix,argv[optInd+1]);
            optInd += 2;
            optCount++;
        }
        else if(strcmp(argv[optInd],"--help")==0)
        {
            system("cat /home/wpornpat/PTime/src/docs/sPV_AIMD.txt");
            exit(1);
        }
        else
        {
            if(!optCount)
            {
                printf("Invalid options!\n");
                exit(1);
            }
            else
                break;
        }
    }

    velNext = alloc1D(3*numPart);

    pathTrajtIn = allocChar1D(MAXLEN);
    pathTrajtOut = allocChar1D(MAXLEN);
    pathTrajOut = allocChar1D(MAXLEN);
    pathForcePrevIn = allocChar1D(MAXLEN);
    pathForceNextIn = allocChar1D(MAXLEN);
    pathForceOut = allocChar1D(MAXLEN);
    pathNWOut = allocChar1D(MAXLEN);
    tTrajtIn = allocChar1D(LEN);
    tTrajtOut = allocChar1D(LEN);
    tForcePrevIn = allocChar1D(LEN);
    tForceNextIn = allocChar1D(LEN);
    tForceOut = allocChar1D(LEN);
    tNWOut = allocChar1D(LEN);

    sprintf(tTrajtIn,"step%d_2.trajt",ind);
    sprintf(tTrajtOut,"step%d_1.trajt",ind+1);
    sprintf(tForceNextIn,"step%d_2.force",ind);
    if(isIni)
        sprintf(tForcePrevIn,"step%d_1.force",ind);
    else
        sprintf(tForcePrevIn,"step%d_2.force",ind-1);
    //sprintf(tNWOut,"step%d_1.nw",ind+1);
    //sprintf(tForceOut,"step%d_1.force",ind+1);

    strcpy(pathTrajtIn,dirPrefix);
    strcpy(pathTrajtOut,dirPrefix);
    strcpy(pathForcePrevIn,dirPrefix);
    strcpy(pathForceNextIn,dirPrefix);
    strcpy(pathForceOut,dirPrefix);
    strcpy(pathNWOut,dirPrefix);
    strcpy(pathTrajOut,dirPrefix);

    strcat(pathTrajtIn,tTrajtIn);
    strcat(pathTrajtOut,tTrajtOut);
    strcat(pathForcePrevIn,tForcePrevIn);
    strcat(pathForceNextIn,tForceNextIn);
    strcat(pathForceOut,tForceOut);
    strcat(pathNWOut,tNWOut);
    strcat(pathTrajOut,"trajectory.traj");

    free(tTrajtIn);
    free(tTrajtOut);
    free(tForcePrevIn);
    free(tForceNextIn);
    free(tNWOut);
    free(tForceOut);

    // READ FROM .TRAJT FILE
    readTrajt(&coord,&ID,&mass,&vel,numPart,pathTrajtIn);

    // READ FORCES FROM INPUT
    readForce(&fPrev,ID,numPart,pathForcePrevIn);
    readForce(&fNext,ID,numPart,pathForceNextIn);

    // PROPAGATION
    // [JUNE 9, 2014] - Everything already in atomic units. 
    sum = 0.0;
    if((ind > 0) && (ind % 3 == 0))
    {
        keave = keprevsum/((double)(ind));
        Tave = 2*keave/(kb*(3*numPart-6));
    }
    for(i = 0; i < numPart; i++)
    {
        massFac = mass[i]*toMe;
        velNext[3*i] = propagateV_au(vel[3*i],fPrev[3*i]/massFac,fNext[3*i]/massFac,dt);
        velNext[3*i+1] = propagateV_au(vel[3*i+1],fPrev[3*i+1]/massFac,fNext[3*i+1]/massFac,dt);
        velNext[3*i+2] = propagateV_au(vel[3*i+2],fPrev[3*i+2]/massFac,fNext[3*i+2]/massFac,dt);
        if((ind > 0) && (ind % 3 == 0))
        {
            velNext[3*i] *= sqrt(TTarget/Tave);
            velNext[3*i+1] *= sqrt(TTarget/Tave);
            velNext[3*i+2] *= sqrt(TTarget/Tave);
        }
        sum += 0.5*massFac*(velNext[3*i]*velNext[3*i]+velNext[3*i+1]*velNext[3*i+1]+velNext[3*i+2]*velNext[3*i+2]);
    }

    // CALCULATE THE AVERAGE TEMPERATURE
    kesum = keprevsum + sum;
    keave = kesum/(ind+1);
    Tave = 2*keave/(kb*(3*numPart-6));
   
    if((strcmp(argv[TOptionIndex],"Y")==0)||(strcmp(argv[TOptionIndex],"y")==0))
    {
        pathTempOut = allocChar1D(MAXLEN);
        strcpy(pathTempOut,dirPrefix);
        strcat(pathTempOut,"temperature.dat");
        writeTemp = fopen(pathTempOut,"a");
        fprintf(writeTemp,"%d %.2lf\n",ind+1,Tave);
        fclose(writeTemp);
        free(pathTempOut);
    }
 
    // OUTPUTS (NO NEED TO WRITE .NW, BECAUSE IT IS ALREADY UPDATED TO THE CURRENT COORDINATES)
    //IDtoName(&name,mass,numPart);
    writeTrajt(coord,ID,mass,velNext,numPart,pathTrajtOut);
    //writeNW(coord,name,numPart,pathNWOut,pathForceOut,dirPrefix);
    writeTrajOne(coord,velNext,mass,ID,ind,numPart,pathTrajOut);

    printf("%.6e",kesum);

    free1D(&coord);
    free1D(&vel);
    free1D(&velNext);
    free1D(&fPrev);
    free1D(&fNext);
    free1D(&mass);
    freeInt1D(&ID);
    free(pathTrajtIn);
    free(pathTrajtOut);
    free(pathForcePrevIn);
    free(pathForceNextIn);
    free(pathNWOut);
    free(pathForceOut);
    free(pathTrajOut);
    free(dirPrefix);

    return 0;
}
