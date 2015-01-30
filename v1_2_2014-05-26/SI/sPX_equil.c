#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fileIO.h"
#include "allocRoutine.h"
#include "mathOps.h"

#define MAXLEN 120
#define LEN 20

int main(int argc, char *argv[])
{
    // ASSUME -N option for numpart, and argv[2] is number of particle
    double *coord,*vel,*mass,*massKg,*force,*coordNext,dt;
    int *ID,numPart,i,ind,optInd,optCount,dirPrefixIndex;
    char **name;
    char *pathTrajtIn,*pathTrajtOut,*pathForceIn,*pathNWOut,*pathForceOut;
    char *tTrajtIn,*tTrajtOut,*tForceIn,*tNWOut,*tForceOut;

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
        else if(strcmp(argv[optInd],"--dt")==0)
        {
            dt = (double)(atof(argv[optInd+1]));
            optInd += 2;
            optCount++;
        }
        else if(strcmp(argv[optInd],"-i")==0)
        {
            ind = atoi(argv[optInd+1]);
            optInd += 2;
            optCount++;
        }
        else if(strcmp(argv[optInd],"-P")==0)
        {
            dirPrefixIndex = optInd+1;
            optInd += 2;
            optCount++;
        }
        else if(strcmp(argv[optInd],"--help")==0)
        {
            system("cat /home/wpornpat/PTime/src/docs/sPX_AIMD.txt");
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

    pathTrajtIn = allocChar1D(MAXLEN);
    pathTrajtOut = allocChar1D(MAXLEN);
    pathForceIn = allocChar1D(MAXLEN);
    pathNWOut = allocChar1D(MAXLEN);
    pathForceOut = allocChar1D(MAXLEN);
    tTrajtIn = allocChar1D(LEN);
    tTrajtOut = allocChar1D(LEN);
    tForceIn = allocChar1D(LEN);
    tNWOut = allocChar1D(LEN);
    tForceOut = allocChar1D(LEN);

    sprintf(tTrajtIn,"step%d_1.trajt",ind);
    sprintf(tTrajtOut,"step%d_2.trajt",ind);
    sprintf(tForceIn,"step%d_1.force",ind);
    sprintf(tNWOut,"step%d_2.nw",ind);
    sprintf(tForceOut,"step%d_2.force",ind);

    strcpy(pathTrajtIn,argv[dirPrefixIndex]);
    strcpy(pathTrajtOut,argv[dirPrefixIndex]);
    strcpy(pathForceIn,argv[dirPrefixIndex]);
    strcpy(pathNWOut,argv[dirPrefixIndex]);
    strcpy(pathForceOut,argv[dirPrefixIndex]);

    strcat(pathTrajtIn,tTrajtIn);
    strcat(pathTrajtOut,tTrajtOut);
    strcat(pathForceIn,tForceIn);
    strcat(pathNWOut,tNWOut);
    strcat(pathForceOut,tForceOut);

    free(tTrajtIn);
    free(tTrajtOut);
    free(tForceIn);
    free(tNWOut);
    free(tForceOut);

    const double auToKg = 1.6605e-27;

    numPart = atoi(argv[2]);
    massKg = alloc1D(numPart);
    coordNext = alloc1D(3*numPart);

    // READ FROM .TRAJT FILE
    readTrajt(&coord,&ID,&mass,&vel,numPart,pathTrajtIn);
    for(i = 0; i < numPart; i++)
        massKg[i] = mass[i]*auToKg;

    // READ FORCE FROM INPUT
    readForce(&force,ID,numPart,pathForceIn);

    // PROPAGATION
    for(i = 0; i < numPart; i++)
    {
        coordNext[3*i] = propagateX(coord[3*i],vel[3*i],force[3*i]/massKg[i],dt);
        coordNext[3*i+1] = propagateX(coord[3*i+1],vel[3*i+1],force[3*i+1]/massKg[i],dt);
        coordNext[3*i+2] = propagateX(coord[3*i+2],vel[3*i+2],force[3*i+2]/massKg[i],dt);
    }

    // OUTPUT TO ANOTHER .TRAJT FILE AND TO .NW FILE
    IDtoName(&name,mass,numPart);
    writeNW(coordNext,name,numPart,pathNWOut,pathForceOut,argv[dirPrefixIndex]);
    writeTrajt(coordNext,ID,mass,vel,numPart,pathTrajtOut);

    free1D(&coord);
    free1D(&coordNext);
    free1D(&mass);
    free1D(&massKg);
    free1D(&vel);
    free1D(&force);
    freeInt1D(&ID);
    freeChar2D(name,numPart);
    free(pathTrajtIn);
    free(pathTrajtOut);
    free(pathForceIn);
    free(pathForceOut);
    free(pathNWOut);
    return 0;
}
