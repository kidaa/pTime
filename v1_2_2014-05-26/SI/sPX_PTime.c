#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fileIO.h"
#include "allocRoutine.h"
#include "mathOps.h"

#define MAXLEN 120
#define LEN 50

int main(int argc, char *argv[])
{
    // ASSUME -N option for numpart, and argv[2] is number of particle
    double *coord,*vel,*mass,*massKg,*force,*coordNext,dt,massFac;
    int *ID,numPart,i,ind,optInd,optCount,dirPrefixIndex;
    char **name;
    char *pathTrajtIn,*pathTrajtOut,*pathForceIn,*pathNWOut,*pathForceOut;
    char *tTrajtIn,*tTrajtOut,*tForceIn,*tNWOut,*tForceOut;

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
        else if(strcmp(argv[optInd],"--ind")==0)
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
        else if(strcmp(argv[optInd],"--dt")==0)
        {
            dt = (double)(atof(argv[optInd+1]));
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

    // ASSUME argv[2] = -N numPart, argv[4] = dir_prefix (-P), argv[6] = -I loop index
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

    freeChar1D(&tTrajtIn);
    freeChar1D(&tTrajtOut);
    freeChar1D(&tForceIn);
    freeChar1D(&tNWOut);
    freeChar1D(&tForceOut);

    const double auToKg = 1.6605e-27;

    coordNext = alloc1D(3*numPart);

    // READ FROM .TRAJT FILE
    readTrajt(&coord,&ID,&mass,&vel,numPart,pathTrajtIn);

    // READ FORCE FROM INPUT
    readForce(&force,ID,numPart,pathForceIn);

    // PROPAGATION
    for(i = 0; i < numPart; i++)
    {
        massFac = mass[i]/auToKg;
        coordNext[3*i] = propagateX(coord[3*i],vel[3*i],force[3*i]/massFac,dt);
        coordNext[3*i+1] = propagateX(coord[3*i+1],vel[3*i+1],force[3*i+1]/massFac,dt);
        coordNext[3*i+2] = propagateX(coord[3*i+2],vel[3*i+2],force[3*i+2]/massFac,dt);
    }

    // OUTPUT TO ANOTHER .TRAJT FILE AND TO .NW FILE
    IDtoName(&name,mass,numPart);
    writeNW(coordNext,name,numPart,pathNWOut,pathForceOut,argv[dirPrefixIndex]);
    writeTrajt(coordNext,ID,mass,vel,numPart,pathTrajtOut);

    free1D(&coord);
    free1D(&coordNext);
    free1D(&mass);
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
