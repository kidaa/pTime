#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fileIO.h"
#include "allocRoutine.h"

#define MAXLEN 120
#define LEN 20

// READ FROM INITIAL XYZ FILE AND OUTPUT NWCHEM INPUT FILE AS WELL AS .TRAJT FORMAT FILE FOR THE INPUT OF SERIALPROPAGATEX

int main(int argc, char *argv[])
{
    const double auToKg = 1.6605e-27;
    double *coord,*mass,*vel;
    int *ID,numPart,i,optInd,optCount;
    char **name,*path,*pathNW,*pathF,*pathTrajt,*mainDir,*xyzName;

    // CHECKING OPTIONS - CALL A FUNCTION IN "fileIO" ROUTINES!!
    optInd = 1;
    optCount = 0;
    while(optInd < argc)
    {
        if(strcmp(argv[optInd],"-P")==0)
        {
            mainDir = allocChar1D(MAXLEN);
            strcpy(mainDir,argv[optInd+1]);
            optInd += 2;
            optCount++;
        }
        else if(strcmp(argv[optInd],"-x")==0)
        {
            xyzName = allocChar1D(LEN);
            strcpy(xyzName,argv[optInd+1]);
            optInd += 2;
            optCount++;
        }
        else if(strcmp(argv[optInd],"--help")==0)
        {
            system("cat /home/wpornpat/PTime/src/docs/readXYZ.txt");
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

    // FIX PATTERN!! argv[2] = dir_prefix (-p), argv[4] = xyz_filename (-x)
    path = allocChar1D(MAXLEN);
    pathNW = allocChar1D(MAXLEN);
    pathF = allocChar1D(MAXLEN);
    pathTrajt = allocChar1D(MAXLEN);

    strcpy(path,mainDir);
    strcpy(pathNW,mainDir);
    strcpy(pathF,mainDir);
    strcpy(pathTrajt,mainDir);

    strcat(path,xyzName);
    freeChar1D(&xyzName);
    strcat(pathNW,"step0_1.nw");
    strcat(pathF,"step0_1.force");
    strcat(pathTrajt,"step0_1.trajt");

    readXYZ(&coord,&mass,&ID,&numPart,path);
    vel = alloc1D(3*numPart);
    setZero1D(vel,3*numPart); 		// IN CASE VELOCITIES ARE NOT READ FROM XYZ FILE

    IDtoName(&name,mass,numPart);

    // WRITE TO .TRAJT FILE AS WELL AS WRITE THE .NW NWCHEM INPUT FILE
    writeNW(coord,name,numPart,pathNW,pathF,mainDir); 
    writeTrajt(coord,ID,mass,vel,numPart,pathTrajt);
    printf("%d",numPart);		// FOR RETURNING NUMPART TO BASH ENV VARIABLE
    free1D(&coord);
    freeChar1D(&mainDir);
    free1D(&mass);
    free1D(&vel);
    freeInt1D(&ID);
    freeChar2D(name,numPart);
    free(path);
    free(pathNW);
    free(pathF);
    free(pathTrajt);

    return 0;
}
