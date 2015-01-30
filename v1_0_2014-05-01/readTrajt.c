#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fileIO.h"
#include "allocRoutine.h"

#define MAXLEN 120
#define LEN 20

// READ FROM INITIAL .TRAJT FILE AND OUTPUT NWCHEM INPUT FILE FOR THE INPUT OF SERIALPROPAGATEX

int main(int argc, char *argv[])
{
    double *coord,*mass,*vel;
    int *ID,numPart,optInd,optCount;
    char **name,*path,*pathNW,*pathF,*mainDir;

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
        else if(strcmp(argv[optInd],"-N")==0)
        {
            numPart = atoi(argv[optInd+1]);
            optInd += 2;
            optCount++;
        }
        else if(strcmp(argv[optInd],"--help")==0)
        {
            system("cat /home/wpornpat/PTime/src/docs/readXYZ.txt");		// UPDATE PATH TO DOCUMENTATION!!
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

    path = allocChar1D(MAXLEN);
    pathNW = allocChar1D(MAXLEN);
    pathF = allocChar1D(MAXLEN);

    strcpy(path,mainDir);
    strcpy(pathNW,mainDir);
    strcpy(pathF,mainDir);

    strcat(path,"step0_1.trajt");
    strcat(pathNW,"step0_1.nw");
    strcat(pathF,"step0_1.force");

    readTrajt(&coord,&ID,&mass,&vel,numPart,path);

    IDtoName(&name,mass,numPart);

    // WRITE TO .TRAJT FILE AS WELL AS WRITE THE .NW NWCHEM INPUT FILE
    writeNW(coord,name,numPart,pathNW,pathF,mainDir); 
    free(mainDir);
    free1D(&coord);
    free1D(&mass);
    free1D(&vel);
    freeInt1D(&ID);
    freeChar2D(name,numPart);
    free(path);
    free(pathNW);
    free(pathF);

    return 0;
}
