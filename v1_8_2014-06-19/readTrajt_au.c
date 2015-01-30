#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fileIO.h"
#include "allocRoutine.h"

#define MAXLEN 120
#define LEN 20

// READ FROM INITIAL .TRAJT FILE AND OUTPUT NWCHEM INPUT FILE FOR THE INPUT OF SERIALPROPAGATEX
// UPDATED: May 26, 2014 --> CONVERSION TO ATOMIC UNITS

int main(int argc, char *argv[])
{
    double *coord,*mass,*vel;
    int i,*ID,numPart,optInd,optCount,dtPathIndex;
    char **name,*path,*pathNW,*pathF,*mainDir;
    char *tPath,*tNW,*tF;

    int ind = 0;

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
        // [JUNE 16, 2014] - Add "-I" option
        else if(strcmp(argv[optInd],"-I")==0)
        {
            ind = atoi(argv[optInd+1]);
            optInd += 2;
            optCount++;
        }
        // [JUNE 29, 2014] - Add "-D" option for /dtemp path prefix
        else if(strcmp(argv[optInd],"-D")==0)
        {
            dtPathIndex = optInd + 1;
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
    tPath = allocChar1D(LEN);
    tNW = allocChar1D(LEN);
    tF = allocChar1D(LEN);

    strcpy(path,mainDir);
    strcpy(pathNW,mainDir);
    strcpy(pathF,mainDir);

    sprintf(tPath,"step%d_1.trajt",ind);
    sprintf(tNW,"step%d_1.nw",ind);
    sprintf(tF,"step%d_1.force",ind);

    strcat(path,tPath);
    strcat(pathNW,tNW);
    strcat(pathF,tF);

    free(tPath);
    free(tNW);
    free(tF);

    readTrajt(&coord,&ID,&mass,&vel,numPart,path);

    // [MAY 28, 2014] - Update coord to atomic unit and invoke writeTrajt function to overwrite the configuration.
    //for(i = 0; i < 3*numPart; i++)
    //{
    //    coord[i] /= 0.529;
    //    vel[i] /= 2.18769e+6;
    //}
    
    //writeTrajt(coord,ID,mass,vel,numPart,path);

    IDtoName(&name,mass,numPart);

    // WRITE TO .TRAJT FILE AS WELL AS WRITE THE .NW NWCHEM INPUT FILE
    writeNW(coord,name,numPart,pathNW,pathF,mainDir,argv[dtPathIndex]); 
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
