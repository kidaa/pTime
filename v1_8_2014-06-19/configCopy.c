#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fileIO.h"
#include "allocRoutine.h"

#define MAXSTRLEN 100

// VERSION 1.1 - Last updated: April 29, 2014
//
// BUG FIXES:
//   - Command line option checking added
//   - Add -P option for users to enter a specific working directory. WILL ADD THE DEFAULT CASE LATER..

// SRC DIRECTORY FOR THIS VERSION:
// "/home/wpornpat/PTime/src/v1_1_2014-04-29/"

int main(int argc, char *argv[])
{
    double *coord,*vel,*mass;
    int *ID,i,j,numPart,IValue,M,optInd,optCount,dirPrefixIndex,dtPathIndex;
    char **name,**trajtFilenames,**nwFilenames,**forceFilenames,**tfn_t,**nfn_t,**ffn_t;
    char *pathTrajtIn,*pathNWOut,*pathForceOut,*tno,*tfo;
    char *tTrajtIn;

    // CHECKING OPTIONS
    // -P, -N, -M, -I
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
        // [JUNE 19, 2014] - Add "-D" option for dtemp
        else if(strcmp(argv[optInd],"-D")==0)
        {
            dtPathIndex = optInd+1;
            optInd += 2;
            optCount++;
        }
        else if(strcmp(argv[optInd],"-I")==0)
        {
            IValue = atoi(argv[optInd+1]);
            optInd += 2;
            optCount++;
        }
        else if(strcmp(argv[optInd],"-M")==0)
        {
            M = atoi(argv[optInd+1]);
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
            system("more /home/wpornpat/PTime/src/v1_1_2014-04-29/docs/configCopy.txt");
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

    // LAST UPDATED: May 1, 2014
    // Read from step{0,M,2M,...}_1.trajt instead of step{1,M+1,2M+1,...} to reflect the actual algorithm in the paper!
    // The initial file is step0_1.trajt

    pathTrajtIn = allocChar1D(MAXSTRLEN);
    tTrajtIn = allocChar1D(MAXSTRLEN);
    sprintf(tTrajtIn,"step%d_1.trajt",(M*IValue));
    strcpy(pathTrajtIn,argv[dirPrefixIndex]);
    strcat(pathTrajtIn,tTrajtIn);
    freeChar1D(&tTrajtIn);

    // LAST UPDATED: May 1, 2014
    // Creath paths for .nw and .force file for the initial pTime config.

    pathNWOut = allocChar1D(MAXSTRLEN);
    tno = allocChar1D(40);
    sprintf(tno,"step%d_1.nw",(M*IValue));
    strcpy(pathNWOut,argv[dirPrefixIndex]);
    strcat(pathNWOut,tno);
    freeChar1D(&tno);

    pathForceOut = allocChar1D(MAXSTRLEN);
    tfo = allocChar1D(40);
    sprintf(tfo,"step%d_1.force",(M*IValue));
    strcpy(pathForceOut,argv[dirPrefixIndex]);
    strcat(pathForceOut,tfo);
    freeChar1D(&tfo);

    // DYNAMIC ARRAYS ALLOCATION
    trajtFilenames = allocChar2D(M,MAXSTRLEN);
    nwFilenames = allocChar2D(M,MAXSTRLEN);
    forceFilenames = allocChar2D(M,MAXSTRLEN);
    tfn_t = allocChar2D(M,40);
    nfn_t = allocChar2D(M,40);
    ffn_t = allocChar2D(M,40);

    // READING FROM THE .TRAJT FILE
    readTrajt(&coord,&ID,&mass,&vel,numPart,pathTrajtIn);
    freeChar1D(&pathTrajtIn);

    // ASSIGN ELEMENT NAMES
    IDtoName(&name,mass,numPart);

    // WRITING COPIES. MAKE SURE ALL THESE DIRECTORIES ARE CREATED IN THE SHELL SCRIPT BEFOREHAND, ONLY AT THE FIRST STEP!!
    for(i = 0; i < M; i++)
    {
        sprintf(tfn_t[i],"copy%d/input%d.trajt",i+1,i+1);
        sprintf(nfn_t[i],"copy%d/input%d.nw",i+1,i+1);
        sprintf(ffn_t[i],"copy%d/force%d.force",i+1,i+1);

        strcpy(trajtFilenames[i],argv[dirPrefixIndex]);
        strcpy(nwFilenames[i],argv[dirPrefixIndex]);
        strcpy(forceFilenames[i],argv[dirPrefixIndex]);

        strcat(trajtFilenames[i],tfn_t[i]);
        strcat(nwFilenames[i],nfn_t[i]);
        strcat(forceFilenames[i],ffn_t[i]);

        writeTrajt(coord,ID,mass,vel,numPart,trajtFilenames[i]);

        sprintf(tfn_t[i],"copy%d/",i+1);
        strcpy(trajtFilenames[i],argv[dirPrefixIndex]);
        strcat(trajtFilenames[i],tfn_t[i]);

        writeNW(coord,name,numPart,nwFilenames[i],forceFilenames[i],trajtFilenames[i],argv[dtPathIndex]);
    }

    // LAST UPDATED: May 1, 2014 --> Writing .nw file to the root directory as an initial config file for PTime
    writeNW(coord,name,numPart,pathNWOut,pathForceOut,argv[dirPrefixIndex],argv[dtPathIndex]);

    freeChar1D(&pathNWOut);
    freeChar1D(&pathForceOut);
    freeChar2D(tfn_t,M);
    freeChar2D(nfn_t,M);
    freeChar2D(ffn_t,M);
    freeChar2D(nwFilenames,M);
    freeChar2D(trajtFilenames,M);
    freeChar2D(forceFilenames,M);
    freeChar2D(name,numPart);
    free1D(&coord);
    free1D(&vel);
    free1D(&mass);
    freeInt1D(&ID);

    return 0;
}
