#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allocRoutine.h"

int main(int argc, char *argv[])
{
    int optInd, optCount, dirPrefixIndex, numNodes, i;
    FILE *confInitial, *confMain, *confFinal;
    char *pathInitial, *pathMain, *pathFinal;

    // CHECKING OPTIONS
    optInd = 1;
    optCount = 0;
    while(optInd < argc)
    {
        if(strcmp(argv[optInd],"-P")==0)
        {
            dirPrefixIndex = optInd+1;
            optInd += 2;
            optCount++;
        }
        else if(strcmp(argv[optInd],"-M")==0)
        {
            numNodes = atoi(argv[optInd+1]);
            optInd += 2;
            optCount++;
        }
        else if(strcmp(argv[optInd],"--help")==0)
        {
            system("cat /home/wpornpat/PTime/src/docs/readXYZ.txt");            // UPDATE PATH TO DOCUMENTATION!!
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

    pathInitial = allocChar1D(120);
    pathMain = allocChar1D(120);
    pathFinal = allocChar1D(120);

    strcpy(pathInitial,argv[dirPrefixIndex]);
    strcpy(pathMain,argv[dirPrefixIndex]);
    strcpy(pathFinal,argv[dirPrefixIndex]);

    strcat(pathInitial,"mpdb.conf");
    strcat(pathMain,"slurmConf.conf");
    strcat(pathFinal,"ex.conf"); 

    confInitial = fopen(pathInitial,"w");
    confMain = fopen(pathMain,"w");
    confFinal = fopen(pathFinal,"w");

    for(i = 0; i < numNodes; i++)
    {
        fprintf(confInitial,"%d   mpdboot -n 1\n",i);
        fprintf(confMain,"%d   /dtemp/scicons/bin/nwchem %scopy%d/input%d.nw\n",i,argv[dirPrefixIndex],i+1,i+1);
        fprintf(confFinal,"%d   mpdallexit\n",i);
    }

    fclose(confInitial);
    fclose(confMain);
    fclose(confFinal);

    free(pathInitial);
    free(pathMain);
    free(pathFinal);

    return 0;
}
