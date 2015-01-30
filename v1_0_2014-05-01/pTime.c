/***********************************************************************************************
 *											       *
 *	      pTime.c - parallel in time algorithm for ab initio molecular dynamics            *
 *                      WEARE GROUP, UNIVERSITY OF CALIFORNIA, SAN DIEGO                       *
 *			    LAST UPDATED: April 29, 2014, 3:38PM PDT                           *
 *											       *
 ***********************************************************************************************
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "mathOps.h"
#include "fileIO.h"
#include "arrayManip.h"
#include "allocRoutine.h"

#define MAXNAMELEN 150
#define MIDNAMELEN 50

// VERSION 1.1
// Change log:
//   - Options checking feature added.
//

int main(int argc, char *argv[])
{
    // NOTE: bigIndX = X-index for 6MN size vector
    // NOTE: midIndY = Y-index for 3MN size vector
    // NOTE: sIndZ = Z-index for 3N size vector
    int M,numPart,histNum,i,j,kValue,IValue,*ID,midIndX,midIndY,midIndZ,midPrevX,midPrevY,midPrevZ,sIndX,sIndY,sIndZ;
    double *Fr,*Fv,*force,*mass,*massKg,**rHist,**aHist,**Hr,*coord,*coordPrev,*velPrev,*forcePrev,*vel,*tcoord,*tvel,*tforce;
    double *x0,*v0,*f0,massFac,dt,tolR,tolV;
    double *rTilde,*vTilde,*XNext;
    char **elemName,*pathTrajtIn,*pathStep0In,*pathForce0In,*pathForceIn,*pathHistIn;									// INPUT PATHS
    char *pathTrajtOutConv,*pathTrajtOutNConv,*pathNWOutConv,*pathNWOutNConv,*pathForceOutConv,*pathForceOutNConv,*pathTrajOut,*pathHist;		// OUTPUT PATHS
    int optInd,optCount,dirPrefixIndex;
    char *tps0in,*tpf0in,*tTrajtIn,*tForceIn;
    char *tTrajtOutConv,*tNWOutConv,*tForceOutConv,*tTrajOut,*tTrajtOutNConv,*tNWOutNConv,*tForceOutNConv,*tHist;
    char *pathChunkOut,*tChunkOut,*pathFOut,*pathFIn;

    const double auToKg = 1.6605e-27;

    // CHECK OPTIONS!! - ADD -P and --dt options
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
        else if(strcmp(argv[optInd],"-M")==0)
        {
            M = atoi(argv[optInd+1]);
            optInd += 2;
            optCount++;
        }
        else if(strcmp(argv[optInd],"-I")==0)
        {
            IValue = atoi(argv[optInd+1]);
            optInd += 2;
            optCount++;
        }
        else if(strcmp(argv[optInd],"-k")==0)
        {
            kValue = atoi(argv[optInd+1]);
            optInd += 2;
            optCount++;
        }
        else if(strcmp(argv[optInd],"-h")==0)
        {
            histNum = atoi(argv[optInd+1]);
            optInd += 2;
            optCount++;
        }
        else if(strcmp(argv[optInd],"--dt")==0)
        {
            dt = (double)(atof(argv[optInd+1]));
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
            system("more /home/wpornpat/PTime/src/v1_1_2014-04-29/docs/pTime.txt");
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

    // ARRAY ALLOCATIONS (FOR GOOD PRACTICE, ALLOCATE ARRAYS RIGHT BEFORE THE USAGE AND DEALLOCATE RIGHT AFTER THE USAGE FOR MEMORY EFFICIENCY!!)
    pathTrajtIn = allocChar1D(MAXNAMELEN);
    pathForceIn = allocChar1D(MAXNAMELEN);
    tTrajtIn = allocChar1D(MIDNAMELEN);
    tForceIn = allocChar1D(MIDNAMELEN);
    pathStep0In = allocChar1D(MAXNAMELEN);
    pathForce0In = allocChar1D(MAXNAMELEN);
    Fr = alloc1D(3*M*numPart);
    Fv = alloc1D(3*M*numPart);
    coord = alloc1D(3*M*numPart);
    vel = alloc1D(3*M*numPart);
    force = alloc1D(3*M*numPart);
    tcoord = alloc1D(3*numPart);
    tvel = alloc1D(3*numPart);
    tforce = alloc1D(3*numPart);
    x0 = alloc1D(3*numPart);
    v0 = alloc1D(3*numPart);
    f0 = alloc1D(3*numPart);
    mass = alloc1D(numPart);
    ID = allocInt1D(numPart);

    // TEMPORARY NUMERICAL FILENAME STRINGS ALLOCATIONS

    // READ FROM .TRAJT FILE USING THE FAMILIAR METHOD (DO NOT ALLOCATE X YET!!)
    // DO IT SERIALLY TO REDUCE HDD I/O OVERHEAD!!
    for(i = 0; i < M; i++)
    {
        // READ THE PREVIOUS STEP OUTSIDE THE K-LEVEL LOOP (STEP 0,M,2M,...). THESE GET READ ONLY ONCE PER EACH PTIME ITERATION!
        // LAST UPDATED: May 1, 2014 --> Change the filename template to reflect the initial file!
        if(i == 0)
        {
            tps0in = allocChar1D(MIDNAMELEN);
            tpf0in = allocChar1D(MIDNAMELEN);

            sprintf(tps0in,"step%d_1.trajt",(M*IValue));
            sprintf(tpf0in,"step%d_1.force",(M*IValue));

            strcpy(pathStep0In,argv[dirPrefixIndex]);
            strcpy(pathForce0In,argv[dirPrefixIndex]);

            strcat(pathStep0In,tps0in);
            strcat(pathForce0In,tpf0in);

            freeChar1D(&tps0in);
            freeChar1D(&tpf0in);

            readTrajtPTime(x0,ID,mass,v0,numPart,pathStep0In,i);
            readForcePTime(f0,ID,numPart,pathForce0In);			// THERE WAS AN ERROR, BUT I HAVE EDITTED IN 'fileIO.c' already on Apr 29, 2014, 4:02PM PDT!!
        }

        sprintf(tTrajtIn,"copy%d/input%d.trajt",i+1,i+1);			// TEST IN A SAMPLE CODE FIRST WHETHER THIS ASSIGNMENT WORKS!!
        sprintf(tForceIn,"copy%d/force%d.force",i+1,i+1);

        strcpy(pathTrajtIn,argv[dirPrefixIndex]);	// CHECK IF 'strcpy' replaces the new string to the old variable after the new round in the loop!!
        strcpy(pathForceIn,argv[dirPrefixIndex]);

        strcat(pathTrajtIn,tTrajtIn);
        strcat(pathForceIn,tForceIn);

        if(i == 0)
        {
            readTrajtPTime(tcoord,ID,mass,tvel,numPart,pathTrajtIn,1);					// MASS AND ID ARE ALREADY READ!! PASS 1 TO PREVENT RE-READING!
            readForcePTime(tforce,ID,numPart,pathForceIn);
        }
        else
        {
            readTrajtPTime(tcoord,ID,mass,tvel,numPart,pathTrajtIn,i);                                  // PASS i > 0 TO PREVENT MASS AND ID RE-READING!!
            readForcePTime(tforce,ID,numPart,pathForceIn);
        }

        // ASSIGN VALUES TO THE MAIN CONFIGURATION ARRAYS (3MN SIZE)
        for(j = 0; j < numPart; j++)
        {
            //bigIndX = 3*(2*i)*numPart+(3*j);			// TO REFER TO THE VELOCITY PART IN BIG VECTORS, USE "bigIndX+3*numPart"
            //bigIndY = 3*(2*i)*numPart+(3*j+1);
            //bigIndZ = 3*(2*i)*numPart+(3*j+2);
            midIndX = 3*i*numPart+(3*j);
            midIndY = 3*i*numPart+(3*j+1);
            midIndZ = 3*i*numPart+(3*j+2);
            sIndX = 3*j;
            sIndY = 3*j+1;
            sIndZ = 3*j+2;
            coord[midIndX] = tcoord[sIndX];
            coord[midIndY] = tcoord[sIndY];
            coord[midIndZ] = tcoord[sIndZ];
            vel[midIndX] = tvel[sIndX];
            vel[midIndY] = tvel[sIndY];
            vel[midIndZ] = tvel[sIndZ];
            force[midIndX] = tforce[sIndX];
            force[midIndY] = tforce[sIndY];
            force[midIndZ] = tforce[sIndZ];

            // FOR EFFICIENCY, ASSIGN VALUES FOR F VECTOR RIGHT HERE!!
            if(i == 0)
            {
                massFac = mass[j]*auToKg;
                Fr[midIndX] = coord[midIndX] - propagateX(x0[sIndX],v0[sIndX],f0[sIndX]/massFac,dt);			// AT THIS POINT, Fr are in angstroms!
                Fr[midIndY] = coord[midIndY] - propagateX(x0[sIndY],v0[sIndY],f0[sIndY]/massFac,dt);
                Fr[midIndZ] = coord[midIndZ] - propagateX(x0[sIndZ],v0[sIndZ],f0[sIndZ]/massFac,dt);
                Fv[midIndX] = vel[midIndX] - propagateV(v0[sIndX],f0[sIndX]/massFac,force[midIndX]/massFac,dt);
                Fv[midIndY] = vel[midIndY] - propagateV(v0[sIndY],f0[sIndY]/massFac,force[midIndY]/massFac,dt);
                Fv[midIndZ] = vel[midIndZ] - propagateV(v0[sIndZ],f0[sIndZ]/massFac,force[midIndZ]/massFac,dt);
            }
            else
            {
                massFac = mass[j]*auToKg;
                midPrevX = 3*(i-1)*numPart+(3*j);
                midPrevY = 3*(i-1)*numPart+(3*j+1);
                midPrevZ = 3*(i-1)*numPart+(3*j+2);
                Fr[midIndX] = coord[midIndX] - propagateX(coord[midPrevX],vel[midPrevX],force[midPrevX]/massFac,dt);	// AT THIS POINT, Fr are in angstroms!
                Fr[midIndY] = coord[midIndY] - propagateX(coord[midPrevY],vel[midPrevY],force[midPrevY]/massFac,dt);
                Fr[midIndZ] = coord[midIndZ] - propagateX(coord[midPrevZ],vel[midPrevZ],force[midPrevZ]/massFac,dt);
                Fv[midIndX] = vel[midIndX] - propagateV(vel[midPrevX],force[midPrevX]/massFac,force[midIndX]/massFac,dt);
                Fv[midIndY] = vel[midIndY] - propagateV(vel[midPrevY],force[midPrevY]/massFac,force[midIndY]/massFac,dt);
                Fv[midIndZ] = vel[midIndZ] - propagateV(vel[midPrevZ],force[midPrevZ]/massFac,force[midIndZ]/massFac,dt);
            }
        }
    }

    // CONVERGENCE CHECK HERE!! IF CONVERGED, EXIT! IF NOT, GO ON (NOT WRITTEN)
    // WE DISPLAY VALUES AT 7 SIG FIGS IN THE FILES. IF WE WANT THAT HIGH RESOLUTION, OUR NORM SHOULD BE AT THE ORDER OF MAGNETUDE LESS THAN -14
    tolR = 1.0e-8;
    tolV = 1.0e-6;

    // FOR DEBUGGING PROCESS ONLY
    printf("At k = %d: %.3e	%.3e\n",kValue,vecNormSq(Fr,3*M*numPart),vecNormSq(Fv,3*M*numPart));

    if((vecNormSq(Fr,3*M*numPart) < tolR) && (vecNormSq(Fv,3*M*numPart) < tolV))
    {
        //if(kValue == 1) CHECK IF IT'S ALREADY CONVERGED IN THE FIRST STEP!! OUTPUT THE FILES IMMEDIATELY, OTHERWISE PROCEED!!
        printf("1");	// CONVERGENCE ACHIEVED!! SEND THIS AS A SIGNAL TO BREAK OUT FROM K-LEVEL LOOP IN THE SCRIPT!!
        freeChar1D(&pathTrajtIn);
        freeChar1D(&pathForceIn);
        freeChar1D(&pathStep0In);
        freeChar1D(&pathForce0In);
        freeChar1D(&tTrajtIn);
        freeChar1D(&tForceIn);
        free1D(&tcoord);
        free1D(&tvel);
        free1D(&tforce);
        free1D(&x0);
        free1D(&v0);
        free1D(&f0);
        free1D(&force);
        free1D(&coord);
        free1D(&vel);
        free1D(&mass);
        freeInt1D(&ID);
        free1D(&Fr);
        free1D(&Fv);
        exit(1);
    }

    // FREE ALLOCATED MEMORY AFTER USAGE FOR EFFICIENCY (THESE MEMORY ARE ONLY USED WHEN READING FILES, AND HAVE NO MORE USES AFTER FILE I/O IS DONE, SO BETTER FREE THEM!!
    freeChar1D(&pathTrajtIn);
    freeChar1D(&pathForceIn);
    freeChar1D(&pathStep0In);
    freeChar1D(&pathForce0In);
    freeChar1D(&tTrajtIn);
    freeChar1D(&tForceIn);
    free1D(&tcoord);
    free1D(&tvel);
    free1D(&tforce);

    // HISTORY STORAGE: At first k, only store the current and the step 0 trajectories.
    // ALLOCATE HISTORY ARRAYS HERE!!
    rHist = alloc2D(histNum,3*M*numPart);
    aHist = alloc2D(histNum,3*M*numPart);
    setZero2D(rHist,histNum,3*M*numPart);
    setZero2D(aHist,histNum,3*M*numPart);
    pathHistIn = allocChar1D(MAXNAMELEN);
    strcpy(pathHistIn,argv[dirPrefixIndex]);
    strcat(pathHistIn,"history.hist");
    storeHist(rHist,aHist,coord,x0,force,f0,mass,auToKg,1.0e-10,M,histNum,numPart,kValue,pathHistIn);
    freeChar1D(&pathHistIn);

    // FREE x0,v0,f0,force
    free1D(&x0);
    free1D(&v0);
    free1D(&f0);
    free1D(&force);

    // ALLOCATE *rTilde and *vTilde and pass these 2 arrays to 'calcTilde' function
    rTilde = alloc1D(3*M*numPart);
    vTilde = alloc1D(3*M*numPart);
    calcTilde(Fr,Fv,rTilde,vTilde,rHist,aHist,M,numPart,histNum,kValue,dt);

    // CALCULATE XNEXT
    XNext = alloc1D(6*M*numPart);
    calcXNext(XNext,rTilde,vTilde,coord,vel,numPart,M);

    // FREE coord AND vel, rTilde and vTilde
    free1D(&coord);
    free1D(&vel);
    free1D(&rTilde);
    free1D(&vTilde);

    // OUTPUT THE RESULTS...
    pathChunkOut = allocChar1D(MAXNAMELEN);
    pathTrajtOutConv = allocChar1D(MAXNAMELEN);
    pathNWOutConv = allocChar1D(MAXNAMELEN);
    pathForceOutConv = allocChar1D(MAXNAMELEN);
    tChunkOut = allocChar1D(MIDNAMELEN);
    tTrajtOutConv = allocChar1D(MIDNAMELEN);
    tNWOutConv = allocChar1D(MIDNAMELEN);
    tForceOutConv = allocChar1D(MIDNAMELEN);

    sprintf(tChunkOut,"chunk%d.traj",IValue+1);
    sprintf(tTrajtOutConv,"step%d_1.trajt",(IValue+1)*M);
    sprintf(tNWOutConv,"step%d_1.nw",(IValue+1)*M);
    sprintf(tForceOutConv,"step%d_1.force",(IValue+1)*M);
    strcpy(pathChunkOut,argv[dirPrefixIndex]);
    strcpy(pathTrajtOutConv,argv[dirPrefixIndex]);
    strcpy(pathNWOutConv,argv[dirPrefixIndex]);
    strcpy(pathForceOutConv,argv[dirPrefixIndex]);
    strcat(pathChunkOut,tChunkOut);
    strcat(pathTrajtOutConv,tTrajtOutConv);
    strcat(pathNWOutConv,tNWOutConv);
    strcat(pathForceOutConv,tForceOutConv);

    tcoord = alloc1D(3*numPart);
    tvel = alloc1D(3*numPart);
    writeTraj(XNext,mass,ID,IValue,M,numPart,pathChunkOut);
    sliceX(XNext,tcoord,tvel,M-1,numPart);
    writeTrajt(tcoord,ID,mass,tvel,numPart,pathTrajtOutConv);
    IDtoName(&elemName,mass,numPart);
    writeNW(tcoord,elemName,numPart,pathNWOutConv,pathForceOutConv,argv[dirPrefixIndex]);

    freeChar1D(&tForceOutConv);
    freeChar1D(&pathForceOutConv);
    freeChar1D(&tNWOutConv);
    freeChar1D(&pathNWOutConv);
    freeChar1D(&tTrajtOutConv);
    freeChar1D(&pathTrajtOutConv);
    freeChar1D(&tChunkOut);
    freeChar1D(&pathChunkOut);
    free1D(&tcoord);
    free1D(&tvel);

    pathTrajtOutNConv = allocChar1D(MAXNAMELEN);
    pathNWOutNConv = allocChar1D(MAXNAMELEN);
    pathForceOutNConv = allocChar1D(MAXNAMELEN);
    pathHist = allocChar1D(MAXNAMELEN);

    tTrajtOutNConv = allocChar1D(MIDNAMELEN);
    tNWOutNConv = allocChar1D(MIDNAMELEN);
    tForceOutNConv = allocChar1D(MIDNAMELEN);
    tHist = allocChar1D(MIDNAMELEN);

    tcoord = alloc1D(3*numPart);
    tvel = alloc1D(3*numPart);

    for(i = 0; i < M; i++)
    {
        strcpy(pathTrajtOutNConv,argv[dirPrefixIndex]);                     // Test in a simple program whether 'strcpy' replaces every time!!
        strcpy(pathNWOutNConv,argv[dirPrefixIndex]);
        strcpy(pathForceOutNConv,argv[dirPrefixIndex]);

        sprintf(tTrajtOutNConv,"copy%d/input%d.trajt",i+1,i+1);
        sprintf(tNWOutNConv,"copy%d/input%d.nw",i+1,i+1);
        sprintf(tForceOutNConv,"copy%d/force%d.force",i+1,i+1);

        strcat(pathTrajtOutNConv,tTrajtOutNConv);
        strcat(pathNWOutNConv,tNWOutNConv);
        strcat(pathForceOutNConv,tForceOutNConv);

        sliceX(XNext,tcoord,tvel,i,numPart);
        writeTrajt(tcoord,ID,mass,tvel,numPart,pathTrajtOutNConv);
        writeNW(tcoord,elemName,numPart,pathNWOutNConv,pathForceOutNConv,argv[dirPrefixIndex]);
    }

    strcpy(pathHist,argv[dirPrefixIndex]);
    sprintf(tHist,"history.hist");
    strcat(pathHist,tHist);
    writeHist(rHist,aHist,histNum,kValue,pathHist,M,numPart);

    free1D(&tcoord);
    free1D(&tvel);
    freeChar2D(elemName,numPart);
    freeChar1D(&pathTrajtOutNConv);
    freeChar1D(&pathNWOutNConv);
    freeChar1D(&pathForceOutNConv);
    freeChar1D(&pathHist);
    freeChar1D(&tHist);
    freeChar1D(&tTrajtOutNConv);
    freeChar1D(&tNWOutNConv);
    freeChar1D(&tForceOutNConv);

    // DEALLOCATE THE ARRAYS THAT IS USED THROUGHOUT THE CODE!!
    free1D(&mass);
    freeInt1D(&ID);
    free2D(rHist,histNum);		// CHECK IF I COULD DEALLOCATE THIS EARLIER!!
    free2D(aHist,histNum);		// CHECK IF I COULD DEALLOCATE THIS EARLIER!!
    free1D(&Fr);
    free1D(&Fv);

    printf("0");		// THIS IS A SIGNAL THAT TELLS THAT THE K-LEVEL ITERATION HAS NOT CONVERGED

    return 0;
}
