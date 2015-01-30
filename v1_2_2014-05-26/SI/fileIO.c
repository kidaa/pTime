#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fileIO.h"
#include "allocRoutine.h"

char elemName[20][3] = {"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca"};
double elemMass[20] = {1.008,4.003,6.941,9.012,10.811,12.011,14.007,15.999,18.998,20.180,22.990,24.305,26.982,28.086,30.974,32.066,35.453,39.948,39.098,40.078};

void readXYZ(double **coord, double **mass, int **index, int *numPart, const char *path)
{
    double *tmp_coord, *tmp_mass, tmp_X, tmp_Y, tmp_Z;
    int tmp_num,count,count2,*ind,i;
    FILE *fp;
    char tmp_name[3],buff[100],*temp,*tok;
    fp = fopen(path,"r");
    count = 1;
    tmp_num = 0;
    while(!feof(fp))
    {
        if(count > 2+tmp_num)
            break;
        if(count == 1)
        {
            fgets(buff,100,fp);
            tmp_num = atoi(buff);
            tmp_coord = alloc1D(3*tmp_num);
            tmp_mass = alloc1D(tmp_num);
            ind = allocInt1D(tmp_num);
        }
        else if(count == 2)
            fgets(buff,100,fp);
        else if((count > 2) && (count <= 2+tmp_num))
        {
            fgets(buff,100,fp);
            temp = buff;
            count2 = 0;
            while((tok = strsep(&temp," ")) != NULL)
            {
                if(count2 == 0)
                {
                    strcpy(tmp_name,tok);
                    count2++;
                    for(i = 0; i < 20; i++)
                    {
                        if(strcmp(tmp_name,elemName[i]) == 0)
                        {
                            tmp_mass[count-3] = elemMass[i];
                            break;
                        }
                    }
                    ind[count-3] = (count-2);
                }
                else if(count2 == 1)
                {
                    tmp_X = (double)(atof(tok));
                    count2++;
                    tmp_coord[3*(count-3)] = tmp_X;
                }
                else if(count2 == 2)
                {
                    tmp_Y = (double)(atof(tok));
                    count2++;
                    tmp_coord[3*(count-3)+1] = tmp_Y;
                }
                else if(count2 == 3)
                {
                    tmp_Z = (double)(atof(tok));
                    count2++;
                    tmp_coord[3*(count-3)+2] = tmp_Z;
                }
                else
                    break;
            }
        }
        count++;
    }

    *coord = tmp_coord;
    *mass = tmp_mass;
    *numPart = tmp_num;
    *index = ind;

    fclose(fp);
}

void writeNW(double *coord, char **name, int numPart, char *pathNW, char *pathF, char *pathPrefix)
{
    FILE *fp;
    int i;
    fp = fopen(pathNW,"w");
    fprintf(fp,"title \"Force calculator\"\n");
    fprintf(fp,"\n");
    fprintf(fp,"charge 0\n");
    fprintf(fp,"\n");
    fprintf(fp,"scratch_dir     %sscratch\n",pathPrefix);
    fprintf(fp,"permanent_dir   %sperm\n",pathPrefix);
    fprintf(fp,"\n");
    fprintf(fp,"geometry units angstrom noautoz noautosym noprint\n");
    for(i = 0; i < numPart; i++)
        fprintf(fp,"   %s      %.7e      %.7e      %.7e\n",name[i],coord[3*i],coord[3*i+1],coord[3*i+2]);
    fprintf(fp,"end\n");
    fprintf(fp,"\n");
    fprintf(fp,"basis noprint\n");
    fprintf(fp,"   *   library   6-31G*\n");            // MAY CHANGE LATER OR ADD BASIS AND THEORY TO THE FUNCTION INPUT ARGUMENTS
    fprintf(fp,"end\n");
    fprintf(fp,"\n");
    fprintf(fp,"python noprint\n");
    fprintf(fp,"   (energy,gradient) = task_gradient('mp2')\n");
    fprintf(fp,"   if(ga_nodeid()==0):\n");
    fprintf(fp,"      fN = 8.2387e-8\n");
    fprintf(fp,"      fileWrite = open('%s','w')\n",pathF);
    fprintf(fp,"      for i in range(len(gradient)/3):\n");
    fprintf(fp,"         print >> fileWrite, str(i+1) + ' ' + str((-1.0)*gradient[3*i]*fN) + ' ' + str((-1.0)*gradient[3*i+1]*fN) + ' ' + str((-1.0)*gradient[3*i+2]*fN)\n");
    fprintf(fp,"      fileWrite.close()\n");
    fprintf(fp,"end\n");
    fprintf(fp,"\n");
    fprintf(fp,"task python");
    fclose(fp);
}

void IDtoName(char ***name, double *massLocal, int numPart)
{
    int i,j;
    char **tmp_name;

    tmp_name = allocChar2D(numPart,3);

    for(i = 0; i < numPart; i++)
    {
        for(j = 0; j < 20; j++)
        {
            if(elemMass[j] == massLocal[i])
            {
                strcpy(tmp_name[i],elemName[j]);
                break;
            }
        }
    }
    *name = tmp_name; 
}

void writeTrajt(double *coord, int *ID, double *mass, double *vel, int numPart, const char *path)
{
    int i;
    FILE *fp;
    const double v_au = 2.18769e+6;
    fp = fopen(path,"w");
    for(i = 0; i < numPart; i++)
        fprintf(fp,"%d %.3lf %.7e %.7e %.7e %.7e %.7e %.7e\n",ID[i],mass[i],coord[3*i],coord[3*i+1],coord[3*i+2],vel[3*i]/v_au,vel[3*i+1]/v_au,vel[3*i+2]/v_au);
    fclose(fp);
}

void readTrajt(double **coord, int **ID, double **mass, double **vel, int numPart, const char *path)
{
    FILE *fp;
    double *temp_coord, *temp_vel, *temp_mass;
    double tcx,tcy,tcz,tvx,tvy,tvz,tmass;
    int *temp_ID, i, c1, c2, tID;
    char buff[150],*temp,*tok;
    fp = fopen(path,"r");
    c1 = 1;
    temp_coord = alloc1D(3*numPart);
    temp_vel = alloc1D(3*numPart);
    temp_mass = alloc1D(numPart);
    temp_ID = allocInt1D(numPart);
    while(!feof(fp))
    {
        if(c1 <= numPart)
        {
            fgets(buff,150,fp);
            temp = buff;
            c2 = 0;
            while((tok = strsep(&temp," ")) != NULL)
            {
                if(c2 == 0)
                {
                    tID = atoi(tok);
                    c2++;
                    temp_ID[c1-1] = tID;
                }
                else if(c2 == 1)
                {
                    tmass = (double)(atof(tok));
                    c2++;
                    temp_mass[c1-1] = tmass;
                }
                else if(c2 == 2)
                {
                    tcx = (double)(atof(tok));
                    c2++;
                    temp_coord[3*(c1-1)] = tcx;
                }
                else if(c2 == 3)
                {
                    tcy = (double)(atof(tok));
                    c2++;
                    temp_coord[3*(c1-1)+1] = tcy;
                }
                else if(c2 == 4)
                {
                    tcz = (double)(atof(tok));
                    c2++;
                    temp_coord[3*(c1-1)+2] = tcz;
                }
                else if(c2 == 5)
                {
                    tvx = (double)(atof(tok));
                    c2++;
                    temp_vel[3*(c1-1)] = tvx;
                }
                else if(c2 == 6)
                {
                    tvy = (double)(atof(tok));
                    c2++;
                    temp_vel[3*(c1-1)+1] = tvy;
                }
                else if(c2 == 7)
                {
                    tvz = (double)(atof(tok));
                    c2++;
                    temp_vel[3*(c1-1)+2] = tvz;
                }
                else
                    break;
            }
        }
        else
            break;
        c1++;
    }
    *coord = temp_coord;
    *vel = temp_vel;
    *mass = temp_mass;
    *ID = temp_ID;
    fclose(fp);
}

void readTrajtPTime(double *coord, int *ID, double *mass, double *vel, int numPart, const char *path, int index)
{
    FILE *fp;
    double tcx,tcy,tcz,tvx,tvy,tvz,tmass;
    int i, c1, c2, tID;
    char buff[150],*temp,*tok;
    fp = fopen(path,"r");
    c1 = 1;
    while(!feof(fp))
    {
        if(c1 <= numPart)
        {
            fgets(buff,150,fp);
            temp = buff;
            c2 = 0;
            while((tok = strsep(&temp," ")) != NULL)
            {
                if(c2 == 0)
                {
                    if(index == 0)
                    {
                        tID = atoi(tok);
                        c2++;
                        ID[c1-1] = tID;
                    }
                    else
                        c2++;
                }
                else if(c2 == 1)
                {
                    if(index == 0)
                    {
                        tmass = (double)(atof(tok));
                        c2++;
                        mass[c1-1] = tmass;
                    }
                    else
                        c2++;
                }
                else if(c2 == 2)
                {
                    tcx = (double)(atof(tok));
                    c2++;
                    coord[3*(c1-1)] = tcx;
                }
                else if(c2 == 3)
                {
                    tcy = (double)(atof(tok));
                    c2++;
                    coord[3*(c1-1)+1] = tcy;
                }
                else if(c2 == 4)
                {
                    tcz = (double)(atof(tok));
                    c2++;
                    coord[3*(c1-1)+2] = tcz;
                }
                else if(c2 == 5)
                {
                    tvx = (double)(atof(tok));
                    c2++;
                    vel[3*(c1-1)] = tvx;
                }
                else if(c2 == 6)
                {
                    tvy = (double)(atof(tok));
                    c2++;
                    vel[3*(c1-1)+1] = tvy;
                }
                else if(c2 == 7)
                {
                    tvz = (double)(atof(tok));
                    c2++;
                    vel[3*(c1-1)+2] = tvz;
                }
                else
                    break;
            }
        }
        else
            break;
        c1++;
    }
    fclose(fp);
}

void readForce(double **force, int *ID, int numPart, const char *path)
{
    FILE *fp;
    double *tmp_force,tfx,tfy,tfz;
    int i,j,tID;
    char buff[120],*temp,*tok;
    fp = fopen(path,"r");
    i = 0;
    tmp_force = alloc1D(3*numPart);
    while(!feof(fp))
    {
        if(i < numPart)
        {
            fgets(buff,120,fp);
            temp = buff;
            j = 0;
            while((tok = strsep(&temp," ")) != NULL)
            {
                if(j == 0)
                {
                    tID = atoi(tok);
                    if(tID != ID[i])
                    {
                        printf("Particle ID doesn't match!\n");
                        exit(1);
                    }
                }
                else if(j == 1)
                {
                    tfx = (double)(atof(tok));
                    tmp_force[3*i] = tfx;
                }
                else if(j == 2)
                {
                    tfy = (double)(atof(tok));
                    tmp_force[3*i+1] = tfy;
                }
                else if(j == 3)
                {
                    tfz = (double)(atof(tok));
                    tmp_force[3*i+2] = tfz;
                }
                else
                    break;
                j++;
            }
        }
        else
            break;
        i++;
    }
    *force = tmp_force;
    fclose(fp);
}

void readForcePTime(double *force, int *ID, int numPart, const char *path)
{
    FILE *fp;
    double tfx,tfy,tfz;
    int i,j,tID;
    char buff[120],*temp,*tok;
    fp = fopen(path,"r");
    i = 0;
    while(!feof(fp))
    {
        if(i < numPart)
        {
            fgets(buff,120,fp);
            temp = buff;
            j = 0;
            while((tok = strsep(&temp," ")) != NULL)
            {
                if(j == 0)
                {
                    tID = atoi(tok);
                    if(tID != ID[i])
                    {
                        printf("Particle ID doesn't match!\n");
                        exit(1);
                    }
                }
                else if(j == 1)
                {
                    tfx = (double)(atof(tok));
                    force[3*i] = tfx;
                }
                else if(j == 2)
                {
                    tfy = (double)(atof(tok));
                    force[3*i+1] = tfy;
                }
                else if(j == 3)
                {
                    tfz = (double)(atof(tok));		// PREVIOUSLY 'tfz' was 'tfy', and this caused an error!!  UPDATED: Apr 29, 2014 4:02PM PDT
                    force[3*i+2] = tfz;
                }
                else
                    break;
                j++;
            }
        }
        else
            break;
        i++;
    }
    fclose(fp);
}

void writeHist(double **rHist, double **aHist, int histNum, int kValue, char *path, int M, int numPart)
{
    // NOTE: rHist and aHist are already in SI units!! No need to further convert the values!
    int i,j,k;
    int midIndX,midIndY,midIndZ;
    FILE *fp;
    fp = fopen(path,"w");
    fprintf(fp,"%d\n",((kValue <= histNum-1)?(kValue+1):histNum));
    for(i = 0; i < ((kValue <= histNum-1)?(kValue+1):histNum); i++)
    {
        for(j = 0; j < M; j++)
        {
            for(k = 0; k < numPart; k++)
            {
                midIndX = 3*j*numPart + (3*k);
                midIndY = 3*j*numPart + (3*k+1);
                midIndZ = 3*j*numPart + (3*k+2);
                fprintf(fp,"%.7e %.7e %.7e %.7e %.7e %.7e\n",rHist[i][midIndX],rHist[i][midIndY],rHist[i][midIndZ],aHist[i][midIndX],aHist[i][midIndY],aHist[i][midIndZ]);
            }
        }
    }
    fclose(fp);
}

void writeTraj(double *XNext, double *mass, int *ID, int IValue, int M, int numPart, char *path)
{
    int step,bigIndX,bigIndY,bigIndZ,bigIndVX,bigIndVY,bigIndVZ;
    int i,j;
    FILE *fp;
    fp = fopen(path,"w");
    for(i = 0; i < M; i++)
    {
        step = IValue*M+i+1;
        for(j = 0; j < numPart; j++)
        {
            bigIndX = 2*i*(3*numPart) + (3*j);
            bigIndY = 2*i*(3*numPart) + (3*j+1);
            bigIndZ = 2*i*(3*numPart) + (3*j+2);
            bigIndVX = bigIndX + 3*numPart;
            bigIndVY = bigIndY + 3*numPart;
            bigIndVZ = bigIndZ + 3*numPart;

            fprintf(fp,"%d %d %.3lf %.7e %.7e %.7e %.7e %.7e %.7e\n",step,ID[j],mass[j],XNext[bigIndX],XNext[bigIndY],XNext[bigIndZ],XNext[bigIndVX],XNext[bigIndVY],XNext[bigIndVZ]);
        }
    }
    fclose(fp);
}

void writeTrajOne(double *x, double *v, double *mass, int *ID, int iVal, int numPart, char *path)
{
    FILE *fp;
    int i;
    fp = fopen(path,"a");
    for(i = 0; i < numPart; i++)
        fprintf(fp,"%d %d %.3lf %.7e %.7e %.7e %.7e %.7e %.7e\n",iVal+1,ID[i],mass[i],x[3*i],x[3*i+1],x[3*i+2],v[3*i],v[3*i+1],v[3*i+2]);
    fclose(fp);
}

void readHist(double **rHist, double **aHist, int M, int numPart, char *path, int histNum)
{
    FILE *fp;
    int i,j,tIter;
    char buff[150],*temp,*tok;

    fp = fopen(path,"r");
    i = 0;

    while(!feof(fp))
    {
        if((i>0)&&(i>tIter*M*numPart))
            break;
        if(i == 0)
        {
            fgets(buff,150,fp);
            tIter = atoi(buff);
        }
        else if(i <= tIter*M*numPart)
        {
            fgets(buff,150,fp);
            temp = buff;
            j = 0;
            while((tok = strsep(&temp," ")) != NULL)
            {
                if(j == 0)									// DON'T FORGET CHECKING MECHANISMS!!
                {
                    if(((i-1)/(M*numPart))+1 < histNum)
                        rHist[((i-1)/(M*numPart))+1][3*((i-1)%(M*numPart))+j] = (double)(atof(tok));	// Added '+1' after (i-1)/(M*numPart) to update rHist & aHist at the same time!
                }
                else if(j == 1)
                {
                    if(((i-1)/(M*numPart))+1 < histNum)
                        rHist[((i-1)/(M*numPart))+1][3*((i-1)%(M*numPart))+j] = (double)(atof(tok));
                }
                else if(j == 2)
                {
                    if(((i-1)/(M*numPart))+1 < histNum)
                        rHist[((i-1)/(M*numPart))+1][3*((i-1)%(M*numPart))+j] = (double)(atof(tok));
                }
                else if(j == 3)
                {
                    if(((i-1)/(M*numPart))+1 < histNum)
                        aHist[((i-1)/(M*numPart))+1][3*((i-1)%(M*numPart))+(j-3)] = (double)(atof(tok));
                }
                else if(j == 4)
                {
                    if(((i-1)/(M*numPart))+1 < histNum)
                        aHist[((i-1)/(M*numPart))+1][3*((i-1)%(M*numPart))+(j-3)] = (double)(atof(tok));
                }
                else if(j == 5)
                {
                    if(((i-1)/(M*numPart))+1 < histNum)
                        aHist[((i-1)/(M*numPart))+1][3*((i-1)%(M*numPart))+(j-3)] = (double)(atof(tok));
                }

                if(j > 5)
                    break;

                j++;
            }
        }
        i++;
    }

    fclose(fp);
}

/*void writeF(double *Fr, double *Fv, int M, int numPart, char *path)
{
    FILE *fp;
    int i;
    fp = fopen(path,"w");
    for(i = 0; i < 3*M*numPart; i++)
        fprintf(fp,"%.7e %.7e\n",Fr[i],Fv[i]);
    fclose(fp);    
}

void readF(double *FrPrev, double *FvPrev, int M, int numPart, char *path)
{
    FILE *fp;
    int i,j;
    char *tmp,*tok,buff[100];

    fp = fopen(path,"r");   
    i = 0;

    while(!feof(fp))
    {
        if(i >= 3*M*numPart)
            break;
        fgets(buff,100,fp);
        tmp = buff;
        j = 0;
        while((tok = strsep(&temp," ")) != NULL)
        {
            if(j == 0)
                FrPrev[i] = (double)(atof(tok));
            else if(j == 1)
                FvPrev[i] = (double)(atof(tok));   
            else
                break;
            j++;
        }
        i++;
    }
    fclose(fp);
}*/

/*int checkOptions(int argCount, char *argVal[], char argIn[])
{
    // RETURN THE CORRESPONDING INDEX OF ARGVAL
    int i,index;
    for(i = 0; i < argCount; i++)
    {
        if(strcmp(argVal[i],argIn)==0)
        {
            index = i+1;
            break;
        }
    }
    return index;
}*/
