#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>

float **fmatrix(int,int);
double **dmatrix(int,int);
char **cmatrix(int,int);
char ***c3tensor(int,int,int);
float ***f3tensor(int,int,int);
void free_cmatrix(char **);
void free_fmatrix(float **);
void free_f3tensor(float ***);

float **fmatrix(int nr,int nc)
{
  int i;
  float **m;

  m=NULL;
  m=(float **) calloc(nr, sizeof(float*));
  if (!m) {printf("allocation failure 1 in matrix()\n");exit(0);}

  m[0]=(float*)calloc(nr*nc, sizeof(float));
  if (!m[0]) {printf("allocation failure 1 in matrix()\n");exit(0);}

  for(i=1;i<nr;i++)
    m[i]=m[i-1]+nc;

  return m;
}

void free_fmatrix(float **m)
{
  free(m[0]); free(m);
}

double **dmatrix(int nr,int nc)
{
  int i;
  double **m;

  m=NULL;
  m=(double **) calloc(nr, sizeof(double*));
  if (!m) {printf("allocation failure 1 in matrix()\n");exit(0);}

  m[0]=(double*)calloc(nr*nc, sizeof(double));
  if (!m[0]) {printf("allocation failure 1 in matrix()\n");exit(0);}

  for(i=1;i<nr;i++)
    m[i]=m[i-1]+nc;

  return m;
}

void free_dmatrix(double **m)
{
  free(m[0]); free(m);
}


char **cmatrix(int nr,int nc)
{
  int i;
  char **m;

  m=NULL;
  m=(char **) calloc(nr, sizeof(char*));
  if (!m) {printf("allocation failure 1 in matrix()\n");exit(0);}

  m[0]=(char*)calloc(nr*nc, sizeof(char));
  if (!m[0]) {printf("allocation failure 1 in matrix()\n");exit(0);}

  for(i=1;i<nr;i++)
    m[i]=m[i-1]+nc;

  return m;
}

void free_cmatrix(char **m)
{
  free(m[0]); free(m);
}


char ***c3tensor(int nr, int nc, int nd)
{
  int i,j;
  char ***t;

  t=(char ***) calloc(nr, sizeof(char**));
  if (!t) {printf("allocation failure 1 in tensor()\n");exit(0);}

  t[0]=(char**)calloc(nr*nc, sizeof(char*));
  if (!t[0]) {printf("allocation failure 1 in tensor()\n");exit(0);}

  t[0][0]=(char*)calloc(nr*nc*nd, sizeof(char));
  if (!t[0][0]) {printf("allocation failure 3 in tensor()\n");exit(0);}

  for(j=1;j<nc;j++) t[0][j]=t[0][j-1]+nd;
  for(i=1;i<nr;i++)
    {
      t[i]=t[i-1]+nc;
      t[i][0]=t[i-1][0]+nc*nd;
      for(j=1;j<nc;j++) t[i][j]=t[i][j-1]+nd;
    }
  return t;
}


float ***f3tensor(int nr, int nc, int nd)
{
  int i,j;
  float ***t;

  t=(float ***) calloc(nr, sizeof(float**));
  if (!t) {printf("allocation failure 1 in tensor()\n");exit(0);}

  t[0]=(float**)calloc(nr*nc, sizeof(float*));
  if (!t[0]) {printf("allocation failure 1 in tensor()\n");exit(0);}

  t[0][0]=(float*)calloc(nr*nc*nd, sizeof(float));
  if (!t[0][0]) {printf("allocation failure 3 in tensor()\n");exit(0);}

  for(j=1;j<nc;j++) t[0][j]=t[0][j-1]+nd;
  for(i=1;i<nr;i++)
    {
      t[i]=t[i-1]+nc;
      t[i][0]=t[i-1][0]+nc*nd;
      for(j=1;j<nc;j++) t[i][j]=t[i][j-1]+nd;
    }
  return t;
}


void free_3tensor(float ***m)
{
  free(m[0][0]);  free(m[0]); free(m);
}




