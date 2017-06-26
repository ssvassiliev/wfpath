#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdbool.h>

#include "../include/fdtimes.h"
#include "time_3d.c"
#include "/home/svassili/src/Lib/memalloc.c"


#include "../include/eikonal.h"
#include "io.c"

#define maxline 200

typedef struct
{
  char pathfile[maxline];
  char speedmap[maxline];
  char timemap[maxline];
  float start_coord[3];
  float end_coord[3];
  float fd_start_coord[3];
  float integration_step;
  float output_step;
  float scale_pot;
  bool do_fdsolve;
  bool do_raytrace;
  bool do_pmf;
  bool invert_init_direct;
} SETUP;


typedef struct vector {
  float x; float y; float z;
} vec3d;

static inline float vnorm(float *vec) { 
  return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]); }

static inline float dist3d(float *p1, float *p2){
  return sqrt((p2[0]-p1[0])*(p2[0]-p1[0])+(p2[1]-p1[1])*(p2[1]-p1[1])+(p2[2]-p1[2])*(p2[2]-p1[2]));}

static inline float vecdot(float *b, float *c)
{return (b[0]*c[0] + b[1]*c[1] + b[2]*c[2]);}


void grad6(float ***GX, float ***GY, float ***GZ, int x, int y, int z) {
  float g1[3],g2[3];
  int inx,ipx,iny,ipy,inz,ipz;
  // Central differences
  inx=ipx=iny=ipy=inz=ipz=0;
  if(x > 0) inx = 1;
  if(y > 0) iny = 1;
  if(z > 0) inz = 1;
  if(x < xsize-1) ipx = 1;
  if(y < ysize-1) ipy = 1;
  if(z < zsize-1) ipz = 1;

  g1[0] = M[x-inx][y][z];
  g2[0] = M[x+ipx][y][z];
  g1[1] = M[x][y-iny][z];
  g2[1] = M[x][y+ipy][z];
  g1[2] = M[x][y][z-inz]; 
  g2[2] = M[x][y][z+ipz];

  GX[x][y][z] = g2[0]-g1[0];
  GY[x][y][z] = g2[1]-g1[1];
  GZ[x][y][z] = g2[2]-g1[2];

  if(inx && ipx) GX[x][y][z] *= 0.5; 
  if(iny && ipy) GY[x][y][z] *= 0.5; 
  if(inz && ipz) GZ[x][y][z] *= 0.5; 
}


void grad3i(float ***GX, float ***GY, float ***GZ, int x, int y, int z) {
// Intermediate differences
  int ipx,ipy,ipz;
  ipx=ipy=ipz=0;
  if(x < xsize-1) ipx = 1;
  if(y < ysize-1) ipy = 1;
  if(z < zsize-1) ipz = 1;

  GX[x][y][z] = M[x+ipx][y][z] - M[x][y][z];
  GY[x][y][z] = M[x][y+ipy][z] - M[x][y][z];
  GZ[x][y][z] = M[x][y][z+ipz] - M[x][y][z];
}


void interpgrad(float *Greturn, float ***GX, float ***GY, float ***GZ, float *pointCoord, int grad_type )
{
  int i, j, k;
  float x, y, z;
  double tmpG[3], norm;

  x = pointCoord[0]; y = pointCoord[1]; z = pointCoord[2];
  if(!grad_type)
    {x -= 0.5; y -= 0.5; z -= 0.5;} // only if intermediate differences gradient (type 0) is used
  i = floor(x); j = floor(y); k = floor(z);
  x -= i; y -= j; z -= k;

  tmpG[0]= \
    GX[i][j][k]*(1-x)*(1-y)*(1-z) +
    GX[i+1][j][k]*x*(1-y)*(1-z) +
    GX[i][j+1][k]*(1-x)*y*(1-z) +
    GX[i][j][k+1]*(1-x)*(1-y)*z +
    GX[i+1][j][k+1]*x*(1-y)*z +
    GX[i][j+1][k+1]*(1-x)*y*z +
    GX[i+1][j+1][k]*x*y*(1-z) +
    GX[i+1][j+1][k+1]*x*y*z;
  tmpG[1]= \
    GY[i][j][k]*(1-x)*(1-y)*(1-z) +
    GY[i+1][j][k]*x*(1-y)*(1-z) +
    GY[i][j+1][k]*(1-x)*y*(1-z) +
    GY[i][j][k+1]*(1-x)*(1-y)*z +
    GY[i+1][j][k+1]*x*(1-y)*z +
    GY[i][j+1][k+1]*(1-x)*y*z +
    GY[i+1][j+1][k]*x*y*(1-z) +
    GY[i+1][j+1][k+1]*x*y*z;
  tmpG[2]= \
    GZ[i][j][k]*(1-x)*(1-y)*(1-z) +
    GZ[i+1][j][k]*x*(1-y)*(1-z) +
    GZ[i][j+1][k]*(1-x)*y*(1-z) +
    GZ[i][j][k+1]*(1-x)*(1-y)*z +
    GZ[i+1][j][k+1]*x*(1-y)*z +
    GZ[i][j+1][k+1]*(1-x)*y*z +
    GZ[i+1][j+1][k]*x*y*(1-z) +
    GZ[i+1][j+1][k+1]*x*y*z;

  norm=sqrt(tmpG[0]*tmpG[0]+tmpG[1]*tmpG[1]+tmpG[2]*tmpG[2]);
  Greturn[0]=tmpG[0]/norm;
  Greturn[1]=tmpG[1]/norm;
  Greturn[2]=tmpG[2]/norm;
} 

float interpval(float ***M, vec3d pointCoord)
{
  int i,j,k;
  float x, y, z;

  x = pointCoord.x; y = pointCoord.y; z = pointCoord.z;
  i = floor(x); j = floor(y); k = floor(z);
  x -= i; y -= j; z -= k;

  return(					\
	 M[i][j][k]*(1-x)*(1-y)*(1-z) +
	 M[i+1][j][k]*x*(1-y)*(1-z) +
	 M[i][j+1][k]*(1-x)*y*(1-z) +
	 M[i][j][k+1]*(1-x)*(1-y)*z +
	 M[i+1][j][k+1]*x*(1-y)*z +
	 M[i][j+1][k+1]*(1-x)*y*z +
	 M[i+1][j+1][k]*x*y*(1-z) +
	 M[i+1][j+1][k+1]*x*y*z); 
}

void rk4step_3d( float ***GX, float ***GY, float ***GZ, float *startPoint, float *nextPoint,  float *Greturn, float stepSize, int grad_type)
{
  float k1[3], k2[3], k3[3], k4[3], k5[3];
  float tempPoint[3];
  float tempnorm;
  float ss;
  float alpha = -0.0;
  
  /*Calculate k1 */
  interpgrad(k1, GX, GY, GZ, startPoint, grad_type);
  tempnorm=vnorm(k1);
  k1[0] = k1[0]*stepSize/tempnorm;
  k1[1] = k1[1]*stepSize/tempnorm;
  k1[2] = k1[2]*stepSize/tempnorm;
  tempPoint[0]=startPoint[0] - k1[0]*0.5;
  tempPoint[1]=startPoint[1] - k1[1]*0.5;
  tempPoint[2]=startPoint[2] - k1[2]*0.5;
  /*Calculate k2 */
  interpgrad(k2, GX, GY, GZ, tempPoint, grad_type);
  tempnorm=vnorm(k2);
  k2[0] = k2[0]*stepSize/tempnorm;
  k2[1] = k2[1]*stepSize/tempnorm;
  k2[2] = k2[2]*stepSize/tempnorm;
  tempPoint[0]=startPoint[0] - k2[0]*0.5;
  tempPoint[1]=startPoint[1] - k2[1]*0.5;
  tempPoint[2]=startPoint[2] - k2[2]*0.5; 
  /*Calculate k3 */
  interpgrad(k3, GX, GY, GZ, tempPoint, grad_type);
  tempnorm=vnorm(k3);
  k3[0] = k3[0]*stepSize/tempnorm;
  k3[1] = k3[1]*stepSize/tempnorm;
  k3[2] = k3[2]*stepSize/tempnorm;
  tempPoint[0]=startPoint[0] - k3[0];
  tempPoint[1]=startPoint[1] - k3[1];
  tempPoint[2]=startPoint[2] - k3[2];
  /*Calculate k4 */  
  interpgrad(k4, GX, GY, GZ, tempPoint, grad_type);
  tempnorm=vnorm(k4);
  k4[0] = k4[0]*stepSize/tempnorm;
  k4[1] = k4[1]*stepSize/tempnorm;
  k4[2] = k4[2]*stepSize/tempnorm;
  /*Calculate final point */
  k5[0] = (k1[0] + k2[0]*2.0 + k3[0]*2.0 + k4[0])/6.0;
  k5[1] = (k1[1] + k2[1]*2.0 + k3[1]*2.0 + k4[1])/6.0;
  k5[2] = (k1[2] + k2[2]*2.0 + k3[2]*2.0 + k4[2])/6.0;
  tempnorm=vnorm(k5);
  k5[0] = k5[0]*stepSize/tempnorm;
  k5[1] = k5[1]*stepSize/tempnorm;
  k5[2] = k5[2]*stepSize/tempnorm;

//  k5[0] = k1[0]*stepSize/tempnorm;
//  k5[1] = k1[1]*stepSize/tempnorm;
//  k5[2] = k1[2]*stepSize/tempnorm;



  nextPoint[0] = startPoint[0] - k5[0];
  nextPoint[1] = startPoint[1] - k5[1];   
  nextPoint[2] = startPoint[2] - k5[2]; 

  Greturn[0]=k5[0]/tempnorm;
  Greturn[1]=k5[1]/tempnorm;
  Greturn[2]=k5[2]/tempnorm;	
}


void inverse_step(float *oldpoint, float *newpoint, float *g, float step)
{
newpoint[0] = oldpoint[0] + g[0]*step;
newpoint[1] = oldpoint[1] + g[1]*step;   
newpoint[2] = oldpoint[2] + g[2]*step; 
}

int gradient_descent_3d(SETUP *setup)
{
  FILE *fp;
  float *TmpG, *OldG, *startCoord, *endCoord, *newCoord, tempnorm;
  float step, outstep, distance, newdistance;
  int i, ii, j, k, c, stride;
  vec3d *path;
  int ntrials, pointcount; 
  float tc1[3], tc2[3], tmpc[3];

  ntrials=200000; 
  startCoord=malloc(3*sizeof(float));
  endCoord=malloc(3*sizeof(float));
  newCoord=malloc(3*sizeof(float));
  path=malloc(sizeof(vec3d));

  startCoord[0] = setup->start_coord[0];
  startCoord[1] = setup->start_coord[1];
  startCoord[2] = setup->start_coord[2];
  endCoord[0] = setup->end_coord[0];
  endCoord[1] = setup->end_coord[1];
  endCoord[2] = setup->end_coord[2];

  printf("\n ** Raytracing module **\n\n"); 
  read_dx_file(setup->timemap); 
  print_stats();
 
  printf("Starting from: %8.3f %8.3f %8.3f\n", startCoord[0], startCoord[1], startCoord[2]);
  printf("Source at:     %8.3f %8.3f %8.3f\n", endCoord[0], endCoord[1], endCoord[2]);
  printf("Gradient descent:\n"); 
  // Allocate memory for 3d data and gradients
  M=f3tensor(xsize,ysize,zsize);
  GX_i=f3tensor(xsize,ysize,zsize);
  GY_i=f3tensor(xsize,ysize,zsize);
  GZ_i=f3tensor(xsize,ysize,zsize);
  GX_c=f3tensor(xsize,ysize,zsize);
  GY_c=f3tensor(xsize,ysize,zsize);
  GZ_c=f3tensor(xsize,ysize,zsize);
  TmpG=malloc(3*sizeof(float));
  OldG=malloc(3*sizeof(float));

  step=setup->integration_step;
  outstep=setup->output_step;
  stride=outstep/step;
  // Convert coordinates to grid units
  startCoord[0] = (startCoord[0]-origin[0])/vnorm(xdelta)+0.5;
  startCoord[1] = (startCoord[1]-origin[1])/vnorm(ydelta)+0.5;
  startCoord[2] = (startCoord[2]-origin[2])/vnorm(zdelta)+0.5;
  endCoord[0] = (endCoord[0]-origin[0])/vnorm(xdelta)+0.5;
  endCoord[1] = (endCoord[1]-origin[1])/vnorm(ydelta)+0.5;
  endCoord[2] = (endCoord[2]-origin[2])/vnorm(zdelta)+0.5;
  // Convert data to 3D array
  ii=0;
  for(k=0;k<zsize;k++)
    for(j=0;j<ysize;j++)
      for(i=0;i<xsize;i++)
	{ M[i][j][k] = data[ii]; ii++;}
  // Compute gradient volume
  for(k=0;k<zsize;k++)
    for(j=0;j<ysize;j++)
      for(i=0;i<xsize;i++)
	{
	  grad6(GX_c, GY_c, GZ_c, i,j,k );
	  grad3i(GX_i, GY_i, GZ_i, i, j, k );
	}

  // rk4 integration
  distance=dist3d(startCoord,endCoord);
  interpgrad(OldG, GX_i, GY_i, GZ_i, startCoord, 0);

  path[0].x = startCoord[0];
  path[0].y = startCoord[1];
  path[0].z = startCoord[2];
  for(pointcount=1;;)
    {
      // save initial coord
      tmpc[0]=startCoord[0];
      tmpc[1]=startCoord[1];
      tmpc[2]=startCoord[2];

      // try intermediate differences
      for(i=0;i<ntrials;i++)
	{
	  rk4step_3d( GX_i, GY_i, GZ_i, startCoord, startCoord, TmpG, step, 0);
	  if(dist3d(startCoord,tmpc) > setup->output_step)break;
	}

    
      if(i == ntrials) // propagation failed
      	{
      	  printf("I.Diff. step failed, trying C.Diff.\n");
	  // restore initial coord
      	  startCoord[0]= tmpc[0];
      	  startCoord[1]= tmpc[1];
      	  startCoord[2]= tmpc[2];
	  // try central differences
      	  for(i=0;i<ntrials;i++)
      	    {
      	      rk4step_3d( GX_c, GY_c, GZ_c, startCoord, startCoord, TmpG, step, 1);
     	      if(dist3d(startCoord,tmpc) > setup->output_step)break;
      	    }
      	}

      newdistance=dist3d(startCoord,endCoord);

      if((i == ntrials) && (newdistance*vnorm(xdelta)  < 1.0))goto normalterm;

      if(i==ntrials) // Both methods failed, terminate
      	{printf("Terminated:  Central Differences gradient step failed %f\n",dist3d(startCoord,tmpc));goto errorterm;}
  
    // check if we are still moving in same direction
      // direction of previous step 
      tc1[0]=path[pointcount-2].x-path[pointcount-1].x;
      tc1[1]=path[pointcount-2].y-path[pointcount-1].y;
      tc1[2]=path[pointcount-2].z-path[pointcount-1].z;
      // direction of current step
      tc2[0]=path[pointcount-1].x-startCoord[0];
      tc2[1]=path[pointcount-1].y-startCoord[1];
      tc2[2]=path[pointcount-1].z-startCoord[2];

      if(vecdot(tc1,tc2) < 0)
	{
	  // if direction changed by > 90, take step in opposite direction
	  startCoord[0]=path[pointcount-1].x+tc2[0];
	  startCoord[1]=path[pointcount-1].y+tc2[1];
	  startCoord[2]=path[pointcount-1].z+tc2[2];
	  printf(" ** WARNING **  Direction changed by > 90 grad, inverting step\n");
	}

      if(setup->invert_init_direct == true && pointcount == 1)
	{
	  startCoord[0]=2*tmpc[0]-startCoord[0];
	  startCoord[1]=2*tmpc[1]-startCoord[1];
	  startCoord[2]=2*tmpc[2]-startCoord[2];
	}
     
      newdistance=dist3d(startCoord,endCoord);
      if((newdistance > distance) && (distance*vnorm(xdelta) < 1.0))
	{pointcount--;goto normalterm;}
      if(newdistance > distance)
	{printf(" ** WARNING ** Delta D = %g in rk4 step\n",(newdistance-distance)*vnorm(xdelta));}

      printf("%i D = %f\n", pointcount, newdistance*vnorm(xdelta));
      distance=newdistance;
      path=realloc(path,sizeof(vec3d)*(pointcount+1));
      if(path==NULL) {printf("Out of memory\n");exit(0);}
      path[pointcount].x = startCoord[0];
      path[pointcount].y = startCoord[1];
      path[pointcount].z = startCoord[2];
      pointcount++; 
    }

  normalterm:
  {printf("Normal termination: wave source reached\n");
    // Save path
    fp=fopen("path.xyz","w");
    fprintf(fp,"%i\ncreated by eikonal\n",pointcount);
    for(ii=0;ii<pointcount;ii++){
      path[ii].x=(path[ii].x-0.5)*vnorm(xdelta)+origin[0];
      path[ii].y=(path[ii].y-0.5)*vnorm(ydelta)+origin[1];
      path[ii].z=(path[ii].z-0.5)*vnorm(zdelta)+origin[2];
      fprintf(fp," O       %f %f %f\n",path[ii].x, path[ii].y, path[ii].z);
    }
    fclose(fp); 
    free_3tensor(M); free(path); 
    free_3tensor(GX_i); free_3tensor(GY_i); free_3tensor(GZ_i);
    free_3tensor(GX_c); free_3tensor(GY_c); free_3tensor(GZ_c);
    return(0);}
 errorterm:
  {printf("Abnormal termination: wave source not reached\n");
    // Save path
    fp=fopen("path.xyz","w");
    fprintf(fp,"%i\n\n",pointcount);
    for(ii=0;ii<pointcount;ii++){
      path[ii].x=(path[ii].x-0.5)*vnorm(xdelta)+origin[0];
      path[ii].y=(path[ii].y-0.5)*vnorm(ydelta)+origin[1];
      path[ii].z=(path[ii].z-0.5)*vnorm(zdelta)+origin[2];
      fprintf(fp," O       %f %f %f\n", path[ii].x, path[ii].y, path[ii].z);
    }
    fclose(fp); 
    free_3tensor(M); free(path); 
    free_3tensor(GX_i); free_3tensor(GY_i); free_3tensor(GZ_i);
    free_3tensor(GX_c); free_3tensor(GY_c); free_3tensor(GZ_c);
    return(1);}
  
}

void calc_pmf(SETUP *setup)
{
  FILE *fp;
  char *line_buf;
  int i, ii, j, k, n; 
  vec3d *path;
  float val;

  printf("\n ** PMF module **\n\n");  
  read_dx_file(setup->speedmap); 
  print_stats();
  printf("Projecting PMF on reaction coordinate\n");  

  line_buf= malloc(maxline*sizeof(char));	
  fp=fopen(setup->pathfile, "rt");
  if(fp==NULL){printf("** E ** xyz file not found\n");exit (1);}
  fgets(line_buf,maxline,fp); 
  sscanf(line_buf,"%i",&n);
  path=malloc(n*sizeof(vec3d));
  fgets(line_buf,maxline,fp);
  for(i=0;i<n;i++){
    fgets(line_buf,maxline,fp);
    sscanf(&line_buf[3],"%f%f%f",&path[i].x,&path[i].y,&path[i].z);

  } 
  free(line_buf);
  fclose(fp);

  M=f3tensor(xsize,ysize,zsize);
  // Convert data to 3D array
  ii=0;
  for(k=0;k<zsize;k++)
    for(j=0;j<ysize;j++)
      for(i=0;i<xsize;i++)
	{ M[i][j][k] = data[ii]; ii++;}
      
  fp=fopen("pmf.dat", "w");
  for(i=0;i<n;i++)
    {
      path[i].x = (path[i].x-origin[0])/vnorm(xdelta);
      path[i].y = (path[i].y-origin[1])/vnorm(ydelta);
      path[i].z = (path[i].z-origin[2])/vnorm(zdelta);

      val=interpval(M,path[i]);
      fprintf(fp,"%i %f\n",i, val*0.6);
    }
  fclose(fp);  

  printf("pmf.dat saved\n\n");
  read_dx_file("times.dx"); 
  print_stats();
  printf("Saving arrival times\n");  

  ii=0;
  for(k=0;k<zsize;k++)
    for(j=0;j<ysize;j++)
      for(i=0;i<xsize;i++)
	{ M[i][j][k] = data[ii]; ii++;}

 fp=fopen("time.dat", "w");
  for(i=0;i<n;i++)
    {
      val=interpval(M,path[i]);
      fprintf(fp,"%i %f\n",i, val);
    }
  fclose(fp);
  free_3tensor(M);
  printf("time.dat saved\n");
  printf("Normal termination\n\n");  
}    


void fdsolve(SETUP *setup)
{
  FILE *fp;
  float xs, ys, zs;
  int i, j, k, ii;
 
  xs = setup->fd_start_coord[0];
  ys = setup->fd_start_coord[1];
  zs = setup->fd_start_coord[2];
      
  printf("\n ** Wavefront propagation module **\n\n");  
  read_dx_file(setup->speedmap); 
  print_stats();

  // Convert starting point from coord to grid units
  xs -= origin[0]; ys -= origin[1]; zs -= origin[2];
  k = (int)(xs*xaxis[0] + ys*xaxis[1] + zs*xaxis[2])/(vnorm(xaxis)*vnorm(xdelta));
  j = (int)(xs*yaxis[0] + ys*yaxis[1] + zs*yaxis[2])/(vnorm(yaxis)*vnorm(ydelta));
  i = (int)(xs*zaxis[0] + ys*zaxis[1] + zs*zaxis[2])/(vnorm(zaxis)*vnorm(zdelta));
  
  // Data exponential and scaling
  for (ii=0; ii < N; ii++) 
    {
      data[ii]=exp((data[ii]-min_E)*setup->scale_pot);
      if(data[i] > 1e8)data[i]=1e8;
    }
  
  timefld=malloc(N*sizeof(float));
  // Propagate wavefront
  time_3d(data, timefld, xsize, ysize, zsize, i, j, k, 0.001, 1);
  write_dx_file("times.dx", timefld);
  printf("Normal termination\n\n");  
}


void read_setup(char *filename, SETUP *setup)
{ 
  FILE *fp;
  char *line_buf,  *var, *value;
  float x,y,z;

  line_buf=malloc(maxline*sizeof(char));	
  var=malloc(maxline*sizeof(char));	
  value=malloc(maxline*sizeof(char));	

  fp=fopen(filename, "rt");
  if(fp==NULL)
    {printf("** E ** setup file not found\n");exit(1);}

  while(1)
    {
      if(fgets(line_buf,maxline,fp)==NULL)
	break;
      sscanf(line_buf,"%s%s",var,value);
      if(!strncasecmp(var,"TimeMap",strlen(var)))     
        strcpy(setup->timemap,value);
      if(!strncasecmp(var,"VolumeMap",strlen(var)))     
	strcpy(setup->speedmap,value);
     if(!strncasecmp(var,"PathFile",strlen(var)))     
	strcpy(setup->pathfile,value);

      if(!strncasecmp(var,"FDSolveStartCoord",strlen(var)))
	{
	  sscanf(&line_buf[strlen(var)],"%f %f %f",&x,&y,&z );
	  setup->fd_start_coord[0]=x;
	  setup->fd_start_coord[1]=y;
	  setup->fd_start_coord[2]=z;
	}
      if(!strncasecmp(var,"RaytraceStartCoord",strlen(var)))
	{
	  sscanf(&line_buf[strlen(var)],"%f %f %f",&x,&y,&z );
	  setup->start_coord[0]=x;
	  setup->start_coord[1]=y;
	  setup->start_coord[2]=z;
	}
      if(!strncasecmp(var,"RaytraceEndCoord",strlen(var)))
	{
	  sscanf(&line_buf[strlen(var)],"%f %f %f",&x,&y,&z );
	  setup->end_coord[0]=x;
	  setup->end_coord[1]=y;
	  setup->end_coord[2]=z;
	}
      if(!strncasecmp(var,"PotentialScaling",strlen(var)))	
	setup->scale_pot=atof(value);
      if(!strncasecmp(var,"IntegrationStep",strlen(var)))	
	setup->integration_step=atof(value);
      if(!strncasecmp(var,"OutputStep",strlen(var)))	
	setup->output_step=atof(value);
      if(!strncasecmp(var,"FDSolve",strlen(var)))	
	if(!strncmp(value,"yes", 3))setup->do_fdsolve=true;
      if(!strncasecmp(var,"Raytrace",strlen(var)))	
	if(!strncmp(value,"yes", 3))setup->do_raytrace=true;
      if(!strncasecmp(var,"CalcPMF",strlen(var)))	
	if(!strncmp(value,"yes", 3))setup->do_pmf=true;
     if(!strncasecmp(var,"InvertInitialDirection",strlen(var)))	
	if(!strncmp(value,"yes", 3))setup->invert_init_direct=true;
    }
  fclose(fp);
}


int main(int argc, char **argv)
{  
  SETUP *setup = malloc(sizeof(SETUP));
  
  read_setup("eikonal.in",setup);
  
  if(setup->do_fdsolve)
    fdsolve(setup);
  
  if (setup->do_raytrace)
    gradient_descent_3d(setup);
  
  
  if(setup->do_pmf)
    calc_pmf(setup);
}


  
