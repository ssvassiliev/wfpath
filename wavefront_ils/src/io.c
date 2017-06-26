// Get a string from a stream, printing any errors that occur
char *dxgets(char *s, int n, FILE *stream) {
  char *returnVal;

  if (feof(stream)) {
    fprintf(stderr, "Unexpected end-of-file.\n");
    return NULL;
  } else if (ferror(stream)) {
    fprintf(stderr, "Error reading file.\n");
    return NULL;
  } else {
    returnVal = fgets(s, n, stream);
    if (returnVal == NULL) {
      fprintf(stderr, "Error reading line (error = %d).\n", errno);
    }
  }
  return returnVal;
}

int read_dx_file (const char *filename) {
  const int LINESIZE=512;
  char inbuf[LINESIZE];
  int i;
  
  printf("Loading DX map from file \"%s\"\n", filename);
  weight=1;
  
  FILE *fd = fopen(filename, "r");
  if (!fd) {
    fprintf(stderr, "Error: Could not open file \"%s\" for reading\n", filename);
    return -1;
  };
  

  if (!dxgets(inbuf, LINESIZE, fd)) return -1;
  if(!strncmp(inbuf,"# weight",8))
    sscanf(&inbuf[9],"%g",&weight);

 /* skip comments */
  do {
    if (!dxgets(inbuf, LINESIZE, fd)) return -1;
  } while (inbuf[0] == '#');
 
        
  /* get the number of grid points */
  if (sscanf(&inbuf[35], "%i %i %i", &xsize, &ysize, &zsize) != 3) {
    fprintf(stderr, "Error reading grid dimensions.\n");
    return -1;
  }

  /* get the cell origin */
  if (dxgets(inbuf, LINESIZE, fd) == NULL) {
    return -1;
  }
  if (sscanf(&inbuf[6], "%f %f %f", origin, origin+1, origin+2) != 3) {
    fprintf(stderr, "Error reading grid origin.\n");
    return -1;
  }

  /* get the cell dimensions */
  if (dxgets(inbuf, LINESIZE, fd) == NULL) return -1;
  if (sscanf(&inbuf[5], "%f %f %f", xdelta, xdelta+1, xdelta+2) != 3) {
    fprintf(stderr, "Error reading cell x-dimension.\n");
    return -1;
  }

  if (dxgets(inbuf, LINESIZE, fd) == NULL) return -1;
  if (sscanf(&inbuf[5], "%f %f %f", ydelta, ydelta+1, ydelta+2) != 3) {
    fprintf(stderr, "Error reading cell y-dimension.\n");
    return -1;
  }

  if (dxgets(inbuf, LINESIZE, fd) == NULL) return -1;
  if (sscanf(&inbuf[5], "%f %f %f", zdelta, zdelta+1, zdelta+2) != 3) {
    fprintf(stderr, "Error reading cell z-dimension.\n");
    return -1;
  }

  /* skip the last two lines of the header*/
  if (dxgets(inbuf, LINESIZE, fd) == NULL) return -1;
  if (dxgets(inbuf, LINESIZE, fd) == NULL) return -1;

  /* Set the unit cell origin and basis vectors */
  for (i=0; i<3; i++) {
    xaxis[i] = xdelta[i] * (xsize-1);
    yaxis[i] = ydelta[i] * (ysize-1);
    zaxis[i] = zdelta[i] * (zsize-1);
  }  
  
  /* Read the values from the file */
  int xysize = xsize*ysize;
  int gridsize = xsize*ysize*zsize;
  int gx, gy, gz; 
  float grid[3];
  int count;
    
  data = malloc(gridsize*sizeof(float));
      
  gx = gy = gz = 0;  
  for (count=0; count < gridsize/3; count++) {
    if (dxgets(inbuf, LINESIZE, fd) == NULL ) return -1;
    
    if (sscanf(inbuf, "%f %f %f", grid, grid+1, grid+2) != 3) {
      fprintf(stderr, "Error reading grid data.\n");
      return -1;
    }
  
    for (i=0; i < 3; i++) { 
      data[gx + gy*xsize + gz*xysize] = grid[i];
      gz++;
      if (gz >= zsize) {
        gz = 0;
        gy++;
        if (gy >= ysize) {
          gy = 0;
          gx++;
        }
      }
    }
    
  }

  // This reads the last data line, if it only contains 1 or 2 voxels
  if (gridsize%3) {
    if (dxgets(inbuf, LINESIZE, fd) == NULL )
      return -1;

    count = sscanf(inbuf, "%f %f %f", grid, grid+1, grid+2);
    if (count != (gridsize%3)) {
      fprintf(stderr, "Error: incorrect number of data points.\n");
      return -1;
    }

    for (i=0; i<count; i++) {
      data[gx + gy*xsize + gz*xysize] = grid[i];
      gz++;
    }
  }
  
  fclose(fd);
  
  printf("%d voxels\n", gridsize);
  
  return 0;
   
}

static inline double voxel_value(int gx, int gy, int gz, float *dt) {
  return dt[gx + gy*xsize + gz*ysize*xsize];
}


int  write_dx_file (const char *filename, float *dt){
  int gridsize = xsize*ysize*zsize;
  int i;
  
  printf(" Writing DX map to file \"%s\"\n", filename);
  
  FILE *fout = fopen(filename, "w");
  if (!fout) {
    fprintf(stderr, "volmap: Error: Cannot open file \"%s\" for writing\n", filename);
    return -1;
  };
    
  fprintf(fout, "# Data calculated by volutil (http://www.ks.uiuc.edu)\n");
  fprintf(fout, "# VMDTAG WEIGHT %g\n", weight);
 
  fprintf(fout, "object 1 class gridpositions counts %d %d %d\n", xsize, ysize, zsize);
  fprintf(fout, "origin %g %g %g\n", origin[0], origin[1], origin[2]);
  fprintf(fout, "delta %g %g %g\n", xdelta[0], xdelta[1], xdelta[2]);
  fprintf(fout, "delta %g %g %g\n", ydelta[0], ydelta[1], ydelta[2]);
  fprintf(fout, "delta %g %g %g\n", zdelta[0], zdelta[1], zdelta[2]);
  fprintf(fout, "object 2 class gridconnections counts %d %d %d\n", xsize, ysize, zsize);
  fprintf(fout, "object 3 class array type double rank 0 items %d data follows\n", gridsize);
  
  // This reverses the ordering from x fastest to z fastest changing variable
  float val1,val2,val3;
  int gx=0, gy=0, gz=-1;
  for (i=0; i < (gridsize/3)*3; i+=3)  {
    if (++gz >= zsize) {
      gz=0;
      if (++gy >= ysize) {gy=0; gx++;}
    }
    val1 = voxel_value(gx,gy,gz,dt);
    if (++gz >= zsize) {
      gz=0;
      if (++gy >= ysize) {gy=0; gx++;}
    }
    val2 = voxel_value(gx,gy,gz,dt);
    if (++gz >= zsize) {
      gz=0;
      if (++gy >= ysize) {gy=0; gx++;}
    }
    val3 = voxel_value(gx,gy,gz,dt);    
    fprintf(fout, "%g %g %g\n", val1, val2, val3);
  }
  for (i=(gridsize/3)*3; i < gridsize; i++) {
    if (++gz >= zsize) {
      gz=0;
      if (++gy >= ysize) {gy=0; gx++;}
    }
    fprintf(fout, "%g ", voxel_value(gx,gy,gz,dt));
  }
  if (gridsize%3) fprintf(fout, "\n");
  
  fprintf(fout, "\n");

  // XXX todo: make sure that "dataname" contains no quotes
  fprintf(fout, "object \"%s\" class field\n", "volutil output");
  
  fclose(fout);
  return 0;
}

void print_stats() {  
  int gx, gy, gz;
  float  E;
 
  float sum_E = 0.; 
  min_E = data[0]; 
  max_E = data[0];  
  
  for (gx=0; gx<xsize; gx++)
  for (gy=0; gy<ysize; gy++) 
  for (gz=0; gz<zsize; gz++) {
    E = data[gx + gy*xsize + gz*xsize*ysize];
    sum_E += E;
    if (E<min_E) min_E = E;
    if (E>max_E) max_E = E;
  }
  
  N = xsize*ysize*zsize;
  
  printf("  WEIGHT:    %g\n", weight);
  printf("  AVERAGE:   %g\n", sum_E/N);
  printf("  MIN:       %g\n", min_E);
  printf("  MAX:       %g\n", max_E);
  printf("\n");
 
}
