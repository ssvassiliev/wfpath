float weight;
float origin [3];
float xdelta[3];            // x axis of unit cell
float ydelta[3];            // y axis of unit cell
float zdelta[3];            // z axis of unit cell
float xaxis[3];             // direction and length for X axis (non-unit)
float yaxis[3];             // direction and length for Y axis (non-unit)
float zaxis[3];             // direction and length for Z axis (non-unit)
int xsize, ysize, zsize, N; // number of samples along each axis
float *data;                // raw data, total of xsize*ysize*zsize voxels
float ***M;                 // 3D format data
float *timefld;             // first arrival time
float ***GX_i, ***GY_i, ***GZ_i;  // gradient intermediate differences
float ***GX_c, ***GY_c, ***GZ_c;  // gradient central differences
float ***GX_26, ***GY_26, ***GZ_26;  // gradient central differences     
float min_E, max_E;         // data min & max



