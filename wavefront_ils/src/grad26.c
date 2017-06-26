void grad26(float ***GX, float ***GY, float ***GZ, int x, int y, int z) {
  //  float g1[3],g2[3],g[3];
#define r2  7.071067812e-1
#define r3  5.773502692e-1

  int inx,ipx,iny,ipy,inz,ipz;
  float g1[3][3][3];
  float g2[3][3][3];
  float gCenterX, gSideX, gCornX;
  float gCenterY, gSideY, gCornY;
  float gCenterZ, gSideZ, gCornZ;

  inx=ipx=iny=ipy=inz=ipz=0;
  if(x > 0) inx = 1;
  if(y > 0) iny = 1;
  if(z > 0) inz = 1;
  if(x < xsize-1) ipx = 1;
  if(y < ysize-1) ipy = 1;
  if(z < zsize-1) ipz = 1;

  // X-center
  g1[1][0][0] = M[x-inx][y][z];
  g2[1][0][0] = M[x+ipx][y][z];
  gCenterX = g2[1][0][0]-g1[1][0][0];
  // X-side
  g1[1][1][0] = M[x-inx][y-iny][z];
  g2[1][1][0] = M[x+ipx][y-iny][z];
  g1[1][2][0] = M[x-inx][y+ipy][z];
  g2[1][2][0] = M[x+ipx][y+ipy][z];
  g1[1][0][1] = M[x-inx][y][z-inz];
  g2[1][0][1] = M[x+ipx][y][z-inz];
  g1[1][0][2] = M[x-inx][y][z+ipz];
  g2[1][0][2] = M[x+ipx][y][z+ipz];
  gSideX = \
    g2[1][1][0] - g1[1][1][0] +			\
    g2[1][2][0] - g1[1][2][0] +			\
    g2[1][0][1] - g1[1][0][1] +			\
    g2[1][0][2] - g1[1][0][2];
  // X-corner
  g1[1][1][0] = M[x-inx][y-iny][z-inz];
  g2[1][1][0] = M[x+ipx][y-iny][z-inz];
  g1[1][2][0] = M[x-inx][y+ipy][z+ipz];
  g2[1][2][0] = M[x+ipx][y+ipy][z+ipz];
  g1[1][0][1] = M[x-inx][y-iny][z-inz];
  g2[1][0][1] = M[x+ipx][y-iny][z-inz];
  g1[1][0][2] = M[x-inx][y+ipy][z+ipz];
  g2[1][0][2] = M[x+ipx][y+ipy][z+ipz];
  gCornX = \
    g2[1][1][0] - g1[1][1][0] +			\
    g2[1][2][0] - g1[1][2][0] +			\
    g2[1][0][1] - g1[1][0][1] +			\
    g2[1][0][2] - g1[1][0][2];
  // Y-center
  g1[0][1][0] = M[x][y-iny][z];
  g2[0][1][0] = M[x][y+ipy][z];
  gCenterY = g2[0][1][0]-g1[0][1][0];
  // Y-side
  g1[1][1][0] = M[x-inx][y-iny][z];
  g2[1][1][0] = M[x-inx][y+ipy][z];
  g1[1][2][0] = M[x+ipx][y-iny][z];
  g2[1][2][0] = M[x+ipx][y+ipy][z];
  g1[1][0][1] = M[x][y-iny][z-inz];
  g2[1][0][1] = M[x][y+ipy][z-inz];
  g1[1][0][2] = M[x][y-iny][z+ipz];
  g2[1][0][2] = M[x][y+ipy][z+ipz];
gSideY = \
    g2[1][1][0] - g1[1][1][0] +			\
    g2[1][2][0] - g1[1][2][0] +			\
    g2[1][0][1] - g1[1][0][1] +			\
    g2[1][0][2] - g1[1][0][2];
 // Y-corner
  g1[1][1][0] = M[x-inx][y-iny][z-inz];
  g2[1][1][0] = M[x-inx][y+ipy][z-inz];
  g1[1][2][0] = M[x+ipx][y-iny][z+ipz];
  g2[1][2][0] = M[x+ipx][y+ipy][z+ipz];
  g1[1][0][1] = M[x-inx][y-iny][z-inz];
  g2[1][0][1] = M[x-inx][y+ipy][z-inz];
  g1[1][0][2] = M[x+ipx][y-iny][z+ipz];
  g2[1][0][2] = M[x+ipx][y+ipy][z+ipz];
gCornY = \
    g2[1][1][0] - g1[1][1][0] +			\
    g2[1][2][0] - g1[1][2][0] +			\
    g2[1][0][1] - g1[1][0][1] +			\
    g2[1][0][2] - g1[1][0][2];

  // Z-center
  g1[0][0][1] = M[x][y][z-inz]; 
  g2[0][0][1] = M[x][y][z+ipz];
  gCenterZ = g2[0][0][1]-g1[0][0][1];
 // Z-side
  g1[1][1][0] = M[x-inx][y][z-inz];
  g2[1][1][0] = M[x-inx][y][z+ipz];
  g1[1][2][0] = M[x+ipx][y][z-inz];
  g2[1][2][0] = M[x+ipx][y][z+ipz];
  g1[1][0][1] = M[x][y-iny][z-inz];
  g2[1][0][1] = M[x][y-iny][z+ipz];
  g1[1][0][2] = M[x][y+ipy][z-inz];
  g2[1][0][2] = M[x][y+ipy][z+ipz];
gSideZ = \
    g2[1][1][0] - g1[1][1][0] +			\
    g2[1][2][0] - g1[1][2][0] +			\
    g2[1][0][1] - g1[1][0][1] +			\
    g2[1][0][2] - g1[1][0][2];
 // Z-corner
  g1[1][1][0] = M[x-inx][y-iny][z-inz];
  g2[1][1][0] = M[x-inx][y-iny][z+ipz];
  g1[1][2][0] = M[x+ipx][y+ipy][z-inz];
  g2[1][2][0] = M[x+ipx][y+ipy][z+ipz];
  g1[1][0][1] = M[x-inx][y-iny][z-inz];
  g2[1][0][1] = M[x-inx][y-iny][z+ipz];
  g1[1][0][2] = M[x+ipx][y+ipy][z-inz];
  g2[1][0][2] = M[x+ipx][y+ipy][z+ipz];
gCornZ = \
    g2[1][1][0] - g1[1][1][0] +			\
    g2[1][2][0] - g1[1][2][0] +			\
    g2[1][0][1] - g1[1][0][1] +			\
    g2[1][0][2] - g1[1][0][2];

  /* if(inx && ipx) */ GX[x][y][z] = (gCenterX + gSideX*0.25*r2 + gCornX*0.25*r3); 
  /* if(iny && ipy) */ GY[x][y][z] = (gCenterY + gSideY*0.25*r2 + gCornY*0.25*r3);
  /* if(inz && ipz) */ GZ[x][y][z] = (gCenterZ + gSideZ*0.25*r2 + gCornZ*0.25*r3);
}


