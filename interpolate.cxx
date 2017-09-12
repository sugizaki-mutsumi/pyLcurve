#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

double interpolate_xyarray(int nv, double *vx, double *vy, double x){
  int idx;
  double y;
  if(x<vx[0]) return vy[0];
  else if(vx[nv-1]<=x) return vy[nv-1];

  for(idx=0; idx<(nv-1); idx++){
    if(vx[idx]<=x && x<vx[idx+1]){
      y = ((vx[idx+1]-x)*vy[idx] + (x-vx[idx])*vy[idx+1])/(vx[idx+1]-vx[idx]);
      break;
    }
  }
  return y;
}


int interpolate_xyearray(int nv, double *vx, double *vy, double *vye, double x, double *y, double *ye){

  int idx;
  double a0, a1;

  if(x<vx[0]) {
    *y  = vy[0];      
    *ye = vye[0];      
    return 1;
  } else if(vx[nv-1]<=x) {
    *y  = vy[nv-1];      
    *ye = vye[nv-1];      
    return 1;
  }

  for(idx=0; idx<(nv-1); idx++){
    if(vx[idx]<=x && x<vx[idx+1]){
      a0 = (vx[idx+1]-x)/(vx[idx+1]-vx[idx]);
      a1 = (x-vx[idx])/(vx[idx+1]-vx[idx]);
      *y  = a0*vy[idx] + a1*vy[idx+1];
      *ye = sqrt(a0*a0*vye[idx]*vye[idx] + a1*a1*vye[idx+1]*vye[idx+1]);
      break;
    }
  }
  return 1;
}


