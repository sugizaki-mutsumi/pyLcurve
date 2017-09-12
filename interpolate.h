#ifndef _INTERPOLATE_H_
#define _INTERPOLATE_H_


double interpolate_xyarray(int nv, double *vx, double *vy, double x);
int interpolate_xyearray(int nv, double *vx, double *vy, double *vye, double x, double *y, double *ye);

#endif /* _INTERPOLATE_H_ */
