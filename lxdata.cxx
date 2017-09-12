#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "fitsio.h"
#include "atFunctions.h"

#include "interpolate.h"
#include "lxdata.h"

//#define GW_LIMIT  5.0

Lxdata::Lxdata(){

  m_nbuf = 0;
  m_mjd_start = 0.0;
  m_mjd_end   = 0.0;

  m_poweridx = 0.0;  
  m_lxlowlimit = 0.0;
  m_lxoffset = 0.0;

  return;
}


Lxdata::~Lxdata(){
  //free_buf();
  return;
}


int Lxdata::free_buf(){
  if (m_lx)      delete m_lx;
  if (m_lxpow)   delete m_lxpow;
  if (m_integlxpow)  delete m_integlxpow;
  if (m_mjd)     delete m_mjd;
  if (m_binsize) delete m_binsize;
  m_nbuf = 0;
  return 1;
}

int Lxdata::init_buf(double mjd_start, double mjd_end, int nbuf){
  int idx;
  double fixbinsize;
  
  m_nbuf = nbuf;
  m_lx      = new double[nbuf]; 
  m_lxe     = new double[nbuf]; 
  m_lxpow   = new double[nbuf]; 
  m_integlxpow = new double[nbuf]; 
  m_mjd     = new double[nbuf]; 
  m_binsize = new double[nbuf]; 
  m_mjd_start = mjd_start;
  m_mjd_end   = mjd_end;
  fixbinsize = (mjd_end-mjd_start)/nbuf;

  for(idx=0; idx<nbuf; idx++){
    m_lx[idx]    = 0.0;
    m_lxe[idx]   = 0.0;
    m_lxpow[idx] = 0.0;
    m_integlxpow[idx] = 0.0;
    m_mjd[idx]   = mjd_start + (mjd_end-mjd_start)/nbuf*idx; 
    m_binsize[idx] = fixbinsize;
  }
  return 1;
}


//int Lxdata::set_lxdata_fromlctxt(char *lcfname, double mjd_start, double mjd_end, int nbuf, double lx_lowlimit, double lx_offset){
//int Lxdata::set_lxdata_fromlctxt(char *lcfname, double mjd_start, double mjd_end, int nbuf){
int Lxdata::read_textfile(char *lcfname, double f2lx){

  int idx, ndat;
  double mjd;
  FILE *fp;
  char line[FILENAME_MAX];


  if(m_nbuf==0){
    printf("Error in Lxdata: data buffer is not initialized!\n");
    return 1;
  }


  const int MAXNROWS = 1000000;
  double *vmjd = new double[MAXNROWS];
  double *vlx  = new double[MAXNROWS];
  double *vlxe = new double[MAXNROWS];


  /// read LC text file
  fp = fopen(lcfname, "r");
  if(fp == NULL){
    fprintf(stderr, "fopen error\n");
    return -1;
  }

  idx = 0;
  while( fgets(line, sizeof(line), fp) != NULL){
    if (line[0]=='#') continue;
    if(sscanf(line, "%lf %lf", &vmjd[idx], &vlx[idx])==2){
      vlxe[idx]=0.0;
      idx++;
    }
    if(MAXNROWS<=idx) {
      break;
    }
  }
  fclose(fp);

  ndat = idx;
  /// convert flux to lumi
  for(idx=0; idx<ndat; idx++){
    vlx[idx]*=f2lx;
    vlxe[idx]*=f2lx;
  }

  set_lxdata(ndat, vmjd, vlx, vlxe);

  delete[] vmjd; 
  delete[] vlx; 
  delete[] vlxe; 

  return 1;
}

int Lxdata::read_fitsfile(char *fname, char *extname, char *timecol, char *ratecol, char *rerrcol, double f2lx){
  
  int status;
  long nrows, idx; 
  int  ncols, ncol_time, ncol_rate, ncol_rerr; //, colnum;
  int  hdunum, hdupos, hdutype, anynul;
  fitsfile *fptr;

  status = 0;

  fits_open_table(&fptr, fname, READONLY, &status);
  fits_movnam_hdu(fptr, BINARY_TBL, extname, 0, &status);
  fits_get_num_rows(fptr, &nrows, &status);
  fits_get_num_cols(fptr, &ncols, &status);
  fits_get_colnum(fptr, CASEINSEN, timecol, &ncol_time, &status);
  fits_get_colnum(fptr, CASEINSEN, ratecol, &ncol_rate, &status);
  fits_get_colnum(fptr, CASEINSEN, rerrcol, &ncol_rerr, &status);
  printf("lcfits file=%s, nrow=%ld, ncols=%d\n", fname, nrows, ncols);
  printf("       ncol_time=%d, ncol_rate=%d, ncol_rerr=%d\n", ncol_time, ncol_rate, ncol_rerr);
  printf("       status=%d\n", status);

  double *vmjd = new double[nrows];
  double *vlx  = new double[nrows];
  double *vlxe = new double[nrows];
  
  fits_read_col(fptr, TDOUBLE, ncol_time, 1, 1, nrows, NULL, vmjd, &anynul, &status);
  fits_read_col(fptr, TDOUBLE, ncol_rate, 1, 1, nrows, NULL, vlx,  &anynul, &status);
  fits_read_col(fptr, TDOUBLE, ncol_rerr, 1, 1, nrows, NULL, vlxe, &anynul, &status);
  
  /// convert flux to lumi
  for(idx=0; idx<nrows; idx++){
    vlx[idx]*=f2lx;
    vlxe[idx]*=f2lx;
  }

  set_lxdata((int)nrows, vmjd, vlx, vlxe);

  fits_close_file(fptr, &status);
  delete[] vmjd; 
  delete[] vlx; 
  delete[] vlxe; 

  return 1;

}


int Lxdata::read_gsclcfile(char *fname, char *extname, double snr_llim, double f2lx){
  
  int status;
  long nrows, irow, idx; 
  int  ncols, ncol_time, ncol_rate, ncol_rerr; //, ncol_dflag; //, colnum;
  int  hdunum, hdupos, hdutype, anynul;
  double mjd, rate, rerr;
  short  dflag;

  fitsfile *fptr;

  status = 0;

  fits_open_table(&fptr, fname, READONLY, &status);
  fits_movnam_hdu(fptr, BINARY_TBL, extname, 0, &status);
  fits_get_num_rows(fptr, &nrows, &status);
  fits_get_num_cols(fptr, &ncols, &status);
  fits_get_colnum(fptr, CASEINSEN, "MJD", &ncol_time, &status);
  fits_get_colnum(fptr, CASEINSEN, "RATE", &ncol_rate, &status);
  fits_get_colnum(fptr, CASEINSEN, "RERR", &ncol_rerr, &status);
  //fits_get_colnum(fptr, CASEINSEN, "DATA_FLAG", &ncol_dflag, &status);
  printf("lcfits file=%s, nrow=%ld, ncols=%d\n", fname, nrows, ncols);
  printf("       ncol_time=%d, ncol_rate=%d, ncol_rerr=%d\n", ncol_time, ncol_rate, ncol_rerr);
  printf("       status=%d\n", status);

  double *vmjd = new double[nrows];
  double *vlx  = new double[nrows];
  double *vlxe = new double[nrows];
  
  idx = 0;
  for(irow=0; irow<nrows; irow++){
    fits_read_col(fptr, TDOUBLE, ncol_time,  irow+1, 1, 1, NULL, &mjd,   &anynul, &status);
    fits_read_col(fptr, TDOUBLE, ncol_rate,  irow+1, 1, 1, NULL, &rate,  &anynul, &status);
    fits_read_col(fptr, TDOUBLE, ncol_rerr,  irow+1, 1, 1, NULL, &rerr,  &anynul, &status);
    //fits_read_col(fptr, TSHORT,  ncol_dflag, irow+1, 1, 1, NULL, &dflag, &anynul, &status);
    if ((rate/rerr)>snr_llim){
      vlx[idx] =rate*f2lx;
      vlxe[idx]=rerr*f2lx;
    } else {
      vlx[idx]=0.0;
      vlxe[idx]=0.0;
    }
    vmjd[idx]=mjd;
    idx++;
  }
  set_lxdata((int)idx, vmjd, vlx, vlxe);

  fits_close_file(fptr, &status);
  delete[] vmjd; 
  delete[] vlx; 

  return 1;

}

int Lxdata::read_batlcfile(char *fname, double snr_llim, double f2lx){
  
  int status;
  long nrows, irow, idx; 
  int  ncols, ncol_time, ncol_rate, ncol_rerr, ncol_dflag; //, colnum;
  int  hdunum, hdupos, hdutype, anynul;
  double mjd, rate, rerr;
  short  dflag;

  fitsfile *fptr;

  status = 0;

  fits_open_table(&fptr, fname, READONLY, &status);
  fits_movnam_hdu(fptr, BINARY_TBL, "RATE", 0, &status);
  fits_get_num_rows(fptr, &nrows, &status);
  fits_get_num_cols(fptr, &ncols, &status);
  fits_get_colnum(fptr, CASEINSEN, "TIME", &ncol_time, &status);
  fits_get_colnum(fptr, CASEINSEN, "RATE", &ncol_rate, &status);
  fits_get_colnum(fptr, CASEINSEN, "ERROR", &ncol_rerr, &status);
  fits_get_colnum(fptr, CASEINSEN, "DATA_FLAG", &ncol_dflag, &status);
  printf("lcfits file=%s, nrow=%ld, ncols=%d\n", fname, nrows, ncols);
  printf("       ncol_time=%d, ncol_rate=%d, ncol_rerr=%d\n", ncol_time, ncol_rate, ncol_rerr);
  printf("       status=%d\n", status);

  double *vmjd = new double[nrows];
  double *vlx  = new double[nrows];
  double *vlxe = new double[nrows];
  
  idx = 0;
  for(irow=0; irow<nrows; irow++){
    fits_read_col(fptr, TDOUBLE, ncol_time,  irow+1, 1, 1, NULL, &mjd,   &anynul, &status);
    fits_read_col(fptr, TDOUBLE, ncol_rate,  irow+1, 1, 1, NULL, &rate,  &anynul, &status);
    fits_read_col(fptr, TDOUBLE, ncol_rerr,  irow+1, 1, 1, NULL, &rerr,  &anynul, &status);
    fits_read_col(fptr, TSHORT,  ncol_dflag, irow+1, 1, 1, NULL, &dflag, &anynul, &status);
    if(dflag==0){
      if ((rate/rerr)>snr_llim){
	vlx[idx] =rate*f2lx;
	vlxe[idx]=rerr*f2lx;
      } else {
	vlx[idx] =0.0;
	vlxe[idx]=0.0;
      }
      vmjd[idx]=mjd;
      idx++;
    }
  }
  set_lxdata((int)idx, vmjd, vlx, vlxe);

  fits_close_file(fptr, &status);
  delete[] vmjd; 
  delete[] vlx; 

  return 1;

}


int Lxdata::set_lxdata(int ndat, double *vmjd, double *vlx, double *vlxe){

  int idx;
  double mjd;
  double lx, lxe;

  if(m_nbuf==0){
    printf("Error in Lxdata: data buffer is not initialized!\n");
    return 1;
  }

  //init_buf(mjd_start, mjd_end, nbuf);
  for(idx=0; idx<m_nbuf; idx++){
    mjd = m_mjd[idx]+m_binsize[idx]/2.0;
    //m_lx[idx] = interpolate_xyarray(ndat, vmjd, vlx, mjd);
    interpolate_xyearray(ndat, vmjd, vlx, vlxe, mjd, &lx, &lxe);
    m_lx[idx] =lx;
    m_lxe[idx]=lxe;
  }
  return 1;
}


/*
int Lxdata::set_lxdata_fromgaus(char *parfname, double mjd_start, double mjd_end, int nbuf){
  int idx;
  double gn, gw, gc;
  FILE *fp;
  char line[FILENAME_MAX];

  init_buf(mjd_start, mjd_end, nbuf);

  /// read parameter file
  fp = fopen(parfname, "r");
  if(fp == NULL){
    fprintf(stderr, "fopen error\n");
    return -1;
  }
  idx=0;
  while( fgets(line, sizeof(line), fp) != NULL){
    if (line[0]=='#') continue;
    if(sscanf(line, "%lf %lf %lf", &gn, &gw, &gc)==3){
      //printf("#%d: gn=%lf, gw=%lf, gc=%lf\n", idx, gn, gw, gc); 
      set_lxdata_fillgaus(gn, gw, gc);
      idx++;
    }
  }
  fclose(fp);

  /// set integrals of Lx^6/7
  //set_integlxpow(0.0, 0.0);
  return 1;
}

int Lxdata::set_lxdata_fillgaus(double gn, double gw, double gc){

  int idx;
  double mjd, z; //lx, 

  for(idx=0; idx<m_nbuf; idx++){
    mjd = m_mjd[idx] + m_binsize[idx]*0.5;
    z = (mjd-gc)/gw;
    //printf("idx=%d, mjd=%lf, z=%lf\n", idx, mjd, z);
    if(fabs(z)<GW_LIMIT){
      m_lx[idx] += gn*exp(-0.5*z*z);
    }
  }
  return 1;
}
*/


int Lxdata::set_integlxpow(double lxlowlimit, double lxoffset, double poweridx){

  int idx;
  double integval;
  
  if ( fabs(m_lxlowlimit-lxlowlimit)<1e-6 && fabs(m_lxoffset-lxoffset)<1e-6 && fabs(m_poweridx-poweridx)<1e-6){
    return 1;
  }
  //printf("poweridx = %lf\n", poweridx);

  integval = 0.0;
  for(idx=0; idx<m_nbuf; idx++){
    if (lxlowlimit<m_lx[idx]){
      //m_lxpow[idx] = pow(m_lx[idx]-lx_offset, m_poweridx);
      m_lxpow[idx] = pow(m_lx[idx]-lxoffset, poweridx);
    } else {
      m_lxpow[idx] = 0.0;
    }
    integval += m_lxpow[idx] * m_binsize[idx];
    m_integlxpow[idx] = integval;
  }

  m_lxlowlimit = lxlowlimit;
  m_lxoffset   = lxoffset;
  m_poweridx   = poweridx;

  return 1;
}

double Lxdata::calc_lx(double mjd){
  int idx;
  double lx, w;
  if( mjd<m_mjd_start || m_mjd_end<=mjd ) {
    return 0.0;
  }
  for(idx=0; idx<(m_nbuf-1); idx++){
    if(m_mjd[idx]<=mjd && mjd<m_mjd[idx+1]){
      ///lx = m_lx[idx]; /// need modify
      w = (mjd-m_mjd[idx])/(m_mjd[idx+1]-m_mjd[idx]);
      lx = (1-w)*m_lx[idx] + w*m_lx[idx+1];
      break;
    }
  }
  return lx;
}

double Lxdata::calc_lxe(double mjd){
  int idx;
  double lxe, w;
  if( mjd<m_mjd_start || m_mjd_end<=mjd ) {
    return 0.0;
  }
  for(idx=0; idx<(m_nbuf-1); idx++){
    if(m_mjd[idx]<=mjd && mjd<m_mjd[idx+1]){
      ///lx = m_lx[idx]; /// need modify
      w = (mjd-m_mjd[idx])/(m_mjd[idx+1]-m_mjd[idx]);
      lxe = sqrt((1-w)*(1-w)*m_lxe[idx]*m_lxe[idx] + w*w*m_lxe[idx+1]*m_lxe[idx+1]);
      break;
    }
  }
  return lxe;
}


//double Lxdata::calc_integlxpow(double mjd_base, double mjd){
double Lxdata::calc_integlxpow(double mjd){
  int idx;
  double integlxpow;
  if( mjd<m_mjd_start) {
    return m_integlxpow[0];
  } else if (m_mjd[m_nbuf-1]<=mjd ) { //(m_mjd_end<=mjd ) {
    return m_integlxpow[m_nbuf-1];
  }
  
  for(idx=0; idx<(m_nbuf-1); idx++){
    if(m_mjd[idx]<=mjd && mjd<m_mjd[idx+1]){
      integlxpow = m_integlxpow[idx]; /// need modify
      break;
    }
  }
  return integlxpow;
}


