#ifndef _lxdata_H_
#define _lxdata_H_

class Lxdata {
 private :
  //double    m_dt;

  int    m_nbuf;
  double m_mjd_start, m_mjd_end;

  double *m_lx;
  double *m_lxe;
  double *m_lxpow;
  double *m_integlxpow;
  double *m_mjd;
  double *m_binsize;
  
  double m_poweridx;
  double m_lxlowlimit;
  double m_lxoffset;


 public:
  Lxdata();
  virtual ~Lxdata();

  int free_buf();
  int init_buf(double mjd_start, double mjd_end, int nbuf);

  int read_textfile(char *fname, double f2lx);
  int read_fitsfile(char *fname, char *extname, char *timecol, char *ratecol, char *rerrcol, double f2lx);

  int read_batlcfile(char *fname, double snr_llim, double f2lx);
  int read_gsclcfile(char *fname, char *extname, double snr_llim, double f2lx);

  int set_lxdata(int ndat, double *vmjd, double *vlx, double *vlxe);

  //int set_lxdata_fromgaus(char *parfname, double mjd_start, double mjd_end, int nbuf);
  //int set_lxdata_fillgaus(double gn, double gw, double gc);
  int set_integlxpow(double lxlowlimit, double lxoffset, double poweridx);

  double get_mjd(int idx){return m_mjd[idx];}
  double get_lx(int idx){return m_lx[idx];}
  double get_lxpow(int idx){return m_lxpow[idx];}
  double get_integlxpow(int idx){return m_integlxpow[idx];}

  double get_poweridx(){return m_poweridx;}
  double get_lxlowlimit(){return m_lxlowlimit;}
  double get_lxoffset(){return m_lxoffset;}
  
  double calc_lx(double mjd);
  double calc_lxe(double mjd);
  double calc_integlxpow(double mjd);
  
  //double calc_deltap_ghlmb(double mjd_base, double mjd);

  //int    set_integ_dt(double dt);
  //double get_dt(){return m_dt;};
};

#endif /* _lxdata_H_ */
