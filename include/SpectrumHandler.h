#ifndef SpectrumHandler_h__
#define SpectrumHandler_h__

#include <myheader/myRoot.h>
#include <myheader/myStandard.h>

int binaryLowBoundSearch(double X,int N,double *val);

class SpectrumHandler{
 public:
  SpectrumHandler();
  
  // set data
  bool setFlux(int n, double *E,double *flux, double *fluxE = 0,bool doIntegral = true);
  bool setRate(int n, double *E,double *rate, double *rateE = 0);


  bool isEmpty(){ return (m_nPoint==0)? false : true; }
  int getN(){return m_nPoint;}
  double bin2Ene(int ibin){ return (ibin<0 || ibin>=m_nPoint)? 0 : m_E[ibin];}

  // get flux, rate histogram, and its integral
  double getFlux(double E, bool err = false); 
  double getFluxErr(double E){ return getFlux(E,true);}
  double getRate(double E, bool err = false);
  double getRateErr(double E){ return getFlux(E,true);}
  double getIntegral(double E);
  double getIntegAll(){return (m_nPoint==0)? 0 : m_integ[0];}
  double getRandomE(double rand); // rand: 0--1 

  double *getEneArray(){return m_E;}
  double *getFluxArray(){return m_flux;}
  double *getFluxErrArray(){return m_fluxErr;}
  double *getRateArray(){return m_rate;}
  double *getRateErrArray(){return m_rateErr;}
  double *getIntegArray(){return m_integ;}

  // 
  void printFlux(std::ofstream &ofs);
  void printRate(std::ofstream &ofs,double factor);

 private:
  void initArrays();
  void convertFluxToRate();
  bool convertRateToFlux();
  void integrateRate();

  int m_nPoint;
  static int const MAXARR = 1024;
  double m_E[MAXARR]; // [nPoint], GeV
  double m_flux[MAXARR]; // [nPoint], 1/GeV/m^2/sr/sec
  double m_fluxErr[MAXARR]; // [nPoint], error of flux
  double m_rate[MAXARR]; // [nPoint], 1/m^2/sr/sec, m_rate[i] = integral flux btwn m_E[i]--m_E[i+1]
  double m_rateErr[MAXARR]; // [nPoint], error of rate
  double m_integ[MAXARR]; // [nPoint], 1/m^2/sr/sec, m_integ[i] = integral of m_rate from maximum E to m_E[i].

};

#endif // SpectrumHandler_h__
