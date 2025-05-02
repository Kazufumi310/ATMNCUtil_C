#ifndef Histgrammer_h__
#define Histgrammer_h__

#include <ContainerClass.h>
#include <atmnc-constants.h>
#include <atmnc-particle-code.h>
#include <SpectrumHandler.h>

class Histogrammer{
public:
  Histogrammer(ATMNCInputInfo input);

  // to get/add Histogram
  void getHistogram(int ikind, int isite, int icosz, int iazim,
		    int &num, double ene[], double nobs[], double err[]);
  void addHistogram(int ikind, int isite, int icosz, int iazim,
		    int num, double ene[], double nobs[], double err[]);
  
  // set/add data
  void setTequiv(double tequiv){m_tequiv = tequiv;}
  void setTequiv(TTree *trin);
  void setDarea(int isite, double darea){ 
    m_Darea[isite] = darea; 
    std::cout<<"reset Darea["<<isite<<"] to "<<darea<<std::endl;
  }
  
  // get the aggregated spectrum 
  void getSpectrum(int ik, int is, int ic, int ia, SpectrumHandler &spec);
  void averageFlux(int ik,int is,int icmin,int icmax,int iamin,int iamax, SpectrumHandler &spec);
  double getTequiv(){return m_tequiv;}
  double getDarea(int isite){return m_Darea[isite];}
  double getDomega(int icosz){return m_Domega[icosz];}
  double getNormFactor(int is,int ic);

  // functions to get histogram's index
  int getI(int ikind,int isite,int icosz,int iazim);
  double getCosz(double ic);
  double getAzim(double ia);
  double getEne(int ie){return m_ene[ie];}
  int getCoszIndex(double theta);
  int getAzimIndex(double phi);
  int getEneIndex(double ene);
  int getNKind(){return m_nkind;}
  int getNSite(){return m_nsite;}
  int getNCosz(){return m_ncosz;}
  int getNAzim(){return m_nazim;}
  int getNBin(){return m_nbin;}
  int getPID(int ik);
  double getFactorOfLastAverage(){return m_factorOfLastAverage;}
  
  //private:
 protected:
  void getDirection(TVector3 const &mom, TVector3 const &pos,
		    double &theta,double &phi);
  int  m_nkind;//default = 6;
  int  m_nsite;//default = 8;
  int  m_ncosz;//default = 40;
  int  m_nazim;//default = 12;
  int  m_nbin;// default = 103;
  double m_cmin; 
  //  double m_cmax;
  double m_emin;
  double m_Letmi;
  int m_pid2ik[ParticleCode::nParType];

  double m_Darea[128];//[nsite]
  double m_Domega[128];//[ncosz]
  double m_tequiv;
  double m_factorOfLastAverage;
  
  double **m_nobs;// [nkind*nsite*ncosz*nazim][nbin]
  double **m_err2; // [nkind*nsite*ncosz*nazim][nbin]
  double *m_ene;  // [nbin]

};

#endif //Histgrammer_h__
