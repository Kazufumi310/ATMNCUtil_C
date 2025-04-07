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

SpectrumHandler::SpectrumHandler(){
  m_nPoint = 0;
  initArrays();
}

void SpectrumHandler::initArrays(){
  for(int i=0;i<MAXARR;i++){
    m_E[i] = m_flux[i] = m_fluxErr[i] = m_rate[i] = m_rateErr[i] = m_integ[i] = 0;
  }
}

bool SpectrumHandler::setFlux(int n, double *E,double *flux,double *fluxE,bool doIntegral){
  m_nPoint = n;
  if(doIntegral) initArrays();
  else{
    for(int i=0;i<MAXARR;i++){
      m_flux[i] = m_fluxErr[i] = 0;
    }
  }
  for(int i=0;i<n;i++){
    if(doIntegral){
      m_E[i] = E[i];
    }else{
      if(m_E[i]!=E[i]){
	std::cout<<"SpectrumHandler::setFlux: warning: energy binning is different from the conventional one: ibin,m_E[ibin],E[ibin]:"<<i<<" "<<m_E[i]<<" "<<E[i]<<std::endl;
      }
    }
    m_flux[i] = flux[i];
    if(i>0 && m_E[i]<=m_E[i-1]){
      std::cout<<"error in setFlux: m_E["<<i<<"]<=m_E["<<i-1<<"]"<<std::endl;
      return false;
    }
  }

  if(fluxE!=0){
    for(int i=0;i<n;i++) m_fluxErr[i] = fluxE[i];
  }

  if(doIntegral){
    convertFluxToRate();
    integrateRate();
  }
  return true;
}

bool SpectrumHandler::setRate(int n, double *E,double *rate, double *rateE){
  m_nPoint = n;
  initArrays();
  for(int i=0;i<n;i++){
    m_E[i] = E[i];
    m_rate[i] = rate[i];
  }

  if(rateE!=0){
    for(int i=0;i<n;i++) m_rateErr[i] = rateE[i];
  }

  if(!convertRateToFlux()) return false;
  integrateRate();
  return true;
}

void SpectrumHandler::convertFluxToRate(){
  for(int i=0;i<m_nPoint;i++){
    if(i==m_nPoint-1){
      double frat = m_flux[i]/m_flux[i-1];
      double erat = m_E[i]/m_E[i-1];
      double A = log(frat)/log(erat);
      m_rate[i] = -m_E[i]*m_flux[i]/(1+A);
    }else{
      double frat = m_flux[i+1]/m_flux[i];
      double erat = m_E[i+1]/m_E[i];
      double A = log(frat)/log(erat);
      //      m_rate[i] =  m_E[i]*m_flux[i]/(1+A)*(frat*erat-1);
      m_rate[i] =  m_E[i]*m_flux[i]/(1+A)*(pow(erat,A+1)-1);
    }
  }
}


bool SpectrumHandler::convertRateToFlux(){
  double const delLogE = (log(m_E[m_nPoint-1])-log(m_E[0]))/(m_nPoint-1);
  double const ree = exp(delLogE);
  for(int i=1;i<m_nPoint;i++){
    // check : {m_E} should be in ascending order, and even-spaced in log scale
    if(m_E[i]<m_E[i-1]){
      std::cout<<"error: E["<<i<<"]<E["<<i-1<<"]"<<std::endl;
      return false;
    }
    double space = log(m_E[i])-log(m_E[i-1]);
    if(fabs(space - delLogE)/delLogE > 1e-3){
      std::cout<<"error: logE["<<i<<"]- logE["<<i-1<<"] != "<<delLogE<<std::endl;
      return false;
    }

    // calculate flux from neighboring 2 bins
    if(m_rate[i-1]<=0 || m_rate[i]<=0 ){
      double nn = std::max(0.,m_rate[i-1])+std::max(0.,m_rate[i]);
      double de = m_E[i]*ree-m_E[i]/ree;
      m_flux[i] = nn/de;
    }else if(m_rate[i-1]==m_rate[i]){
      m_flux[i] = m_rate[i]/(m_E[i]*delLogE);
    }else{
      m_flux[i] = (m_rate[i-1]*m_rate[i])/(m_rate[i]-m_rate[i-1])/m_E[i]
	*log(m_rate[i]/m_rate[i-1])/delLogE;
    }
    if(TMath::IsNaN(m_flux[i]) || !TMath::Finite(m_flux[i])){
      if(TMath::IsNaN(m_flux[i]))
	std::cout<<"flux NaN found"<<std::endl;
      else
	std::cout<<"flux Inf found"<<std::endl;
      std::cout<<"i,rate[i-1],rate[i]:"<<i<<" "<<m_rate[i-1]<<" "<<m_rate[i]<<std::endl;
      std::cout<<"E[i],delLogE:"<<m_E[i]<<" "<<delLogE<<std::endl;
      exit(1);
    }
    // for 0th bin
    if(m_rate[0]*m_rate[1] > 0){
      m_flux[0] = m_flux[1]*(m_rate[0]/m_rate[1])*ree;
    }else{
      m_flux[0] = 0;
    }
    
    // calcurate error
    if(m_rateErr[i-1]==0 && m_rateErr[i]==0 ) m_fluxErr[i] = 0;
    else if(m_rate[i-1]<=0||m_rate[i]<=0){
      //      double nn = std::max(0.,m_rate[i])+std::max(0.,m_rate[i-1]);
      double err;
      if(m_rate[i-1]>0) err = m_rateErr[i-1];
      else err = m_rateErr[i];
      double de = m_E[i]*ree-m_E[i]/ree;
      m_fluxErr[i] = err/de;
    }else if(m_rate[i-1]==m_rate[i]){ // alpha = 1
      m_fluxErr[i] = 1./(m_E[i]*delLogE)*m_rateErr[i];
    }else{
      double X = m_rate[i]/m_rate[i-1];
      double dFdE[2]; //d(flux)/dE
      dFdE[0] = 1./(m_E[i]*delLogE)*(pow(X,2)/pow(X-1,2)*log(X)-X/(X-1));
      dFdE[1] = 1./(m_E[i]*delLogE)*(       1/pow(X-1,2)*log(X)-1/(X-1));
      // 210414 found BUG! not m_rateErr[1], but m_rateErr[i] !!!!
      //      m_fluxErr[i] = pow(dFdE[0],2)*pow(m_rateErr[i-1],2)+pow(dFdE[1],2)*pow(m_rateErr[1],2);
      m_fluxErr[i] = pow(dFdE[0],2)*pow(m_rateErr[i-1],2)+pow(dFdE[1],2)*pow(m_rateErr[i],2);
      if(m_fluxErr[i]<0){
	std::cout<<"err < 0"<<std::endl;
	std::cout<<"E0,E1,N0,N1,dFdE0,dFdE1:"<<m_E[i-1]<<" "<<m_E[i]<<" "<<m_rate[i-1]<<" "<<m_rate[i]<<" "<<dFdE[0]<<" "<<dFdE[1]<<std::endl;
	m_fluxErr[i] = 1./(m_E[i]*delLogE)*m_rateErr[i];
      }else{
	m_fluxErr[i] = sqrt(m_fluxErr[i]);
      }
    }
  }
  // for 0th bin
  if(m_rate[0] > 0){
    m_fluxErr[0] = 1/sqrt(m_rate[0]);
  }else{
    m_fluxErr[0] = m_fluxErr[1];
  }
  return true;
}


void SpectrumHandler::integrateRate(){
  for(int i=m_nPoint-1;i>=0;i--){
    if(i==m_nPoint-1){
      m_integ[i] = m_rate[i];
    }else{
      m_integ[i] = m_integ[i+1]+m_rate[i];
    }
  }
}


double SpectrumHandler::getFlux(double E, bool err){
  /// should add err calculaiton
  int iE = binaryLowBoundSearch(E,m_nPoint,m_E);
  if(iE<0) return 0;
  if(iE==m_nPoint-1){
    if(E > m_E[iE]) return 0;
    return  (err) ? m_fluxErr[iE] : m_flux[iE];
  }
  if(m_flux[iE]==0 || m_flux[iE+1]==0)
    return  (err) ? m_fluxErr[iE] : m_flux[iE];
    
  double frat = m_flux[iE+1]/m_flux[iE];
  double erat = m_E[iE+1]/m_E[iE];
  double A = log(frat)/log(erat);
  double f = m_flux[iE]*pow(E/m_E[iE],A);
  if(!err) return f;
  double dfdfi = f/m_flux[iE] * (1-log(E/m_E[iE])/log(erat));
  double dfdfj = f/m_flux[iE+1] * log(E/m_E[iE])/log(erat);
  return sqrt(pow(dfdfi*m_fluxErr[iE],2) 
	      + pow(dfdfj*m_fluxErr[iE+1],2));
}


double SpectrumHandler::getIntegral(double E){
  int iE = binaryLowBoundSearch(E,m_nPoint,m_E);
  if(iE<0) return 0;
  if(iE==m_nPoint-1){
    if(E > m_E[iE]) return 0;
    return  m_integ[iE];
  }

  double frat = m_integ[iE+1]/m_integ[iE];
  double erat = m_E[iE+1]/m_E[iE];
  double A = log(frat)/log(erat);
  double f = m_integ[iE]*pow(E/m_E[iE],A);
  return f;
}

double SpectrumHandler::getRate(double E,bool err){
  int iE = binaryLowBoundSearch(E,m_nPoint,m_E);
  if(iE<0) return 0;
  if(iE==m_nPoint-1){
    if(E > m_E[iE]) return 0;
  }
  return  (err) ? m_rateErr[iE] : m_rate[iE];
}

double SpectrumHandler::getRandomE(double rand){ // rand: 0--1 
  double Imax = m_integ[0];
  double I = Imax*rand;
  int iE = binaryLowBoundSearch(I,m_nPoint,m_integ);
  
  double E;
  if(iE == m_nPoint-1){
    double y1 = I / Imax;
    double frat = m_flux[iE]/m_flux[iE-1];
    double erat = m_E[iE]/m_E[iE-1];
    double A = 1./(1.+log(frat)/log(erat));
    E = pow(y1,A); 
  }else{
    double frat = m_integ[iE+1]/m_integ[iE];
    double erat = m_E[iE+1]/m_E[iE];
    double A = log(frat)/log(erat);
    E = m_E[iE]*pow(I/m_integ[iE],1/A);
  }
  return E;
}


int binaryLowBoundSearch(double X,int N,double *val){
  if(N==0){
    std::cout<<"error in binaryLowBoundSearch: N==0"<<std::endl;
    exit(1);
  }
  bool ascending = (val[0]<val[N-1]) ? true : false;
  if(ascending){
    if(X < val[0]) return -1;
    if(X >= val[N-1]) return N-1;
  }else{
    if(X > val[0]) return -1;
    if(X <= val[N-1]) return N-1;
  }
  
  int left=0,right=N-1;
  int mid;
  while( left+1 < right ){
    mid = (left+right)/2;
    if(X==val[mid]){
      left = mid;
      break;
    }
    if(ascending){
      if(X>val[mid]){
	left = mid;
      }
      if(X<val[mid]){
	right = mid;
      }
    }else{
      if(X<val[mid]){
	left = mid;
      }
      if(X>val[mid]){
	right = mid;
      }
    }
  }
  return left;
}

void SpectrumHandler::printFlux(std::ofstream &ofs){
  ofs<<"Ek"<<"\t"<<"flux"<<"\t"<<"err"<<std::endl;
  for(int i=0;i<m_nPoint;i++){
    ofs<<m_E[i]<<"\t"<<m_flux[i]<<"\t"<<m_fluxErr[i]<<std::endl;
  }
}

void SpectrumHandler::printRate(std::ofstream &ofs, double factor){
  ofs<<"Ek"<<"\t"<<"rate"<<"\t"<<"err"<<std::endl;
  for(int i=0;i<m_nPoint;i++){
    ofs<<m_E[i]<<"\t"<<m_rate[i]*factor<<"\t"<<m_rateErr[i]*factor<<std::endl;
  }
}
