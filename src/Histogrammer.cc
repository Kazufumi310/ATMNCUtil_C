#include "Histogrammer.h"


Histogrammer::Histogrammer(ATMNCInputInfo input)
{
  m_factorOfLastAverage = 1.;
  
  m_nkind = input.nkind;//default = 6;
  m_nsite = input.nsite;//default = 8;
  m_ncosz = input.nczo;//default = 40;
  m_nazim = input.nazo;//default = 12;
  m_nbin = 103;// default = 103; variable bin number is not yet implemented.
  std::cout<<"Histogrammer: # of histograms:"<<std::endl;
  int nall = m_nkind*m_nsite*m_ncosz*m_nazim;
  std::cout<<"nkind x nsite x ncos x nazim = "<<m_nkind<<" x "<<m_nsite<<" x "<<m_ncosz<<" x "<<m_nazim<<" = "<<nall<<std::endl;
  m_nobs = new double*[nall];
  m_err2 = new double*[nall];
  for(int i=0;i<nall;i++){
    m_nobs[i] = new double[m_nbin];
    m_err2[i] = new double[m_nbin];
    for(int j=0;j<m_nbin;j++){
      m_nobs[i][j] = 0;
      m_err2[i][j] = 0;
    }
  }
  
  m_ene =  new double[m_nbin+1];
  m_emin = input.etmi;
  m_Letmi = int(20*log10(m_emin) + 100.5) - m_nbin + 1;  // +1 is needed to adjust C to Frotran
  for(int ibin=0;ibin<m_nbin+1;ibin++){
    //    m_ene[ibin] = pow(10,(ibin+Letmi)/20.);
    m_ene[ibin] = pow(10,(ibin + m_Letmi)/20.);
  }
  
  m_cmin = input.comi;
  //  m_cmax = 1; 
  
  m_tequiv = 1;
  
  for(int i=0;i<ParticleCode::nParType;i++) m_pid2ik[i] = -1; 
  for(int ik=0;ik<m_nkind;ik++){
    m_pid2ik[input.obsp[ik]] = ik;
  }
  std::cout<<"# of detector: "<<m_nsite<<std::endl;
  for(int is=0;is<m_nsite;is++){
    double r = input.siteParam[is][3]+Constants::rEarth;
    double lat = input.siteParam[is][0];
    double deg = input.siteParam[is][2];
    if(deg>0.01 && lat<=90.1){
      double theta = deg/TMath::RadToDeg();
      m_Darea[is] = 2*TMath::Pi()*pow(r,2)*(1-cos(theta));
    }else{
      double rinj = input.hinj + Constants::rEarth;
      double czpr = cos(input.rsmp/TMath::RadToDeg());
      m_Darea[is] = 2*TMath::Pi()*pow(rinj,2)*(1-czpr);
    }
    std::cout<<"isite, Darea:"<<is<<" "<<m_Darea[is]<<std::endl;
  }
  
  for(int ic=0;ic<m_ncosz;ic++){
    double c1 = getCosz(ic);
    double c2 = getCosz(ic+1);
    m_Domega[ic] = (c1*c2 >0)
      ? TMath::Pi()*fabs(c1*c1 - c2*c2)/m_nazim
      : TMath::Pi()*(c1*c1 + c2*c2)/m_nazim;
    //    std::cout<<"ic,c1,c2,Domega:"<<ic<<" "<<c1<<" "<<c2<<" "<<m_Domega[ic]<<std::endl;
  }
}


void Histogrammer::getHistogram(int ikind, int isite, int icosz, int iazim,
				int &num, double ene[], double nobs[],
				double err[]){
  num = m_nbin;
  int I = getI(ikind,isite,icosz,iazim);
  
  for(int i=0;i<m_nbin;i++){
    ene[i] = m_ene[i];
    nobs[i] = m_nobs[I][i];
    err[i] = sqrt(m_err2[I][i]);
  }
}

void Histogrammer::addHistogram(int ikind, int isite, int icosz, int iazim,
				int num, double ene[], double nobs[],
				double err[]){
  if(num != m_nbin){
    std::cout<<"error in addHistogram: The histogram with different bin division is attempted to add."<<std::endl;
    exit(1);
  }
  int I = getI(ikind,isite,icosz,iazim);
  for(int i=0;i<m_nbin;i++){
    if(fabs(ene[i]- m_ene[i])/m_ene[i]>1e-3){
    }
    m_nobs[I][i] += nobs[i];
    m_err2[I][i] += pow(err[i],2);
  }
}

void Histogrammer::setTequiv(TTree* trinfo){
  int nentry = trinfo->GetEntries();
  double tequiv = 0;
  ATMNCStatInfo output;
  output.setBranchAddress(trinfo);
  for(int ientry=0;ientry<nentry;ientry++){
    trinfo->GetEntry(ientry);
    tequiv += output.tequiv;
  }
  setTequiv(tequiv);
}


void Histogrammer::addEvent(HitInfo &info, double weight){
  for(int ihit = 0; ihit< info.nhit; ihit++){
    int pid = info.khit[ihit];
    int ik = m_pid2ik[pid];

    TVector3 mom(info.phit[ihit]);
    TVector3 pos(info.rhit[ihit]);
    int is = info.site[ihit] - 1; // 1 in fortran is 0 in c++

    if(ik==-1){
      std::cout<<"error! pid:"<<pid<<std::endl; 
      exit(0);
    }
    
    double theta,phi;
    getDirection(mom,pos,theta,phi);
    int ic = getCoszIndex(theta);
    if( ic>=m_ncosz || ic<0 ){
      continue;
    }
    int ia = getAzimIndex(phi);
    if( ia>=m_nazim || ia<0 ){
      continue;
    }
    
    int id = getI(ik,is,ic,ia);
    
    double kineVal = mom.Mag();
    if(pid == ParticleCode::kpro || pid == ParticleCode::kprobar || pid == ParticleCode::kneut || pid == ParticleCode::kneutbar ){
      double m = ParticleCode::mass[pid];
      kineVal = sqrt(mom.Mag2()+pow(m,2)) - m;
    }
    int ibin = getEneIndex(kineVal);

    if(ibin<m_nbin && ibin >= 0){
      double w = (weight<0)?1:weight;
      m_nobs[id][ibin] += w;
      m_err2[id][ibin] += w*w;
    }
  }
}




double Histogrammer::getNormFactor(int is,int ic){
  double facInv = m_tequiv*m_Darea[is]*m_Domega[ic];
  if(facInv == 0){
    static bool alreadyWarn=false;
    if(!alreadyWarn){
      alreadyWarn=true;
      std::cout<<"Histogrammer::getNormaFactor: not a proper normalization factor"<<std::endl
	       <<"is, ic: "<<is<<" "<<ic<<std::endl
	       <<"  tequiv, Darea, Domega:"<<m_tequiv<<" "<<m_Darea[is]<<" "<<m_Domega[ic]<<std::endl
	       <<"the factor is set to 1."<<std::endl;
    }
    return 1;
  }
  return 1./facInv;
}


void Histogrammer::getSpectrum(int ik,int is,int ic,int ia, SpectrumHandler &spec){
  
  //  double denomi = 0;
  //  for(int ic=icmin;ic<=icmax;ic++){
  //    double dom = 0;
  //    for(int ia=iamin;ia<=iamax;ia++){
  int I = getI(ik,is,ic,ia);

  double *num = new double[m_nbin];
  double *nume = new double[m_nbin];
  for(int ibin=0;ibin<m_nbin;ibin++){
    num[ibin] = m_nobs[I][ibin];
    nume[ibin] = sqrt(m_err2[I][ibin]);
  }
  // change #observed to rate
  //  double factor = 1./(m_tequiv*m_Darea[is]*m_Domega[ic]);
  double factor = getNormFactor(is,ic);
  for(int ibin=0;ibin<m_nbin;ibin++){
    num[ibin]  *= factor;
    nume[ibin] *= factor;
  }

  // rate was set to SpectrumHandler. flux is calculated inside it.
  spec.setRate(m_nbin,m_ene,num,nume);
  
  delete [] num;
  delete [] nume;
}


void Histogrammer::averageFlux(int ik,int is,int icmin,int icmax,int iamin,int iamax, SpectrumHandler &spec){
  double *flux = 0, *fluxe = 0;
  double *num = 0, *nume = 0;
  flux = new double[m_nbin];
  fluxe = new double[m_nbin];
  num = new double[m_nbin];
  nume = new double[m_nbin];
  for(int i=0;i<m_nbin;i++){
    flux[i] = fluxe[i] = 0;
  }
  if(icmin<0||iamin<0||icmax>=m_ncosz||iamax>=m_nazim){
    std::cout<<"Histogrammer::averageFlux() : ERROR: invalid range"<<std::endl;
    std::cout<<"required icmin,icmax,iamin,iamax:"<<icmin<<" "<<icmax<<" "<<iamin<<" "<<iamax<<std::endl;
    exit(1);
  }
  
  double denomi = 0;
  for(int ic=icmin;ic<=icmax;ic++){
    for(int i=0;i<m_nbin;i++){
      num[i] = nume[i] = 0;
    }
    // averaging over azimuthal angle
    double dom = 0;
    for(int ia=iamin;ia<=iamax;ia++){
      int I = getI(ik,is,ic,ia);
      dom += m_Domega[ic];
      for(int ibin=0;ibin<m_nbin;ibin++){
        num[ibin] += m_nobs[I][ibin];
        nume[ibin] += m_err2[I][ibin];
	if(TMath::IsNaN(num[ibin]) || TMath::IsNaN(nume[ibin])){
	  std::cout<<"Histogrammer::averageFlux: there is NaN : ibin,ia,num,nume,m_nobs, m_err2: "<<ibin<<" "<<ia<<" "<<num[ibin]<<" "<<nume[ibin]<<" "<<m_nobs[I][ibin]<<" "<< m_err2[I][ibin]<<std::endl;
	  exit(1);
	}
      }
    }
    for(int ibin=0;ibin<m_nbin;ibin++){
      nume[ibin]=sqrt(nume[ibin]);
    }

    // change #observed to rate
    double factor = 1./(m_tequiv*m_Darea[is]*dom);
    m_factorOfLastAverage = factor;
    //double factor = 1.;
    for(int ibin=0;ibin<m_nbin;ibin++){
      num[ibin]  *= factor;
      nume[ibin] *= factor;

      if(TMath::IsNaN(num[ibin]) || TMath::IsNaN(nume[ibin])){
	std::cout<<"Histogrammer::averageFlux: there is NaN : ibin,num,nume,factor: "<<ibin<<" "<<num[ibin]<<" "<<nume[ibin]<<" "<<factor<<std::endl;
	exit(1);
      }
    }
    // rate was converted to flux by using SpectrumHandler
    spec.setRate(m_nbin,m_ene,num,nume);
    
    // averaging the flux over zenith angle
    double integFactor = getCosz(ic) - getCosz(ic+1);
    denomi += integFactor;
    
    //    for(int ibin=1;ibin<m_nbin;ibin++){
    for(int ibin=0;ibin<m_nbin;ibin++){ // 260516 modified, 
      double E = spec.bin2Ene(ibin);
      flux[ibin] += spec.getFlux(E) * integFactor;
      fluxe[ibin] += pow(spec.getFluxErr(E) * integFactor,2);
    }
  }
  
  //  for(int ibin=1;ibin<m_nbin;ibin++){
  for(int ibin=0;ibin<m_nbin;ibin++){ // 260516 modified
    flux[ibin] /= denomi;
    fluxe[ibin] = sqrt(fluxe[ibin]);
    fluxe[ibin] /= denomi;
    
  }
  // 260516 comment out
  /*
    flux[0] = flux[1];
    fluxe[0] = fluxe[1];
  */
  //  spec.setFlux(m_nbin,m_ene,flux,fluxe);
  spec.setFlux(m_nbin,m_ene,flux,fluxe,false);
  
  
  delete flux;
  delete fluxe;
  delete num;
  delete nume;
}



//////////////////////////
/////////////////////////
// private functions
/////////////////////////


int Histogrammer::getI(int ik,int is,int ic,int ia){
  return ia + m_nazim*(ic + m_ncosz*(is + m_nsite*ik));
}

int Histogrammer::getEneIndex(double ene){
  //  double Letmi = int(20*log10(m_emin) + 100.5) - m_nbin;
  int iE = int(20*log10(ene) + 100) - 100 - (m_Letmi + 1);
  return iE;
}

double Histogrammer::getCosz(double ic){
  double const d = (1-m_cmin)/m_ncosz;
  return 1 - ic*d;
}
double Histogrammer::getAzim(double ia){
  double const d = 2*TMath::Pi()/m_nazim;
  return ia*d;
}
int Histogrammer::getCoszIndex(double theta){
  double const d = (1-m_cmin)/m_ncosz;
  return int((1-cos(theta))/d);
}
int Histogrammer::getAzimIndex(double phi){
  if(phi<0) phi += 2*TMath::Pi();
  double const d = 2*TMath::Pi()/m_nazim;
  return int(phi/d);
}
int Histogrammer::getPID(int ik){
  for(int i=1;i<ParticleCode::nParType;i++){
    if(m_pid2ik[i] == ik){
      return i;
    }
  }
  return -1;
}


void Histogrammer::getDirection(TVector3 const &mom, TVector3 const &pos,
				double &thetaDir,double &phiDir){
  double rphi = pos.Phi();
  double rtheta = pos.Theta();
  double ct = cos(rtheta);
  double st = sin(rtheta);
  double cp = cos(rphi);
  double sp = sin(rphi);
  
  double rotX[3] = {-ct*cp,-ct*sp,st};
  double rotY[3] = {sp,-cp,0};
  double rotZ[3] = {-st*cp,-st*sp,-ct};

  double dirX = rotX[0]*mom[0]+rotX[1]*mom[1]+rotX[2]*mom[2];
  double dirY = rotY[0]*mom[0]+rotY[1]*mom[1]+rotY[2]*mom[2];
  double dirZ = rotZ[0]*mom[0]+rotZ[1]*mom[1]+rotZ[2]*mom[2];
  TVector3 dir(dirX,dirY,dirZ);
  thetaDir = dir.Theta();
  phiDir = dir.Phi();
}
