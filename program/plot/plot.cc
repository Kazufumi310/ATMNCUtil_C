#include <myheader/myRoot.h>
#include <myheader/myStandard.h>
#include "atmnc-constants.h"
#include "atmnc-particle-code.h"
#include "ContainerClass.h"
#include "Histogrammer.h"

void plotAverage(Histogrammer &hist,TFile &fout,TCanvas* canvas);
void plotAverageAzim(Histogrammer &hist,TFile &fout,TCanvas* canvas);
void plotVsTheta(Histogrammer &hist,TFile &fout,TCanvas* canvas);
void plotVsPhi(Histogrammer &hist,TFile &fout,TCanvas* canvas);
void writeEach(Histogrammer &hister,TFile &fout);

class ArgHandler{
public:
  ArgHandler();
  bool getArg(int &argc,char** argv);
  enum{
    vsE_average=1,vsE_averageAzim=2,vsTheta=4,vsPhi=8
  };
  
  bool isTargetPID(int pid);
  bool isTargetSite(int site);
  double getFactor(){return m_eFactor;}
  void getRange(double &emin,double &emax){emin = pow(10,m_eRange[0]);emax = pow(10,m_eRange[1]);}
  int getPlotFlag(int targetFlag);
  std::string getInputFile(){return m_inputFile;}
  std::string getOutputFile(){return m_outputFile;}
  void showUsage();

  bool isRatePlot(){return false;} // legacy
  
private:
  std::vector<int> m_targetPID;
  std::vector<int> m_targetSite;
  double m_eFactor;
  double m_eRange[2];
  int m_plotFlag;
  std::string m_inputFile;
  std::string m_outputFile;
};

ArgHandler::ArgHandler(){
  m_eFactor = 0;
  m_eRange[0] = -1;
  m_eRange[1] = 4;
  m_outputFile = "out.root";
  m_plotFlag = 13;
}

int ArgHandler::getPlotFlag(int targetFlag){
  // 1: vsE, average, 2: vs theat, 4: vs phi
  if(targetFlag==0) return targetFlag;
  if( (m_plotFlag & targetFlag)!=0 ) return 1;
  return 0;
}

bool ArgHandler::isTargetSite(int site){
  if(m_targetSite.size()!=0){
    bool find = false;
    for(int i=0;i<int(m_targetSite.size());i++){
      if(m_targetSite[i]==site){
	find=true;
	break;
      }
    }
    if(!find) return false;
  }
  return true;
}

void ArgHandler::showUsage(){
  std::cout<<"usage: bin/plot -i inputfile (-XXX option)"<<std::endl;
  std::cout<<"option:"<<std::endl;
  std::cout<<"\t -s : detector site ex) -s 0,1,2 (Note that the site ID begins from 0) "<<std::endl;
  std::cout<<"\t -k : PID ex) -k 19,20"<<std::endl;
  std::cout<<"\t -e : energy range ex) -e -1,4 (=10^-1 to 10^4)"<<std::endl;
  std::cout<<"\t -f : energy factor ex) -f 3 (=E^3 is multiplyed)"<<std::endl;
  std::cout<<"\t -i : input file name"<<std::endl;
  std::cout<<"\t -o : output file name"<<std::endl;
}

bool ArgHandler::isTargetPID(int pid){
  if(m_targetPID.size()!=0){
    bool find = false;
    for(int i=0;i<int(m_targetPID.size());i++){
      if(pid==m_targetPID[i]){
	find=true;
	break;
      }
    }
    if(!find) return false;
  }
  return true;
} 

bool ArgHandler::getArg(int &argc, char** argv){
  std::cout<<"command: ";
  for(int i=0;i<argc;i++){
    std::cout<<argv[i]<<" ";
  }
  std::cout<<std::endl;
  for(int iarg=1;iarg<argc;iarg++){
    std::string arg = argv[iarg];
    if(arg[0]=='-'){
      if(arg=="-h"){ // help
	showUsage();
	return false;
      }else{
	iarg++;
	std::string opt = argv[iarg];
	for(int i=0;i<opt.length();i++){
	  if(opt[i]==',') opt[i]=' ';
	}
	std::stringstream ss(opt);
      
	if(arg=="-s"){ // site 
	  int is;
	  while(ss>>is){
	    m_targetSite.push_back(is);
	  }
	  std::cout<<"target site:";
	  for(int i=0;i<int(m_targetSite.size());i++){
	    std::cout<<" "<<m_targetSite[i];
	  }
	  std::cout<<std::endl;
	}else if(arg=="-k"){ // PID of observed particle
	  int ik;
	  while(ss>>ik){
	    m_targetPID.push_back(ik);
	  }
	  std::cout<<"target particle:";
	  for(int i=0;i<int(m_targetPID.size());i++){
	    std::cout<<" "<<m_targetPID[i];
	  }
	  std::cout<<std::endl;
	}else if(arg=="-e"){ // energy range
	  double emin,emax;
	  if(! (ss>>emin>>emax) ){
	    std::cout<<"strange argument for option -e: "<<opt<<std::endl;
	    exit(1);
	  }else{
	    m_eRange[0] = emin;
	    m_eRange[1] = emax;
	    std::cout<<"energy range : "<<pow(10,m_eRange[0])<<" "<<pow(10,m_eRange[1])<<std::endl;
	  }
	}else if(arg=="-F"){
	  ss>>m_plotFlag;
	}else if(arg=="-f"){
	  ss>>m_eFactor;
	}else if(arg=="-i"){ // input file
	  ss>>m_inputFile;
	}else if(arg=="-o"){ // output file
	  ss>>m_outputFile;
	}else{
	  std::cout<<"argument error! "<<iarg<<"-th argument"<<std::endl;
	  showUsage();
	  return false;
	}
      }
    }else{
      std::cout<<"argument error! "<<iarg<<"-th argument"<<std::endl;
      showUsage();
      return false;
    }
  }

  if(m_inputFile.empty()){
    std::cout<<"argument error! set input file name with -i option"<<std::endl;
    return false;
  }

  
  std::cout<<"plotFlag:"
	   <<"vsE(average all):"<<getPlotFlag(vsE_average)<<" "
    	   <<"vsE(average phi):"<<getPlotFlag(vsE_averageAzim)<<" "
	   <<"vsTheta:"<<getPlotFlag(vsTheta)<<std::endl<<" "
	   <<"vsPhi:"<<getPlotFlag(vsPhi)<<std::endl;
  std::cout<<"additional factor: E^("<<m_eFactor<<")"<<std::endl;
  
  return true;
}

TH1D* s_histView = 0;
std::string outputPDF;
ArgHandler gArg;

////////////////
// main
///////////////
int main(int argc, char** argv){
  if(!gArg.getArg(argc,argv)) return 0;
  std::string inputRoot = gArg.getInputFile();
  std::string outputRoot = gArg.getOutputFile();
  outputPDF = outputRoot+".pdf";
  
  std::cout<<"input : "<<inputRoot<<std::endl;
  std::cout<<"output: "<<outputRoot<<std::endl;

  TFile fout(outputRoot.c_str(),"recreate");
  TCanvas *canvas = new TCanvas("canvas","");
  canvas->Print((outputPDF+"[").c_str());

  TChain chhist("tr");
  TChain chinfo("trinfo");
  chinfo.Add(inputRoot.c_str());
  chhist.Add(inputRoot.c_str());
  std::cout<<"read "<<inputRoot<<std::endl;
  
  ATMNCInputInfo inputInfo;
  inputInfo.setBranchAddress(&chinfo);
  ATMNCStatInfo statInfo;
  statInfo.setBranchAddress(&chinfo);
  ATMNCHistInfo histInfo;
  histInfo.setBranchAddress(&chhist);

  chinfo.GetEntry(0);
  Histogrammer hister(inputInfo);
  hister.setTequiv(statInfo.tequiv);
  
  int const nEntry = chhist.GetEntries();
  std::cout<<"nEntry:"<<nEntry<<std::endl;
  for(int ientry=0;ientry<nEntry;ientry++){
    //    if(ientry % int(nEntry/100.) == 0) std::cout<<ientry / int(nEntry/100.)<<"%"<<std::endl;
    chhist.GetEntry(ientry);
    for(int ik=0;ik<histInfo.nkind;ik++){
      hister.addHistogram(ik, histInfo.isite, histInfo.icosz, histInfo.iazim,
			  histInfo.nbin, histInfo.ebin, histInfo.nobs[ik], histInfo.nobsErr[ik]);
    }
  }

  s_histView = new TH1D("histView",";E[GeV];# events",hister.getNBin(),log10(hister.getEne(0)),log10(hister.getEne(hister.getNBin())));

  if(gArg.getPlotFlag(ArgHandler::vsE_average)){
    plotAverage(hister,fout,canvas);
    plotAverageAzim(hister,fout,canvas);
  }
  if(gArg.getPlotFlag(ArgHandler::vsTheta)){
    plotVsTheta(hister,fout,canvas);
  }
  if(gArg.getPlotFlag(ArgHandler::vsPhi)){
    plotVsPhi(hister,fout,canvas);
  }

  writeEach(hister,fout);
  
  canvas->Print((outputPDF+"]").c_str());
  fout.Close();
}


////////////////
//end of main
///////////////


void plotAverage(Histogrammer &hist,TFile &fout, TCanvas *canvas){
  std::cout<<"plot averaged flux ..."<<std::endl;
  int const nazim = hist.getNAzim();
  int const ncosz = hist.getNCosz();
  int const nsite = hist.getNSite();
  int const nkind = hist.getNKind();
  fout.cd();

  std::string yaxisTitle = (gArg.getFactor()==0)
    ? "flux [1/s/sr/m^{2}/GeV]"
    : Form("flux #times E^{%g} [GeV^{%g}/s/sr/m^{2}]",gArg.getFactor(),gArg.getFactor()-1);

  for(int is=0;is<nsite;is++){
    if(!gArg.isTargetSite(is)) continue;

    TGraph *gra[32]={0};
    double emin,emax;
    gArg.getRange(emin,emax);
    double mini = 1e30,maxi=-1e30;
    for(int ik=0;ik<nkind;ik++){
      int pid = hist.getPID(ik);
      if(!gArg.isTargetPID(pid))continue;
      
      s_histView->Reset();
      s_histView->SetName(Form("hAverage%d_%d",is,ik));
      s_histView->SetTitle(Form("# hits obserbed, average of all direction, detID=%d, PID=%d",is,pid));
      
      SpectrumHandler spec;
      hist.averageFlux(ik,is,0,ncosz-1,0,nazim-1,spec);
      int nn = spec.getN();
      double *E = spec.getEneArray();
      double *rate = spec.getRateArray();
      double *rateErr = spec.getRateErrArray();
      double scale = (gArg.isRatePlot()) ? 1 : hist.getFactorOfLastAverage();
      
      double flux[1024],fluxErr[1024];
      for(int i=0;i<nn;i++){
	double factor = pow(E[i],gArg.getFactor());
	flux[i] = factor*spec.getFlux(E[i]);
	fluxErr[i] = factor*spec.getFluxErr(E[i]);
	if(E[i] > emin && E[i] < emax){
	  if(flux[i] > 0 && mini>flux[i]) mini = flux[i];
	  if(flux[i] > 0 && maxi<flux[i]) maxi = flux[i];
	}
	s_histView->SetBinContent(i+1,rate[i]/scale);
	s_histView->SetBinContent(i+1,rateErr[i]/scale);
      }
      gra[ik] = new TGraphErrors(nn,E,flux,0,fluxErr);
      gra[ik]->SetName(Form("gAverage%d_%d",is,ik));
      std::string title = Form("flux, average of all direction, detID=%d, PID=%d;E[GeV];",is,pid)+yaxisTitle;
      gra[ik]->SetTitle(title.c_str());
      //    gra[ik]->SetLineColor(color[ik]);
      //    gra[ik]->SetMarkerColor(color[ik]);
      gra[ik]->SetLineColor(ik+1);
      gra[ik]->SetMarkerColor(ik+1);
      gra[ik]->SetMarkerStyle(20);
      gra[ik]->SetMarkerSize(0.7);
      
      s_histView->Write();
    }
    
    //    std::cout<<"mini,maxi:"<<mini<<" "<<maxi<<std::endl;
    mini = log10(mini);
    mini = (mini<0)? int(mini)-2. : int(mini)-1.;
    mini = pow(10,mini);
    maxi = log10(maxi);
    maxi = (maxi<0)? int(maxi)+1. : int(maxi)+2.;
    maxi = pow(10,maxi);
    std::cout<<"mini,maxi:"<<mini<<" "<<maxi<<std::endl;
    TH2F frame("frame",";E[GeV];flux #times E^{3} [GeV^{2}/s/sr/m^{2}]",100,emin,emax,100,mini,maxi);
    frame.GetYaxis()->SetTitle(yaxisTitle.c_str());
    frame.SetStats(0);
    
    canvas->SetLogx();
    canvas->SetLogy();
    for(int ik=0;ik<nkind;ik++){
      if(gra[ik]==0) continue;
      frame.SetTitle(Form("average all direction,site=%d,  PID=%d",is,hist.getPID(ik)));
      frame.Draw();
      gra[ik]->Draw("p");
      if( gArg.getPlotFlag(ArgHandler::vsE_average) ){
	canvas->Print(outputPDF.c_str());
      }
    }
    frame.SetTitle(Form("average all direction, site=%d",is));
    frame.Draw();
    for(int ik=0;ik<nkind;ik++){
      if(gra[ik]==0) continue;
      gra[ik]->Draw("p");
    }
    if( gArg.getPlotFlag(ArgHandler::vsE_average) ){
      canvas->Print(outputPDF.c_str());
    }
    
    canvas->SetLogx(0);
    canvas->SetLogy(0);
    
    for(int ik=0;ik<nkind;ik++){
      if(gra[ik]==0) continue;
      //      gra[ik]->SetName(Form("gVsE_pid%d_site%d",hist.getPID(ik),is));
      gra[ik]->Write();
    }
    
    for(int ik=0;ik<nkind;ik++){
      if(gra[ik]!=0)  delete gra[ik];
    }
  }
}

void plotAverageAzim(Histogrammer &hist,TFile &fout, TCanvas *canvas){
  std::cout<<"plot fluxes averaged around phi ..."<<std::endl;
  int const nazim = hist.getNAzim();
  int const ncosz = hist.getNCosz();
  int const nsite = hist.getNSite();
  int const nkind = hist.getNKind();


  fout.cd();


  std::string yaxisTitle = (gArg.getFactor()==0)
    ? "flux [1/s/sr/m^{2}/GeV]"
    : Form("flux #times E^{%g} [GeV^{%g}/s/sr/m^{2}]",gArg.getFactor(),gArg.getFactor()-1);

  for(int is=0;is<nsite;is++){
    if(!gArg.isTargetSite(is)) continue;

    TGraph *gra[32]={0};
    double emin,emax;
    gArg.getRange(emin,emax);

    for(int ic=0;ic<ncosz;ic++){
      double mini = 1e30,maxi=-1e30;
      for(int ik=0;ik<nkind;ik++){
	int pid = hist.getPID(ik);
	//	std::cout<<"ik,pid:"<<ik<<" "<<pid<<std::endl;
	if(!gArg.isTargetPID(pid))continue;

	s_histView->Reset();
	s_histView->SetName(Form("hAveAzim%d_%d_%d",is,ik,ic));
	s_histView->SetTitle(Form("# hits obserbed, average over azimuthal angle, detID=%d, PID=%d, cosz: %g-%g",is,pid,hist.getCosz(ic),hist.getCosz(ic+1)));
	
	SpectrumHandler spec;
	hist.averageFlux(ik,is,ic,ic,0,nazim-1,spec);
	int nn = spec.getN();
	double *E = spec.getEneArray();
	double *rate = spec.getRateArray();
	double *rateErr = spec.getRateErrArray();
	double scale = (gArg.isRatePlot()) ? 1 : hist.getFactorOfLastAverage();

	double flux[1024],fluxErr[1024];
	for(int i=0;i<nn;i++){
	  double factor = pow(E[i],gArg.getFactor());
	  flux[i] = factor*spec.getFlux(E[i]);
	  fluxErr[i] = factor*spec.getFluxErr(E[i]);
	  if(E[i] > emin && E[i] < emax){
	    if(flux[i] > 0 && mini>flux[i]) mini = flux[i];
	    if(flux[i] > 0 && maxi<flux[i]) maxi = flux[i];
	  }

	  s_histView->SetBinContent(i+1,rate[i]/scale);
	  s_histView->SetBinError(i+1,rateErr[i]/scale);
	}
	canvas->SetLogx(0);
	canvas->SetLogy();
  
	s_histView->Write();
	s_histView->Draw("e");
	//	Plotter::getPlotter()->logAxisLabel(0,s_histView);
	if( gArg.getPlotFlag(ArgHandler::vsE_averageAzim) ){
	  canvas->Print(outputPDF.c_str());
	}
	
	gra[ik] = new TGraphErrors(nn,E,flux,0,fluxErr);
	//    gra[ik]->SetLineColor(color[ik]);
	//    gra[ik]->SetMarkerColor(color[ik]);
	gra[ik]->SetName(Form("gAveAzim%d_%d_%d",is,ik,ic));
	std::string title = Form("flux, average over azimuthal angle, detID=%d, PID=%d, cosz: %g-%g; E[GeV];",is,pid,hist.getCosz(ic),hist.getCosz(ic+1)) + yaxisTitle;
	gra[ik]->SetTitle(title.c_str());

	gra[ik]->SetLineColor(ik+1);
	gra[ik]->SetMarkerColor(ik+1);
	gra[ik]->SetMarkerStyle(20);
	gra[ik]->SetMarkerSize(0.7);
	//	s_histView->Write();
      }
      
      std::cout<<"mini,maxi:"<<mini<<" "<<maxi<<std::endl;
      mini = log10(mini);
      mini = (mini<0)? int(mini)-2. : int(mini)-1.;
      mini = pow(10,mini);
      maxi = log10(maxi);
      maxi = (maxi<0)? int(maxi)+1. : int(maxi)+2.;
      maxi = pow(10,maxi);
      std::cout<<"mini,maxi:"<<mini<<" "<<maxi<<std::endl;
      TH2F frame("frame",";E[GeV];flux #times E^{3} [GeV^{2}/s/sr/m^{2}]",100,emin,emax,100,mini,maxi);

      frame.GetYaxis()->SetTitle(yaxisTitle.c_str());
      frame.SetStats(0);
      
      canvas->SetLogx();
      canvas->SetLogy();
      for(int ik=0;ik<nkind;ik++){
	if(gra[ik]==0) continue;
	frame.SetTitle(Form("cosz= %g -- %g, site=%d,PID=%d",hist.getCosz(ic),hist.getCosz(ic+1),is,hist.getPID(ik)));
	frame.Draw();
	gra[ik]->Draw("p");
	if( gArg.getPlotFlag(ArgHandler::vsE_averageAzim) ){
	  canvas->Print(outputPDF.c_str());
	}
      }
      frame.SetTitle(Form("cosz= %g -- %g,site=%d",hist.getCosz(ic),hist.getCosz(ic+1),is));
      frame.Draw();
      for(int ik=0;ik<nkind;ik++){
	if(gra[ik]==0) continue;
	gra[ik]->Draw("p");
      }
      if( gArg.getPlotFlag(ArgHandler::vsE_averageAzim) ){
	canvas->Print(outputPDF.c_str());
      }
      
      canvas->SetLogx(0);
      canvas->SetLogy(0);
      
      for(int ik=0;ik<nkind;ik++){
	if(gra[ik]==0) continue;
	//	gra[ik]->SetName(Form("gVsE_pid%d_site%d",hist.getPID(ik),is));
	gra[ik]->Write();
      }
      
      for(int ik=0;ik<nkind;ik++){
	if(gra[ik]!=0)  delete gra[ik];
      }
    }
  }
}


void plotVsTheta(Histogrammer &hist,TFile &fout, TCanvas *canvas){
  std::cout<<"plot flux vs theta ..."<<std::endl;
  double const targetE = 3.2; // GeV
  int const nazim = hist.getNAzim();
  int const ncosz = hist.getNCosz();
  int const nsite = hist.getNSite();
  int const nkind = hist.getNKind();

  fout.cd();
  for(int is=0;is<nsite;is++){
    if(!gArg.isTargetSite(is)) continue;
    TGraph *gra[32]={0};
    double mini=1e9,maxi=0;
    for(int ik=0;ik<nkind;ik++){
      int pid = hist.getPID(ik);
      if(!gArg.isTargetPID(pid))continue;
      std::cout<<"ik:"<<ik<<std::endl;
      
      int nn = 0;
      double cost[1024],flux[1024],fluxErr[1024];
      for(int ic=0;ic<ncosz;ic++){
	SpectrumHandler spec;
	hist.averageFlux(ik,is,ic,ic,0,nazim-1,spec);
	cost[nn] = hist.getCosz(ic+0.5);
	flux[nn] = spec.getFlux(targetE);
	fluxErr[nn] = spec.getFluxErr(targetE);
	if(flux[nn] > 0 && mini>flux[nn]) mini = flux[nn];
	if(flux[nn] > 0 && maxi<flux[nn]) maxi = flux[nn];

	nn++;
      }
      gra[ik] = new TGraphErrors(nn,cost,flux,0,fluxErr);
      //    gra[ik]->SetLineColor(color[ik]);
      //    gra[ik]->SetMarkerColor(color[ik]);
      gra[ik]->SetLineColor(ik+1);
      gra[ik]->SetMarkerColor(ik+1);
      gra[ik]->SetMarkerStyle(20);
      gra[ik]->SetMarkerSize(0.7);
    }
    std::cout<<"mini,maxi:"<<mini<<" "<<maxi<<std::endl;
    mini = log10(mini);
    mini = (mini<0)? int(mini)-2. : int(mini)-1.;
    mini = pow(10,mini);
    maxi = log10(maxi);
    maxi = (maxi<0)? int(maxi)+1. : int(maxi)+2.;
    maxi = pow(10,maxi);
    std::cout<<"mini,maxi:"<<mini<<" "<<maxi<<std::endl;
  
    TH2F frame("frame",";cos#theta;flux [1/GeV/s/sr/m^{2}]",100,-1,1,100,mini,maxi);
    frame.SetStats(0);

    canvas->SetLogy();
    for(int ik=0;ik<nkind;ik++){
      if(gra[ik]==0) continue;
      frame.SetTitle(Form("site=%d,PID=%d",is,hist.getPID(ik)));
      frame.Draw();
      gra[ik]->Draw("p");
      if( gArg.getPlotFlag(ArgHandler::vsTheta) ){
	canvas->Print(outputPDF.c_str());
      }
    }
  
    frame.SetTitle("");
    frame.Draw();
    for(int ik=0;ik<nkind;ik++){
      if(gra[ik]==0) continue;
      gra[ik]->Draw("p");
    }
    if( gArg.getPlotFlag(ArgHandler::vsTheta) ){
      canvas->Print(outputPDF.c_str());
    }

    canvas->SetLogx(0);
    canvas->SetLogy(0);

    for(int ik=0;ik<nkind;ik++){
      if(gra[ik]==0) continue;
      gra[ik]->SetName(Form("gVsTheta_pid%d_site%d",hist.getPID(ik),is));
      gra[ik]->Write();
    }
  
    for(int ik=0;ik<nkind;ik++){
      if(gra[ik]!=0) delete gra[ik];
    }
  }
}


void plotVsPhi(Histogrammer &hist,TFile &fout, TCanvas *canvas){
  std::cout<<"plot flux vs phi ..."<<std::endl;
  double const targetE = 3.2; // GeV
  int const nazim = hist.getNAzim();
  int const ncosz = hist.getNCosz();
  int const nsite = hist.getNSite();
  int const nkind = hist.getNKind();

  int const nCoszRegion = 5;
  double const dCosz = 0.4;
  int const color[nCoszRegion]={1,4,3,7,2};
  
  canvas->SetLogx(0);
  canvas->SetLogy();
  fout.cd();
  for(int is=0;is<nsite;is++){
    if(!gArg.isTargetSite(is)) continue;
    for(int ik=0;ik<nkind;ik++){
      int pid = hist.getPID(ik);
      if(!gArg.isTargetPID(pid))continue;
      TGraph *gra[nCoszRegion];
      double mini=1e9,maxi=0;
      for(int ic=0;ic<nCoszRegion;ic++){
	int icmin = hist.getCoszIndex(acos(1-dCosz*ic));
	int icmax = hist.getCoszIndex(acos(1-dCosz*(ic+1)))-1;
	double azim[1024],flux[1024],fluxErr[1024];
	int nn = 0;
	for(int ia=0;ia<nazim;ia++){
	  SpectrumHandler spec;
	  hist.averageFlux(ik,is,icmin,icmax,ia,ia,spec);
	  azim[nn] = hist.getAzim(ia+0.5)*TMath::RadToDeg();
	  flux[nn] = spec.getFlux(targetE);
	  fluxErr[nn] = spec.getFluxErr(targetE);
	  if(flux[nn] > 0 && mini>flux[nn]) mini = flux[nn];
	  if(flux[nn] > 0 && maxi<flux[nn]) maxi = flux[nn];
	  nn++;
	}
	gra[ic] = new TGraphErrors(nn,azim,flux,0,fluxErr);
	gra[ic]->SetLineColor(color[ic]);
	gra[ic]->SetMarkerColor(color[ic]);
	gra[ic]->SetMarkerStyle(20);
	gra[ic]->SetMarkerSize(0.7);
      }
      std::cout<<"mini,maxi:"<<mini<<" "<<maxi<<std::endl;
      mini = log10(mini);
      mini = (mini<0)? int(mini)-2. : int(mini)-1.;
      mini = pow(10,mini);
      maxi = log10(maxi);
      maxi = (maxi<0)? int(maxi)+1. : int(maxi)+2.;
      maxi = pow(10,maxi);
      std::cout<<"mini,maxi:"<<mini<<" "<<maxi<<std::endl;
      
      TH2F frame("frame",";azim. angle [deg];flux [1/GeV/s/sr/m^{2}]",100,0,360,100,mini,maxi);
      frame.SetStats(0);
      
      for(int ic=0;ic<nCoszRegion;ic++){
	frame.SetTitle(Form("site=%d, PID=%d, cos#theta=%g-%g",is,hist.getPID(ik),1-dCosz*ic,1-dCosz*(ic+1)));
	frame.Draw();
	gra[ic]->Draw("lp");
	if( gArg.getPlotFlag(ArgHandler::vsPhi) ){
	  canvas->Print(outputPDF.c_str());
	}
      }
      
      frame.Draw();
      for(int ic=0;ic<nCoszRegion;ic++){
	gra[ic]->Draw("pl");
      }
      if( gArg.getPlotFlag(ArgHandler::vsPhi) ){
	canvas->Print(outputPDF.c_str());
      }
      for(int ic=0;ic<nCoszRegion;ic++){
	gra[ic]->SetName(Form("gVsPhi_kind%d_site%d_cos%d",ik,is,ic));
	gra[ic]->Write();
      }
      
      for(int ic=0;ic<nCoszRegion;ic++) delete gra[ic];
    }  
  } 
  canvas->SetLogx(0);
  canvas->SetLogy(0);
}


void writeEach(Histogrammer &hister,TFile &fout){
  std::cout<<"write graphs of all bins in TTree..."<<std::endl;
  fout.cd();

  TTree *trout = new TTree("histtr","each histogram and flux");
  TH1D* hist = new TH1D("histEach",";log10(E[GeV]);# events",hister.getNBin(),log10(hister.getEne(0)),log10(hister.getEne(hister.getNBin())));
  TGraphErrors* gra = 0;
  //  TH1D* hflux = new TH1D("hfluxView",";log10(E[GeV]);flux ",nbin,log10(emin),log10(emax));
  int isite,icosz,iazim,ikind,pid;
  double cmin,amin;
  trout->Branch("hist",hist);
  trout->Branch("flux",&gra);
  trout->Branch("ikind",&ikind,"ikind/I");
  trout->Branch("pid",&pid,"pid/I");
  trout->Branch("isite",&isite,"isite/I");
  trout->Branch("icosz",&icosz,"icosz/I");
  trout->Branch("cmin",&cmin,"cmin/D");
  trout->Branch("iazim",&iazim,"iazim/I");
  trout->Branch("amin",&amin,"amin/D");

  
  std::string yaxisTitle = (gArg.getFactor()==0)
    ? "flux [1/s/sr/m^{2}/GeV]"
    : Form("flux #times E^{%g} [GeV^{%g}/s/sr/m^{2}]",gArg.getFactor(),gArg.getFactor()-1);

  int const nazim = hister.getNAzim();
  int const ncosz = hister.getNCosz();
  int const nsite = hister.getNSite();
  int const nkind = hister.getNKind();

  for(int is=0;is<nsite;is++){
    isite = is;
    if(!gArg.isTargetSite(is)) continue;
    double emin,emax;
    gArg.getRange(emin,emax);
    for(int ik=0;ik<nkind;ik++){    
      ikind = ik;
      pid = hister.getPID(ik);
      for(int ic=0;ic<ncosz;ic++){
	icosz = ic;
	cmin = hister.getCosz(ic);
	for(int ia=0;ia<nazim;ia++){
	  iazim = ia;
	  amin = hister.getAzim(ia);
	  hist->Reset();
	  hist->SetName(Form("hEach%d_%d_%d_%d",is,ik,ic,ia));
	  hist->SetTitle(Form("# hits obserbed, detID=%d, PID=%d, cosz: %g-%g, azim: %g-%g",is,pid,cmin,hister.getCosz(ic+1),amin,hister.getAzim(ia+1)));
	
	  SpectrumHandler spec;
	  hister.averageFlux(ik,is,ic,ic,ia,ia,spec);
	  int nn = spec.getN();
	  double *E = spec.getEneArray();
	  double *rate = spec.getRateArray();
	  double *rateErr = spec.getRateErrArray();
	  double scale = (gArg.isRatePlot()) ? 1 : hister.getFactorOfLastAverage();
	  
	  double flux[1024],fluxErr[1024];
	  double mini = 1e9,maxi = 0;
	  for(int i=0;i<nn;i++){
	    double factor = pow(E[i],gArg.getFactor());
	    flux[i] = factor*spec.getFlux(E[i]);
	    fluxErr[i] = factor*spec.getFluxErr(E[i]);
	    if(E[i] > emin && E[i] < emax){
	      if(flux[i] > 0 && mini>flux[i]) mini = flux[i];
	      if(flux[i] > 0 && maxi<flux[i]) maxi = flux[i];
	    }
	    
	    hist->SetBinContent(i+1,rate[i]/scale);
	    hist->SetBinError(i+1,rateErr[i]/scale);
	  }
	  
	  gra = new TGraphErrors(nn,E,flux,0,fluxErr);
	  gra->SetName(Form("gEach%d_%d_%d_%d",is,ik,ic,ia));
	  std::string title = Form("flux, detID=%d, PID=%d, cosz: %g-%g, azim: %g-%g; E[GeV];",is,pid,hister.getCosz(ic),hister.getCosz(ic+1),hister.getAzim(ia),hister.getAzim(ia+1)) + yaxisTitle;
	  gra->SetTitle(title.c_str());
	  gra->SetMarkerStyle(20);
	  gra->SetMarkerSize(0.7);
	  
	  trout->Fill();
	  
	  delete gra;
	}
      }
    }
  }
  trout->Write();
  delete hist;
  delete trout;
}

